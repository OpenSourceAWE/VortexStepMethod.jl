"""
    NeuralFoilSolver

NeuralFoil backend. Evaluates the deformed section's Kulfan parameters through the
network (vectorized over angle) and reconstructs the surface pressure from the
predicted edge-velocity ratios. Fast, differentiable, and guarded by
`analysis_confidence` (carried into `SectionSolution.confidence`).

# Fields
- `model_size`: network size (default "large").
- `n_crit`: critical amplification factor (default 9.0).
- `xtr_upper`, `xtr_lower`: forced transition locations (default 1.0).
- `weights_dir`: override for the weights directory (default package data).
"""
@with_kw struct NeuralFoilSolver <: AbstractAirfoilSolver
    model_size::String = "large"
    n_crit::Float64 = 9.0
    xtr_upper::Float64 = 1.0
    xtr_lower::Float64 = 1.0
    weights_dir::Union{Nothing,String} = nothing
end

"""
    analyze_sweep(solver::NeuralFoilSolver, def, alpha_range, Re) -> Vector{SectionSolution}

Evaluate all angles (radians) in one vectorized NeuralFoil call on the deformed
Kulfan parameters, then assemble a per-node `cp` on the deformed contour nodes
(`def.x`, `def.y`) from the predicted edge velocities
([`neuralfoil_contour_solution`](@ref)). `cf` is a flat-plate closure
([`flat_plate_cf`](@ref)), NeuralFoil not exposing skin friction.
"""
function analyze_sweep(solver::NeuralFoilSolver, def::DeformedSection, alpha_range, Re)
    res = neuralfoil_section(def.kulfan, rad2deg.(collect(alpha_range)), Re;
        model_size=solver.model_size, weights_dir=solver.weights_dir,
        n_crit=solver.n_crit, xtr_upper=solver.xtr_upper, xtr_lower=solver.xtr_lower)
    x, y = collect(float.(def.x)), collect(float.(def.y))
    le = argmin(x)
    arc = contour_arc(x, y, le)
    cf = [flat_plate_cf(clamp(xk, 0.0, 1.0), Re) for xk in x]
    return [neuralfoil_contour_solution(alpha_range[i], res, i, x, y, le, arc, cf)
            for i in eachindex(alpha_range)]
end

"""
    contour_arc(x, y, le) -> Vector

Signed arc length of every contour node from the leading edge node `le`: positive
along the upper surface toward the trailing edge, negative along the lower. The
natural coordinate across a blunt nose, where a small chordwise step is a long one
along the skin.
"""
function contour_arc(x, y, le)
    arc = zeros(Float64, length(x))
    for k in (le - 1):-1:1
        arc[k] = arc[k + 1] + hypot(x[k] - x[k + 1], y[k] - y[k + 1])
    end
    for k in (le + 1):length(x)
        arc[k] = arc[k - 1] - hypot(x[k] - x[k - 1], y[k] - y[k - 1])
    end
    return arc
end

"""
    arc_at_chord(x, arc, indices, fractions) -> Vector

Signed arc length at each chord fraction, read off the contour nodes `indices` — the
one surface the fractions belong to.
"""
function arc_at_chord(x, arc, indices, fractions)
    nodes = collect(indices)
    order = sortperm(x[nodes])
    xs, as = x[nodes][order], arc[nodes][order]
    return [begin
        j = clamp(searchsortedfirst(xs, f), 2, length(xs))
        gap = max(xs[j] - xs[j - 1], eps())
        as[j - 1] + (f - xs[j - 1]) / gap * (as[j] - as[j - 1])
    end for f in fractions]
end

"""
    trailing_edge_speed(res, i, upper_arc, lower_arc, upper_te, lower_te) -> Float64

One edge speed for both surfaces at the trailing edge: each surface's last two
stations extrapolated to its own trailing-edge arc length, the two magnitudes
averaged, so a sharp edge carries a single pressure (Kutta).
"""
function trailing_edge_speed(res, i, upper_arc, lower_arc, upper_te, lower_te)
    reach(a1, a2, v1, v2, target) =
        v2 + (target - a2) / (a2 - a1 + eps()) * (v2 - v1)
    n = length(upper_arc)
    n < 2 && return abs(res.ue_upper[end, i])
    upper = reach(upper_arc[n - 1], upper_arc[n],
                  res.ue_upper[n - 1, i], res.ue_upper[n, i], upper_te)
    lower = reach(lower_arc[n - 1], lower_arc[n],
                  res.ue_lower[n - 1, i], res.ue_lower[n, i], lower_te)
    return (abs(upper) + abs(lower)) / 2
end

"""
    velocity_knots(res, i, x, arc, le) -> (knots, values)

The whole contour's edge velocity as one curve of (signed arc length, signed
`ue/u∞`), strictly increasing in arc length.

NeuralFoil reports nothing over the first and last `1/2N` of chord. Signing the lower
surface negative joins both surfaces into one curve through the stagnation point, so
the unsampled nose is interpolated between the innermost station on each side rather
than extrapolated off the end of one. The stagnation point enters as its own knot,
where `ue` interpolated linearly between those two stations reaches zero — the
stagnation-point-flow result, `ue` being linear in arc length there.
"""
function velocity_knots(res, i, x, arc, le)
    upper_arc = arc_at_chord(x, arc, 1:le, res.x)
    lower_arc = arc_at_chord(x, arc, le:length(x), res.x)
    speed = trailing_edge_speed(res, i, upper_arc, lower_arc, arc[1], arc[end])
    lower_first, upper_first = abs(res.ue_lower[1, i]), abs(res.ue_upper[1, i])
    stagnation = lower_arc[1] + lower_first / max(lower_first + upper_first, eps()) *
                 (upper_arc[1] - lower_arc[1])
    knots = [arc[end]; reverse(lower_arc); stagnation; upper_arc; arc[1]]
    values = [-speed; reverse(-abs.(res.ue_lower[:, i])); 0.0;
              abs.(res.ue_upper[:, i]); speed]
    order = sortperm(knots)
    knots, values = knots[order], values[order]
    keep = [1; [k for k in 2:length(knots) if knots[k] > knots[k - 1] + 1e-9]]
    return knots[keep], values[keep]
end

"""
    neuralfoil_contour_solution(alpha, res, i, x, y, le, arc, cf) -> SectionSolution

Assemble a full-contour [`SectionSolution`](@ref) for case `i` of a NeuralFoil sweep
`res`: build the one-curve edge velocity ([`velocity_knots`](@ref)), interpolate it
with a shape-preserving monotone cubic in arc length, and square it into
`Cp = 1 - ue²` at every contour node. `ue` is the interpolated quantity because it is
linear in arc length through a stagnation point, where `Cp` is quadratic.
"""
function neuralfoil_contour_solution(alpha, res, i, x, y, le, arc, cf)
    knots, values = velocity_knots(res, i, x, arc, le)
    speed = interpolate(knots, values, FritschButlandMonotonicInterpolation())
    cp = [1 - speed(clamp(a, first(knots), last(knots)))^2 for a in arc]
    return SectionSolution(alpha, res.cl[i], res.cd[i], res.cm[i],
                           res.confidence[i], x, y, cp, cf)
end

"""
    analyze_section(solver::NeuralFoilSolver, def, alpha, Re) -> SectionSolution

Single-angle convenience; delegates to [`analyze_sweep`](@ref).
"""
function analyze_section(solver::NeuralFoilSolver, def::DeformedSection, alpha, Re)
    return analyze_sweep(solver, def, [alpha], Re)[1]
end

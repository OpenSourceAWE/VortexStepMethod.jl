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
contour_arc(x, y, le) = contour_arc!(zeros(Float64, length(x)), x, y, le)

"""
    contour_arc!(arc, x, y, le) -> arc

[`contour_arc`](@ref) written into `arc`, which is resized to the contour. A live
pressure refresh walks panel after panel through one buffer this way.
"""
function contour_arc!(arc::Vector{Float64}, x, y, le)
    resize!(arc, length(x))
    arc[le] = 0.0
    for k in (le - 1):-1:1
        arc[k] = arc[k + 1] + hypot(x[k] - x[k + 1], y[k] - y[k + 1])
    end
    for k in (le + 1):length(x)
        arc[k] = arc[k - 1] - hypot(x[k] - x[k - 1], y[k] - y[k - 1])
    end
    return arc
end

"""
    ContourPressureScratch()

The buffers a surface-pressure reconstruction reuses, so [`contour_pressure!`](@ref)
allocates nothing per panel or per frame. Each grows to whatever contour and station
count it is first handed and stays that size.
"""
struct ContourPressureScratch
    "One surface's node chord fractions, ascending."
    surface_x::Vector{Float64}
    "The same nodes' signed arc length."
    surface_arc::Vector{Float64}
    "Arc length of the network's stations on the upper surface."
    upper_arc::Vector{Float64}
    "Arc length of the network's stations on the lower surface."
    lower_arc::Vector{Float64}
    "Arc length of every knot of the joined edge-velocity curve."
    knots::Vector{Float64}
    "Signed edge velocity at those knots."
    values::Vector{Float64}
end

ContourPressureScratch() = ContourPressureScratch(Float64[], Float64[], Float64[],
                                                  Float64[], Float64[], Float64[])

"""
    arc_at_chord(x, arc, indices, fractions) -> Vector

Signed arc length at each chord fraction, read off the contour nodes `indices` — the
one surface the fractions belong to.
"""
function arc_at_chord(x, arc, indices, fractions)
    scratch = ContourPressureScratch()
    return arc_at_chord!(zeros(length(fractions)), scratch, x, arc, indices, fractions)
end

"""
    arc_at_chord!(out, scratch, x, arc, indices, fractions) -> out

[`arc_at_chord`](@ref) written into `out`, taking the sorted copy of the surface from
`scratch` rather than allocating one per call.
"""
function arc_at_chord!(out::Vector{Float64}, scratch::ContourPressureScratch, x, arc,
                       indices, fractions)
    xs, as = sorted_surface!(scratch, x, arc, indices)
    resize!(out, length(fractions))
    @inbounds for (i, f) in enumerate(fractions)
        j = clamp(searchsortedfirst(xs, f), 2, length(xs))
        gap = max(xs[j] - xs[j - 1], eps())
        out[i] = as[j - 1] + (f - xs[j - 1]) / gap * (as[j] - as[j - 1])
    end
    return out
end

"""
    sorted_surface!(scratch, x, arc, indices) -> (surface_x, surface_arc)

One surface's nodes as chord fraction and signed arc length, ascending in chord, in the
scratch's own storage. Sorted by insertion, which a surface listed from one end to the
other is either already in or exactly reversed from, and which needs no scratch of its
own.
"""
function sorted_surface!(scratch::ContourPressureScratch, x, arc, indices)
    surface_x, surface_arc = scratch.surface_x, scratch.surface_arc
    n = length(indices)
    resize!(surface_x, n)
    resize!(surface_arc, n)
    @inbounds for (k, node) in enumerate(indices)
        surface_x[k] = x[node]
        surface_arc[k] = arc[node]
    end
    sort_pairs!(surface_x, surface_arc)
    return (surface_x, surface_arc)
end

"""
    sort_pairs!(sorted, carried)

Insertion-sort `sorted` ascending, carrying `carried` along with it. Stable, in place,
and linear on input that is already ordered, which both of its callers hand it.
"""
function sort_pairs!(sorted::Vector{Float64}, carried::Vector{Float64})
    @inbounds for k in 2:length(sorted)
        key, along = sorted[k], carried[k]
        j = k - 1
        while j >= 1 && sorted[j] > key
            sorted[j + 1], carried[j + 1] = sorted[j], carried[j]
            j -= 1
        end
        sorted[j + 1], carried[j + 1] = key, along
    end
    return nothing
end

"""
    trailing_edge_speed(ue_upper, ue_lower, upper_arc, lower_arc, upper_te, lower_te)
        -> Float64

One edge speed for both surfaces at the trailing edge: each surface's last two
stations extrapolated to its own trailing-edge arc length, the two magnitudes
averaged, so a sharp edge carries a single pressure (Kutta).
"""
function trailing_edge_speed(ue_upper, ue_lower, upper_arc, lower_arc,
                             upper_te, lower_te)
    reach(a1, a2, v1, v2, target) =
        v2 + (target - a2) / (a2 - a1 + eps()) * (v2 - v1)
    n = length(upper_arc)
    n < 2 && return abs(ue_upper[end])
    upper = reach(upper_arc[n - 1], upper_arc[n],
                  ue_upper[n - 1], ue_upper[n], upper_te)
    lower = reach(lower_arc[n - 1], lower_arc[n],
                  ue_lower[n - 1], ue_lower[n], lower_te)
    return (abs(upper) + abs(lower)) / 2
end

"""
    velocity_knots(station_x, ue_upper, ue_lower, x, arc, le) -> (knots, values)

The whole contour's edge velocity as one curve of (signed arc length, signed
`ue/u∞`), strictly increasing in arc length.

NeuralFoil reports nothing over the first and last `1/2N` of chord. Signing the lower
surface negative joins both surfaces into one curve through the stagnation point, so
the unsampled nose is interpolated between the innermost station on each side rather
than extrapolated off the end of one. The stagnation point enters as its own knot,
where `ue` interpolated linearly between those two stations reaches zero — the
stagnation-point-flow result, `ue` being linear in arc length there.
"""
velocity_knots(station_x, ue_upper, ue_lower, x, arc, le) =
    velocity_knots!(ContourPressureScratch(), station_x, ue_upper, ue_lower, x, arc, le)

"""
    velocity_knots!(scratch, station_x, ue_upper, ue_lower, x, arc, le)
        -> (knots, values)

[`velocity_knots`](@ref) assembled in `scratch`. The two vectors returned are the
scratch's own and are overwritten by the next call.
"""
function velocity_knots!(scratch::ContourPressureScratch, station_x, ue_upper, ue_lower,
                         x, arc, le)
    n_stations = length(station_x)
    upper_arc = arc_at_chord!(scratch.upper_arc, scratch, x, arc, 1:le, station_x)
    lower_arc = arc_at_chord!(scratch.lower_arc, scratch, x, arc, le:length(x),
                              station_x)
    speed = trailing_edge_speed(ue_upper, ue_lower, upper_arc, lower_arc,
                                arc[1], arc[end])
    lower_first, upper_first = abs(ue_lower[1]), abs(ue_upper[1])
    stagnation = lower_arc[1] + lower_first / max(lower_first + upper_first, eps()) *
                 (upper_arc[1] - lower_arc[1])
    knots = resize!(scratch.knots, 2 * n_stations + 3)
    values = resize!(scratch.values, 2 * n_stations + 3)
    knots[1], values[1] = arc[end], -speed
    for k in 1:n_stations
        knots[1 + k], values[1 + k] = lower_arc[n_stations + 1 - k],
                                      -abs(ue_lower[n_stations + 1 - k])
    end
    knots[n_stations + 2], values[n_stations + 2] = stagnation, 0.0
    for k in 1:n_stations
        knots[n_stations + 2 + k], values[n_stations + 2 + k] = upper_arc[k],
                                                                abs(ue_upper[k])
    end
    knots[end], values[end] = arc[1], speed
    sort_pairs!(knots, values)
    kept, previous = 1, knots[1]
    for k in 2:length(knots)
        current = knots[k]
        if current > previous + 1e-9
            kept += 1
            knots[kept], values[kept] = current, values[k]
        end
        previous = current
    end
    return (resize!(knots, kept), resize!(values, kept))
end

"""
    contour_pressure(station_x, ue_upper, ue_lower, x, arc, le) -> Vector{Float64}

`Cp` at every node of a closed contour `x` with signed arc length `arc` and nose at
`le`, from one case's edge-velocity ratios at NeuralFoil's `station_x`. Builds the
one-curve edge velocity ([`velocity_knots`](@ref)), interpolates it with a
shape-preserving monotone cubic in arc length, and squares it into `Cp = 1 - ue²`.
"""
contour_pressure(station_x, ue_upper, ue_lower, x, arc, le) =
    contour_pressure!(zeros(length(x)), ContourPressureScratch(), station_x, ue_upper,
                      ue_lower, x, arc, le)

"""
    contour_pressure!(cp, scratch, station_x, ue_upper, ue_lower, x, arc, le) -> cp

[`contour_pressure`](@ref) written into `cp`, with `scratch` carrying the edge-velocity
curve. Only the monotone interpolation itself still allocates, once per contour.
"""
function contour_pressure!(cp, scratch::ContourPressureScratch, station_x, ue_upper,
                           ue_lower, x, arc, le)
    knots, values = velocity_knots!(scratch, station_x, ue_upper, ue_lower, x, arc, le)
    speed = interpolate(knots, values, FritschButlandMonotonicInterpolation())
    @inbounds for k in eachindex(cp, arc)
        cp[k] = 1 - speed(clamp(arc[k], first(knots), last(knots)))^2
    end
    return cp
end

"""
    neuralfoil_contour_solution(alpha, res, i, x, y, le, arc, cf) -> SectionSolution

Assemble a full-contour [`SectionSolution`](@ref) for case `i` of a NeuralFoil sweep
`res`: reconstruct `Cp` on the contour nodes from the predicted edge velocities
([`contour_pressure`](@ref)), carrying the precomputed `cf`.
"""
function neuralfoil_contour_solution(alpha, res, i, x, y, le, arc, cf)
    cp = contour_pressure(res.x, view(res.ue_upper, :, i),
                          view(res.ue_lower, :, i), x, arc, le)
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

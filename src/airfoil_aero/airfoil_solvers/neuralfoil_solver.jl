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
Kulfan parameters, then assemble a single per-node `cp` on the deformed contour nodes
(`def.x`, `def.y`) from NeuralFoil's separate upper/lower surface pressures. `cf` is a
flat-plate closure ([`flat_plate_cf`](@ref)), NeuralFoil not exposing skin friction.
"""
function analyze_sweep(solver::NeuralFoilSolver, def::DeformedSection, alpha_range, Re)
    res = neuralfoil_section(def.kulfan, rad2deg.(collect(alpha_range)), Re;
        model_size=solver.model_size, weights_dir=solver.weights_dir,
        n_crit=solver.n_crit, xtr_upper=solver.xtr_upper, xtr_lower=solver.xtr_lower)
    x, y = collect(float.(def.x)), collect(float.(def.y))
    le = argmin(x)
    cf = [flat_plate_cf(clamp(xk, 0.0, 1.0), Re) for xk in x]
    return [neuralfoil_contour_solution(alpha_range[i], res, i, x, y, le, cf)
            for i in eachindex(alpha_range)]
end

"""
    neuralfoil_contour_solution(alpha, res, i, x, y, le, cf) -> SectionSolution

Assemble a full-contour [`SectionSolution`](@ref) for case `i` of a NeuralFoil sweep
`res`: interpolate `res.cp_upper`/`res.cp_lower` (at `res.x`) onto the contour nodes
`(x, y)` split at the leading edge `le`, carrying the precomputed `cf`.
"""
function neuralfoil_contour_solution(alpha, res, i, x, y, le, cf)
    up = linear_interpolation(res.x, res.cp_upper[:, i]; extrapolation_bc=Line())
    lo = linear_interpolation(res.x, res.cp_lower[:, i]; extrapolation_bc=Line())
    cp = [(k <= le ? up : lo)(clamp(x[k], 0.0, 1.0)) for k in eachindex(x)]
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

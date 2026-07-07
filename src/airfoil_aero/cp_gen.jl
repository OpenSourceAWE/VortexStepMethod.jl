"""
    generate_cp_polar(solver, base; alpha_range, delta_range, n_chord,
                      reynolds_number, crease_frac=0.9, thickness_frac=1.0,
                      flip_thickness_neg=true, remove_nan=true) -> CpPolar

Build a chord-slice `Cp(alpha, delta)` table for a base airfoil (Kulfan
parameters) using the given [`AbstractAirfoilSolver`](@ref). For each `delta` the
section is deformed and refit ([`deform_section`](@ref)), analysed over
`alpha_range`, and each surface's `Cp(x/c)` resampled onto the `n_chord` chord
slices. Non-converged points stay `NaN` and, when `remove_nan`, are filled per
slice with [`interpolate_matrix_nans!`](@ref). `flip_thickness_neg` flips the
through-thickness pivot for negative `delta` (see [`deform_section`](@ref)).
"""
function generate_cp_polar(solver::AbstractAirfoilSolver, base::KulfanParameters;
        alpha_range, delta_range, n_chord, reynolds_number,
        crease_frac=0.9, thickness_frac=1.0, flip_thickness_neg=true, remove_nan=true)
    chord_x = [(i - 0.5) / n_chord for i in 1:n_chord]
    n_alpha, n_delta = length(alpha_range), length(delta_range)
    cp_upper = fill(NaN, n_chord, n_alpha, n_delta)
    cp_lower = fill(NaN, n_chord, n_alpha, n_delta)

    for (jd, delta) in enumerate(delta_range)
        def = deform_section(base, delta; crease_frac, thickness_frac,
                             flip_thickness_neg)
        sols = analyze_sweep(solver, def, alpha_range, reynolds_number)
        for (ia, sol) in enumerate(sols)
            isempty(sol.cp_upper) && continue
            cp_upper[:, ia, jd] = resample_x(sol.x_upper, sol.cp_upper, chord_x)
            cp_lower[:, ia, jd] = resample_x(sol.x_lower, sol.cp_lower, chord_x)
        end
    end

    if remove_nan
        fill_slice_nans!(cp_upper)
        fill_slice_nans!(cp_lower)
    end
    data = CpData(n_chord, chord_x, collect(float.(alpha_range)),
                  collect(float.(delta_range)), cp_upper, cp_lower)
    return CpPolar(data)
end

"""
    generate_cp_polar(solver, x::Vector, y::Vector; kwargs...) -> CpPolar

Convenience: fit base Kulfan parameters from coordinates (via `EnvelopeFit`) first.
"""
function generate_cp_polar(solver::AbstractAirfoilSolver, x::Vector, y::Vector; kwargs...)
    return generate_cp_polar(solver, fit_kulfan_parameters(x, y, EnvelopeFit()); kwargs...)
end

"""
    fill_slice_nans!(grid)

Fill `NaN`s in each chord slice's `(alpha, delta)` matrix with
[`interpolate_matrix_nans!`](@ref); all-`NaN` slices are left untouched.
"""
function fill_slice_nans!(grid)
    for i in axes(grid, 1)
        slice = grid[i, :, :]
        all(isnan, slice) && continue
        interpolate_matrix_nans!(slice; prn=false)
        grid[i, :, :] = slice
    end
    return nothing
end

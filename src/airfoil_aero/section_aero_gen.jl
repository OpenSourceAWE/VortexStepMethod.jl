"""
    generate_section_aero(solver, base; alpha_range, delta_range, reynolds_number,
                          crease_frac=0.9, remove_nan=true) -> SectionAero

Build a native-resolution surface table (closed contour + `Cp` + `cf` per node over the
`(alpha, delta)` grid) for a base airfoil (Kulfan parameters) using the given
[`AbstractAirfoilSolver`](@ref). For each `delta` the section is deformed
([`deform_section`](@ref)) and swept over `alpha_range` (radians); each converged angle
contributes its full-contour `(x, y, cp, cf)`. Non-converged points stay `NaN` and,
when `remove_nan`, are filled per node with `interpolate_matrix_nans!`.
"""
function generate_section_aero(solver::AbstractAirfoilSolver, base::KulfanParameters;
        alpha_range, delta_range, reynolds_number, crease_frac=0.9, remove_nan=true)
    x0, y0 = kulfan_to_coordinates(base)
    n_alpha, n_delta = length(alpha_range), length(delta_range)
    sols = Vector{Vector{SectionSolution}}(undef, n_delta)
    for (jd, delta) in enumerate(delta_range)
        def = deform_section(x0, y0, delta; crease_frac)
        sols[jd] = analyze_sweep(solver, def, alpha_range, reynolds_number)
    end

    n_node = 0
    for col in sols, s in col
        if !isempty(s.cp)
            n_node = length(s.cp)
            break
        end
    end
    n_node == 0 && throw(ErrorException("generate_section_aero: no converged angle"))

    x = fill(NaN, n_node, n_delta)
    y = fill(NaN, n_node, n_delta)
    cp = fill(NaN, n_node, n_alpha, n_delta)
    cf = fill(NaN, n_node, n_alpha, n_delta)
    for jd in 1:n_delta, (ia, s) in enumerate(sols[jd])
        (isempty(s.cp) || length(s.cp) != n_node) && continue
        cp[:, ia, jd] .= s.cp
        cf[:, ia, jd] .= s.cf
        if any(isnan, view(x, :, jd))
            x[:, jd] .= s.x
            y[:, jd] .= s.y
        end
    end

    if remove_nan
        for i in 1:n_node
            fill_node_nans!(cp, i)
            fill_node_nans!(cf, i)
        end
    end
    return SectionAero(collect(float.(alpha_range)), collect(float.(delta_range)),
                       x, y, cp, cf)
end

"""
    generate_section_aero(solver, x::Vector, y::Vector; kwargs...) -> SectionAero

Convenience: [`shrink_wrap`](@ref) the coordinates and fit base Kulfan parameters
([`LeastSquaresFit`](@ref)) first.
"""
function generate_section_aero(solver::AbstractAirfoilSolver, x::Vector, y::Vector;
        kwargs...)
    xw, yw = shrink_wrap(x, y, ShrinkWrap())
    return generate_section_aero(solver, fit_kulfan_parameters(xw, yw, LeastSquaresFit());
                                 kwargs...)
end

"""
    fill_node_nans!(grid, i)

Fill `NaN`s in node `i`'s `(alpha, delta)` matrix `grid[i, :, :]` with
`interpolate_matrix_nans!`; an all-`NaN` node is left untouched.
"""
function fill_node_nans!(grid, i)
    m = grid[i, :, :]
    all(isnan, m) && return grid
    interpolate_matrix_nans!(m; prn=false)
    grid[i, :, :] = m
    return grid
end

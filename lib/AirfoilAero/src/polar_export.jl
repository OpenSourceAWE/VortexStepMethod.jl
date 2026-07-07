"""
    generate_aero_matrices(solver, base; alpha_range, delta_range, Re,
                           crease_frac=0.75, remove_nan=true) -> (cl, cd, cm)

Build `(alpha × delta)` coefficient matrices for a base airfoil (`KulfanParameters`).
Each `delta` deflects the trailing edge ([`deform_section`](@ref)); the deflected
shape is then swept over `alpha_range` (radians) with `solver` — any
[`AbstractAirfoilSolver`](@ref), so this works identically for XFoil and NeuralFoil.
`Re` is the Reynolds number. With `remove_nan` the (non-converged) `NaN` entries are
interpolated away.
"""
function generate_aero_matrices(solver::AbstractAirfoilSolver, base::KulfanParameters;
        alpha_range, delta_range, Re, crease_frac=0.75, remove_nan=true)
    na, nd = length(alpha_range), length(delta_range)
    cl = fill(NaN, na, nd)
    cd = fill(NaN, na, nd)
    cm = fill(NaN, na, nd)
    alphas = collect(Float64, alpha_range)
    for (j, delta) in enumerate(delta_range)
        def = deform_section(base, delta; crease_frac)
        sols = analyze_sweep(solver, def, alphas, Re)
        for (i, s) in enumerate(sols)
            cl[i, j] = s.cl
            cd[i, j] = s.cd
            cm[i, j] = s.cm
        end
    end
    if remove_nan
        interpolate_matrix_nans!(cl)
        interpolate_matrix_nans!(cd)
        interpolate_matrix_nans!(cm)
    end
    return cl, cd, cm
end

"""
    create_2d_polars(; dat_path, cl_polar_path, cd_polar_path, cm_polar_path, wind_vel,
                     area, width, crease_frac, alpha_range, delta_range,
                     solver=XFoilSolver(), remove_nan=true)

Generate 2D airfoil-section `(alpha, delta)` `POLAR_MATRICES` CSVs for an airfoil
`.dat`, deflecting the trailing edge over `delta_range` and sweeping `alpha_range`
(both radians) with `solver`. Pass `solver=NeuralFoilSolver()` to use NeuralFoil
instead of XFoil — both emit the same three CSV files, so VSM loads them identically.

The Reynolds number is `wind_vel * (area / width) / ν` (kinematic viscosity of air).
`crease_frac` is the chordwise hinge location (0–1). Writes lift, drag, and moment
coefficient matrices to the three output paths.
"""
function create_2d_polars(; dat_path, cl_polar_path, cd_polar_path, cm_polar_path,
        wind_vel, area, width, crease_frac, alpha_range, delta_range,
        solver::AbstractAirfoilSolver=XFoilSolver(), remove_nan=true)
    x, y = read_dat_coordinates(dat_path)
    isempty(x) && error("No valid coordinates found in $dat_path")
    base = fit_kulfan_parameters(x, y)
    Re = wind_vel * (area / width) / KINEMATIC_VISCOSITY
    @info "Generating polars with $(nameof(typeof(solver))) at Re=$(round(Re))."
    cl, cd, cm = generate_aero_matrices(solver, base;
        alpha_range, delta_range, Re, crease_frac, remove_nan)
    write_aero_matrix(cl_polar_path, cl, alpha_range, delta_range, "C_l")
    write_aero_matrix(cd_polar_path, cd, alpha_range, delta_range, "C_d")
    write_aero_matrix(cm_polar_path, cm, alpha_range, delta_range, "C_m")
    return nothing
end

"""
    read_dat_coordinates(path) -> (x, y)

Read airfoil coordinates from a Selig-format `.dat` file, skipping header and
comment lines.
"""
function read_dat_coordinates(path::String)
    x = Float64[]
    y = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || !(isdigit(s[1]) || s[1] == '-' || s[1] == '.')) && continue
        parts = split(s)
        length(parts) >= 2 || continue
        xp = tryparse(Float64, parts[1])
        yp = tryparse(Float64, parts[2])
        (xp === nothing || yp === nothing) && continue
        push!(x, xp)
        push!(y, yp)
    end
    return x, y
end

"""
    write_dat(filepath, name, x, y) -> filepath

Write airfoil coordinates to a Selig-format `.dat` file.
"""
function write_dat(filepath::String, name::String, x::Vector, y::Vector)
    open(filepath, "w") do io
        println(io, name)
        for i in eachindex(x)
            println(io, @sprintf("%.5f %.5f", x[i], y[i]))
        end
    end
    return filepath
end

"""
    write_polar_csv(filepath, result::NeuralFoilResult)

Write a NeuralFoil result to a `POLAR_VECTORS` CSV (`alpha, Cd, Cs, Cl, Cm`).
"""
function write_polar_csv(filepath::String, result::NeuralFoilResult)
    n = length(result.alpha)
    open(filepath, "w") do io
        println(io, "alpha,Cd,Cs,Cl,Cm")
        for i in 1:n
            row = (result.alpha[i], result.CD[i], 0.0, result.CL[i], result.CM[i])
            println(io, join((@sprintf("%.16g", v) for v in row), ","))
        end
    end
    return filepath
end

"""
    generate_polar_from_coordinates(x, y, output_path; Re, alpha_range=-180:1:180,
                                    solver=NeuralFoilSolver(), delta_range=nothing,
                                    crease_frac=0.75, dat_prefix=nothing)

Sweep `solver` over the airfoil coordinates `(x, y)` and write the polar CSV. XFoil
uses the coordinates directly; NeuralFoil fits [`LeastSquaresFit`](@ref) Kulfan
parameters ([`deform_section`](@ref)). Wrap a raw or open single-membrane slice with
[`shrink_wrap`](@ref) before calling this. Pass a [`NeuralFoilSolver`](@ref) or
[`XFoilSolver`](@ref) to pick the backend. With `delta_range === nothing` the sweep is
over `alpha_range` only and written as a `POLAR_VECTORS` CSV (returns the
`Vector{SectionSolution}`); pass a `delta_range` of trailing-edge deflections to sweep
`(alpha, delta)` and write a long-format `POLAR_MATRICES` CSV (returns the `(cl, cd,
cm)` matrices). Both angle ranges are in degrees. `crease_frac` is the chordwise hinge
location (0–1) about which each `delta_range` deflection pivots. With `dat_prefix` set,
each deflected shape is also written to `{dat_prefix}_d{deg}.dat` (the deflection in
degrees).
"""
function generate_polar_from_coordinates(x::Vector, y::Vector, output_path::String;
                                         Re::Real, alpha_range=-180:1:180,
                                         solver::AbstractAirfoilSolver=NeuralFoilSolver(),
                                         delta_range=nothing, crease_frac=0.75,
                                         dat_prefix=nothing)
    delta_tag(delta_rad) = begin
        deg = round(rad2deg(delta_rad); digits=1)
        "d" * (deg == round(deg) ? string(Int(round(deg))) : string(deg))
    end
    alphas = deg2rad.(collect(Float64, alpha_range))
    if delta_range === nothing
        def = deform_section(x, y, 0.0)
        sols = analyze_sweep(solver, def, alphas, Re)
        write_polar_csv(output_path, sols)
        return sols
    end
    deltas = deg2rad.(collect(Float64, delta_range))
    on_deform = dat_prefix === nothing ? nothing :
        (d, xd, yd) -> write_dat("$(dat_prefix)_$(delta_tag(d)).dat",
                                 "deflection", xd, yd)
    cl, cd, cm = generate_aero_matrices(solver, x, y;
        alpha_range=alphas, delta_range=deltas, Re, crease_frac, on_deform)
    write_polar_matrix_csv(output_path, alphas, deltas, cl, cd, cm)
    return (cl, cd, cm)
end

"""
    resolve_airfoil(type, info, out_dir, id; Re, alpha_range) -> (new_type, new_info)

Resolve one awesIO `wing_airfoils` entry to a core-loadable form. `breukels_regression`
`(t, kappa)` → `poly` coeffs (via [`lei_poly_coeffs`](@ref)); `neuralfoil`
`(dat_file_path, …)` → a `polars` CSV (via [`generate_polar_from_dat`](@ref)) written
under `out_dir`; `polars`/`poly`/`inviscid` pass through. `masure_regression` is not
yet supported. `info` file paths should already be absolute.
"""
function resolve_airfoil(type::AbstractString, info::AbstractDict, out_dir, id;
                         Re, alpha_range)
    if type == "breukels_regression"
        cl, cd, cm = lei_poly_coeffs(Float64(info["t"]), Float64(info["kappa"]))
        return "poly", Dict{String,Any}("cl_coeffs" => cl, "cd_coeffs" => cd,
                                        "cm_coeffs" => cm)
    elseif type == "neuralfoil"
        csv = joinpath(out_dir, "$(id).csv")
        solver = NeuralFoilSolver(
            model_size=String(get(info, "model_size", "large")),
            n_crit=Float64(get(info, "n_crit", 9.0)))
        generate_polar_from_dat(String(info["dat_file_path"]), csv; Re,
            alpha_range=collect(alpha_range), solver)
        return "polars", Dict{String,Any}("csv_file_path" => csv)
    elseif type in ("polars", "poly", "inviscid")
        return String(type), info
    elseif type == "masure_regression"
        error("masure_regression airfoils are not yet supported by AirfoilAero.")
    else
        error("Unknown airfoil type: $type")
    end
end

"""
    generate_polar_from_dat(dat_path, output_path; Re, kwargs...)

Read a `.dat` airfoil and generate its polar CSV via
`generate_polar_from_coordinates`. This is AirfoilAero's `.dat` → CSV entry
point.
"""
function generate_polar_from_dat(dat_path::String, output_path::String;
                                 Re::Real, kwargs...)
    isfile(dat_path) || error("DAT file not found: $dat_path")
    x, y = read_dat_coordinates(dat_path)
    isempty(x) && error("No valid coordinates found in $dat_path")
    return generate_polar_from_coordinates(x, y, output_path; Re, kwargs...)
end

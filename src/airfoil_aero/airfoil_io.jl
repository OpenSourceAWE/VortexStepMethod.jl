"""
    write_dat(filepath, name, x, y) -> filepath

Write airfoil coordinates to a Selig-format `.dat` file.
"""
function write_dat(filepath::String, name::String, x::Vector, y::Vector)
    open(filepath, "w") do io
        println(io, name)
        for i in eachindex(x)
            println(io, @sprintf("%.8f %.8f", x[i], y[i]))
        end
    end
    return filepath
end

"""
    write_polar(filepath, result::NeuralFoilResult)

Write a NeuralFoil result to a `POLAR_VECTORS` table (`alpha, Cd, Cs, Cl, Cm`), CSV or
Arrow as the suffix of `filepath` says (see [`write_node_rows`](@ref
VortexStepMethod.write_node_rows)).
"""
function write_polar(filepath::String, result::NeuralFoilResult)
    values = [result.CD zero(result.CD) result.CL result.CM]
    return write_node_rows(filepath, deg2rad.(result.alpha), nothing, values;
                           columns=["Cd", "Cs", "Cl", "Cm"])
end

"""
    generate_polar_from_coordinates(x, y, output_path; Re, alpha_range=-180:1:180,
                                    solver=NeuralFoilSolver(), delta_range=nothing,
                                    crease_frac=0.75, dat_prefix=nothing,
                                    wrap_method=ShrinkWrap())

Sweep `solver` over the airfoil coordinates `(x, y)` and write the polar table, CSV or
Arrow as the suffix of `output_path` says. XFoil uses the coordinates directly;
NeuralFoil fits [`LeastSquaresFit`](@ref) Kulfan parameters ([`deform_section`](@ref)).
Wrap a raw or open single-membrane slice with [`shrink_wrap`](@ref) before calling this
and pass that [`ShrinkWrap`](@ref) as `wrap_method`. Pass a [`NeuralFoilSolver`](@ref)
or [`XFoilSolver`](@ref) to pick the backend. With `delta_range === nothing` the sweep
is over `alpha_range` only and written as a `POLAR_VECTORS` table (returns the
`Vector{SectionSolution}`); pass a `delta_range` of trailing-edge deflections to sweep
`(alpha, delta)` and write a long-format `POLAR_MATRICES` table (returns the `(cl, cd,
cm)` matrices). Both angle ranges are in degrees. `crease_frac` is the chordwise hinge
location (0–1) about which each `delta_range` deflection pivots. With `dat_prefix` set,
each deflected shape is also written to `{dat_prefix}_{delta_suffix(δ)}.dat`.
"""
function generate_polar_from_coordinates(x::Vector, y::Vector, output_path::String;
                                         Re::Real, alpha_range=-180:1:180,
                                         solver::AbstractAirfoilSolver=NeuralFoilSolver(),
                                         delta_range=nothing, crease_frac=0.75,
                                         dat_prefix=nothing,
                                         wrap_method::ShrinkWrap=ShrinkWrap())
    alphas = deg2rad.(collect(Float64, alpha_range))
    if delta_range === nothing
        def = deform_section(x, y, 0.0; wrap_method)
        sols = analyze_sweep(solver, def, alphas, Re)
        write_polar(output_path, sols)
        return sols
    end
    deltas = deg2rad.(collect(Float64, delta_range))
    on_deform = dat_prefix === nothing ? nothing :
        (d, xd, yd) -> write_dat("$(dat_prefix)_$(delta_suffix(d)).dat",
                                 "deflection", xd, yd)
    cl, cd, cm = generate_aero_matrices(solver, x, y;
        alpha_range=alphas, delta_range=deltas, Re, crease_frac, on_deform, wrap_method)
    write_polar_matrix(output_path, alphas, deltas, cl, cd, cm)
    return (cl, cd, cm)
end

"""
    resolve_airfoil(type, info, out_dir, id; Re, alpha_range, table_format=:csv)
        -> (new_type, new_info)

Resolve one awesIO `wing_airfoils` entry to a core-loadable form. `breukels_regression`
`(t, kappa)` → `poly` coeffs (via [`lei_poly_coeffs`](@ref)); `neuralfoil`
`(dat_file_path, …)` → a `polars` table `{id}.{table_format}` (`:csv` or `:arrow`, via
[`generate_polar_from_dat`](@ref)) written under `out_dir`; `polars`/`poly`/`inviscid`
pass through. `masure_regression` is not yet supported. `info` file paths should
already be absolute.
"""
function resolve_airfoil(type::AbstractString, info::AbstractDict, out_dir, id;
                         Re, alpha_range, table_format::Symbol=:csv)
    if type == "breukels_regression"
        cl, cd, cm = lei_poly_coeffs(Float64(info["t"]), Float64(info["kappa"]))
        return "poly", Dict{String,Any}("cl_coeffs" => cl, "cd_coeffs" => cd,
                                        "cm_coeffs" => cm)
    elseif type == "neuralfoil"
        polar = joinpath(out_dir, "$(id).$(table_format)")
        solver = NeuralFoilSolver(
            model_size=String(get(info, "model_size", "large")),
            n_crit=Float64(get(info, "n_crit", 9.0)))
        generate_polar_from_dat(String(info["dat_file_path"]), polar; Re,
            alpha_range=collect(alpha_range), solver)
        return "polars", Dict{String,Any}("polar_file_path" => polar)
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

Read a `.dat` airfoil and generate its polar table via
`generate_polar_from_coordinates`. This is AirfoilAero's `.dat` → polar entry point.
"""
function generate_polar_from_dat(dat_path::String, output_path::String;
                                 Re::Real, kwargs...)
    isfile(dat_path) || error("DAT file not found: $dat_path")
    x, y = read_dat_coordinates(dat_path)
    isempty(x) && error("No valid coordinates found in $dat_path")
    return generate_polar_from_coordinates(x, y, output_path; Re, kwargs...)
end

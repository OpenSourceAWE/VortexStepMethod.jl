"""
Adapter that converts a 3D wing `.obj` mesh into the package's native YAML
geometry route: per-section airfoil `.dat` files, NeuralFoil polar CSVs, and a
geometry YAML referencing them. After conversion everything follows the standard
`Wing(geometry_file)` path.
"""

"""
    airfoils_from_yaml(geometry_file) -> Vector of (; id, x, y, x_raw, y_raw)

Read each airfoil's coordinates from the `.dat` files referenced by a YAML
geometry's `wing_airfoils`. `x, y` come from `dat_file` (the fitted airfoil);
`x_raw, y_raw` come from `raw_dat_file` if present (the raw sliced points), and
are empty otherwise. Relative paths resolve against the geometry file's directory.
"""
function airfoils_from_yaml(geometry_file::String)
    data = YAML.load_file(geometry_file)
    airfoils = data["wing_airfoils"]
    headers = airfoils["headers"]
    base = dirname(geometry_file)
    resolve(rel) = isabspath(rel) ? rel : joinpath(base, rel)
    out = NamedTuple[]
    for row in airfoils["data"]
        d = Dict(zip(headers, row))
        dat_rel = get(d["info_dict"], "dat_file", "")
        (isempty(dat_rel) || !isfile(resolve(dat_rel))) && continue
        x, y = read_dat_coordinates(resolve(dat_rel))
        x_raw, y_raw = Float64[], Float64[]
        raw_rel = get(d["info_dict"], "raw_dat_file", "")
        if !isempty(raw_rel) && isfile(resolve(raw_rel))
            x_raw, y_raw = read_dat_coordinates(resolve(raw_rel))
        end
        push!(out, (id=Int(d["airfoil_id"]), x=x, y=y, x_raw=x_raw, y_raw=y_raw))
    end
    return out
end

"""
    obj_to_yaml(obj_path, output_dir; n_sections, Re, alpha_range=-180:1:180,
                aero_solver=NeuralFoilSolver(...), model_size="xlarge",
                weights_dir=nothing, n_crit=9.0, fit_method=LeastSquaresFit(),
                spanwise_direction=[0.0, 1.0, 0.0], verbose=true)

Convert a 3D wing `.obj` mesh to the native YAML geometry route.

Stations are placed at equal leading-edge arc-length intervals and sliced
perpendicular to the local span (see [`perpendicular_sections`](@ref)), which
keeps the airfoil undistorted near curved tips; each shape is then fitted to
Kulfan parameters and evaluated with `aero_solver`.

`aero_solver` selects the 2D-airfoil backend: [`NeuralFoilSolver`](@ref) (default,
fast) or [`XFoilSolver`](@ref) (viscous panel code). The default solver is built
from `model_size`, `n_crit`, and `weights_dir`; pass `aero_solver=XFoilSolver()`
to use XFoil instead. Each section's polar is written as `POLAR_VECTORS`.

`fit_method` selects the Kulfan fit (see [`fit_kulfan_parameters`](@ref)). Pass
an [`EnvelopeFit`](@ref) to wrap each section tightly *around* the slice points;
this is robust to the noisy interior-structure points of a ram-air kite slice
(ribs, spars) that otherwise pull the default [`LeastSquaresFit`](@ref) inward.

A slice whose fitted airfoil is implausibly thick relative to the others (e.g. a
near-vanishing wingtip slice that fits to a blob) is flagged as degenerate: a
section is degenerate when its fitted `max|y|` exceeds `max_thickness_ratio` times
the median over all sections. With `reuse_valid_airfoils=true` (default) each
degenerate section reuses the nearest valid section's airfoil shape and polar,
while keeping its own leading- and trailing-edge positions, so the station is
still placed but carries a sane airfoil; a warning lists which sections were
reused. With `reuse_valid_airfoils=false` the degenerate fit is written as is.

# Output
Writes into `output_dir`, indexed by source-airfoil id `j` (degenerate sections
share a neighbour's id):
- `airfoils/{j}.dat`     — fitted Kulfan airfoil coordinates (matches the polar)
- `airfoils/{j}_raw.dat` — raw sliced section points (the fit is wrapped around these)
- `polars/{j}.csv`    — NeuralFoil polar (alpha, Cd, Cs, Cl, Cm)
- `geometry.yaml`     — `wing_sections` + `wing_airfoils` referencing the above

# Returns
- Path to the written `geometry.yaml`.
"""
function obj_to_yaml(obj_path::String, output_dir::String;
                     n_sections::Int, Re::Real,
                     alpha_range=-180:1:180,
                     model_size::String="xlarge", weights_dir=nothing, n_crit=9.0,
                     aero_solver::AbstractAirfoilSolver=NeuralFoilSolver(;
                         model_size, n_crit, weights_dir),
                     fit_method::KulfanFitMethod=LeastSquaresFit(),
                     reuse_valid_airfoils::Bool=true, max_thickness_ratio::Real=2.0,
                     spanwise_direction=[0.0, 1.0, 0.0], rotation=I,
                     verbose::Bool=true)
    (!endswith(obj_path, ".obj")) && (obj_path *= ".obj")
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) &&
        throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0]"))

    airfoil_dir = joinpath(output_dir, "airfoils")
    polar_dir = joinpath(output_dir, "polars")
    mkpath(airfoil_dir)
    mkpath(polar_dir)

    vertices, faces = read_faces(obj_path)
    alphas = collect(Float64, alpha_range)
    aero_solver isa NeuralFoilSolver &&
        load_neuralfoil_model(aero_solver.model_size; weights_dir=aero_solver.weights_dir)

    raw_sections = perpendicular_sections(vertices, faces, n_sections; rotation)
    isempty(raw_sections) && error("No valid sections sliced from $obj_path")

    stations = NamedTuple[]
    for (i, sec) in enumerate(raw_sections)
        verbose && print("  Section $i (y=$(round(sec.LE_point[2], digits=3)))... ")
        params = fit_kulfan_parameters(sec.x_airfoil, sec.y_airfoil, fit_method)
        x_fit, y_fit = kulfan_to_coordinates(params)
        thickness = maximum(abs, y_fit)
        verbose && println("fitted (thickness=$(round(thickness, digits=3)))")
        push!(stations, (; sec.LE_point, sec.TE_point, xa=sec.x_airfoil,
                         ya=sec.y_airfoil, params, x_fit, y_fit, thickness))
    end

    n = length(stations)
    thickness_limit = max_thickness_ratio * median(s.thickness for s in stations)
    degenerate = [s.thickness > thickness_limit for s in stations]
    valid = [k for k in 1:n if !degenerate[k]]
    isempty(valid) && error("All sliced sections are degenerate in $obj_path")
    ids = [degenerate[k] && reuse_valid_airfoils ?
           valid[argmin(abs.(valid .- k))] : k for k in 1:n]
    reused = [k => ids[k] for k in 1:n if degenerate[k] && ids[k] != k]
    isempty(reused) || @warn "Degenerate section fits reused the nearest valid " *
        "airfoil: " * join(["section $k → airfoil $j" for (k, j) in reused], ", ")

    section_rows = Vector{Any}[]
    airfoil_rows = Vector{Any}[]
    for j in unique(ids)
        s = stations[j]
        dat_rel = joinpath("airfoils", "$j.dat")
        raw_rel = joinpath("airfoils", "$(j)_raw.dat")
        csv_rel = joinpath("polars", "$j.csv")
        dat_abs = joinpath(output_dir, dat_rel)
        write_dat(dat_abs, "section_$j", s.x_fit, s.y_fit)
        write_dat(joinpath(output_dir, raw_rel), "section_$(j)_raw", s.xa, s.ya)
        def = DeformedSection(s.params, s.x_fit, s.y_fit)
        sols = analyze_sweep(aero_solver, def, deg2rad.(alphas), Float64(Re))
        write_polar_csv(joinpath(output_dir, csv_rel), sols)
        push!(airfoil_rows, Any[j, "polar_vectors",
                                Dict("dat_file" => dat_rel, "raw_dat_file" => raw_rel,
                                     "csv_file_path" => csv_rel)])
        cl_max = maximum((sol.cl for sol in sols if !isnan(sol.cl)); init=NaN)
        verbose && println("  Airfoil $j: CL_max=$(round(cl_max, digits=2))")
    end
    for k in 1:n
        s = stations[k]
        push!(section_rows, Any[ids[k], s.LE_point[1], s.LE_point[2], s.LE_point[3],
                                s.TE_point[1], s.TE_point[2], s.TE_point[3]])
    end

    yaml_path = joinpath(output_dir, "geometry.yaml")
    write_geometry_yaml(yaml_path, section_rows, airfoil_rows)
    verbose && @info "Wrote geometry to $yaml_path ($(length(section_rows)) sections)"
    return yaml_path
end

"""
    obj_to_matrix_yaml(obj_path, output_dir; n_sections, wind_vel, alpha_range,
                       delta_range, foil_dat=nothing, solver=NeuralFoilSolver(),
                       crease_frac=0.2, remove_nan=true, verbose=true) -> yaml_path

Convert an `.obj` wing mesh into a standard geometry YAML backed by one shared
`(alpha, delta)` `POLAR_MATRICES`. Section leading/trailing edges come from
perpendicular mesh slices ([`perpendicular_sections`](@ref)); the shared airfoil's
cl/cd/cm matrices are generated once with [`create_2d_polars`](@ref) and referenced by
every section. Load the result with `Wing(yaml_path; n_panels)`. This replaces the
old live `ObjWing` route with convert-then-load.

The shared airfoil is, by default, the mesh's median-chord slice (a representative
interior section, not the midspan). Pass `foil_dat` to override it with an explicit
`.dat`. `solver` selects the 2D backend: [`NeuralFoilSolver`](@ref) (default) or
[`XFoilSolver`](@ref).
"""
function obj_to_matrix_yaml(obj_path::String, output_dir::String;
        n_sections::Int, wind_vel::Real, alpha_range, delta_range,
        foil_dat=nothing, solver::AbstractAirfoilSolver=NeuralFoilSolver(),
        crease_frac=0.2, remove_nan=true, rotation=I, verbose=true)
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    polar_dir = joinpath(output_dir, "polars")
    airfoil_dir = joinpath(output_dir, "airfoils")
    mkpath(polar_dir)
    mkpath(airfoil_dir)

    vertices, faces = read_faces(obj_path)
    secs = perpendicular_sections(vertices, faces, n_sections; rotation)
    isempty(secs) && error("No valid sections sliced from $obj_path")

    mean_chord = sum(norm(s.TE_point .- s.LE_point) for s in secs) / length(secs)

    if foil_dat === nothing
        chords = [norm(s.TE_point .- s.LE_point) for s in secs]
        rep = secs[argmin(abs.(chords .- median(chords)))]
        x, y = rep.x_airfoil, rep.y_airfoil
    else
        isfile(foil_dat) || error("DAT file not found: $foil_dat")
        x, y = read_dat_coordinates(foil_dat)
    end
    dat_abs = joinpath(airfoil_dir, "1.dat")
    write_dat(dat_abs, "airfoil_1", x, y)

    cl_path = joinpath(polar_dir, "cl.csv")
    cd_path = joinpath(polar_dir, "cd.csv")
    cm_path = joinpath(polar_dir, "cm.csv")
    create_2d_polars(; dat_path=dat_abs, cl_polar_path=cl_path, cd_polar_path=cd_path,
        cm_polar_path=cm_path, wind_vel=Float64(wind_vel), area=mean_chord, width=1.0,
        crease_frac, alpha_range, delta_range, solver, remove_nan)

    section_rows = [Any[1, s.LE_point[1], s.LE_point[2], s.LE_point[3],
                        s.TE_point[1], s.TE_point[2], s.TE_point[3]] for s in secs]
    airfoil_rows = [Any[1, "polar_matrices", Dict(
        "dat_file" => joinpath("airfoils", "1.dat"),
        "cl_file_path" => joinpath("polars", "cl.csv"),
        "cd_file_path" => joinpath("polars", "cd.csv"),
        "cm_file_path" => joinpath("polars", "cm.csv"))]]
    yaml_path = joinpath(output_dir, "geometry.yaml")
    write_geometry_yaml(yaml_path, section_rows, airfoil_rows)
    verbose && @info "Wrote $yaml_path ($(length(section_rows)) sections, POLAR_MATRICES)"
    return yaml_path
end

"""
    resolve_aero_geometry(yaml_in, out_dir; verbose=true) -> yaml_out

Read an awesIO-style geometry YAML and resolve every `wing_airfoils` entry to a
core-loadable form via [`resolve_airfoil`](@ref) (`breukels_regression` → `poly`,
`neuralfoil` → `polars` CSV, others pass through), writing generated CSVs and a
resolved `geometry.yaml` under `out_dir`. `wing_sections` (incl. any `VUP` up-vectors)
pass through unchanged. Load the result with `Wing(yaml_out)`.
"""
function resolve_aero_geometry(yaml_in::String, out_dir::String; verbose=true)
    mkpath(out_dir)
    polar_dir = joinpath(out_dir, "polars")
    mkpath(polar_dir)
    data = YAML.load_file(yaml_in)
    wa = data["wing_airfoils"]
    headers = wa["headers"]
    col(name) = findfirst(==(name), headers)
    ti, ii, idi = col("type"), col("info_dict"), col("airfoil_id")
    aspec = get(wa, "alpha_range", [-180, 180, 1.0])
    Re = Float64(get(wa, "reynolds", 1.0e6))
    alpha_range = Float64(aspec[1]):Float64(aspec[3]):Float64(aspec[2])
    base = dirname(abspath(yaml_in))
    for row in wa["data"]
        info = Dict{String,Any}(row[ii])
        for k in ("dat_file_path", "csv_file_path", "cl_file_path",
                  "cd_file_path", "cm_file_path")
            if haskey(info, k) && !isabspath(String(info[k]))
                info[k] = abspath(joinpath(base, String(info[k])))
            end
        end
        new_type, new_info = resolve_airfoil(String(row[ti]), info, polar_dir,
                                             row[idi]; Re, alpha_range)
        row[ti] = new_type
        row[ii] = new_info
    end
    yaml_out = joinpath(out_dir, "geometry.yaml")
    write_yaml(yaml_out, data)
    verbose && @info "Resolved geometry -> $yaml_out"
    return yaml_out
end

"""
    write_geometry_yaml(path, section_rows, airfoil_rows)

Write a geometry YAML via [`write_yaml`], one line per section/airfoil row.
`section_rows` are `[airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]`; `airfoil_rows`
are `[airfoil_id, type, info_dict]` where `info_dict` holds `dat_file`,
`csv_file_path`, and optionally `raw_dat_file`.
"""
function write_geometry_yaml(path::String, section_rows, airfoil_rows)
    data = Dict(
        "wing_sections" => Dict(
            "headers" => ["airfoil_id", "LE_x", "LE_y", "LE_z",
                          "TE_x", "TE_y", "TE_z"],
            "data" => section_rows),
        "wing_airfoils" => Dict(
            "headers" => ["airfoil_id", "type", "info_dict"],
            "data" => [[id, type, info] for (id, type, info) in airfoil_rows]))
    return write_yaml(path, data)
end

yaml_scalar(x::Bool) = string(x)
yaml_scalar(x::Integer) = string(x)
yaml_scalar(x::Real) = string(round(Float64(x); digits=3))
yaml_scalar(x::AbstractString) = "\"$x\""
yaml_scalar(::Nothing) = "null"
yaml_scalar(x) = string(x)

yaml_sorted(x::AbstractDict) = sort!(collect(x); by = p -> string(first(p)))
yaml_isscalar(x) = x isa Union{Real,AbstractString,Bool,Nothing}
yaml_scalardict(x) = x isa AbstractDict && all(yaml_isscalar, values(x))
yaml_flowable(x) = yaml_isscalar(x) || yaml_scalardict(x) ||
    (x isa AbstractVector && all(e -> yaml_isscalar(e) || yaml_scalardict(e), x))

yaml_flow(x) =
    yaml_isscalar(x) ? yaml_scalar(x) :
    x isa AbstractDict ?
        "{" * join(("$k: $(yaml_flow(v))" for (k, v) in yaml_sorted(x)), ", ") * "}" :
    x isa AbstractVector ? "[" * join(map(yaml_flow, x), ", ") * "]" :
    yaml_scalar(x)

function yaml_emit(io, x, indent)
    pad = "  "^indent
    if x isa AbstractDict
        for (k, v) in yaml_sorted(x)
            yaml_flowable(v) ? println(io, pad, k, ": ", yaml_flow(v)) :
                (println(io, pad, k, ":"); yaml_emit(io, v, indent + 1))
        end
    elseif x isa AbstractVector
        for el in x
            yaml_flowable(el) ? println(io, pad, "- ", yaml_flow(el)) :
                (println(io, pad, "-"); yaml_emit(io, el, indent + 1))
        end
    else
        println(io, pad, yaml_scalar(x))
    end
end

"""
    write_yaml(path, data)

Write `data` (nested `Dict`s, vectors, scalars) to `path` as YAML with every float
rounded to millimetre precision (3 decimals) and every leaf list on one line. The
single writer for all generated geometry YAMLs, so their formatting stays consistent.
Mapping keys are emitted in sorted order.
"""
function write_yaml(path::String, data)
    open(path, "w") do io
        yaml_emit(io, data, 0)
    end
    return path
end

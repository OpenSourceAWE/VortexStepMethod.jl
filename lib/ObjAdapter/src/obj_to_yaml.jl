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
                model_size="xlarge", weights_dir=nothing, n_crit=9.0,
                fit_method=LeastSquaresFit(),
                spanwise_direction=[0.0, 1.0, 0.0], verbose=true)

Convert a 3D wing `.obj` mesh to the native YAML geometry route.

Stations are placed at equal leading-edge arc-length intervals and sliced
perpendicular to the local span (see [`perpendicular_sections`](@ref)), which
keeps the airfoil undistorted near curved tips; each shape is then fitted to
Kulfan parameters and evaluated with NeuralFoil.

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
                     model_size::String="xlarge", weights_dir=nothing,
                     n_crit=9.0, fit_method::KulfanFitMethod=LeastSquaresFit(),
                     reuse_valid_airfoils::Bool=true, max_thickness_ratio::Real=2.0,
                     spanwise_direction=[0.0, 1.0, 0.0],
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
    load_neuralfoil_model(model_size; weights_dir)

    raw_sections = perpendicular_sections(vertices, faces, n_sections)
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
        result = generate_polar_from_dat(dat_abs, joinpath(output_dir, csv_rel);
            Re=Float64(Re), alpha_range=alphas, model_size, weights_dir, n_crit)
        push!(airfoil_rows, Any[j, "polar_vectors",
                                Dict("dat_file" => dat_rel, "raw_dat_file" => raw_rel,
                                     "csv_file_path" => csv_rel)])
        verbose && println("  Airfoil $j: CL_max=$(round(maximum(result.CL), digits=2))")
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
    obj_to_matrix_yaml(obj_path, foil_dat, output_dir; n_sections, wind_vel,
                       alpha_range, delta_range, crease_frac=0.2, remove_nan=true,
                       verbose=true) -> yaml_path

Convert an `.obj` wing mesh plus a single airfoil `.dat` into a standard geometry YAML
backed by one shared XFoil `(alpha, delta)` `POLAR_MATRICES`. Section leading/trailing
edges come from perpendicular mesh slices ([`perpendicular_sections`](@ref)); the
airfoil's cl/cd/cm matrices are generated once with [`create_polars`](@ref) and
referenced by every section. Load the result with `Wing(yaml_path; n_panels)`. This
replaces the old live `ObjWing` route with convert-then-load.
"""
function obj_to_matrix_yaml(obj_path::String, foil_dat::String, output_dir::String;
        n_sections::Int, wind_vel::Real, alpha_range, delta_range,
        crease_frac=0.2, remove_nan=true, verbose=true)
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    isfile(foil_dat) || error("DAT file not found: $foil_dat")
    polar_dir = joinpath(output_dir, "polars")
    airfoil_dir = joinpath(output_dir, "airfoils")
    mkpath(polar_dir)
    mkpath(airfoil_dir)

    vertices, faces = read_faces(obj_path)
    secs = perpendicular_sections(vertices, faces, n_sections)
    isempty(secs) && error("No valid sections sliced from $obj_path")

    mean_chord = sum(norm(s.TE_point .- s.LE_point) for s in secs) / length(secs)

    cl_path = joinpath(polar_dir, "cl.csv")
    cd_path = joinpath(polar_dir, "cd.csv")
    cm_path = joinpath(polar_dir, "cm.csv")
    x, y = read_dat_coordinates(foil_dat)
    write_dat(joinpath(output_dir, "airfoils", "1.dat"), "airfoil_1", x, y)
    create_polars(; dat_path=foil_dat, cl_polar_path=cl_path, cd_polar_path=cd_path,
        cm_polar_path=cm_path, wind_vel=Float64(wind_vel), area=mean_chord, width=1.0,
        crease_frac, alpha_range, delta_range, remove_nan)

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
    YAML.write_file(yaml_out, data)
    verbose && @info "Resolved geometry -> $yaml_out"
    return yaml_out
end

"""
    write_geometry_yaml(path, section_rows, airfoil_rows)

Write a geometry YAML in compact flow style: one-line headers and one line per
data row. `section_rows` are `[airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]`;
`airfoil_rows` are `[airfoil_id, type, info_dict]` where `info_dict` holds
`dat_file`, `csv_file_path`, and optionally `raw_dat_file`.
"""
function write_geometry_yaml(path::String, section_rows, airfoil_rows)
    open(path, "w") do io
        println(io, "wing_sections:")
        println(io, "  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]")
        println(io, "  data:")
        for row in section_rows
            vals = [v isa Integer ? string(v) : string(round(v; digits=3)) for v in row]
            println(io, "    - [", join(vals, ", "), "]")
        end
        println(io, "wing_airfoils:")
        println(io, "  headers: [airfoil_id, type, info_dict]")
        println(io, "  data:")
        for (id, type, info) in airfoil_rows
            keys_ordered = ["dat_file", "raw_dat_file", "csv_file_path",
                            "cl_file_path", "cd_file_path", "cm_file_path"]
            parts = ["$k: \"$(info[k])\"" for k in keys_ordered if haskey(info, k)]
            inline = "{" * join(parts, ", ") * "}"
            println(io, "    - [", id, ", ", type, ", ", inline, "]")
        end
    end
    return path
end

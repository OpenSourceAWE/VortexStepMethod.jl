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
                aero_solver=NeuralFoilSolver(), wrap_method=ShrinkWrap(),
                spanwise_direction=[0.0, 1.0, 0.0], crease_frac=0.75, verbose=true)

Convert a 3D wing `.obj` mesh to the native YAML geometry route.

Stations are placed at equal leading-edge arc-length intervals and sliced
perpendicular to the local span (see [`perpendicular_sections`](@ref)), which
keeps the airfoil undistorted near curved tips; each shape is then shrink-wrapped
into a clean airfoil and evaluated with `aero_solver`.

`aero_solver` selects the 2D-airfoil backend: [`NeuralFoilSolver`](@ref) (default,
fast) or [`XFoilSolver`](@ref) (viscous panel code); pass `aero_solver=XFoilSolver()`
to use XFoil instead. Each section's polar is written as `POLAR_VECTORS`.

`crease_frac` is the chordwise hinge location (0–1) about which each `delta_range`
trailing-edge deflection pivots.

With `force=false` (default) an existing `geometry.yaml` in `output_dir` is reused;
`force=true` regenerates it (e.g. after changing `delta_range` or the mesh).

`wrap_method` ([`ShrinkWrap`](@ref)) wraps each slice's point cloud into a clean closed
airfoil, robust both to the noisy interior-structure points (ribs, spars) of a ram-air
kite slice that otherwise pull a plain least-squares fit inward, and to the complex,
harsh corners of an LEI kite slice; NeuralFoil then fits Kulfan parameters to it and
XFoil uses it directly.

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
- `airfoils/{j}.dat`     — shrink-wrapped airfoil coordinates (matches the polar)
- `airfoils/{j}_raw.dat` — raw sliced section points (the wrap encloses these)
- `airfoils/{j}_d{deg}.dat` — each trailing-edge-deflected shape (with a `delta_range`)
- `polars/{j}.csv`    — NeuralFoil polar (alpha, Cd, Cs, Cl, Cm)
- `geometry.yaml`     — `wing_sections` + `wing_airfoils` referencing the above

# Returns
- Path to the written `geometry.yaml`.
"""
function obj_to_yaml(obj_path::String, output_dir::String;
                     n_sections::Int, Re::Real,
                     alpha_range=-180:1:180, delta_range=nothing,
                     aero_solver::AbstractAirfoilSolver=NeuralFoilSolver(),
                     wrap_method::ShrinkWrap=ShrinkWrap(),
                     reuse_valid_airfoils::Bool=true, max_thickness_ratio::Real=2.0,
                     spanwise_direction=[0.0, 1.0, 0.0], rotation=I,
                     wingtip_distance=0.05, crease_frac=0.75, force::Bool=false,
                     verbose::Bool=true)
    (!endswith(obj_path, ".obj")) && (obj_path *= ".obj")
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) &&
        throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0]"))

    yaml_path = joinpath(output_dir, "geometry.yaml")
    if !force && isfile(yaml_path)
        verbose && @info "Reusing existing geometry (force=true to regenerate)" yaml_path
        return yaml_path
    end

    vertices, faces = read_faces(obj_path)

    raw_sections = perpendicular_sections(vertices, faces, n_sections;
                                          rotation, wingtip_distance)
    isempty(raw_sections) && error("No valid sections sliced from $obj_path")

    stations = NamedTuple[]
    for (i, sec) in enumerate(raw_sections)
        verbose && print("  Section $i (y=$(round(sec.LE_point[2], digits=3)))... ")
        x_fit, y_fit = shrink_wrap(sec.x_airfoil, sec.y_airfoil, wrap_method)
        thickness = maximum(abs, y_fit)
        verbose && println("wrapped (thickness=$(round(thickness, digits=3)))")
        push!(stations, (; sec.LE_point, sec.TE_point, xa=sec.x_airfoil,
                         ya=sec.y_airfoil, x_fit, y_fit, thickness))
    end

    n = length(stations)
    thickness_limit = max_thickness_ratio * median(s.thickness for s in stations)
    degenerate = [s.thickness > thickness_limit for s in stations]
    valid = [k for k in 1:n if !degenerate[k]]
    isempty(valid) && error("All sliced sections are degenerate in $obj_path")
    ids = [degenerate[k] && reuse_valid_airfoils ?
           valid[argmin(abs.(valid .- k))] : k for k in 1:n]
    airfoils = [(; id = j, x_fit = stations[j].x_fit, y_fit = stations[j].y_fit,
                 x_raw = stations[j].xa, y_raw = stations[j].ya) for j in unique(ids)]
    airfoil_rows, ok = generate_airfoils(airfoils, output_dir; Re, alpha_range,
        delta_range, aero_solver, reuse_valid_airfoils, crease_frac, verbose)
    isempty(ok) && error("No section produced a valid polar in $obj_path")

    # Each section uses its nearest airfoil that actually produced a polar — covering
    # both too-thick (degenerate) fits and sections the solver could not converge.
    final_ids = [ok[argmin(abs.(ok .- k))] for k in 1:n]
    reused = [k => final_ids[k] for k in 1:n if final_ids[k] != k]
    isempty(reused) || @warn "Reused the nearest valid airfoil for: " *
        join(["section $k → airfoil $j" for (k, j) in reused], ", ")
    section_rows = Vector{Any}[]
    for k in 1:n
        s = stations[k]
        push!(section_rows, Any[final_ids[k], s.LE_point[1], s.LE_point[2], s.LE_point[3],
                                s.TE_point[1], s.TE_point[2], s.TE_point[3]])
    end

    sort!(section_rows; by = row -> row[3])          # clean spanwise order (by LE_y)
    sort!(airfoil_rows; by = row -> row[1])          # airfoils by id
    write_geometry_yaml(yaml_path, section_rows, airfoil_rows)
    verbose && @info "Wrote geometry to $yaml_path ($(length(section_rows)) sections)"
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

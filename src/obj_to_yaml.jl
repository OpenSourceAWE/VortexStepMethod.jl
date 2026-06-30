"""
Adapter that converts a 3D wing `.obj` mesh into the package's native YAML
geometry route: per-section airfoil `.dat` files, NeuralFoil polar CSVs, and a
geometry YAML referencing them. After conversion everything follows the standard
`Wing(geometry_file)` path.
"""

using Printf: @sprintf

"""
    read_dat_coordinates(path) -> (x, y)

Read airfoil coordinates from a `.dat` file, skipping header and comment lines.
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
    slice_obj_section(vertices, faces, y)

Slice a mesh at spanwise position `y` and return
`(LE_point, TE_point, x_airfoil, y_airfoil)`:
the un-normalized 3D leading- and trailing-edge points (min-x and max-x of the
slice) and the normalized Selig airfoil coordinates. Returns `nothing` if the
slice is degenerate.
"""
function slice_obj_section(vertices, faces, y)
    segments = slice_mesh_at_y(vertices, faces, Float64(y))
    isempty(segments) && return nothing

    contour = order_segments_to_contour(segments)
    length(contour) < 5 && return nothing

    xs = [p[1] for p in contour]
    le_idx = argmin(xs)
    te_idx = argmax(xs)
    LE_point = [contour[le_idx][1], Float64(y), contour[le_idx][2]]
    TE_point = [contour[te_idx][1], Float64(y), contour[te_idx][2]]

    x_airfoil, y_airfoil = contour_to_airfoil(contour)
    isempty(x_airfoil) && return nothing
    x_airfoil, y_airfoil = reorder_airfoil_selig(x_airfoil, y_airfoil)

    return LE_point, TE_point, x_airfoil, y_airfoil
end

"""
    write_dat(filepath, name, x, y)

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
    obj_to_yaml(obj_path, output_dir; n_sections, Re, alpha_range=-180:1:180,
                model_size="xlarge", weights_dir=nothing, n_crit=9.0,
                fit_method=LeastSquaresFit(),
                spanwise_direction=[0.0, 1.0, 0.0], verbose=true)

Convert a 3D wing `.obj` mesh to the native YAML geometry route.

Slices the mesh at `n_sections` spanwise stations; for each station the
leading-/trailing-edge points come from the slice itself (min-x / max-x), and
the airfoil shape is fitted to Kulfan parameters and evaluated with NeuralFoil.

`fit_method` selects the Kulfan fit (see [`fit_kulfan_parameters`](@ref)). Pass
an [`EnvelopeFit`](@ref) to wrap each section tightly *around* the slice points;
this is robust to the noisy interior-structure points of a ram-air kite slice
(ribs, spars) that otherwise pull the default [`LeastSquaresFit`](@ref) inward.

# Output
Writes into `output_dir`:
- `airfoils/{i}.dat`     — fitted Kulfan airfoil coordinates (matches the polar)
- `airfoils/{i}_raw.dat` — raw sliced section points (the fit is wrapped around these)
- `polars/{i}.csv`    — NeuralFoil polar (alpha, Cd, Cs, Cl, Cm)
- `geometry.yaml`     — `wing_sections` + `wing_airfoils` referencing the above

# Returns
- Path to the written `geometry.yaml`.
"""
function obj_to_yaml(obj_path::String, output_dir::String;
                     n_sections::Int, Re::Real,
                     alpha_range=-180:1:180,
                     model_size::String="xlarge", weights_dir=nothing,
                     n_crit=9.0, fit_method::KulfanFitMethod=LeastSquaresFit(),
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
    y_coords = [v[2] for v in vertices]
    y_min, y_max = minimum(y_coords), maximum(y_coords)
    delta = (y_max - y_min) / n_sections
    y_positions = [y_min + delta / 2 + (i - 1) * delta for i in 1:n_sections]

    alphas = collect(Float64, alpha_range)
    load_neuralfoil_model(model_size; weights_dir)

    section_rows = Vector{Any}[]
    airfoil_rows = Vector{Any}[]

    for (i, y) in enumerate(y_positions)
        verbose && print("  Section $i (y=$(round(y, digits=3)))... ")
        slice = slice_obj_section(vertices, faces, y)
        if slice === nothing
            verbose && println("degenerate, skipped")
            continue
        end
        LE_point, TE_point, xa, ya = slice

        dat_rel = joinpath("airfoils", "$i.dat")
        raw_rel = joinpath("airfoils", "$(i)_raw.dat")
        csv_rel = joinpath("polars", "$i.csv")

        params = fit_kulfan_parameters(xa, ya, fit_method)
        x_fit, y_fit = kulfan_to_coordinates(params)
        x_raw, y_raw = inset_airfoil(xa, ya, fit_clearance(fit_method))
        write_dat(joinpath(output_dir, dat_rel), "section_$i", x_fit, y_fit)
        write_dat(joinpath(output_dir, raw_rel), "section_$(i)_raw", x_raw, y_raw)

        result = neuralfoil_aero(params, alphas, Float64(Re);
                                 model_size, weights_dir, n_crit)
        write_polar_csv(joinpath(output_dir, csv_rel), result)

        push!(section_rows, Any[i, LE_point[1], LE_point[2], LE_point[3],
                                TE_point[1], TE_point[2], TE_point[3]])
        push!(airfoil_rows, Any[i, "polar_vectors",
                                Dict("dat_file" => dat_rel, "raw_dat_file" => raw_rel,
                                     "csv_file_path" => csv_rel)])
        verbose && println("done (CL_max=$(round(maximum(result.CL), digits=2)))")
    end

    isempty(section_rows) && error("No valid sections sliced from $obj_path")

    yaml_path = joinpath(output_dir, "geometry.yaml")
    write_geometry_yaml(yaml_path, section_rows, airfoil_rows)
    verbose && @info "Wrote geometry to $yaml_path ($(length(section_rows)) sections)"
    return yaml_path
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
            keys_ordered = ["dat_file", "raw_dat_file", "csv_file_path"]
            parts = ["$k: \"$(info[k])\"" for k in keys_ordered if haskey(info, k)]
            inline = "{" * join(parts, ", ") * "}"
            println(io, "    - [", id, ", ", type, ", ", inline, "]")
        end
    end
    return path
end

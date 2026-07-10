# Data structures for YAML wing geometry
@with_kw struct WingAirfoilInfo
    csv_file_path::String
    dat_file::String = ""
    cp_file_path::String = ""
    cl_file_path::String = ""
    cd_file_path::String = ""
    cm_file_path::String = ""
end

@with_kw struct WingSectionData
    airfoil_id::Int64
    LE_x::Float64
    LE_y::Float64
    LE_z::Float64
    TE_x::Float64
    TE_y::Float64
    TE_z::Float64
end

@with_kw struct WingAirfoilData
    airfoil_id::Int64
    type::String
    info_dict::WingAirfoilInfo
end

"""
    load_polar_data(csv_file_path::String) -> Tuple{Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}, Symbol}

Load aerodynamic polar data from a CSV file using only `readlines`.

The CSV file must contain a header row with columns for `alpha`, `cl`, `cd`, and `cm` (case-insensitive, order arbitrary). Each subsequent row should contain numeric values for these columns.

# Arguments
- `csv_file_path::String`: Path to the CSV file containing polar data.

# Returns
- A tuple `(aero_data, model_type)` where:
    - `aero_data`: A tuple of vectors `(alpha, cl, cd, cm)` if the file is valid, or `nothing` if invalid or missing.
        - `alpha`: Angle of attack in degrees (converted to radians internally).
        - `cl`: Lift coefficient.
        - `cd`: Drag coefficient.
        - `cm`: Moment coefficient.
    - `model_type`: `POLAR_VECTORS` if data is loaded, or `INVISCID` if not.

# Behavior
- If the file is missing, empty, or invalid, a warning is issued and `(nothing, INVISCID)` is returned.

# Example
```julia
# Create a YAML-based wing from configuration file
wing = Wing(
    "path/to/wing_config.yaml";
    n_panels=40,
    n_groups=4
)
```
"""
function assemble_polar_matrix(dv::Dict{String,Vector{Float64}})
    alphas = sort(unique(dv["alpha"]))
    deltas = sort(unique(dv["delta"]))
    ai = Dict(a => i for (i, a) in enumerate(alphas))
    di = Dict(d => j for (j, d) in enumerate(deltas))
    cl = fill(NaN, length(alphas), length(deltas))
    cd = fill(NaN, length(alphas), length(deltas))
    cm = fill(NaN, length(alphas), length(deltas))
    for k in eachindex(dv["alpha"])
        i, j = ai[dv["alpha"][k]], di[dv["delta"][k]]
        cl[i, j], cd[i, j], cm[i, j] = dv["cl"][k], dv["cd"][k], dv["cm"][k]
    end
    return (deg2rad.(alphas), deg2rad.(deltas), cl, cd, cm)
end

function load_polar_data(csv_file_path::String)
    # Return early for empty path
    if isempty(csv_file_path)
        @warn "Empty CSV file path provided"
        return (nothing, INVISCID)
    end

    # Check if file exists
    if !isfile(csv_file_path)
        @warn "CSV file not found: $csv_file_path"
        return (nothing, INVISCID)
    end

    try
        # Read all lines from the file
        lines = readlines(csv_file_path)
        
        # Check if file is empty
        if isempty(lines)
            @warn "CSV file is empty: $csv_file_path"
            return (nothing, INVISCID)
        end

        # Parse header - make case insensitive
        header_line = strip(lines[1])
        if isempty(header_line)
            @warn "CSV file has empty header: $csv_file_path"
            return (nothing, INVISCID)
        end

        # Split header and normalize to lowercase. A `delta` column marks the
        # long-format `(alpha, delta)` polar, loaded as POLAR_MATRICES.
        header_parts = map(strip ∘ lowercase, split(header_line, ','))
        has_delta = "delta" in header_parts
        wanted = has_delta ? ["alpha", "delta", "cl", "cd", "cm"] :
                             ["alpha", "cl", "cd", "cm"]
        col_indices = Dict{String, Int}()
        for (i, col_name) in enumerate(header_parts)
            col_name in wanted && (col_indices[col_name] = i)
        end

        missing_cols = setdiff(wanted, keys(col_indices))
        if !isempty(missing_cols)
            @warn "CSV file missing required columns: $(join(missing_cols, ", ")) in $csv_file_path"
            return (nothing, INVISCID)
        end

        # Parse data rows
        data_vectors = Dict(col => Float64[] for col in wanted)
        for line_num in eachindex(lines)
            line_num == firstindex(lines) && continue
            line = strip(lines[line_num])
            isempty(line) && continue  # Skip empty lines

            parts = split(line, ',')
            if length(parts) != length(header_parts)
                @warn "Line $line_num has incorrect number of columns in $csv_file_path"
                return (nothing, INVISCID)  # Strict: any malformed data rejects the file
            end

            try
                for col in wanted
                    push!(data_vectors[col], parse(Float64, strip(parts[col_indices[col]])))
                end
            catch e
                @warn "Failed to parse line $line_num in $csv_file_path: $e"
                return (nothing, INVISCID)  # Strict: any parsing error rejects the file
            end
        end

        # Check if we got any valid data
        if isempty(data_vectors["alpha"])
            @warn "No valid data rows found in $csv_file_path"
            return (nothing, INVISCID)
        end

        has_delta && return (assemble_polar_matrix(data_vectors), POLAR_MATRICES)

        # Convert alpha from degrees to radians and create tuple
        alpha_rad = deg2rad.(data_vectors["alpha"])
        aero_data = (alpha_rad, data_vectors["cl"], data_vectors["cd"], data_vectors["cm"])

        return (aero_data, POLAR_VECTORS)

    catch e
        @warn "Error reading CSV file $csv_file_path: $e"
        return (nothing, INVISCID)
    end
end

"""
    load_matrix_polar_data(cl_path, cd_path, cm_path) -> (aero_data, POLAR_MATRICES)

Read the three `(alpha × delta)` coefficient matrices (see [`read_aero_matrix`](@ref))
and assemble the `POLAR_MATRICES` `aero_data = (alpha, delta, cl, cd, cm)`.
"""
function load_matrix_polar_data(cl_path::String, cd_path::String, cm_path::String)
    cl_matrix, alpha, delta = read_aero_matrix(cl_path)
    cd_matrix, _, _ = read_aero_matrix(cd_path)
    cm_matrix, _, _ = read_aero_matrix(cm_path)
    return (collect(alpha), collect(delta), cl_matrix, cd_matrix, cm_matrix),
        POLAR_MATRICES
end

"""
    Wing(geometry_file::String; n_panels=20, spanwise_distribution=LINEAR,
         spanwise_direction=[0.0, 1.0, 0.0], remove_nan=true,
         use_prior_polar=false, billowing_percentage=0.0, prn=false)

Constructs a `Wing` object from a YAML geometry file.

# Arguments
- `geometry_file::String`: Path to a YAML file describing the wing geometry and airfoils.

# Keyword Arguments
- `n_panels::Int`: Number of spanwise panels (default: 20).
- `spanwise_distribution`: Spanwise panel distribution type (default: `LINEAR`).
- `spanwise_direction::Vector{Float64}`: Direction of the spanwise axis (default: `[0.0, 1.0, 0.0]`). Must be the global Y axis.
- `remove_nan::Bool`: Remove NaN values from the geometry (default: `true`).
- `use_prior_polar::Bool`: Reuse prior refined/panel polar mapping on geometry updates (default: `false`).
- `billowing_percentage::Float64`: TE billow as percentage of arc length (default: `0.0`).
- `prn::Bool`: Print informational messages during construction (default: `false`).

# Returns
- `Wing`: A fully constructed and initialized `Wing` object.

# Description
This function reads a YAML configuration file to define the geometry and airfoil data for a multi-section wing.
Each section and corresponding airfoil is parsed from the YAML file, polar data is loaded, and each section is added
to the wing. The geometry logic currently assumes the spanwise direction is `[0.0, 1.0, 0.0]` (aligned with the global Y axis).
The number of unrefined sections is automatically inferred from the sections in the geometry file.

# Errors
- Throws an `ArgumentError` if `spanwise_direction` is not `[0.0, 1.0, 0.0]`.

# Example
```julia
wing = Wing("wing_geometry.yaml"; n_panels=30, prn=true)
```
"""
function Wing(
    geometry_file::String;
    n_panels=20,
    spanwise_distribution=LINEAR,
    spanwise_direction=[0.0, 1.0, 0.0],
    remove_nan=true,
    use_prior_polar=false,
    billowing_percentage=0.0,
    prn=false,
    sort_sections=true
)

    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))

    prn && @info "Reading YAML wing configuration from $geometry_file"
    
    # Load YAML file following Uwe's suggestion
    data = YAML.load_file(geometry_file)

    # Convert wing sections
    wing_sections_data = data["wing_sections"]
    sections = WingSectionData[]
    for row in wing_sections_data["data"]
        section_dict = Dict(zip(wing_sections_data["headers"], row))
        push!(sections, WingSectionData(
            airfoil_id = section_dict["airfoil_id"],
            LE_x = section_dict["LE_x"],
            LE_y = section_dict["LE_y"], 
            LE_z = section_dict["LE_z"],
            TE_x = section_dict["TE_x"],
            TE_y = section_dict["TE_y"],
            TE_z = section_dict["TE_z"]
        ))
    end
    
    # Convert wing airfoils
    wing_airfoils_data = data["wing_airfoils"]
    airfoils = WingAirfoilData[]
    airfoil_poly_map = Dict{Int64, NTuple{3, Vector{Float64}}}()
    airfoil_type_map = Dict{Int64, String}()
    for row in wing_airfoils_data["data"]
        airfoil_dict = Dict(zip(wing_airfoils_data["headers"], row))
        airfoil_type_map[airfoil_dict["airfoil_id"]] = String(airfoil_dict["type"])
        info = airfoil_dict["info_dict"]
        if haskey(info, "cl_coeffs")
            airfoil_poly_map[airfoil_dict["airfoil_id"]] = (
                Float64.(info["cl_coeffs"]), Float64.(info["cd_coeffs"]),
                Float64.(info["cm_coeffs"]))
        end
        push!(airfoils, WingAirfoilData(
            airfoil_id = airfoil_dict["airfoil_id"],
            type = airfoil_dict["type"],
            info_dict = WingAirfoilInfo(
                csv_file_path = get(airfoil_dict["info_dict"], "csv_file_path", ""),
                dat_file = get(airfoil_dict["info_dict"], "dat_file", ""),
                cp_file_path = get(airfoil_dict["info_dict"], "cp_file_path", ""),
                cl_file_path = get(airfoil_dict["info_dict"], "cl_file_path", ""),
                cd_file_path = get(airfoil_dict["info_dict"], "cd_file_path", ""),
                cm_file_path = get(airfoil_dict["info_dict"], "cm_file_path", ""))
        ))
    end

    # Create CSV file mapping from airfoils
    airfoil_csv_map = Dict{Int64, String}()
    airfoil_cp_map = Dict{Int64, String}()
    airfoil_matrix_map = Dict{Int64, NTuple{3, String}}()
    for airfoil in airfoils
        if !isempty(airfoil.info_dict.csv_file_path)
            airfoil_csv_map[airfoil.airfoil_id] = airfoil.info_dict.csv_file_path
        end
        if !isempty(airfoil.info_dict.cp_file_path)
            airfoil_cp_map[airfoil.airfoil_id] = airfoil.info_dict.cp_file_path
        end
        if !isempty(airfoil.info_dict.cl_file_path)
            airfoil_matrix_map[airfoil.airfoil_id] = (airfoil.info_dict.cl_file_path,
                airfoil.info_dict.cd_file_path, airfoil.info_dict.cm_file_path)
        end
    end
    
    # n_unrefined_sections is set automatically as sections are added
    wing = Wing(n_panels;
        spanwise_distribution=spanwise_distribution,
        spanwise_direction=MVec3(spanwise_direction),
        remove_nan=remove_nan,
        use_prior_polar=use_prior_polar,
        billowing_percentage=Float64(billowing_percentage)
    )
    
    # Parse sections and populate wing
    for section in sections
        # Get coordinates directly from struct fields
        le_coord = [section.LE_x, section.LE_y, section.LE_z]
        te_coord = [section.TE_x, section.TE_y, section.TE_z]

        base_dir = dirname(geometry_file)
        resolve(p) = (!isempty(p) && !isabspath(p)) ? joinpath(base_dir, p) : p

        # Core accepts only resolved forms (poly/polars/inviscid); rich types resolve in AirfoilAero.
        airfoil_type = get(airfoil_type_map, section.airfoil_id, "")
        if haskey(airfoil_poly_map, section.airfoil_id)
            aero_data = airfoil_poly_map[section.airfoil_id]
            aero_model = POLY
        elseif haskey(airfoil_matrix_map, section.airfoil_id)
            cl_p, cd_p, cm_p = airfoil_matrix_map[section.airfoil_id]
            aero_data, aero_model = load_matrix_polar_data(
                resolve(cl_p), resolve(cd_p), resolve(cm_p))
        elseif airfoil_type in ("neuralfoil", "breukels_regression", "masure_regression")
            throw(ArgumentError("airfoil_id $(section.airfoil_id) has unresolved " *
                "type \"$airfoil_type\"; resolve it with AirfoilAero " *
                "(resolve_aero_geometry) into polars/poly before loading."))
        else
            csv_file_path = resolve(get(airfoil_csv_map, section.airfoil_id, ""))
            aero_data, aero_model = load_polar_data(csv_file_path)
        end

        cp_file_path = resolve(get(airfoil_cp_map, section.airfoil_id, ""))
        cp_data = isempty(cp_file_path) ? nothing : read_cp_data(cp_file_path)

        prn && println("Section airfoil_id $(section.airfoil_id): Using $aero_model model")

        add_section!(wing, le_coord, te_coord, aero_model, aero_data, cp_data)
    end

    refine!(wing; sort_sections)
    return wing
end

"""
    Wing(settings::VSMSettings)

Create a wing model from VSM settings configuration.

This constructor is a convenience wrapper that extracts wing configuration
from VSMSettings and creates a Wing using either:
- YAML geometry file (geometry_file field), or
- OBJ + DAT files (obj_file and dat_file fields)

The constructor automatically determines which path to use based on which
fields are populated in the settings.

# Arguments
- `settings`: VSMSettings object containing wing configuration

# Returns
A fully initialized `Wing` instance ready for aerodynamic simulation.

# Example
```julia
# Using YAML geometry
settings = VSMSettings("path/to/vsm_settings.yaml")
wing = Wing(settings)

# Settings can specify either:
# - geometry_file: "path/to/wing.yaml"  # YAML-based
# - obj_file + dat_file                  # OBJ-based
```
"""
function Wing(settings::VSMSettings; sort_sections::Bool=true)
    wing_settings = settings.wings[1]

    # Check which geometry format to use
    has_yaml = !isempty(wing_settings.geometry_file)
    has_obj = !isempty(wing_settings.obj_file)
    has_dat = !isempty(wing_settings.dat_file)

    if has_yaml && (has_obj || has_dat)
        throw(ArgumentError(
            "Cannot specify both geometry_file and obj_file/dat_file"
        ))
    end

    if has_obj && !has_dat
        throw(ArgumentError(
            "obj_file requires dat_file to be specified"
        ))
    end

    if has_dat && !has_obj
        throw(ArgumentError(
            "dat_file requires obj_file to be specified"
        ))
    end

    if has_yaml
        # Use YAML geometry constructor
        Wing(wing_settings.geometry_file;
            n_panels=wing_settings.n_panels,
            spanwise_distribution=wing_settings.spanwise_panel_distribution,
            remove_nan=wing_settings.remove_nan,
            use_prior_polar=wing_settings.use_prior_polar,
            billowing_percentage=wing_settings.billowing_percentage,
            sort_sections
        )
    elseif has_obj && has_dat
        throw(ArgumentError(
            "OBJ/DAT geometry is handled by the ObjAdapter package: convert to a " *
            "standard YAML with ObjAdapter first, then load it via geometry_file."
        ))
    else
        throw(ArgumentError(
            "WingSettings must specify either geometry_file or " *
            "both obj_file and dat_file"
        ))
    end
end

"""
    ObjWing(obj_path[, dat_path]; n_panels, Re, alpha_range, delta_range,
            n_sections, spanwise_direction, aero_solver, remake, output_dir, verbose) → Wing

Convenience constructor retained for backward compatibility. Converts an OBJ mesh
to a YAML wing geometry via [`ObjAdapter.obj_to_yaml`](@ref) and returns a [`Wing`](@ref).

`dat_path` is accepted but ignored — airfoil shapes are extracted directly from the
OBJ geometry. Use `aero_solver=AirfoilAero.XFoilSolver()` to reproduce old XFoil-based
polars; the default is `AirfoilAero.NeuralFoilSolver()`.

By default (`remake=false`) an existing `geometry.yaml` in `output_dir` is reused,
skipping the expensive polar generation. Set `remake=true` to force regeneration.
"""
function ObjWing(obj_path, dat_path=nothing;
                 n_panels::Int=56,
                 Re::Real=1e6,
                 alpha_range=deg2rad.(-5:1:20),
                 delta_range=deg2rad.(-5:1:20),
                 n_sections::Union{Nothing, Int}=nothing,
                 spanwise_direction=[0.0, 1.0, 0.0],
                 spanwise_distribution=UNCHANGED,
                 remove_nan::Bool=true,
                 aero_solver=AirfoilAero.NeuralFoilSolver(),
                 remake::Bool=false,
                 output_dir::String=mktempdir(),
                 verbose::Bool=false)
    yaml_path = joinpath(output_dir, "geometry.yaml")
    if remake || !isfile(yaml_path)
        n_sec = isnothing(n_sections) ? n_panels + 1 : n_sections
        yaml_path = ObjAdapter.obj_to_yaml(obj_path, output_dir;
                                           n_sections=n_sec,
                                           Re,
                                           alpha_range,
                                           delta_range,
                                           aero_solver,
                                           spanwise_direction,
                                           verbose)
    end
    return Wing(yaml_path; n_panels, spanwise_distribution, spanwise_direction, remove_nan)
end

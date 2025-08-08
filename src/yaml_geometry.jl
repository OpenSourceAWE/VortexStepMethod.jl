# Data structures for YAML wing geometry
@with_kw struct WingAirfoilInfo
    csv_file_path::String
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

        # Split header and normalize to lowercase
        header_parts = map(strip ∘ lowercase, split(header_line, ','))
        
        # Find column indices for required columns
        required_cols = ["alpha", "cl", "cd", "cm"]
        col_indices = Dict{String, Int}()
        
        for (i, col_name) in enumerate(header_parts)
            if col_name in required_cols
                col_indices[col_name] = i
            end
        end
        
        # Check if all required columns are present
        missing_cols = setdiff(required_cols, keys(col_indices))
        if !isempty(missing_cols)
            @warn "CSV file missing required columns: $(join(missing_cols, ", ")) in $csv_file_path"
            return (nothing, INVISCID)
        end

        # Parse data rows
        data_vectors = Dict{String, Vector{Float64}}()
        for col in required_cols
            data_vectors[col] = Float64[]
        end

        for line_num in 2:length(lines)
            line = strip(lines[line_num])
            isempty(line) && continue  # Skip empty lines
            
            parts = split(line, ',')
            if length(parts) != length(header_parts)
                @warn "Line $line_num has incorrect number of columns in $csv_file_path"
                return (nothing, INVISCID)  # Strict: any malformed data rejects the file
            end

            try
                for col in required_cols
                    value_str = strip(parts[col_indices[col]])
                    value = parse(Float64, value_str)
                    push!(data_vectors[col], value)
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
    # Module: yaml_geometry.jl

    This module provides functionality for parsing and handling geometry data from YAML files.
    It is intended for use within the VortexStepMethod.jl package to facilitate the import,
    validation, and manipulation of geometric configurations specified in YAML format.

    Functions and types in this module allow users to:
    - Load geometry definitions from YAML files.
    - Convert YAML data into Julia-native structures.
    - Validate the integrity and consistency of imported geometry data.

    Typical use cases include reading wing or blade geometry for aerodynamic simulations.

    # Example
    ```julia
    geometry = load_yaml_geometry("geometry.yaml")
    ```

    # See Also
    - [`YAML`](@ref)
    - [`VortexStepMethod`](@ref)
"""
function Wing(
    geometry_file::String;
    n_panels=20,
    n_groups=1,
    spanwise_distribution=LINEAR,
    spanwise_direction=[0.0, 1.0, 0.0],
    remove_nan=true,
    prn=false
)
    !(n_panels % n_groups == 0) && throw(ArgumentError("Number of panels should be divisible by number of groups"))
    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))
    
    prn && @info "Reading YAML wing configuration from $geometry_file"
    
    # Load YAML file following Uwe's suggestion
    data = YAML.load_file(geometry_file)

    # Convert YAML data to our struct format
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
    for row in wing_airfoils_data["data"]
        airfoil_dict = Dict(zip(wing_airfoils_data["headers"], row))
        push!(airfoils, WingAirfoilData(
            airfoil_id = airfoil_dict["airfoil_id"],
            type = airfoil_dict["type"],
            info_dict = WingAirfoilInfo(csv_file_path = get(airfoil_dict["info_dict"], "csv_file_path", ""))
        ))
    end
    
    # Create CSV file mapping from airfoils
    airfoil_csv_map = Dict{Int64, String}()
    for airfoil in airfoils
        if !isempty(airfoil.info_dict.csv_file_path)
            airfoil_csv_map[airfoil.airfoil_id] = airfoil.info_dict.csv_file_path
        end
    end
    
    # Create Wing using the standard constructor
    wing = Wing(n_panels; 
        n_groups=n_groups, 
        spanwise_distribution=spanwise_distribution,
        spanwise_direction=MVec3(spanwise_direction), 
        remove_nan=remove_nan
    )
    
    # Parse sections and populate wing
    for section in sections
        # Get coordinates directly from struct fields
        le_coord = [section.LE_x, section.LE_y, section.LE_z]
        te_coord = [section.TE_x, section.TE_y, section.TE_z]
        
        # Load polar data and create section
        csv_file_path = get(airfoil_csv_map, section.airfoil_id, "")
        if !isempty(csv_file_path) && !isabspath(csv_file_path)
            # NOTE: The spanwise direction is currently restricted to [0.0, 1.0, 0.0] (the global Y axis).
            # This is required because downstream geometry and panel generation code assumes the spanwise axis is aligned with Y.
            # If you need to support arbitrary spanwise directions, refactor the geometry logic accordingly.
            csv_file_path = joinpath(dirname(geometry_file), csv_file_path)
        end
        aero_data, aero_model = load_polar_data(csv_file_path)
        
        prn && println("Section airfoil_id $(section.airfoil_id): Using $aero_model model")
        
        add_section!(wing, le_coord, te_coord, aero_model, aero_data)
    end
    
    # Initialize the wing after adding all sections
    reinit!(wing)
    
    return wing
end

"""
    Wing(settings::VSMSettings)

Create a wing model from VSM settings configuration.

This constructor is a convenience wrapper that extracts wing configuration 
from VSMSettings and creates a Wing using the YAML geometry file path and 
parameters specified in the settings.

# Arguments
- `settings`: VSMSettings object containing wing configuration

# Returns
A fully initialized `Wing` instance ready for aerodynamic simulation.

# Example
```julia
# Load settings and create wing in one step
settings = VSMSettings("path/to/vsm_settings.yaml")
wing = Wing(settings)
```
"""
function Wing(settings::VSMSettings)
    Wing(settings.wings[1].geometry_file; 
        n_panels=settings.wings[1].n_panels, 
        n_groups=settings.wings[1].n_groups,
        spanwise_distribution=settings.wings[1].spanwise_panel_distribution
    )
end

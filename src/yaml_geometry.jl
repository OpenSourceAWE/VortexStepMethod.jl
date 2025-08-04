@with_kw mutable struct AirfoilSettings
    airfoil_id::Int64
    type::String
    info_dict::Dict{String, Any}
end

"""
    YamlWing <: AbstractWing

A wing model created from YAML configuration files with CSV polar data.

## Core Features
- Wing geometry defined through YAML section coordinates
- Aerodynamic properties based on CSV polar data files
- Support for multiple airfoil types per wing
- Configurable panel distribution and geometry parameters

## Notable Fields
- `n_panels::Int16`: Number of panels in aerodynamic mesh
- `n_groups::Int16`: Number of control groups
- `spanwise_distribution::PanelDistribution`: Panel distribution type
- `sections::Vector{Section}`: Wing cross-sections with aerodynamic data
- `airfoil_settings::Vector{AirfoilSettings}`: Airfoil configuration data
- `alpha_range::Vector{Float64}`: Angle of attack range for polars
- `reynolds::Float64`: Reynolds number for aerodynamic analysis

See constructor `YamlWing(yaml_path; kwargs...)` for usage details.
"""
mutable struct YamlWing <: AbstractWing
    n_panels::Int16
    n_groups::Int16
    spanwise_distribution::PanelDistribution
    panel_props::PanelProperties
    spanwise_direction::MVec3
    sections::Vector{Section}
    refined_sections::Vector{Section}
    remove_nan::Bool
    
    # Additional fields for YamlWing
    airfoil_settings::Vector{AirfoilSettings}
    alpha_range::Vector{Float64}
    reynolds::Float64
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}
end

"""
    load_polar_data(csv_file_path)

Load polar data from CSV file at the given path.
Expected CSV format: columns named alpha, cl, cd, cm (case insensitive)

Returns (aero_data, aero_model) tuple.
"""
function load_polar_data(csv_file_path::String)
    if isempty(csv_file_path) || !isfile(csv_file_path)
        @warn "Polar file not found or empty path: $csv_file_path. Using INVISCID model instead."
        return nothing, INVISCID
    end
    
    try
        # Dynamically load CSV and DataFrames if not already loaded
        if !isdefined(Main, :CSV)
            @eval Main using CSV
        end
        if !isdefined(Main, :DataFrames)
            @eval Main using DataFrames
        end
        
        # Read CSV file
        df = Main.DataFrame(Main.CSV.File(csv_file_path))
        
        # Get column names and find required columns (case insensitive)
        col_names = names(df)
        col_names_lower = lowercase.(col_names)
        
        # Find required columns by matching lowercase names
        alpha_col = findfirst(x -> x == "alpha", col_names_lower)
        cl_col = findfirst(x -> x == "cl", col_names_lower)
        cd_col = findfirst(x -> x == "cd", col_names_lower)
        cm_col = findfirst(x -> x == "cm", col_names_lower)
        
        # Check if all required columns are found
        if isnothing(alpha_col) || isnothing(cl_col) || isnothing(cd_col) || isnothing(cm_col)
            missing_cols = String[]
            isnothing(alpha_col) && push!(missing_cols, "alpha")
            isnothing(cl_col) && push!(missing_cols, "cl")
            isnothing(cd_col) && push!(missing_cols, "cd")
            isnothing(cm_col) && push!(missing_cols, "cm")
            @warn "Missing columns $missing_cols in $csv_file_path. Using INVISCID model instead."
            return nothing, INVISCID
        end
        
        # Extract data using actual column names
        alpha_data = df[!, col_names[alpha_col]]
        cl_data = df[!, col_names[cl_col]]
        cd_data = df[!, col_names[cd_col]]
        cm_data = df[!, col_names[cm_col]]
        
        # Create aero_data tuple
        aero_data = (alpha_data, cl_data, cd_data, cm_data)
        return aero_data, POLAR_VECTORS
        
    catch e
        @warn "Error reading polar file $csv_file_path: $e. Using INVISCID model instead."
        return nothing, INVISCID
    end
end

"""
    YamlWing(yaml_path; kwargs...)

Create a wing model from YAML configuration file with CSV polar data.

This constructor builds a complete aerodynamic model by:
1. Loading wing geometry and airfoil configurations from YAML file
2. Loading aerodynamic polars from CSV files specified in the YAML
3. Creating wing sections with appropriate aerodynamic models
4. Setting up panel distribution and geometric properties

# Arguments
- `yaml_path`: Path to YAML file containing wing geometry and airfoil specifications

# Keyword Arguments
- `n_panels=20`: Number of aerodynamic panels across wingspan
- `n_groups=4`: Number of control groups
- `spanwise_distribution=LINEAR`: Panel distribution type
- `spanwise_direction=[0.0, 1.0, 0.0]`: Spanwise direction vector
- `remove_nan=true`: Interpolate NaN values in aerodynamic data
- `prn=true`: Print info messages during construction

# Returns
A fully initialized `YamlWing` instance ready for aerodynamic simulation.

# Example YAML format
```yaml
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 10.0, 0.0, 1.0, 10.0, 0.0]
    - [2, 0.0, -10.0, 0.0, 1.0, -10.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 31, 0.5]
  reynolds: 1e6
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "polars/1.csv"}]
    - [2, polars, {csv_file_path: "polars/2.csv"}]
```

# Example
```julia
# Create a YAML-based wing from configuration file
wing = YamlWing(
    "path/to/wing_config.yaml";
    n_panels=40,
    n_groups=4
)
```
"""
function YamlWing(
    yaml_path;
    n_panels=20,
    n_groups=4,
    spanwise_distribution=LINEAR,
    spanwise_direction=[0.0, 1.0, 0.0],
    remove_nan=true,
    prn=true
)
    !(n_panels % n_groups == 0) && throw(ArgumentError("Number of panels should be divisible by number of groups"))
    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))
    
    prn && @info "Reading YAML wing configuration from $yaml_path"
    
    # Read and parse YAML file
    yaml_data = YAML.load_file(yaml_path)
    
    # Parse airfoils and create CSV file mapping
    airfoil_csv_map = Dict{Int64, String}()
    airfoil_settings = AirfoilSettings[]
    
    for row in yaml_data["wing_airfoils"]["data"]
        airfoil_dict = Dict(zip(yaml_data["wing_airfoils"]["headers"], row))
        airfoil_id, airfoil_type, info_dict = airfoil_dict["airfoil_id"], airfoil_dict["type"], airfoil_dict["info_dict"]
        
        push!(airfoil_settings, AirfoilSettings(airfoil_id, airfoil_type, info_dict))
        haskey(info_dict, "csv_file_path") && (airfoil_csv_map[airfoil_id] = info_dict["csv_file_path"])
    end
    
    # Parse sections and create wing sections
    sections = Section[]
    refined_sections = Section[]
    
    for row in yaml_data["wing_sections"]["data"]
        section_dict = Dict(zip(yaml_data["wing_sections"]["headers"], row))
        airfoil_id = section_dict["airfoil_id"]
        
        # Get coordinates
        le_coord = [section_dict["LE_x"], section_dict["LE_y"], section_dict["LE_z"]]
        te_coord = [section_dict["TE_x"], section_dict["TE_y"], section_dict["TE_z"]]
        
        # Load polar data and create section
        csv_file_path = get(airfoil_csv_map, airfoil_id, "")
        if !isempty(csv_file_path) && !isabspath(csv_file_path)
            # Make relative paths relative to YAML file directory
            csv_file_path = joinpath(dirname(yaml_path), csv_file_path)
        end
        aero_data, aero_model = load_polar_data(csv_file_path)
        
        prn && println("Section airfoil_id $airfoil_id: Using $aero_model model")
        
        section = Section(le_coord, te_coord, aero_model, aero_data)
        push!(sections, section)
        push!(refined_sections, Section(le_coord, te_coord, aero_model, aero_data))
    end
    
    # Extract configuration
    alpha_range = get(yaml_data["wing_airfoils"], "alpha_range", [-10.0, 31.0, 0.5])
    reynolds = get(yaml_data["wing_airfoils"], "reynolds", 1e6)
    
    # Create YamlWing
    YamlWing(
        n_panels, n_groups, spanwise_distribution, PanelProperties{n_panels}(), 
        MVec3(spanwise_direction), sections, refined_sections, remove_nan,
        airfoil_settings, alpha_range, reynolds, [LazyBufferCache()]
    )
end

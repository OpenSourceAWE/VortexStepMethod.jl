

"""
    load_polar_data(csv_file_path)

Load polar data from a CSV file using only `readlines`.
Expected format: header row and columns alpha, cl, cd, cm (case-insensitive, order arbitrary).
- alpha: Angle of attack in degrees
- cl: Lift coefficient
- cd: Drag coefficient
- cm: Moment coefficient

Returns (aero_data, POLAR_VECTORS) or (nothing, INVISCID) if missing/invalid.
"""
function load_polar_data(csv_file_path::String)
    if isempty(csv_file_path) || !isfile(csv_file_path)
        @warn "Polar file not found or empty path: $csv_file_path. Using INVISCID model instead."
        return nothing, INVISCID
    end
    try
        lines = readlines(csv_file_path)
        isempty(lines) && throw(ArgumentError("File is empty"))
        header = lowercase.(split(strip(lines[1]), ','))

        # Create a mapping from column name to index (1-based for Julia arrays)
        col_map = Dict(name => i for (i, name) in enumerate(header))
        required = ["alpha", "cl", "cd", "cm"]
        missing = filter(x -> !haskey(col_map, x), required)
        if !isempty(missing)
            @warn "Missing columns $missing in $csv_file_path. Using INVISCID model instead."
            return nothing, INVISCID
        end

        # Preallocate arrays
        n = length(lines) - 1
        alpha = Vector{Float64}(undef, n)
        cl    = Vector{Float64}(undef, n)
        cd    = Vector{Float64}(undef, n)
        cm    = Vector{Float64}(undef, n)



        # Read each data line
        for (i, line) in enumerate(lines[2:end])
            fields = split(strip(line), ',')
            alpha[i] = parse(Float64, fields[col_map["alpha"]])
            cl[i]    = parse(Float64, fields[col_map["cl"]])
            cd[i]    = parse(Float64, fields[col_map["cd"]])
            cm[i]    = parse(Float64, fields[col_map["cm"]])
        end

        # Converting to radians
        alpha = deg2rad.(alpha)

        return (alpha, cl, cd, cm), POLAR_VECTORS
    catch e
        @warn "Error reading polar file $csv_file_path: $e. Using INVISCID model instead."
        return nothing, INVISCID
    end
end


"""
    YamlWing(geometry_file; kwargs...)

Create a wing model from YAML configuration file with CSV polar data.

This constructor builds a complete aerodynamic model by:
1. Loading wing geometry and airfoil configurations from YAML file
2. Loading aerodynamic polars from CSV files specified in the YAML
3. Creating wing sections with appropriate aerodynamic models
4. Setting up panel distribution and geometric properties

# Arguments
- `geometry_file`: Path to YAML file containing wing geometry and airfoil specifications

# Keyword Arguments
- `n_panels=20`: Number of aerodynamic panels across wingspan
- `n_groups=4`: Number of control groups
- `spanwise_distribution=LINEAR`: Panel distribution type
- `spanwise_direction=[0.0, 1.0, 0.0]`: Spanwise direction vector
- `remove_nan=true`: Interpolate NaN values in aerodynamic data
- `prn=true`: Print info messages during construction

# Returns
A fully initialized `Wing` instance ready for aerodynamic simulation.

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
    geometry_file;
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
    
    # Read and parse YAML file
    yaml_data = YAML.load_file(geometry_file)
    
    # Parse airfoils and create CSV file mapping
    airfoil_csv_map = Dict{Int64, String}()
    
    for row in yaml_data["wing_airfoils"]["data"]
        airfoil_dict = Dict(zip(yaml_data["wing_airfoils"]["headers"], row))
        airfoil_id, airfoil_type, info_dict = airfoil_dict["airfoil_id"], airfoil_dict["type"], airfoil_dict["info_dict"]
        
        haskey(info_dict, "csv_file_path") && (airfoil_csv_map[airfoil_id] = info_dict["csv_file_path"])
    end
    
    # Create Wing using the standard constructor
    wing = Wing(n_panels; 
        n_groups=n_groups, 
        spanwise_distribution=spanwise_distribution,
        spanwise_direction=MVec3(spanwise_direction), 
        remove_nan=remove_nan
    )
    
    # Parse sections and populate wing
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
            csv_file_path = joinpath(dirname(geometry_file), csv_file_path)
        end
        aero_data, aero_model = load_polar_data(csv_file_path)
        
        prn && println("Section airfoil_id $airfoil_id: Using $aero_model model")
        
        add_section!(wing, le_coord, te_coord, aero_model, aero_data)
    end
    
    # Initialize the wing after adding all sections
    reinit!(wing)
    
    return wing
end

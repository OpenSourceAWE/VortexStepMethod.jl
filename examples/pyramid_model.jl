using LinearAlgebra
using VortexStepMethod
using ControlPlots
using YAML
using CSV
using DataFrames


function read_yaml(yaml_file_path)
    # Read the YAML file to get the sections
    yaml_data = YAML.load_file(yaml_file_path)
    
    # Extract headers and data
    headers = yaml_data["wing_sections"]["headers"]
    data_rows = yaml_data["wing_sections"]["data"]
    
    # Convert YAML data to the expected format
    sections = []
    for row in data_rows
        # Create a dictionary mapping headers to values
        section_dict = Dict(zip(headers, row))
        
        # Keep airfoil_id as integer for CSV filename lookup
        airfoil_id = section_dict["airfoil_id"]
        
        # Extract coordinates
        le_coord = [section_dict["LE_x"], section_dict["LE_y"], section_dict["LE_z"]]
        te_coord = [section_dict["TE_x"], section_dict["TE_y"], section_dict["TE_z"]]
        
        section = (
            airfoil_id = airfoil_id,
            le_coord = le_coord,
            te_coord = te_coord
        )
        push!(sections, section)
    end

    return sections
end

function load_polar_data(airfoil_id, polars_dir)
    """
    Load polar data from CSV file for given airfoil_id.
    Expected CSV format: columns named alpha, cl, cd, cm
    """
    csv_filename = "$(airfoil_id).csv"
    csv_filepath = joinpath(polars_dir, csv_filename)
    
    if !isfile(csv_filepath)
        @warn "Polar file not found: $csv_filepath. Using INVISCID model instead."
        return nothing, INVISCID
    end
    
    try
        df = CSV.read(csv_filepath, DataFrame)
        
        # Verify required columns exist
        required_cols = ["alpha", "cl", "cd", "cm"]
        missing_cols = setdiff(required_cols, names(df))
        if !isempty(missing_cols)
            @warn "Missing columns $missing_cols in $csv_filepath. Using INVISCID model instead."
            return nothing, INVISCID
        end
        
        # Create aero_data tuple
        aero_data = (df.alpha, df.cl, df.cd, df.cm)
        return aero_data, POLAR_VECTORS
        
    catch e
        @warn "Error reading polar file $csv_filepath: $e. Using INVISCID model instead."
        return nothing, INVISCID
    end
end


function instantiate_from_yaml(yaml_file_path, polars_dir; n_panels=20, spanwise_distribution=LINEAR)
    """
    Create a wing geometry from YAML file with polar data loaded from CSV files.
    
    # Arguments
    - `yaml_file_path`: Path to YAML file containing wing section definitions
    - `polars_dir`: Path to directory containing polar CSV files (named <airfoil_id>.csv)
    - `n_panels`: Number of panels for the wing (default: 20)
    - `spanwise_distribution`: Panel distribution type (default: LINEAR)
    
    # Returns
    - `BodyAerodynamics`: Initialized aerodynamics object containing the wing
    
    # Example
    ```julia
    yaml_path = "pyramid.yaml"
    polars_dir = "/path/to/polars"
    body_aero = instantiate_from_yaml(yaml_path, polars_dir; n_panels=30)
    ```
    """
    # Create wing with specified panel distribution
    wing = Wing(n_panels, spanwise_distribution=spanwise_distribution)

    # Read sections from YAML file
    sections = read_yaml(yaml_file_path)

    # Add each section to the wing
    for section in sections
        # Load polar data for this airfoil
        aero_data, aero_model = load_polar_data(section.airfoil_id, polars_dir)
        
        println("Section airfoil_id $(section.airfoil_id): Using $(aero_model) model")
        
        add_section!(wing, 
            section.le_coord,     # Leading edge coordinate
            section.te_coord,     # Trailing edge coordinate
            aero_model,          # Aerodynamic model (POLAR_VECTORS or INVISCID)
            aero_data)           # Aerodynamic data tuple or nothing
    end
    
    # Initialize and return aerodynamics object
    body_aero = BodyAerodynamics([wing])
    return body_aero
end

# getting the project_dir
project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root
yaml_file_path = joinpath(project_dir, "data", "pyramid_model", "pyramid.yaml")
polars_dir = joinpath(project_dir, "data", "pyramid_model", "polars")

sections = read_yaml(yaml_file_path)
println("Read sections from YAML file: $(length(sections)) sections")
println("Sections: $(sections)")

body_aero = instantiate_from_yaml(yaml_file_path, polars_dir; n_panels=20)
println("Created wing geometry with $(length(body_aero.panels)) panels")

# Plotting Geometrys
PLOT = true
USE_TEX = false
PLOT && plot_geometry(
    body_aero,
    "";
    data_type=".svg",
    save_path="",
    is_save=false,
    is_show=true,
    view_elevation=15,
    view_azimuth=-120,
    use_tex=USE_TEX
)
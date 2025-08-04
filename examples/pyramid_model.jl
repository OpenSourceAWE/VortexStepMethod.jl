using LinearAlgebra
using VortexStepMethod
using ControlPlots

# settings
project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root
yaml_geometry_path = joinpath(project_dir, "data", "pyramid_model", "geometry.yaml")

# Create wing geometry from YAML (following RamAirWing pattern)
wing = YamlWing(yaml_geometry_path; n_panels=20, spanwise_distribution=LINEAR)
body_aero = BodyAerodynamics([wing])

# Plotting Geometry
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
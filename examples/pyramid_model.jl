using Pkg
if Base.active_project() \!= joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using LinearAlgebra
using VortexStepMethod

# Load VSM vsm_settings from YAML configuration file
vsm_settings = VSMSettings("pyramid_model/vsm_settings.yaml")

# Create wing, body_aero, and solver objects using vsm_settings
wing = Wing(vsm_settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
solver = Solver(body_aero, vsm_settings)

# Set flight conditions from settings
set_va!(body_aero, vsm_settings)

# Extract values for plotting (optional - for reference)
wind_speed = vsm_settings.condition.wind_speed
angle_of_attack_deg = vsm_settings.condition.alpha
sideslip_deg = vsm_settings.condition.beta
yaw_rate = vsm_settings.condition.yaw_rate

# Run the solver
results = VortexStepMethod.solve(solver, body_aero; log=true)

# Using plotting modules, to create more comprehensive plots
PLOT = true
SAVE_ALL = false
USE_TEX = false
OUTPUT_DIR = joinpath(dirname(@__DIR__), "output")

# Plotting polars
PLOT && plot_polars(
    [solver],
    [body_aero],
    ["VSM Pyramid Model"],
    angle_range=range(-5, 25, length=30),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="$(wing.n_panels)_panels_$(wing.spanwise_distribution)_pyramid_model",
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    use_tex=USE_TEX
)

# Plotting geometry
PLOT && plot_geometry(
    body_aero,
    "Pyramid model geometry";
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    view_elevation=15,
    view_azimuth=-120,
    use_tex=USE_TEX
)

# Plotting spanwise distributions
body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

PLOT && plot_distribution(
    [body_y_coordinates],
    [results],
    ["VSM"];
    title="pyramid_spanwise_distributions_alpha_$(round(angle_of_attack_deg, digits=1))_delta_$(round(sideslip_deg, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(wind_speed, digits=1))",
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    use_tex=USE_TEX
)

nothing
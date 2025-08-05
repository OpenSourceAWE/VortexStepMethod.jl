using LinearAlgebra
using VortexStepMethod
using ControlPlots

project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root

# Load VSM settings from YAML configuration file
settings = vs("TUDELFT_V3_KITE/vsm_settings.yaml")

# Extract flight conditions from settings
wind_speed = settings.condition.wind_speed
angle_of_attack_deg = settings.condition.alpha
sideslip_deg = settings.condition.beta
yaw_rate = settings.condition.yaw_rate

literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","Python_VSM_Rey_5e5_n36_CFD_PCHIP_polars_alpha_sweep.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","WindTunnel_Re_5e5_Poland2025_alpha_sweep_beta_0.csv"),
    ]
labels= [
    "Julia VSM 2D CFD PCHIP",
    "CFD RANS Re=5e5", 
    "CFD RANS Re=10e5 (With Struts)", 
    "Python VSM 2D CFD PCHIP Re=5e5",
    "Wind Tunnel Re=5e5 (With Struts)"
    ]

# Load VSM settings from YAML configuration file
settings = vs("TUDELFT_V3_KITE/vsm_settings.yaml")

# Create wing, body_aero, and solver objects using settings
wing = YamlWing(settings.wings[1].geometry_file; 
    n_panels=settings.wings[1].n_panels, 
    n_groups=settings.wings[1].n_groups,
    spanwise_distribution=settings.wings[1].spanwise_panel_distribution
)
body_aero = BodyAerodynamics([wing])
solver = Solver(body_aero;
    aerodynamic_model_type=settings.solver_settings.aerodynamic_model_type,
    density=settings.solver_settings.density,
    max_iterations=settings.solver_settings.max_iterations,
    rtol=settings.solver_settings.rtol,
    relaxation_factor=settings.solver_settings.relaxation_factor,
    core_radius_fraction=settings.solver_settings.core_radius_fraction,
)
# Set flight conditions
α = deg2rad(settings.condition.alpha)
β = deg2rad(settings.condition.beta)
va = settings.condition.wind_speed * [
    cos(α)*cos(β),  # X_b (forward)
    sin(β),         # Y_b (right)
    sin(α)*cos(β)   # Z_b (down)
]
set_va!(body_aero, va)
# Run the solver
results = VortexStepMethod.solve(solver, body_aero; log=true)



# Using plotting modules, to create more comprehensive plots
PLOT = true
USE_TEX = false

# Plotting polars
PLOT && plot_polars(
    [solver],
    [body_aero],
    labels,
    literature_path_list=literature_paths,
    angle_range=range(-5, 25, length=30),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="$(wing.n_panels)_panels_$(wing.spanwise_distribution)_from_yaml_settings",
    data_type=".pdf",
    is_save=false,
    is_show=true,
    use_tex=USE_TEX
)


# Plotting geometry
results = VortexStepMethod.solve(solver, body_aero; log=true)
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


# Plotting spanwise distributions
body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

PLOT && plot_distribution(
    [body_y_coordinates],
    [results],
    ["VSM"];
    title="CAD_spanwise_distributions_alpha_$(round(angle_of_attack_deg, digits=1))_delta_$(round(sideslip_deg, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(wind_speed, digits=1))",
    data_type=".pdf",
    is_save=false,
    is_show=true,
    use_tex=USE_TEX
)


nothing
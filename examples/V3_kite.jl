using LinearAlgebra
using VortexStepMethod
using ControlPlots

# --- User-specified parameters ---
wind_speed = 3.05        # [m/s]
angle_of_attack_deg = 5.0   # [deg]
sideslip_deg = 0.0
yaw_rate = 0.0             # [rad/s] (not used in static analysis)

project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root
yaml_geometry_path = joinpath(project_dir, "data", "TUDELFT_V3_KITE", "geometry_CAD_CFD_polars_pchip_fitted.yaml")
literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","Python_VSM_Rey_5e5_n36_CFD_PCHIP_polars_alpha_sweep.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","WindTunnel_Re_5e5_Poland2025_alpha_sweep_beta_0.csv"),
    ]
labels= [
    "Julia VSM 2D CFD PCHIP",
    "CFD RANS Re=5e5", 
    "CFD RANS Re=10e5", 
    "Python VSM 2D CFD PCHIP Re=5e5",
    "Wind Tunnel Re=5e5"
    ]

# Create wing, body_aero, and solver objects
wing = YamlWing(yaml_geometry_path; n_panels=36, spanwise_distribution=LINEAR)
body_aero = BodyAerodynamics([wing])
solver = Solver(body_aero;
    aerodynamic_model_type=VSM,
    is_with_artificial_damping=false,
    solver_type=LOOP,
    density=1.225,
    max_iterations=5000,
    rtol = 1e-6,
    relaxation_factor=0.01,
    core_radius_fraction=1e-20,
)

# Using plotting modules, to create more comprehensive plots
PLOT = true
USE_TEX = false

# Plotting polars
PLOT && plot_polars(
    [solver],
    [body_aero],
    labels,
    literature_path_list=literature_paths,
    angle_range=range(1, 20, length=20),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="$(wing.n_panels)_distribution_$(wing.spanwise_distribution)",
    data_type=".pdf",
    is_save=false,
    is_show=true,
    use_tex=USE_TEX
)


# Plotting geometry
#--- Set apparent wind vector in body axes ---
α = deg2rad(angle_of_attack_deg)
β = deg2rad(sideslip_deg)
va = wind_speed * [
    cos(α)*cos(β),  # X_b (forward)
    sin(β),         # Y_b (right)
    sin(α)*cos(β)   # Z_b (down)
]
set_va!(body_aero, va)

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
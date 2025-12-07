using LinearAlgebra
using VortexStepMethod
using GLMakie

PLOT = true
USE_TEX = false
DEFORM = false

project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root#
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
settings = VSMSettings("TUDELFT_V3_KITE/vsm_settings.yaml")

# Create wing, body_aero, and solver objects using settings
wing = Wing(settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)

if DEFORM
    VortexStepMethod.unrefined_deform!(
        wing,
        deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
        deg2rad.(range(0, 0, length=wing.n_unrefined_sections));
        smooth=true
    )
    VortexStepMethod.reinit!(body_aero; init_aero=false)
end

solver = Solver(body_aero, settings)

# Set flight conditions from settings
set_va!(body_aero, settings)

# Extract values for plotting (optional - for reference)
wind_speed = settings.condition.wind_speed
angle_of_attack_deg = settings.condition.alpha
sideslip_deg = settings.condition.beta
yaw_rate = settings.condition.yaw_rate

# Using plotting modules, to create more comprehensive plots
PLOT = true
USE_TEX = false

# Solve and plot combined analysis
results = VortexStepMethod.solve(solver, body_aero; log=true)
PLOT && plot_combined_analysis(
    solver,
    body_aero,
    results;
    solver_label="VSM",
    literature_path_list=literature_paths,
    angle_range=range(-5, 25, length=30),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="TU Delft V3 Kite",
    is_show=true,
    use_tex=USE_TEX
)


nothing

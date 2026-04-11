using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
using VortexStepMethod
using GLMakie

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
USE_TEX = false

# Plotting combined analysis
PLOT && plot_combined_analysis(
    solver,
    body_aero,
    results;
    solver_label="VSM",
    angle_range=range(-5, 25, length=30),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="Pyramid Model",
    is_show=true,
    use_tex=USE_TEX
)

nothing
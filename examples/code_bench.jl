using Revise
using LinearAlgebra
using ControlPlots
using VortexStepMethod
using StaticArrays
using BenchmarkTools

# Step 1: Define wing parameters
n_panels = 20          # Number of panels
span = 20.0            # Wing span [m]
chord = 1.0            # Chord length [m]
v_a = 20.0            # Magnitude of inflow velocity [m/s]
density = 1.225        # Air density [kg/m³]
alpha_deg = 30.0       # Angle of attack [degrees]
alpha = deg2rad(alpha_deg)

# Step 2: Create wing geometry with linear panel distribution
wing = Wing(n_panels, spanwise_panel_distribution="linear")

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],   # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    "inviscid")
add_section!(wing, 
    [0.0, -span/2, 0.0],  # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    "inviscid")

# Step 3: Initialize aerodynamics
wa = WingAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
set_va!(wa, (vel_app, 0.0))  # Second parameter is yaw rate

# Step 4: Initialize solvers for both LLT and VSM methods
vsm_solver = Solver(aerodynamic_model_type="VSM")

# Benchmark setup
velocity_induced = @MVector zeros(3)  # StaticArraysCore.MVector{3, Float64}
tempvel = @MVector zeros(3)          # StaticArraysCore.MVector{3, Float64}
panel = wa.panels[1]
ep = @MVector [0.25, 9.5, 0.0]      # StaticArraysCore.MVector{3, Float64}
evaluation_point_on_bound = true     # Bool
va_norm = 20.0                       # Float64
va_unit = @MVector [0.8660254037844387, 0.0, 0.4999999999999999]  # StaticArraysCore.MVector{3, Float64}
core_radius_fraction = 1.0e-20       # Float64
gamma = 1.0                          # Float64

# Create work vectors tuple of MVector{3, Float64}
work_vectors = ntuple(_ -> @MVector(zeros(3)), 10)  # NTuple{10, StaticArraysCore.MVector{3, Float64}}

@btime VortexStepMethod.calculate_velocity_induced_single_ring_semiinfinite!(
    $velocity_induced,
    $tempvel,
    $panel.filaments,
    $ep,
    $evaluation_point_on_bound,
    $va_norm,
    $va_unit,
    $gamma,
    $core_radius_fraction,
    $work_vectors
)

@btime VortexStepMethod.calculate_velocity_induced_single_ring_semiinfinite!(
    velocity_induced,
    tempvel,
    panel.filaments,
    ep,
    evaluation_point_on_bound,
    va_norm,
    va_unit,
    gamma,
    core_radius_fraction,
    work_vectors
)

U_2D = @MVector zeros(3)  # StaticArraysCore.MVector{3, Float64}

@btime VortexStepMethod.calculate_velocity_induced_bound_2D!(
    $U_2D,
    $panel, 
    $ep,
    $work_vectors
)

model = "VSM"
va_norm_array = zeros(wa.n_panels)
va_unit_array = zeros(wa.n_panels, 3)
@btime VortexStepMethod.calculate_AIC_matrices!(
    $wa, model,
    $core_radius_fraction,
    $va_norm_array, 
    $va_unit_array)

nothing
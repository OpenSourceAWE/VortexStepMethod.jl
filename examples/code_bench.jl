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
wing = Wing(n_panels, spanwise_panel_distribution=:linear)

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],   # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    :inviscid)
add_section!(wing, 
    [0.0, -span/2, 0.0],  # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    :inviscid)

# Step 3: Initialize aerodynamics
body_aero = BodyAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
set_va!(body_aero, (vel_app, 0.0))  # Second parameter is yaw rate

# Step 4: Initialize solvers for both LLT and VSM methods
vsm_solver = Solver(aerodynamic_model_type=:VSM)

# Benchmark setup
velocity_induced = @MVector zeros(3)  # StaticArraysCore.MVector{3, Float64}
tempvel = @MVector zeros(3)          # StaticArraysCore.MVector{3, Float64}
panel = body_aero.panels[1]
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

model = :VSM
n_panels = length(body_aero.panels)
va_norm_array = zeros(n_panels)
va_unit_array = zeros(n_panels, 3)
@btime VortexStepMethod.calculate_AIC_matrices!(
    $body_aero, model,
    $core_radius_fraction,
    $va_norm_array, 
    $va_unit_array)

n_panels = length(body_aero.panels)
gamma_new = zeros(n_panels)
va_array = zeros(n_panels, 3)
chord_array = zeros(n_panels)
x_airf_array = zeros(n_panels, 3)
y_airf_array = zeros(n_panels, 3)
z_airf_array = zeros(n_panels, 3)
relaxation_factor = 0.5

# Fill arrays with data from body_aero
for (i, panel) in enumerate(body_aero.panels)
    va_array[i, :] .= panel.va
    chord_array[i] = panel.chord
    x_airf_array[i, :] .= panel.x_airf
    y_airf_array[i, :] .= panel.y_airf
    z_airf_array[i, :] .= panel.z_airf
end

# Benchmark gamma_loop
@btime VortexStepMethod.gamma_loop(
    $vsm_solver,
    $body_aero,
    $gamma_new,
    $va_array,
    $chord_array,
    $x_airf_array,
    $y_airf_array,
    $z_airf_array,
    $body_aero.panels,
    $relaxation_factor;
    log = false
)

nothing
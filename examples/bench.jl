using LinearAlgebra
using ControlPlots
using VortexStepMethod
using Interpolations
using BenchmarkTools

plot = true

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

# Test data generation
n_angles = 5
n_flaps = 5

# Define angle ranges
alphas = range(-deg2rad(10), deg2rad(10), n_angles)  # AoA range
d_trailing_edge_angles = range(-deg2rad(30), deg2rad(30), n_flaps)  # Flap deflection range

# Initialize coefficient matrices
cl_matrix = zeros(n_angles, n_flaps)
cd_matrix = zeros(n_angles, n_flaps)
cm_matrix = zeros(n_angles, n_flaps)

# Fill matrices with realistic aerodynamic data
for (i, α) in enumerate(alphas)
    for (j, δ) in enumerate(d_trailing_edge_angles)
        cl_matrix[i,j] = 2π * α + π/2 * δ
        cd_matrix[i,j] = 0.01 + 0.05 * α^2 + 0.03 * δ^2
        cm_matrix[i,j] = -0.1 * α - 0.2 * δ
    end
end

cl_interp = extrapolate(scale(interpolate(cl_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
cd_interp = extrapolate(scale(interpolate(cd_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
cm_interp = extrapolate(scale(interpolate(cm_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)

add_section!(wing, 
    [0.0, span/2, 0.0],   # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    (:interpolations, (cl_interp, cd_interp, cm_interp)))
add_section!(wing, 
    [0.0, -span/2, 0.0],  # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    (:interpolations, (cl_interp, cd_interp, cm_interp)))

# Step 3: Initialize aerodynamics
wa = BodyAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
set_va!(wa, (vel_app, 0.0))  # Second parameter is yaw rate

# Step 4: Initialize solvers for both LLT and VSM methods
llt_solver = Solver(aerodynamic_model_type=:LLT)
vsm_solver = Solver(aerodynamic_model_type=:VSM)

# Step 5: Solve using both methods
# @btime results_llt = solve($vsm_solver, $wa; log=false)
# @btime results_vsm = solve($vsm_solver, $wa; log=false)
@time results_vsm = solve(vsm_solver, wa; log=false)
# time Python: 32.0 ms Ryzen 7950x
# time Julia:   0.6 ms Ryzen 7950x
#               0.8 ms laptop, performance mode, battery 

nothing
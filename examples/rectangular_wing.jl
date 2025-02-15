using Revise
using LinearAlgebra
using Plots
using VortexStepMethod

# Step 1: Define wing parameters
n_panels = 20            # Number of panels
span = 20.0             # Wing span [m]
chord = 1.0             # Chord length [m]
Umag = 20.0            # Magnitude of inflow velocity [m/s]
density = 1.225        # Air density [kg/m³]
alpha_deg = 30.0       # Angle of attack [degrees]
alpha = deg2rad(alpha_deg)

# Step 2: Create wing geometry with linear panel distribution
wing = Wing(n_panels, spanwise_panel_distribution="linear")

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],   # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    ["inviscid"])
add_section!(wing, 
    [0.0, -span/2, 0.0],  # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    ["inviscid"])

# Step 3: Initialize aerodynamics
wa = WingAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* Umag
set_va!(wa, (vel_app, 0.0))  # Second parameter is yaw rate

# Step 4: Initialize solvers for both LLT and VSM methods
llt_solver = Solver(aerodynamic_model_type="LLT")
vsm_solver = Solver(aerodynamic_model_type="VSM")

# Step 5: Solve using both methods
@time results_llt = solve(llt_solver, wa)
@time results_llt = solve(llt_solver, wa)
@time results_vsm = solve(vsm_solver, wa)
@time results_vsm = solve(vsm_solver, wa)

# Print results comparison
println("\nLifting Line Theory Results:")
println("CL = $(round(results_llt["cl"], digits=4))")
println("CD = $(round(results_llt["cd"], digits=4))")
println("\nVortex Step Method Results:")
println("CL = $(round(results_vsm["cl"], digits=4))")
println("CD = $(round(results_vsm["cd"], digits=4))")
println("Projected area = $(round(results_vsm["projected_area"], digits=4)) m²")

# Step 6: Plot geometry
plot_geometry(
      wa,
      "rectangular_wing_geometry";
      data_type=".pdf",
      save_path=".",
      is_save=false,
      is_show=true,
)

# # Step 7: Plot spanwise distributions
# y_coordinates = [panel.aerodynamic_center[2] for panel in wa.panels]

# plot_distribution(
#     [y_coordinates, y_coordinates],
#     [results_vsm, results_llt],
#     ["VSM", "LLT"],
#     title="Spanwise Distributions"
# )

# # Step 8: Plot polar curves
# angle_range = range(0, 20, 20)
# plot_polars(
#     [llt_solver, vsm_solver],
#     [wa, wa],
#     ["LLT", "VSM"],
#     angle_range=angle_range,
#     angle_type="angle_of_attack",
#     Umag=Umag,
#     title="Rectangular Wing Polars"
# )

# Save plots if needed
# savefig("geometry.pdf")
# savefig("distributions.pdf")
# savefig("polars.pdf")
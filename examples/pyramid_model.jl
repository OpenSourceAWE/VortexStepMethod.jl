using LinearAlgebra
using VortexStepMethod
using ControlPlots

# --- User-specified parameters ---
wind_speed = 10.0           # [m/s]
angle_of_attack_deg = 5.0   # [deg]
sideslip_deg = 0.0
yaw_rate = 0.0             # [rad/s] (not used in static analysis)
air_density = 1.225         # [kg/m^3] Air density at sea level
project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root
yaml_geometry_path = joinpath(project_dir, "data", "pyramid_model", "geometry.yaml")

println("=== Static Aerodynamic Analysis (VSM) ===")
println("Wind speed: $wind_speed m/s")
println("Angle of attack: $angle_of_attack_deg deg")
println("Sideslip: $sideslip_deg deg")

# Create wing geometry from YAML (following RamAirWing pattern)
wing = YamlWing(yaml_geometry_path; n_panels=20, spanwise_distribution=LINEAR)
body_aero = BodyAerodynamics([wing])

# --- Set apparent wind vector in body axes ---
α = deg2rad(angle_of_attack_deg)
β = deg2rad(sideslip_deg)
va = wind_speed * [
    cos(α)*cos(β),  # X_b (forward)
    sin(β),         # Y_b (right)
    sin(α)*cos(β)   # Z_b (down)
]
set_va!(body_aero, va)

# # --- Create and run VSM solver ---
# solver = Solver(body_aero;
#     aerodynamic_model_type=VSM,
#     is_with_artificial_damping=false,
#     rtol=1e-5,
#     solver_type=NONLIN
# )

# # One can run the solver seperately, given back a Dict{String, Any}
# solver = Solver(body_aero; density=air_density)
# solve!(solver, body_aero)
# forces = solver.sol.force
# println("\nType of solver: ", typeof(solver), " Type of forces: ", typeof(forces))

# # Better to use the VortexStepMethod.solve method, returns a dictionary with the results
# results = VortexStepMethod.solve(solver, body_aero; log=true)
# println("Type of results: ", typeof(results))


# # --- Run over an alpha sweep ---
# # using Plots  # Add this at the top if you want to plot
# alphas = 1:10
# f_xs = Float64[]   # Fx over alpha sweep
# f_ys = Float64[]   # Fy over alpha sweep
# f_zs = Float64[]   # Fz over alpha sweep
# m_xs = Float64[]   # Mx over alpha sweep
# m_ys = Float64[]   # My over alpha sweep
# m_zs = Float64[]   # Mz over alpha sweep

# for alpha in alphas
#     local_α = deg2rad(alpha)
#     local_β = deg2rad(sideslip_deg)
#     local_va = wind_speed * [
#         cos(local_α)*cos(local_β),
#         sin(local_β),
#         sin(local_α)*cos(local_β)
#     ]
#     set_va!(body_aero, local_va)

#     solve!(solver, body_aero)
#     results_i = VortexStepMethod.solve(solver, body_aero; log=true)

#     push!(f_xs, results_i["Fx"])
#     push!(f_ys, results_i["Fy"])
#     push!(f_zs, results_i["Fz"])
#     push!(m_xs, results_i["Mx"])
#     push!(m_ys, results_i["My"])
#     push!(m_zs, results_i["Mz"])

# end


# # --- Plot results: true multi-row plotx style as in simple_lin_model.jl ---
# p = plotx(
#     alphas,
#     f_xs,
#     f_ys,
#     f_zs,
#     m_xs,
#     m_ys,
#     m_zs;
#     ylabels=["Fx [N]", "Fy [N]", "Fz [N]", "Mx [Nm]", "My [Nm]", "Mz [Nm]"],
#     xlabel="alpha [deg]",
#     fig="Aerodynamic Forces and Moments"
# )
# display(p)


# Using plotting modules, to create more comprehensive plots
PLOT = true
USE_TEX = false


# Plotting geometry
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

# # Plotting polars
# PLOT && plot_polars(
#     [solver],
#     [body_aero],
#     ["VSM from Ram Air Kite OBJ and DAT file"],
#     angle_range=range(0, 20, length=20),
#     angle_type="angle_of_attack",
#     angle_of_attack=5,
#     side_slip=0,
#     v_a=10,
#     title="pyramid $(wing.n_panels)_distribution_$(wing.spanwise_distribution)",
#     data_type=".pdf",
#     is_save=false,
#     is_show=true,
#     use_tex=USE_TEX
# )

# # Plotting spanwise distributions
# body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

# PLOT && plot_distribution(
#     [body_y_coordinates],
#     [results],
#     ["VSM"];
#     title="CAD_spanwise_distributions_alpha_$(round(angle_of_attack_deg, digits=1))_delta_$(round(sideslip_deg, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(wind_speed, digits=1))",
#     data_type=".pdf",
#     is_save=false,
#     is_show=true,
#     use_tex=USE_TEX
# )


nothing
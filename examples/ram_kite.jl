using Revise
using VortexStepMethod
using CSV
using DataFrames
using LinearAlgebra
using Plots

# Create wing geometry
wing = KiteWing("data/HL5_ram_air_kite_body.obj")

alphas = -10:1:25  # Range of angles from -10 to 25 degrees
polars = zeros(length(alphas), 4)  # Matrix for [alpha, CD, CL, CM]
for (i, alpha) in enumerate(alphas)
    # Simplified aerodynamic coefficients
    cd = 0.015 + 0.015 * abs(alpha/10)^1.5  # Drag increases with angle
    cl = 0.1 * alpha + 0.01 * alpha^2 * exp(-alpha/20)  # Lift with stall behavior
    cm = -0.02 * alpha  # Linear pitching moment
    
    polars[i, :] .= [alpha, cd, cl, cm]
end

for gamma in range(wing.gamma_tip - wing.gamma_tip/10, -wing.gamma_tip + wing.gamma_tip/10, 20)
    add_section!(wing, gamma, ("polar_data", polars))
end
wing_aero = WingAerodynamics([wing])

# Create solvers
VSM = Solver(
    aerodynamic_model_type="VSM",
    is_with_artificial_damping=false
)
VSM_with_stall_correction = Solver(
    aerodynamic_model_type="VSM",
    is_with_artificial_damping=true
)

# Setting velocity conditions
Umag = 15.0
aoa = 17.0
side_slip = 0.0
yaw_rate = 0.0
aoa_rad = deg2rad(aoa)
vel_app = [
    cos(aoa_rad) * cos(side_slip),
    sin(side_slip),
    sin(aoa_rad)
] * Umag
wing_aero.va = vel_app

# # Plotting geometry
# plot_geometry(
#     wing_aero,
#     "";
#     data_type=".svg",
#     save_path="",
#     is_save=false,
#     is_show=true,
#     view_elevation=15,
#     view_azimuth=-120
# )

# Solving and plotting distributions
results = solve(VSM, wing_aero)
@time results_with_stall = solve(VSM_with_stall_correction, wing_aero)
@time results_with_stall = solve(VSM_with_stall_correction, wing_aero)

# CAD_y_coordinates = [panel.aerodynamic_center[2] for panel in wing_aero.panels]

# plot_distribution(
#     [CAD_y_coordinates, CAD_y_coordinates],
#     [results, results_with_stall],
#     ["VSM", "VSM with stall correction"];
#     title="CAD_spanwise_distributions_alpha_$(round(aoa, digits=1))_beta_$(round(side_slip, digits=1))_yaw_$(round(yaw_rate, digits=1))_Umag_$(round(Umag, digits=1))",
#     data_type=".pdf",
#     save_path=joinpath(save_folder, "spanwise_distributions"),
#     is_save=false,
#     is_show=true
# )

# # Plotting polar
# save_path = joinpath(root_dir, "results", "TUD_V3_LEI_KITE")
# path_cfd_lebesque = joinpath(
#     root_dir,
#     "data",
#     "TUDELFT_V3_LEI_KITE",
#     "literature_results",
#     "V3_CL_CD_RANS_Lebesque_2024_Rey_300e4.csv"
# )

# plot_polars(
#     [VSM, VSM_with_stall_correction],
#     [wing_aero, wing_aero],
#     [
#         "VSM CAD 19ribs",
#         "VSM CAD 19ribs , with stall correction",
#         "CFD_Lebesque Rey 30e5"
#     ];
#     literature_path_list=[path_cfd_lebesque],
#     angle_range=range(0, 25, length=25),
#     angle_type="angle_of_attack",
#     angle_of_attack=0,
#     side_slip=0,
#     yaw_rate=0,
#     Umag=10,
#     title="tutorial_testing_stall_model_n_panels_$(n_panels)_distribution_$(spanwise_panel_distribution)",
#     data_type=".pdf",
#     save_path=joinpath(save_folder, "polars"),
#     is_save=true,
#     is_show=true
# )
using VortexStepMethod
using CSV
using DataFrames
using LinearAlgebra

# Create wing geometry
wing = KiteWing("data/HL5_ram_air_kite_body.obj", "data/HL5_ram_air_kite_foil.dat")

# for gamma in range(wing.gamma_tip - wing.gamma_tip/10, -wing.gamma_tip + wing.gamma_tip/10, 20)
#     add_section!(wing, gamma, ("dat_file", "data/centre_line_with_profile.dat"))
# end
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
v_a = 15.0
aoa = 15.0
side_slip = 0.0
yaw_rate = 0.0
aoa_rad = deg2rad(aoa)
vel_app = [
    cos(aoa_rad) * cos(side_slip),
    sin(side_slip),
    sin(aoa_rad)
] * v_a
wing_aero.va = vel_app

# Plotting geometry
plot_geometry(
    wing_aero,
    "";
    data_type=".svg",
    save_path="",
    is_save=false,
    is_show=true,
    view_elevation=15,
    view_azimuth=-120
)

# Solving and plotting distributions
@time results = solve(VSM, wing_aero)
@time results_with_stall = solve(VSM_with_stall_correction, wing_aero)

CAD_y_coordinates = [panel.aerodynamic_center[2] for panel in wing_aero.panels]

plot_distribution(
    [CAD_y_coordinates, CAD_y_coordinates],
    [results, results_with_stall],
    ["VSM", "VSM with stall correction"];
    title="CAD_spanwise_distributions_alpha_$(round(aoa, digits=1))_beta_$(round(side_slip, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(v_a, digits=1))",
    data_type=".pdf",
    is_save=false,
    is_show=true
)

plot_polars(
    [VSM, VSM_with_stall_correction],
    [wing_aero, wing_aero],
    [
        "VSM CAD 19ribs",
        "VSM CAD 19ribs , with stall correction",
    ];
    angle_range=range(0, 25, length=25),
    angle_type="angle_of_attack",
    angle_of_attack=0,
    side_slip=0,
    yaw_rate=0,
    v_a=10,
    title="tutorial_testing_stall_model_n_panels_$(wing.n_panels)_distribution_$(wing.spanwise_panel_distribution)",
    data_type=".pdf",
    is_save=false,
    is_show=true
)
nothing
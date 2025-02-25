using VortexStepMethod
using LinearAlgebra
using Pkg

if !("CSV" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
end
using CSV
using DataFrames

plot = false

# Create wing geometry
wing = KiteWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat")
body_aero = BodyAerodynamics([wing])

# Create solvers
VSM = Solver(
    aerodynamic_model_type=:VSM,
    is_with_artificial_damping=false
)
VSM_with_stall_correction = Solver(
    aerodynamic_model_type=:VSM,
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
body_aero.va = vel_app

# Plotting geometry
plot && plot_geometry(
    body_aero,
    "";
    data_type=".svg",
    save_path="",
    is_save=false,
    is_show=true,
    view_elevation=15,
    view_azimuth=-120
)

# Solving and plotting distributions
results = solve(VSM, body_aero)
@time results = solve(VSM, body_aero)

CAD_y_coordinates = [panel.aerodynamic_center[2] for panel in body_aero.panels]

plot && plot_distribution(
    [CAD_y_coordinates],
    [results],
    ["VSM"];
    title="CAD_spanwise_distributions_alpha_$(round(aoa, digits=1))_beta_$(round(side_slip, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(v_a, digits=1))",
    data_type=".pdf",
    is_save=false,
    is_show=true
)

plot && plot_polars(
    [VSM],
    [body_aero],
    [
        "VSM from Ram Air Kite OBJ and DAT file",
    ];
    angle_range=range(0, 20, length=20),
    angle_type="angle_of_attack",
    angle_of_attack=0,
    side_slip=0,
    v_a=10,
    title="ram_kite_panels_$(wing.n_panels)_distribution_$(wing.spanwise_panel_distribution)",
    data_type=".pdf",
    is_save=false,
    is_show=true
)
nothing
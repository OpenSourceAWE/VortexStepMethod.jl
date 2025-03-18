using ControlPlots
using VortexStepMethod
using LinearAlgebra

PLOT = true
USE_TEX = false
DEFORM = false

# Create wing geometry
wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat")
body_aero = BodyAerodynamics([wing];)

if DEFORM
    # Linear interpolation of alpha from 10° at one tip to 0° at the other
    n_panels = wing.n_panels
    theta_start = deg2rad(10)
    theta_end = deg2rad(0)
    delta_start = deg2rad(10)
    delta_end = deg2rad(0)
    theta_dist = [theta_start - i * (theta_start - theta_end)/(n_panels-1) for i in 0:(n_panels-1)]
    delta_dist = [delta_start - i * (delta_start - delta_end)/(n_panels-1) for i in 0:(n_panels-1)]
    @time VortexStepMethod.deform!(wing, theta_dist, delta_dist)
    @time VortexStepMethod.init!(body_aero)
end

# Create solvers
P = length(body_aero.panels)
vsm_solver = Solver{P}(
    aerodynamic_model_type=VSM,
    is_with_artificial_damping=false
)

# Setting velocity conditions
v_a = 15.0
aoa = 10.0
side_slip = 0.0
yaw_rate = 0.0
aoa_rad = deg2rad(aoa)
vel_app = [
    cos(aoa_rad) * cos(side_slip),
    sin(side_slip),
    sin(aoa_rad)
] * v_a
set_va!(body_aero, vel_app)

# Plotting polar data
PLOT && plot_polar_data(body_aero)

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

# Solving and plotting distributions
results = solve(vsm_solver, body_aero; log=true)
@time results = solve(vsm_solver, body_aero; log=true)

body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

PLOT && plot_distribution(
    [body_y_coordinates],
    [results],
    ["VSM"];
    title="CAD_spanwise_distributions_alpha_$(round(aoa, digits=1))_delta_$(round(side_slip, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(v_a, digits=1))",
    data_type=".pdf",
    is_save=false,
    is_show=true,
    use_tex=USE_TEX
)

PLOT && plot_polars(
    [vsm_solver],
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
    is_show=true,
    use_tex=USE_TEX
)
nothing
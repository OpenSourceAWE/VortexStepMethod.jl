using ControlPlots
using VortexStepMethod
using LinearAlgebra
using Pkg

if !("CSV" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
end
using CSV
using DataFrames

PLOT = true
USE_TEX = false
DEFORM = false

# Create wing geometry
wing = KiteWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat")
body_aero = BodyAerodynamics([wing];)

if DEFORM
    alpha = [deg2rad(10), 0]
    beta = [deg2rad(10), 0]
    @time VortexStepMethod.deform!(wing, alpha, beta; width=1.0)
    @time VortexStepMethod.init!(body_aero)
end

# Create solvers
vsm_solver = Solver(
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
    title="CAD_spanwise_distributions_alpha_$(round(aoa, digits=1))_beta_$(round(side_slip, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(v_a, digits=1))",
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
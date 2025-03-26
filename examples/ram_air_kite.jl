using ControlPlots
using VortexStepMethod
using LinearAlgebra

PLOT = false
USE_TEX = false
DEFORM = false
LINEARIZE = true

# Create wing geometry
wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat")
body_aero = BodyAerodynamics([wing];)
println("First init")
@time VortexStepMethod.init!(body_aero)

if DEFORM
    # Linear interpolation of alpha from 10° at one tip to 0° at the other
    println("Deform")
    @time VortexStepMethod.smooth_deform!(wing, deg2rad.([10,20,10,0]), deg2rad.([-10,0,-10,0]))
    println("Deform init")
    @time VortexStepMethod.init!(body_aero; init_aero=false)
end

# Create solvers
solver = Solver(body_aero;
    aerodynamic_model_type=VSM,
    is_with_artificial_damping=false,
    rtol=1e-5,
    solver_type=NONLIN
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

if LINEARIZE
    println("Linearize")
    jac, res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        wing, 
        [zeros(4); vel_app; zeros(3)]; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        moment_frac=0.1)
    @time jac, res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        wing, 
        [zeros(4); vel_app; zeros(3)]; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        moment_frac=0.1)
end

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
println("Solve")
results = VortexStepMethod.solve(solver, body_aero; log=true)
@time results = solve(solver, body_aero; log=true)

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
    [solver],
    [body_aero],
    [
        "VSM from Ram Air Kite OBJ and DAT file",
    ];
    angle_range=range(0, 20, length=20),
    angle_type="angle_of_attack",
    angle_of_attack=0,
    side_slip=0,
    v_a=10,
    title="ram_kite_panels_$(wing.n_panels)_distribution_$(wing.spanwise_distribution)",
    data_type=".pdf",
    is_save=false,
    is_show=true,
    use_tex=USE_TEX
)
nothing
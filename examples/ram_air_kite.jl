using GLMakie
using VortexStepMethod
using LinearAlgebra

PLOT = true
PRN = true
USE_TEX = false
DEFORM = true
LINEARIZE = true

# Create wing geometry
wing = ObjWing(
    joinpath("data", "ram_air_kite", "ram_air_kite_body.obj"),
    joinpath("data", "ram_air_kite", "ram_air_kite_foil.dat");
    n_unrefined_sections=4,
    prn=PRN
)
body_aero = BodyAerodynamics([wing];)
println("First init")
@time VortexStepMethod.reinit!(body_aero)

if DEFORM
    # Linear interpolation of alpha from 10° at one tip to 0° at the other
    println("Deform")
    @time VortexStepMethod.unrefined_deform!(wing, deg2rad.([10,20,10,0]), deg2rad.([-10,0,-10,0]); smooth=true)
    println("Deform init")
    @time VortexStepMethod.reinit!(body_aero; init_aero=false)
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
        [zeros(4); vel_app; zeros(3)];
        theta_idxs=1:4,
        va_idxs=5:7,
        omega_idxs=8:10,
        moment_frac=0.1)
    @time jac, res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        [zeros(4); vel_app; zeros(3)]; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        moment_frac=0.1)
end

# Solving
println("Solve")
results = VortexStepMethod.solve(solver, body_aero; log=true)
@time results = solve(solver, body_aero; log=true)

body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

if PLOT
    plot_combined_analysis(
        solver,
        body_aero,
        results;
        solver_label="VSM",
        angle_range=range(0, 20, length=20),
        angle_type="angle_of_attack",
        angle_of_attack=aoa,
        side_slip=side_slip,
        v_a=v_a,
        title="Ram Air Kite (α=$(aoa)°, β=$(side_slip)°, v=$(v_a) m/s)",
        view_elevation=15,
        view_azimuth=-120,
        is_show=true,
        use_tex=USE_TEX
    )
end
nothing

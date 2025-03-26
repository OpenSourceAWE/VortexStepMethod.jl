using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!
using LinearAlgebra
using Test
using Logging

if !@isdefined ram_wing
    cp("data/ram_air_kite_body.obj", "/tmp/ram_air_kite_body.obj"; force=true)
    cp("data/ram_air_kite_foil.dat", "/tmp/ram_air_kite_foil.dat"; force=true)
    ram_wing = RamAirWing("/tmp/ram_air_kite_body.obj", "/tmp/ram_air_kite_foil.dat"; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1))
end

@testset "Nonlinear vs Linear" begin
    va = [15.0, 0.0, 0.0]
    theta = zeros(4)
    delta = zeros(4)
    omega = zeros(3)

    body_aero = BodyAerodynamics([ram_wing]; va)
    solver = Solver(body_aero;
        aerodynamic_model_type=VSM,
        is_with_artificial_damping=false,
        solver_type=NONLIN
    )

    jac, res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        wing, 
        [theta; va; omega; delta]; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        delta_idxs=11:14,
        moment_frac=0.1
    )
    results = VortexStepMethod.solve!(solver, body_aero; log=true)

    

end
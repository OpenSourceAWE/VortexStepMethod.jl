using VortexStepMethod
using LinearAlgebra
using Pkg

if !("CSV" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
end
using CSV
using DataFrames
using BenchmarkTools

plot = true

# Create wing geometry
wing = KiteWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat")
body_aero = BodyAerodynamics([wing])

# Create solvers
VSM = Solver(
    aerodynamic_model_type=:VSM,
    is_with_artificial_damping=false
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

# Solving and plotting distributions
results = solve(VSM, body_aero)
@btime results = solve($VSM, $body_aero)

nothing
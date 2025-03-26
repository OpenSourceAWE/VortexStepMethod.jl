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

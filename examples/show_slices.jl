using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using VortexStepMethod
using ObjAdapter

# Change this to the OBJ file you want to slice
OBJ_PATH = joinpath("data", "TUDELFT_V3_KITE", "V3_25.obj")
# OBJ_PATH = joinpath("data", "ram_air_kite", "ram_air_kite_body.obj")

N_SLICES = 10

println("Slicing: $OBJ_PATH")
# Hover a slice to inspect its 2D airfoil (raw points + Kulfan/envelope fits).
fig = plot_slices_3d(OBJ_PATH; n_slices=N_SLICES, rotation=:auto)
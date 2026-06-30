using Pkg
if Base.active_project() \!= joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using MakieControlPlots
using VortexStepMethod
using LinearAlgebra

PLOT = true
PRN = false
USE_TEX = false
DEFORM = true
v_a = 15.0

# Create wing with XFoil polars
println("Creating XFoil wing...")
wing_xfoil = ObjWing(
    joinpath("data", "ram_air_kite", "ram_air_kite_body.obj"),
    joinpath("data", "ram_air_kite", "ram_air_kite_foil.dat");
    n_unrefined_sections=10,
    prn=PRN
)
body_xfoil = BodyAerodynamics([wing_xfoil])
VortexStepMethod.reinit!(body_xfoil)
solver_xfoil = Solver(body_xfoil; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=NONLIN)

# Create wing with NeuralFoil polars
println("Creating NeuralFoil wing...")
wing_nf = ObjWing(
    joinpath("data", "ram_air_kite", "ram_air_kite_body.obj"),
    joinpath("data", "ram_air_kite", "ram_air_kite_foil.dat");
    n_unrefined_sections=10,
    polars_dir=joinpath("data", "ram_air_kite", "polars_neuralfoil"),
    prn=PRN
)
body_nf = BodyAerodynamics([wing_nf])
VortexStepMethod.reinit!(body_nf)
solver_nf = Solver(body_nf; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=NONLIN)

# Compare using plot_polars
if PLOT
    println("Computing and plotting polars...")
    plot_polars(
        [solver_xfoil, solver_nf],
        [body_xfoil, body_nf],
        ["XFoil", "NeuralFoil"];
        angle_range=range(-5, 25, length=31),
        v_a=v_a,
        title="Ram Air Kite: XFoil vs NeuralFoil",
        is_save=false,
        use_tex=USE_TEX
    )
end
nothing
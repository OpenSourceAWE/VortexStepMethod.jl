using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: XFoilSolver, NeuralFoilSolver
using LinearAlgebra

PLOT = true
USE_TEX = false
v_a = 15.0

# Convert-then-load: the .obj mesh is sliced, a representative section is picked,
# and one shared (alpha, delta) POLAR_MATRICES is generated with the chosen solver.
obj_path = joinpath("data", "ram_air_kite", "ram_air_kite_body.obj")
alpha_range = deg2rad.(-8:2:26)
delta_range = deg2rad.(-5:5:5)

function matrix_wing(solver; n_panels=40, n_sections=10)
    yaml = obj_to_matrix_yaml(obj_path, mktempdir();
        n_sections, wind_vel=v_a, alpha_range, delta_range, solver, verbose=false)
    return Wing(yaml; n_panels)
end

println("Creating XFoil wing...")
wing_xfoil = matrix_wing(XFoilSolver())
body_xfoil = BodyAerodynamics([wing_xfoil])
solver_xfoil = Solver(body_xfoil; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=NONLIN)

println("Creating NeuralFoil wing...")
wing_nf = matrix_wing(NeuralFoilSolver())
body_nf = BodyAerodynamics([wing_nf])
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

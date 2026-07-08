using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: XFoilSolver, NeuralFoilSolver, LeastSquaresFit, EnvelopeFit
using LinearAlgebra

PLOT = true
USE_TEX = false
v_a = 15.0
RE = 1e6

# Shared boundary-layer transition settings — both backends use the e^N model.
N_CRIT    = 9.0   # e^N transition criticality (lower = dirtier / earlier)
XTR_UPPER = 0.05  # forced upper-surface transition (chord fraction)
XTR_LOWER = 0.05  # forced lower-surface transition

# XFoil-only settings.
XF_NPAN     = 160  # panel count (XFoil PANE default)
XF_MAX_ITER = 100  # viscous iterations per angle
XF_MACH     = 0.0  # Mach number (0 = incompressible)
XF_SOLVER = XFoilSolver(npan=XF_NPAN, max_iter=XF_MAX_ITER, ncrit=N_CRIT,
                        xtrip=(XTR_UPPER, XTR_LOWER), mach=XF_MACH)

# NeuralFoil-only settings.
NF_MODEL_SIZE = "large"  # network size: xxsmall … xxxlarge (accuracy vs speed)
NF_SOLVER = NeuralFoilSolver(model_size=NF_MODEL_SIZE, n_crit=N_CRIT,
                             xtr_upper=XTR_UPPER, xtr_lower=XTR_LOWER)

# Convert-then-load: the .obj mesh is sliced per section, each section fitted and swept
# over (alpha, delta) with the chosen solver, written as long-format POLAR_MATRICES.
# The result is cached under data/ram_air_kite/<subdir>; delete it to regenerate.
obj_path = joinpath("data", "ram_air_kite", "ram_air_kite_body.obj")
# alpha_range = -8:2:26
alpha_range = -10:5:30
delta_range = -0:1
N_SECTIONS = 50

# Inset the tip stations 10 cm along the leading edge so they avoid the near-zero-chord
# wingtips, which slice to degenerate airfoils that XFoil cannot analyse.
WINGTIP_DISTANCE = 0.1

# FIT_METHOD = LeastSquaresFit()
FIT_METHOD = EnvelopeFit()

# 3D slice diagnostic: mesh + LE/TE curves + section contours and their fits.
# Hover a slice to inspect its 2D airfoil (raw points + the fit).
fig_slices_3d = plot_slices_3d(obj_path; n_slices=N_SECTIONS, fit_method=FIT_METHOD,
                               wingtip_distance=WINGTIP_DISTANCE)
save("ram_air_slices_3d.png", fig_slices_3d)

function matrix_wing(solver, subdir; n_panels=50, n_sections=N_SECTIONS)
    out_dir = joinpath("data", "ram_air_kite", subdir)
    yaml = joinpath(out_dir, "geometry.yaml")
    # if !isfile(yaml)
        obj_to_yaml(obj_path, out_dir; n_sections, Re=RE, alpha_range, delta_range,
            aero_solver=solver, wingtip_distance=WINGTIP_DISTANCE, verbose=false)
    # end
    return Wing(yaml; n_panels)
end

println("Creating XFoil wing...")
wing_xfoil = matrix_wing(XF_SOLVER, "polars_xfoil")
body_xfoil = BodyAerodynamics([wing_xfoil])
solver_xfoil = Solver(body_xfoil; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=NONLIN)

println("Creating NeuralFoil wing...")
wing_nf = matrix_wing(NF_SOLVER, "polars_neuralfoil")
body_nf = BodyAerodynamics([wing_nf])
solver_nf = Solver(body_nf; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=NONLIN)

# Compare using plot_polars
if PLOT
    println("Computing and plotting polars...")
    fig = plot_polars(
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

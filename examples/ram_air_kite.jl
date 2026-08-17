using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using MakieControlPlots
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: XFoilSolver, NeuralFoilSolver, ShrinkWrap
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

# VSM solver stability settings.
RELAXATION         = 0.03   # iteration relaxation factor
ARTIFICIAL_DAMPING = false  # smooth-circulation stabiliser for difficult cases

# Convert-then-load: the .obj mesh is sliced per section, each section fitted and swept
# over (alpha, delta) with the chosen solver, written as long-format POLAR_MATRICES.
# The result is cached under data/ram_air_kite/<subdir>; delete it to regenerate.
obj_path = joinpath("data", "ram_air_kite", "ram_air_kite_body.obj")
alpha_range = -8:2:26
# alpha_range = 0:1
delta_range = 0:1
N_SECTIONS = 10

# At 0.0 XFoil only converges above 16 deg on the near-tip slice this mesh still has.
WINGTIP_DISTANCE = 0.1

# Shrink-wrap each raw slice into a clean closed airfoil: the distance-field wrap
# hugs the cloud at `clearance` (floored at one grid `cell_size`) with a round nose.
WRAP = ShrinkWrap(clearance=0.0)

# 3D slice diagnostic (live preview): mesh + LE/TE curves + section contours, their
# wraps and Kulfan fits. A nonzero `delta` (deg) also overlays the deflected
# re-wrapped airfoil (purple) — the deform_section geometry the solvers consume.
# Hover a slice to inspect its 2D airfoil.
fig_slices_3d = plot_slices_3d(obj_path; n_slices=N_SECTIONS, wrap_method=WRAP,
                               wingtip_distance=WINGTIP_DISTANCE, delta=15.0)
GLMakie.save("ram_air_slices_3d.png", fig_slices_3d)

function matrix_wing(solver, subdir; n_panels=50, n_sections=N_SECTIONS)
    out_dir = joinpath("data", "ram_air_kite", subdir)
    yaml = joinpath(out_dir, "geometry.yaml")
    if !isfile(yaml)
        obj_to_yaml(obj_path, out_dir; n_sections, Re=RE, alpha_range, delta_range,
            aero_solver=solver, wrap_method=WRAP, wingtip_distance=WINGTIP_DISTANCE,
            verbose=false)
    end
    return Wing(yaml; n_panels)
end

println("Creating XFoil wing...")
wing_xfoil = matrix_wing(XF_SOLVER, "polars_xfoil")

# Audit plot: same diagnostic, but read from the generated output directory — the
# written .dat airfoils the polars were actually computed from (delta must be one
# of delta_range; nothing is re-sliced or re-wrapped).
fig_audit = plot_slices_3d(joinpath("data", "ram_air_kite", "polars_xfoil");
                           delta=1.0, obj_path=obj_path)
GLMakie.save("ram_air_slices_audit.png", fig_audit)
body_xfoil = BodyAerodynamics([wing_xfoil])
solver_xfoil = Solver(body_xfoil; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=LOOP,
                      relaxation_factor=RELAXATION, is_with_artificial_damping=ARTIFICIAL_DAMPING)

println("Creating NeuralFoil wing...")
wing_nf = matrix_wing(NF_SOLVER, "polars_neuralfoil")
body_nf = BodyAerodynamics([wing_nf])
solver_nf = Solver(body_nf; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=LOOP,
                   relaxation_factor=RELAXATION, is_with_artificial_damping=ARTIFICIAL_DAMPING)

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

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using MakieControlPlots
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: ShrinkWrap, NeuralFoilSolver, XFoilSolver
using LinearAlgebra

# Configuration
OBJ_PATH = joinpath("data", "TUDELFT_V3_KITE", "V3_25.obj")
POLARS_DIR = joinpath("data", "TUDELFT_V3_KITE", "polars_neuralfoil")
N_SLICES = 50
RE = 1e6

# Boundary-layer transition — shared by both backends (e^N model). A kite canopy is
# effectively tripped near the LE, so transition is forced at 5% chord on both surfaces.
N_CRIT    = 4.0   # e^N transition criticality (lower = dirtier / earlier)
XTR_UPPER = 0.05  # forced upper-surface transition (chord fraction)
XTR_LOWER = 0.05  # forced lower-surface transition

# NeuralFoil settings (the backend used here). DELTA_RANGE=nothing gives a plain
# POLAR_VECTORS sweep; set e.g. `-5:5:5` (deg) for an (alpha, delta) matrix.
NF_MODEL_SIZE = "large"   # network size: xxsmall … xxxlarge (accuracy vs speed)
DELTA_RANGE   = nothing
NF_SOLVER = NeuralFoilSolver(model_size=NF_MODEL_SIZE, n_crit=N_CRIT,
                             xtr_upper=XTR_UPPER, xtr_lower=XTR_LOWER)

# XFoil settings — exposed for completeness. XFoil does not converge on the V3's thin,
# tube-nosed membrane airfoils, so NeuralFoil is used, but the knobs are here to swap in.
XF_NPAN     = 160  # panel count (XFoil PANE default)
XF_MAX_ITER = 100  # viscous iterations per angle
XF_MACH     = 0.0  # Mach number (0 = incompressible)
XF_SOLVER = XFoilSolver(npan=XF_NPAN, max_iter=XF_MAX_ITER, ncrit=N_CRIT,
                        xtrip=(XTR_UPPER, XTR_LOWER), mach=XF_MACH)

# VSM solver stability settings.
RELAXATION         = 0.03   # iteration relaxation factor
ARTIFICIAL_DAMPING = false   # smooth-circulation stabiliser for difficult cases

# V3_25.obj is already in slicer convention (x=chord, y=span, z=up).
ROTATION = I

# Marched leading-edge stations across the span (finer => smoother edge trace).
N_BINS = 100

# Shrink-wrap each raw slice into a clean closed airfoil: the distance-field wrap
# hugs the cloud at `clearance` (floored at one grid `cell_size`); a single-skin
# canopy becomes a thin capsule. `min_concave_radius` sets the fillet bridging the
# concave tube-canopy junction: 0.2 smooths the neck over entirely; the default
# (0.02) traces it, which NeuralFoil also handles fine.
WRAP = ShrinkWrap(clearance=0.0, min_concave_radius=0.2)

# 3D slice diagnostic: mesh + LE/TE curves + section contours and their wraps.
# Hover a slice to inspect its 2D airfoil (raw points + the shrink-wrap).
fig_slices_3d = plot_slices_3d(OBJ_PATH; n_slices=N_SLICES, rotation=ROTATION,
                               n_bins=N_BINS, wrap_method=WRAP)
GLMakie.save("V3_slices_3d.png", fig_slices_3d)

# `obj_to_yaml` generates the whole aero geometry from the OBJ: shape (.dat), polar CSV,
# and Cp/cf tables plus a geometry.yaml. XFoil (`aero_solver=XFoilSolver(...)`) also works
# but doesn't converge on the V3's concave slices (the sharp corner behind the leading
# edge), so NeuralFoil is used here.
println("Generating V3 aero geometry from OBJ (obj_to_yaml)...")
gen_dir = joinpath("data", "TUDELFT_V3_KITE", "generated_neuralfoil")
nf_yaml = obj_to_yaml(OBJ_PATH, gen_dir; n_sections=N_SLICES, Re=RE,
    rotation=ROTATION, wrap_method=WRAP, aero_solver=NF_SOLVER, verbose=true)

# Flight conditions
v_a = 10.0
angle_range = range(-5, 25, length=31)

# Load settings and create wing with CFD polars
println("\nCreating wing with CFD polars...")
settings_cfd = VSMSettings("TUDELFT_V3_KITE/vsm_settings.yaml")
settings_cfd.solver_settings.relaxation_factor = RELAXATION
settings_cfd.solver_settings.artificial_damping = ARTIFICIAL_DAMPING
wing_cfd = Wing(settings_cfd)
refine!(wing_cfd)
body_cfd = BodyAerodynamics([wing_cfd])
VortexStepMethod.reinit!(body_cfd)
solver_cfd = Solver(body_cfd, settings_cfd)

println("Creating wing with NeuralFoil polars...")
wing_nf = Wing(nf_yaml; n_panels=50, spanwise_distribution=LINEAR)
refine!(wing_nf)
body_nf = BodyAerodynamics([wing_nf])
VortexStepMethod.reinit!(body_nf)
solver_nf = Solver(body_nf, settings_cfd)

# Compare CFD-polar and NeuralFoil-polar wings against published references
# (Poland 2025 RANS CFD and wind tunnel). `plot_polars` sweeps each solver over the
# angle range and overlays the literature CSVs; labels cover solvers then references.
lit_dir = joinpath("data", "TUDELFT_V3_KITE", "literature_results")
literature_paths = [
    joinpath(lit_dir, "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(lit_dir, "windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv"),
]
fig = plot_polars(
    [solver_cfd, solver_nf],
    [body_cfd, body_nf],
    ["CFD polars", "NeuralFoil polars", "RANS CFD (Poland 2025)",
     "Wind tunnel (Poland 2025)"];
    literature_path_list=literature_paths,
    angle_range,
    v_a,
    title="TU Delft V3 Kite: CFD vs NeuralFoil (Re=$RE)",
    is_save=false,
)
nothing

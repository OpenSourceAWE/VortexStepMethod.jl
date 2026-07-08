using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: EnvelopeFit, NeuralFoilSolver, XFoilSolver, LeastSquaresFit
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

# Reorient the mesh to the slicer's chord/span/up convention from its extents.
ROTATION = :auto

# Marched leading-edge stations across the span (finer => smoother edge trace).
N_BINS = 100

n_grid = 100
tightness_dist = ones(n_grid)
tightness_dist[1:6] .= 10.0
tightness_dist[90:100] .= 10.0
# tightness_dist[5:35] .= 0.5

FIT_METHOD = EnvelopeFit(min_distance=0.0, tightness=1.0,
                         tightness_dist=tightness_dist)
# FIT_METHOD = LeastSquaresFit()

# 3D slice diagnostic: mesh + LE/TE curves + section contours and their fits.
# Hover a slice to inspect its 2D airfoil (raw points + the envelope fit).
fig_slices_3d = plot_slices_3d(OBJ_PATH; n_slices=N_SLICES, rotation=ROTATION,
                               n_bins=N_BINS, fit_method=FIT_METHOD)
save("V3_slices_3d.png", fig_slices_3d)

# Generate per-slice NeuralFoil polars from the envelope-fitted sections. XFoil is also
# supported (pass `solver=XFoilSolver(...)`) but its panel method does not converge on
# the V3's thin, tube-nosed membrane airfoils, so NeuralFoil is used here.
println("Generating NeuralFoil polars from OBJ mesh...")
generate_section_polars(OBJ_PATH, POLARS_DIR; n_slices=N_SLICES, Re=RE,
    rotation=ROTATION, n_bins=N_BINS, fit_method=FIT_METHOD,
    solver=NF_SOLVER, delta_range=DELTA_RANGE, verbose=true)

# Flight conditions
v_a = 10.0
angle_range = range(-5, 25, length=31)

# Load settings and create wing with CFD polars
println("\nCreating wing with CFD polars...")
settings_cfd = VSMSettings("TUDELFT_V3_KITE/vsm_settings.yaml")
wing_cfd = Wing(settings_cfd)
refine!(wing_cfd)
body_cfd = BodyAerodynamics([wing_cfd])
VortexStepMethod.reinit!(body_cfd)
solver_cfd = Solver(body_cfd, settings_cfd)

# Build a wing whose every airfoil points at a generated POLAR_VECTORS CSV under
# `polars_subdir`. Repoints a copy of the awesIO geometry and loads it as a normal
# uniform POLAR_VECTORS wing — no manual section surgery needed.
using YAML
yaml_path = joinpath("data", "TUDELFT_V3_KITE", "aero_geometry.yaml")
function build_polar_wing(polars_subdir, out_name)
    geom = YAML.load_file(yaml_path)
    for (i, airfoil) in enumerate(geom["wing_airfoils"]["data"])
        airfoil[2] = "polars"
        airfoil[3] = Dict("csv_file_path" => "$polars_subdir/$i.csv")
    end
    out_yaml = joinpath("data", "TUDELFT_V3_KITE", out_name)
    write_yaml(out_yaml, geom)
    wing = Wing(out_yaml; n_panels=50, spanwise_distribution=LINEAR)
    refine!(wing)
    body = BodyAerodynamics([wing])
    VortexStepMethod.reinit!(body)
    return Solver(body, settings_cfd), body
end

println("Creating wing with NeuralFoil polars...")
solver_nf, body_nf = build_polar_wing("polars_neuralfoil", "aero_geometry_neuralfoil.yaml")

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

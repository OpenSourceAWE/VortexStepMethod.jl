using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using MakieControlPlots
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: ShrinkWrap, NeuralFoilSolver

PLOT = true
SAVE_ALL = false
USE_TEX = false
DEFORM = false
NEURALFOIL = true
# Rolling-ball radius of the shrink wrap; fillets the concave tube-canopy junction.
MIN_CONCAVE_RADIUS = 0.4
OUTPUT_DIR = joinpath(dirname(@__DIR__), "output")
REFERENCE_POINT = [0.422646, 0.0, 9.3667]

project_dir = dirname(@__DIR__)
literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results",
        "CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results",
        "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results",
        "python_alpha_sweep.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results",
        "windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv"),
]
beta_literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results",
        "windtunnel_beta_sweep_alpha_07_4_Poland_2025_Rey_5e5.csv"),
]

settings = VSMSettings(joinpath(project_dir, "data", "TUDELFT_V3_KITE",
    "vsm_settings.yaml"); data_prefix=false)
settings.wings[1].geometry_file = joinpath(project_dir,
    settings.wings[1].geometry_file)
wing = Wing(settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)
solver = Solver(body_aero, settings)
solver.reference_point .= REFERENCE_POINT

if DEFORM
    VortexStepMethod.unrefined_deform!(
        wing,
        deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
        deg2rad.(range(0, 0, length=wing.n_unrefined_sections));
        smooth=true
    )
    VortexStepMethod.reinit!(body_aero; init_aero=false)
end

# Second sweep on generated polars: slice V3_25.obj, shrink-wrap every section into a
# closed airfoil and sweep it with NeuralFoil, so the same kite flies on polars derived
# from its own CAD surface instead of the checked-in CFD tables.
if NEURALFOIL
    obj_file = joinpath(project_dir, "data", "TUDELFT_V3_KITE", "V3_25.obj")
    generated_dir = joinpath(project_dir, "data", "TUDELFT_V3_KITE",
        "generated_neuralfoil")
    nf_yaml = obj_to_yaml(obj_file, generated_dir;
        n_sections=settings.wings[1].n_panels, Re=1e6, force=false,
        aero_solver=NeuralFoilSolver(model_size="large", n_crit=4.0,
            xtr_upper=0.05, xtr_lower=0.05),
        wrap_method=ShrinkWrap(clearance=0.0,
            min_concave_radius=MIN_CONCAVE_RADIUS),
    )
    settings_nf = deepcopy(settings)
    settings_nf.wings[1].geometry_file = nf_yaml
    wing_nf = Wing(settings_nf)
    refine!(wing_nf)
    body_nf = BodyAerodynamics([wing_nf])
    VortexStepMethod.reinit!(body_nf)
    solver_nf = Solver(body_nf, settings_nf)
    solver_nf.reference_point .= REFERENCE_POINT

    # Reading the generated directory instead of the OBJ shows the airfoils the polar
    # pipeline actually analysed. Hover a slice to inspect its 2D fit.
    PLOT && plot_slices_3d(generated_dir; obj_path=obj_file)
end

solvers = NEURALFOIL ? [solver, solver_nf] : [solver]
bodies = NEURALFOIL ? [body_aero, body_nf] : [body_aero]
solver_labels = NEURALFOIL ? ["VSM Julia CFD", "VSM Julia NeuralFoil"] :
                ["VSM Julia CFD"]
labels = [solver_labels;
    ["CFD Re=5e5",
     "CFD Re=10e5", #with struts
     "VSM Python Re=5e5",
     "WindTunnel Re=5e5"]] #with struts
beta_labels = [solver_labels; ["Wind Tunnel Re=5e5 beta sweep alpha=7.4"]]

wind_speed = settings.condition.wind_speed
angle_of_attack_deg = settings.condition.alpha
sideslip_deg = settings.condition.beta
yaw_rate = settings.condition.yaw_rate

set_va!(body_aero, settings)
results = VortexStepMethod.solve(solver, body_aero; log=true)

PLOT && plot_polars(
    solvers,
    bodies,
    labels,
    literature_path_list=literature_paths,
    angle_range=range(-5, 25, length=31),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="$(wing.n_panels)_panels_$(wing.spanwise_distribution)_from_yaml_settings",
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    use_tex=USE_TEX,
    show_moments=false,
    cl_over_cd=true
)

# Plotting geometry
PLOT && plot_geometry(
    body_aero,
    "V3 kite geometry";
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    view_elevation=15,
    view_azimuth=-120,
    use_tex=USE_TEX
)

# Plotting spanwise distributions
body_y_coordinates = [panel.aero_center[2]
                      for panel in body_aero.panels]

PLOT && plot_distribution(
    [body_y_coordinates],
    [results],
    ["VSM"];
    title="CAD_spanwise_distributions_alpha_$(round(angle_of_attack_deg, digits=1))_delta_$(round(sideslip_deg, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(wind_speed, digits=1))",
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    use_tex=USE_TEX
)

# --- Beta sweep ---
PLOT && plot_polars(
    solvers,
    bodies,
    beta_labels;
    literature_path_list=beta_literature_paths,
    angle_range=range(0, 12, step=1),
    angle_type="side_slip",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="beta sweep",
    show_moments=true,
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    use_tex=USE_TEX
)

nothing

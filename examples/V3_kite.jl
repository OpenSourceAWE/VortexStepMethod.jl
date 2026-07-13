using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using LinearAlgebra
using VortexStepMethod

PLOT = true
SAVE_ALL = false
USE_TEX = false
DEFORM = false
OUTPUT_DIR = joinpath(dirname(@__DIR__), "output")

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
labels = [
    "VSM Julia Re=5e5",
    "CFD Re=5e5",
    "CFD Re=10e5", #with struts
    "VSM Python Re=5e5",
    "WindTunnel Re=5e5" #with struts
]
beta_literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results",
        "windtunnel_beta_sweep_alpha_07_4_Poland_2025_Rey_5e5.csv"),
]
beta_labels = [
    labels[1],
    "Wind Tunnel Re=5e5 beta sweep alpha=7.4",
]

# Load YAML settings directly
settings_path = joinpath(
    project_dir, "data", "TUDELFT_V3_KITE", "vsm_settings.yaml")
settings_data = VortexStepMethod.YAML.load_file(settings_path)
condition_cfg = settings_data["condition"]
wing_cfg = settings_data["wings"][1]
solver_cfg = settings_data["solver_settings"]

# Create wing, body_aero, and solver objects using settings
wing = Wing(
    joinpath(project_dir, wing_cfg["geometry_file"]);
    n_panels=wing_cfg["n_panels"],
    spanwise_distribution=getproperty(
        VortexStepMethod,
        Symbol(wing_cfg["spanwise_panel_distribution"])),
    spanwise_direction=Float64.(wing_cfg["spanwise_direction"]),
    remove_nan=wing_cfg["remove_nan"],
)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)

if DEFORM
    VortexStepMethod.unrefined_deform!(
        wing,
        deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
        deg2rad.(range(0, 0, length=wing.n_unrefined_sections));
        smooth=true
    )
    VortexStepMethod.reinit!(body_aero; init_aero=false)
end

# Construct Solver using keyword arguments from solver settings
solver = Solver(body_aero;
    solver_type=(solver_cfg["solver_type"] == "NONLIN" ? NONLIN : LOOP),
    aerodynamic_model_type=getproperty(
        VortexStepMethod,
        Symbol(solver_cfg["aerodynamic_model_type"])),
    density=solver_cfg["density"],
    max_iterations=solver_cfg["max_iterations"],
    rtol=solver_cfg["rtol"],
    tol_reference_error=solver_cfg["tol_reference_error"],
    relaxation_factor=solver_cfg["relaxation_factor"],
    is_with_artificial_damping=solver_cfg["artificial_damping"],
    artificial_damping=(k2=solver_cfg["k2"], k4=solver_cfg["k4"]),
    type_initial_gamma_distribution=getproperty(
        VortexStepMethod,
        Symbol(solver_cfg["type_initial_gamma_distribution"])),
    use_gamma_prev=get(solver_cfg, "use_gamma_prev",
        get(solver_cfg, "use_gamme_prev", true)),
    core_radius_fraction=solver_cfg["core_radius_fraction"],
    mu=solver_cfg["mu"],
    is_only_f_and_gamma_output=get(
        solver_cfg, "calc_only_f_and_gamma", false),
    correct_aoa=get(solver_cfg, "correct_aoa", false),
    reference_point=get(solver_cfg, "reference_point",
        [0.422646, 0.0, 9.3667]),
)

# Extract values for plotting
wind_speed = condition_cfg["wind_speed"]
angle_of_attack_deg = condition_cfg["alpha"]
sideslip_deg = condition_cfg["beta"]
yaw_rate = condition_cfg["yaw_rate"]

# Set flight conditions from settings
α0 = deg2rad(angle_of_attack_deg)
β0 = deg2rad(sideslip_deg)
set_va!(body_aero,
    wind_speed .* [cos(α0) * cos(β0), sin(β0), sin(α0) * cos(β0)])

# Solve
results = VortexStepMethod.solve(solver, body_aero; log=true)

# Plotting polars with moment coefficients
PLOT && plot_polars(
    [solver],
    [body_aero],
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
    show_moments=true
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

# --- Dual solver comparison: NONLIN vs LOOP ---
solver_cfg["solver_type"] = "LOOP"
solver_loop = Solver(body_aero;
    solver_type=LOOP,
    aerodynamic_model_type=getproperty(
        VortexStepMethod,
        Symbol(solver_cfg["aerodynamic_model_type"])),
    density=solver_cfg["density"],
    max_iterations=solver_cfg["max_iterations"],
    rtol=solver_cfg["rtol"],
    tol_reference_error=solver_cfg["tol_reference_error"],
    relaxation_factor=solver_cfg["relaxation_factor"],
    is_with_artificial_damping=solver_cfg["artificial_damping"],
    artificial_damping=(k2=solver_cfg["k2"], k4=solver_cfg["k4"]),
    type_initial_gamma_distribution=getproperty(
        VortexStepMethod,
        Symbol(solver_cfg["type_initial_gamma_distribution"])),
    use_gamma_prev=get(solver_cfg, "use_gamma_prev",
        get(solver_cfg, "use_gamme_prev", true)),
    core_radius_fraction=solver_cfg["core_radius_fraction"],
    mu=solver_cfg["mu"],
    is_only_f_and_gamma_output=get(
        solver_cfg, "calc_only_f_and_gamma", false),
    correct_aoa=get(solver_cfg, "correct_aoa", false),
    reference_point=get(solver_cfg, "reference_point",
        [0.422646, 0.0, 9.3667]),
)

PLOT && plot_polars(
    [solver_loop],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 20, step=1),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="LOOP solver",
    show_moments=true,
    save_path=OUTPUT_DIR,
    is_save=false || SAVE_ALL,
    is_show=true,
    use_tex=USE_TEX
)

# --- Beta sweep ---
PLOT && plot_polars(
    [solver_loop],
    [body_aero],
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
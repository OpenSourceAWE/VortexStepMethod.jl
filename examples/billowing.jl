using LinearAlgebra
using VortexStepMethod

PLOT = true
SAVE_ALL = false
USE_TEX = false
OUTPUT_DIR = joinpath(dirname(@__DIR__), "output")

# Data paths (all within this repo)
vsm_src_path = something(pathof(VortexStepMethod), @__FILE__)
project_dir = normpath(joinpath(dirname(vsm_src_path), ".."))
v3_dir = joinpath(project_dir, "data", "TUDELFT_V3_KITE")
polar_dir = joinpath(v3_dir, "polars_CFD_NF_combined")

# Literature results
lit_dir = joinpath(v3_dir, "literature_results")
literature_paths = [
    joinpath(lit_dir,
        "CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv"),
    joinpath(lit_dir,
        "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(lit_dir, "python_alpha_sweep.csv"),
    joinpath(lit_dir,
        "windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv"),
]

# Load solver settings (coarse: 54 panels, matches 10-section geometry)
settings_data = VortexStepMethod.YAML.load_file(
    joinpath(v3_dir, "vsm_settings_coarse.yaml"))
condition_cfg = settings_data["condition"]
solver_cfg = settings_data["solver_settings"]
wing_cfg = settings_data["wings"][1]
n_panels = wing_cfg["n_panels"]

BILLOWING_PCT = get(wing_cfg, "billowing_percentage", 0.0)

labels = [
    "VSM flat",
    "VSM billowing $(BILLOWING_PCT)%",
    "CFD Re=5e5",
    "CFD Re=10e5",
    "VSM Python Re=5e5",
    "WindTunnel Re=5e5",
]

# Load coarse geometry (10 structural rib sections)
geom_data = VortexStepMethod.YAML.load_file(
    joinpath(v3_dir, "aero_geometry_coarse_discretisation.yaml"))
section_headers = geom_data["wing_sections"]["headers"]
section_rows = geom_data["wing_sections"]["data"]

function build_wing(; distribution=SPLIT_PROVIDED,
                      billowing_percentage=0.0)
    wing = Wing(n_panels;
        spanwise_distribution=distribution,
        billowing_percentage=billowing_percentage)
    for row in section_rows
        d = Dict(zip(section_headers, row))
        le = [d["LE_x"], d["LE_y"], d["LE_z"]]
        te = [d["TE_x"], d["TE_y"], d["TE_z"]]
        csv_path = joinpath(polar_dir, "$(d["airfoil_id"]).csv")
        aero_data, aero_model = load_polar_data(csv_path)
        add_section!(wing, le, te, aero_model, aero_data)
    end
    refine!(wing)
    return wing
end

# --- Wing without billowing ---
wing_flat = build_wing()
body_aero_flat = BodyAerodynamics([wing_flat])
VortexStepMethod.reinit!(body_aero_flat)

# --- Wing with billowing ---
wing_bill = build_wing(distribution=BILLOWING,
                       billowing_percentage=BILLOWING_PCT)
body_aero_bill = BodyAerodynamics([wing_bill])
VortexStepMethod.reinit!(body_aero_bill)

# --- Build solvers ---
function make_solver(body_aero)
    Solver(body_aero;
        solver_type=(solver_cfg["solver_type"] == "NONLIN" ?
            NONLIN : LOOP),
        aerodynamic_model_type=getproperty(
            VortexStepMethod,
            Symbol(solver_cfg["aerodynamic_model_type"])),
        density=solver_cfg["density"],
        max_iterations=solver_cfg["max_iterations"],
        rtol=solver_cfg["rtol"],
        tol_reference_error=solver_cfg["tol_reference_error"],
        relaxation_factor=solver_cfg["relaxation_factor"],
        is_with_artificial_damping=solver_cfg["artificial_damping"],
        artificial_damping=(
            k2=solver_cfg["k2"], k4=solver_cfg["k4"]),
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
end

solver_flat = make_solver(body_aero_flat)
solver_bill = make_solver(body_aero_bill)

# --- Set flight conditions ---
wind_speed = condition_cfg["wind_speed"]
angle_of_attack_deg = condition_cfg["alpha"]
sideslip_deg = condition_cfg["beta"]

α0 = deg2rad(angle_of_attack_deg)
β0 = deg2rad(sideslip_deg)
va = wind_speed .* [cos(α0) * cos(β0), sin(β0), sin(α0) * cos(β0)]
set_va!(body_aero_flat, va)
set_va!(body_aero_bill, va)

# --- Solve and compare ---
results_flat = VortexStepMethod.solve(
    solver_flat, body_aero_flat; log=true)
results_bill = VortexStepMethod.solve(
    solver_bill, body_aero_bill; log=true)

println("\nFlat wing: CL=$(round(results_flat["cl"]; digits=4)), " *
        "CD=$(round(results_flat["cd"]; digits=4))")
println("Billowed:  CL=$(round(results_bill["cl"]; digits=4)), " *
        "CD=$(round(results_bill["cd"]; digits=4))")

if PLOT
    # Plot geometry (flat wing)
    plot_geometry(
        body_aero_flat,
        "Flat wing geometry";
        save_path=OUTPUT_DIR,
        is_save=false || SAVE_ALL,
        is_show=true,
        use_tex=USE_TEX
    )

    # Plot spanwise distributions
    y_flat = [panel.aero_center[2]
              for panel in body_aero_flat.panels]
    y_bill = [panel.aero_center[2]
              for panel in body_aero_bill.panels]
    plot_distribution(
        [y_flat, y_bill],
        [results_flat, results_bill],
        ["VSM flat", "VSM billowing"];
        title="Billowing comparison distributions",
        save_path=OUTPUT_DIR,
        is_save=false || SAVE_ALL,
        is_show=true,
        use_tex=USE_TEX
    )

    # Plot polars comparison
    plot_polars(
        [solver_flat, solver_bill],
        [body_aero_flat, body_aero_bill],
        labels;
        literature_path_list=literature_paths,
        angle_range=range(-5, 25, length=31),
        angle_type="angle_of_attack",
        angle_of_attack=angle_of_attack_deg,
        side_slip=sideslip_deg,
        v_a=wind_speed,
        title="V3 Kite flat vs billowing $(BILLOWING_PCT)%",
        save_path=OUTPUT_DIR,
        is_save=false || SAVE_ALL,
        is_show=true,
        use_tex=USE_TEX,
        show_moments=true
    )
end

nothing

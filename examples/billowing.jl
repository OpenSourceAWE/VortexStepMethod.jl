using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using LinearAlgebra
using GLMakie
using MakieControlPlots
using VortexStepMethod

PLOT = true
SAVE_ALL = true
USE_TEX = true
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

# Coarse settings: 54 panels on the 10-section geometry
settings = VSMSettings(joinpath(v3_dir, "vsm_settings_coarse.yaml"); data_prefix=false)
settings.wings[1].geometry_file = joinpath(project_dir, settings.wings[1].geometry_file)
n_panels = settings.wings[1].n_panels

BILLOWING_PCT = settings.wings[1].billowing_percentage

labels = [
    "VSM flat",
    "VSM billowing $(BILLOWING_PCT)%",
    "CFD Re=5e5",
    "CFD Re=10e5",
    "VSM Python Re=5e5",
    "WindTunnel Re=5e5",
]

geom_data = VortexStepMethod.YAML.load_file(settings.wings[1].geometry_file)
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
solver_flat = Solver(settings; reference_point=[0.422646, 0.0, 9.3667])
solver_bill = Solver(settings; reference_point=[0.422646, 0.0, 9.3667])

# --- Set flight conditions ---
va = settings.condition.va
angle_of_attack_deg = 10.0
sideslip_deg = settings.condition.beta

α0 = deg2rad(angle_of_attack_deg)
β0 = deg2rad(sideslip_deg)
va_vec = apparent_wind(α0, β0, va)
set_va!(body_aero_flat, va_vec)
set_va!(body_aero_bill, va_vec)

# --- Solve and compare ---
results_flat = solve!(
    solver_flat, body_aero_flat; log=true)
results_bill = solve!(
    solver_bill, body_aero_bill; log=true)

println("\nFlat wing: CL=$(round(results_flat.cl; digits=4)), " *
        "CD=$(round(results_flat.cd; digits=4))")
println("Billowed:  CL=$(round(results_bill.cl; digits=4)), " *
        "CD=$(round(results_bill.cd; digits=4))")

if PLOT
    # Plot geometry (flat wing)
    plot_geometry(
        body_aero_bill,
        "Billowing wing geometry";
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
        va=va,
        title="V3 Kite flat vs billowing $(BILLOWING_PCT)%",
        save_path=OUTPUT_DIR,
        is_save=false || SAVE_ALL,
        is_show=true,
        use_tex=USE_TEX,
        show_moments=true
    )
end

nothing
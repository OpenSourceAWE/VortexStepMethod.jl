# Renders docs/TU_Delft_V3_Kite.png, the V3 kite figure shown in the README.
using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using CairoMakie
using MakieControlPlots
CairoMakie.activate!()
using VortexStepMethod

project_dir = dirname(@__DIR__)
literature_dir = joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results")
literature_paths = joinpath.(literature_dir, [
    "CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv",
    "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv",
    "python_alpha_sweep.csv",
    "windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv",
])
labels = ["VSM Julia", "CFD Re=5e5", "CFD Re=10e5", "VSM Python Re=5e5",
    "Wind tunnel Re=5e5"]

settings = VSMSettings(joinpath(project_dir, "data", "TUDELFT_V3_KITE",
    "vsm_settings.yaml"); data_prefix=false)
settings.wings[1].geometry_file = joinpath(project_dir,
    settings.wings[1].geometry_file)
wing = Wing(settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)
solver = Solver(settings)
solver.reference_point .= [0.422646, 0.0, 9.3667]
set_va!(body_aero, settings)
results = VortexStepMethod.solve(solver, body_aero)

plot_combined_analysis(solver, body_aero, results;
    labels,
    literature_path_list=literature_paths,
    angle_range=range(-5, 25, length=31),
    va=settings.condition.va,
    angle_of_attack_for_spanwise_distribution=10.0,
    title="TU Delft V3 Kite",
    is_show=false,
    is_save=true,
    save_path=joinpath(project_dir, "docs"),
    data_type=".png",
)

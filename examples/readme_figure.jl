# Renders docs/TU_Delft_V3_Kite.png, the V3 kite figure shown in the README.
using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using CairoMakie
using MakieControlPlots
CairoMakie.activate!()
using VortexStepMethod

include("V3_kite_setup.jl")
labels = ["VSM Julia", "CFD Re=5e5", "CFD Re=10e5", "VSM Python Re=5e5",
    "Wind tunnel Re=5e5"]

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

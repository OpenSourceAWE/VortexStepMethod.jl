# TU Delft V3 kite from vsm_settings.yaml, shared by V3_kite.jl and readme_figure.jl.
using VortexStepMethod

REFERENCE_POINT = [0.422646, 0.0, 9.3667]

project_dir = dirname(@__DIR__)
literature_dir = joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results")
literature_paths = joinpath.(literature_dir, [
    "CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv",
    "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv",
    "python_alpha_sweep.csv",
    "windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv",
])

settings = VSMSettings(joinpath(project_dir, "data", "TUDELFT_V3_KITE",
    "vsm_settings.yaml"); data_prefix=false)
settings.wings[1].geometry_file = joinpath(project_dir,
    settings.wings[1].geometry_file)
wing = Wing(settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)
solver = Solver(settings)
solver.reference_point .= REFERENCE_POINT

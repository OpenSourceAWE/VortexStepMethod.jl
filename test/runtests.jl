using Test, VortexStepMethod

# Make paths robust (avoid cd(".."))
cd(@__DIR__)  # ensure we're in test/ no matter how tests are launched
include("test_data_utils.jl")

# ControlPlots must run in a separate single-threaded process
const _plot_controlplots = "plot-controlplots" in ARGS

# Filter special args from pattern matching
const test_patterns = filter(a -> a != "plot-controlplots", ARGS)

println("Running tests...")
if _plot_controlplots
    println("Running plotting tests with ControlPlots backend")
elseif !isempty(test_patterns)
    println("Filtering tests matching: ", test_patterns)
end

# Helper to check if a test file matches any pattern
function should_run_test(test_path::String)
    isempty(test_patterns) && return true
    for pattern in test_patterns
        # Match directory (e.g., "solver") or specific file (e.g., "test_unrefined_dist")
        if occursin(pattern, test_path)
            return true
        end
    end
    return false
end

# keep your env check as-is...
const build_is_production_build_env_name = "BUILD_IS_PRODUCTION_BUILD"
const build_is_production_build = let v = get(ENV, build_is_production_build_env_name, "true")
    v == "true" || v == "false" || error("unknown value for $build_is_production_build_env_name: $v")
    v == "true"
end::Bool

@testset verbose = true "Testing VortexStepMethod..." begin
    if _plot_controlplots
        include("plotting/test_plotting.jl")
    else
        if build_is_production_build && should_run_test("bench")
            include("bench.jl")
        end
        should_run_test("body_aerodynamics/test_body_aerodynamics.jl") && include("body_aerodynamics/test_body_aerodynamics.jl")
        should_run_test("body_aerodynamics/test_results.jl") && include("body_aerodynamics/test_results.jl")
        should_run_test("test_refinement_validation.jl") && include("test_refinement_validation.jl")
        should_run_test("filament/test_bound_filament.jl") && include("filament/test_bound_filament.jl")
        should_run_test("filament/test_semi_infinite_filament.jl") && include("filament/test_semi_infinite_filament.jl")
        should_run_test("panel/test_panel.jl") && include("panel/test_panel.jl")
        should_run_test("plotting/test_plotting.jl") && include("plotting/test_plotting.jl")
        should_run_test("polars/test_polars.jl") && include("polars/test_polars.jl")
        should_run_test("ram_geometry/test_kite_geometry.jl") && include("ram_geometry/test_kite_geometry.jl")
        should_run_test("settings/test_settings.jl") && include("settings/test_settings.jl")
        should_run_test("solver/test_solver.jl") && include("solver/test_solver.jl")
        should_run_test("solver/test_unrefined_dist.jl") && include("solver/test_unrefined_dist.jl")
        should_run_test("VortexStepMethod/test_VortexStepMethod.jl") && include("VortexStepMethod/test_VortexStepMethod.jl")
        should_run_test("wake/test_wake.jl") && include("wake/test_wake.jl")
        should_run_test("wing_geometry/test_wing_geometry.jl") && include("wing_geometry/test_wing_geometry.jl")
        should_run_test("yaml_geometry/test_yaml_geometry.jl") && include("yaml_geometry/test_yaml_geometry.jl")
        should_run_test("Aqua.jl") && include("Aqua.jl")
    end
end

nothing

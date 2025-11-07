using Test, VortexStepMethod

# Make paths robust (avoid cd(".."))
cd(@__DIR__)  # ensure we're in test/ no matter how tests are launched
include("TestSupport.jl")
using .TestSupport

println("Running tests...")

# keep your env check as-is...
const build_is_production_build_env_name = "BUILD_IS_PRODUCTION_BUILD"
const build_is_production_build = let v = get(ENV, build_is_production_build_env_name, "true")
    v == "true" || v == "false" || error("unknown value for $build_is_production_build_env_name: $v")
    v == "true"
end::Bool

@testset verbose = true "Testing VortexStepMethod..." begin
    if build_is_production_build
        include("bench.jl")
    end
    include("body_aerodynamics/test_body_aerodynamics.jl")
    include("body_aerodynamics/test_results.jl")
    include("filament/test_bound_filament.jl")
    include("filament/test_semi_infinite_filament.jl")
    include("panel/test_panel.jl")
    include("plotting/test_plotting.jl")
    include("polars/test_polars.jl")
    include("ram_geometry/test_kite_geometry.jl")
    include("settings/test_settings.jl")
    include("solver/test_solver.jl")
    include("solver/test_group_coefficients.jl")
    include("VortexStepMethod/test_VortexStepMethod.jl")
    include("wake/test_wake.jl")
    include("wing_geometry/test_wing_geometry.jl")
    include("yaml_geometry/test_yaml_geometry.jl")
    include("Aqua.jl")
end

nothing

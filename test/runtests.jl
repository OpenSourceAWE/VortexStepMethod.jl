using Test, VortexStepMethod

# TODO fix allocation tests
# and: https://github.com/JuliaArrays/FixedSizeArrays.jl/blob/d17373edf4144a672cd80a062bf24d017f01e42f/.github/workflows/UnitTests.yml

# Check if the compilation options allow maximum performance.
const build_is_production_build_env_name = "BUILD_IS_PRODUCTION_BUILD"
const build_is_production_build = let v = get(ENV, build_is_production_build_env_name, "true")
    if v ∉ ("false", "true")
        error("unknown value for environment variable $build_is_production_build_env_name: $v")
    end
    if v == "true"
        true
    else
        false
    end
end::Bool

cd("..")
println("Running tests...")

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
    include("VortexStepMethod/test_VortexStepMethod.jl")
    include("wake/test_wake.jl")
    include("wing_geometry/test_wing_geometry.jl")
    include("yaml_geometry/test_yaml_geometry.jl")
    include("Aqua.jl")

end
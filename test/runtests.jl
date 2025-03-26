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

    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
    cp("data/ram_air_kite_body.obj", body_path; force=true)
    cp("data/ram_air_kite_foil.dat", foil_path; force=true)
    ram_wing = RamAirWing(body_path, foil_path; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1))

    if build_is_production_build
        include("bench.jl")
    end
    include("test_bound_filament.jl")
    include("test_panel.jl")
    include("test_semi_infinite_filament.jl")
    include("test_body_aerodynamics.jl")
    include("test_results.jl")
    include("test_kite_geometry.jl")
    include("test_wing_geometry.jl")
    include("test_plotting.jl")
    include("aqua.jl")
end
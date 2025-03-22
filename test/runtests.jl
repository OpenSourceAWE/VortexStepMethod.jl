using Test

# TODO fix allocation tests
# see: https://github.com/JuliaArrays/FixedSizeArrays.jl/blob/d17373edf4144a672cd80a062bf24d017f01e42f/test/runtests.jl#L6-L17
# and: https://github.com/JuliaArrays/FixedSizeArrays.jl/blob/d17373edf4144a672cd80a062bf24d017f01e42f/.github/workflows/UnitTests.yml
cd("..")
println("Running tests...")
@testset verbose = true "Testing VortexStepMethod..." begin
    include("bench.jl")
    include("test_bound_filament.jl")
    include("test_panel.jl")
    include("test_semi_infinite_filament.jl")
    include("test_body_aerodynamics.jl")
    include("test_kite_geometry.jl")
    include("test_wing_geometry.jl")
    include("test_plotting.jl")
    include("aqua.jl")
end
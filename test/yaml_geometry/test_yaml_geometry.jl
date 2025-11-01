using Test

@testset "YAML Geometry Tests" begin
    # Include specific test files for better organization
    include("test_load_polar_data.jl")
    include("test_wing_constructor.jl")
    include("test_yaml_wing_deformation.jl")
end

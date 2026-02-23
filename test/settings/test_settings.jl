using VortexStepMethod
using Test

@testset "Test settings.jl" begin
    # Use the absolute path to test data to avoid path concatenation issues
    project_root = dirname(dirname(@__DIR__))
    settings_file = joinpath(project_root, "data", "ram_air_kite", "vsm_settings_dual.yaml")
    vss = VSMSettings(settings_file)
    @test vss isa VSMSettings
    @test vss.solver_settings isa SolverSettings
    @test vss.wings isa Vector{WingSettings}
    @test length(vss.wings) == 2
    io = IOBuffer(repr(vss))
    @test countlines(io) > 0
end
nothing

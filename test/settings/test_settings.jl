using Pkg
if ! ("Test" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using VortexStepMethod
using Test

@testset "Test settings.jl" begin
    # Change to project root directory for the test
    project_root = dirname(dirname(@__DIR__))
    cd(project_root) do
        vss = VSMSettings("data/ram_air_kite/vsm_settings_dual.yaml")
        @test vss isa VSMSettings
        @test vss.solver_settings isa SolverSettings
        @test vss.wings isa Vector{WingSettings}
        @test length(vss.wings) == 2
        io = IOBuffer(repr(vss))
        @test countlines(io) == 40  # Updated to match new output format
    end
end
nothing

using Pkg
if ! ("Test" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using VortexStepMethod
using Test

@testset "Test settings.jl" begin
    vss = vs("vsm_settings_dual.yaml")
    @test vss isa VSMSettings
    @test vss.solver_settings isa SolverSettings
    @test vss.wings isa Vector{WingSettings}
    @test length(vss.wings) == 2
end
nothing

using Aqua, VortexStepMethod, Test
@testset "Aqua.jl" begin
    Aqua.test_all(
      VortexStepMethod;
      stale_deps=(ignore=[:Xfoil, :Timers, :PyCall],),
      deps_compat=(ignore=[:PyCall],),
      persistent_tasks=false
    )
end

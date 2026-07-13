using Aqua, VortexStepMethod, Test
@testset "Aqua.jl" begin
    Aqua.test_all(
      VortexStepMethod;
      stale_deps=(ignore=[:Xfoil, :Timers],),
      persistent_tasks=false
    )
end

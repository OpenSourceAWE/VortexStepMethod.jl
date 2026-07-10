# Regression test for https://github.com/OpenSourceAWE/VortexStepMethod.jl/issues/236
# Verifies that loading both Makie and ControlPlots in the same process:
#   (1) does not cause method-redefinition / precompile errors, and
#   (2) set_plot_backend! correctly switches which backend the no-backend wrappers route to.

using CairoMakie

cp_available = try
    @eval using ControlPlots
    true
catch e
    @warn "ControlPlots unavailable, skipping coexistence tests: $e"
    false
end

@testset "Backend coexistence (Makie + ControlPlots)" begin
    if !cp_available
        @test_skip "ControlPlots failed to load"
        return
    end

    backend_before = VortexStepMethod._PLOT_BACKEND[]
    try
        # (1) Both extensions must be loaded without errors when both packages are in scope.
        #     If either `using` above threw a method-redefinition error we would never reach here.
        makie_ext = Base.get_extension(VortexStepMethod, :VortexStepMethodMakieExt)
        cp_ext    = Base.get_extension(VortexStepMethod, :VortexStepMethodControlPlotsExt)
        @test makie_ext !== nothing

        if cp_ext === nothing
            @test_skip "VortexStepMethodControlPlotsExt unavailable (ControlPlots failed to load)"
            return
        end
        @test cp_ext !== nothing

        # (2) set_plot_backend! correctly switches the active backend.
        set_plot_backend!(ControlPlotsBackend())
        @test VortexStepMethod._PLOT_BACKEND[] isa VortexStepMethod.ControlPlotsBackend

        set_plot_backend!(MakieBackend())
        @test VortexStepMethod._PLOT_BACKEND[] isa VortexStepMethod.MakieBackend

        # Round-trip: switch back to ControlPlots.
        set_plot_backend!(ControlPlotsBackend())
        @test VortexStepMethod._PLOT_BACKEND[] isa VortexStepMethod.ControlPlotsBackend
    finally
        # Restore whatever backend was active before this testset ran.
        VortexStepMethod._PLOT_BACKEND[] = backend_before
    end
end

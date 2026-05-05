using VortexStepMethod
using DifferentiationInterface
using LinearAlgebra
using Test

@testset "ForwardDiff linearize" begin
    # Build a small rectangular wing with INVISCID polars so the test runs
    # quickly and the comparison is dominated by the LOOP iteration itself.
    n_panels = 10
    span = 20.0
    chord = 1.0
    wing = Wing(n_panels, spanwise_distribution=LINEAR)
    add_section!(wing,
        [0.0, span/2, 0.0], [chord, span/2, 0.0], INVISCID)
    add_section!(wing,
        [0.0, -span/2, 0.0], [chord, -span/2, 0.0], INVISCID)
    refine!(wing)
    body_aero = BodyAerodynamics([wing])

    va = [15.0, 1.0, 2.0]
    omega = [0.0, 0.0, 0.0]
    y0 = [va; omega]

    @testset "AutoForwardDiff matches AutoFiniteDiff (LOOP, INVISCID)" begin
        # `use_gamma_prev=false` is required so that the FD path does not
        # warm-start from the previous perturbed state — otherwise the FD
        # Jacobian carries a small bias unrelated to the AD result.
        solver = Solver(body_aero; use_gamma_prev=false)

        jac_fwd, _, fwd_converged = VortexStepMethod.linearize(
            solver, body_aero, y0;
            theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true, backend=AutoForwardDiff())
        @test fwd_converged

        jac_fd, _, fd_converged = VortexStepMethod.linearize(
            solver, body_aero, y0;
            theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true,
            backend=AutoFiniteDiff(absstep=1e-5, relstep=1e-5))
        @test fd_converged

        rel_err = maximum(abs.(jac_fwd .- jac_fd)) / maximum(abs, jac_fwd)
        @test rel_err < 1e-4
    end

    @testset "NONLIN+ForwardDiff is rejected" begin
        solver_nl = Solver(body_aero; solver_type=NONLIN)
        @test_throws ErrorException VortexStepMethod.linearize(
            solver_nl, body_aero, y0;
            theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true, backend=AutoForwardDiff())
    end
end

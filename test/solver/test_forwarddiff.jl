using VortexStepMethod
using DifferentiationInterface
using LinearAlgebra
using Test

relative_error(jac, reference) = maximum(abs.(jac .- reference)) / maximum(abs, reference)

@testset "ForwardDiff linearize" begin
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

    va_vec = [15.0, 1.0, 2.0]
    omega = [0.0, 0.0, 0.0]
    y0 = [va_vec; omega]

    turns = ((omega, zeros(3)), ([0.0, 0.0, 0.2], [0.5, 4.0, 0.0]))
    @testset "ForwardDiff matches FiniteDiff about $reference_point (LOOP, INVISCID)" for
            (omega_op, reference_point) in turns
        pivot_body = BodyAerodynamics([wing])
        set_va!(pivot_body, va_vec, omega_op; reference_point)
        y_op = [va_vec; omega_op]
        solver = Solver(wing.n_panels, wing.n_unrefined_sections;
            use_gamma_prev=false,
            type_initial_gamma_distribution=ELLIPTIC)

        jac_fwd, _, fwd_converged = VortexStepMethod.linearize(
            solver, pivot_body, y_op;
            theta_idxs=nothing, va_vec_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true, backend=AutoForwardDiff())
        @test fwd_converged

        jac_fd, _, fd_converged = VortexStepMethod.linearize(
            solver, pivot_body, y_op;
            theta_idxs=nothing, va_vec_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true,
            backend=AutoFiniteDiff(absstep=1e-5, relstep=1e-5))
        @test fd_converged

        @info "INVISCID jacobian norms" norm_fwd=norm(jac_fwd) norm_fd=norm(jac_fd)
        @test relative_error(jac_fd, jac_fwd) < 1e-3
    end

    @testset "NONLIN+ForwardDiff is rejected" begin
        solver_nl = Solver(wing.n_panels, wing.n_unrefined_sections; solver_type=NONLIN)
        @test_throws ErrorException VortexStepMethod.linearize(
            solver_nl, body_aero, y0;
            theta_idxs=nothing, va_vec_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true, backend=AutoForwardDiff())
    end

    @testset "AutoForwardDiff matches AutoFiniteDiff (LOOP, POLAR_MATRICES)" begin
        # At this operating point the LOOP solve converges the mid-span
        # panels' local alpha to ~6.5-9.1deg and the tip panels' to
        # ~14.9-15.0deg. The default 5deg-spaced alpha_range=-5:5:15 puts the
        # tip cluster within ulp-scale distance of the 15deg knot on some
        # platforms (Windows, Julia 1.12), so AutoForwardDiff and the FD
        # reference land on opposite sides of that non-differentiable kink
        # and disagree by several percent (#360). Offset the grid so no knot
        # is near either cluster (margin >0.8deg here vs <0.04deg before).
        ram_wing = ram_air_matrix_wing(; n_panels=8, n_sections=4,
            alpha_range=deg2rad.(-2:5:23),
            delta_range=deg2rad.(-3:3:3),
        )
        ram_body = BodyAerodynamics([ram_wing])
        ram_solver = Solver(ram_wing.n_panels, ram_wing.n_unrefined_sections;
            aerodynamic_model_type=VSM,
            rtol=1e-11,
            solver_type=LOOP,
            use_gamma_prev=false,
        )

        va = 15.0
        aoa_rad = deg2rad(7.5)
        y_op = [zeros(4);
                [cos(aoa_rad), 0.0, sin(aoa_rad)] * va;
                zeros(3)]

        jac_fwd, _, conv_fwd = VortexStepMethod.linearize(
            ram_solver, ram_body, y_op;
            theta_idxs=1:4, va_vec_idxs=5:7, omega_idxs=8:10,
            aero_coeffs=true, backend=AutoForwardDiff())
        @test conv_fwd

        jac_fd, _, conv_fd = VortexStepMethod.linearize(
            ram_solver, ram_body, y_op;
            theta_idxs=1:4, va_vec_idxs=5:7, omega_idxs=8:10,
            aero_coeffs=true, backend=nothing)
        @test conv_fd

        @info "POLAR_MATRICES jacobian norms" norm_fwd=norm(jac_fwd) norm_fd=norm(jac_fd)
        @test relative_error(jac_fd, jac_fwd) < 1e-4

        @testset "the Jacobian does not depend on the ForwardDiff chunk size" begin
            jac_chunk5, _, _ = VortexStepMethod.linearize(
                ram_solver, ram_body, y_op;
                theta_idxs=1:4, va_vec_idxs=5:7, omega_idxs=8:10,
                aero_coeffs=true, backend=AutoForwardDiff(chunksize=5))
            @test relative_error(jac_chunk5, jac_fwd) < 1e-12
        end
    end
end

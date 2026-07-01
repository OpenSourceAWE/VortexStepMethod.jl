using VortexStepMethod
using DifferentiationInterface
using LinearAlgebra
using Test

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

    va = [15.0, 1.0, 2.0]
    omega = [0.0, 0.0, 0.0]
    y0 = [va; omega]

    @testset "AutoForwardDiff matches AutoFiniteDiff (LOOP, INVISCID)" begin
        solver = Solver(body_aero;
            use_gamma_prev=false,
            type_initial_gamma_distribution=ELLIPTIC)

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

        @info "INVISCID linearize jacobian norms" norm_fwd=norm(jac_fwd) norm_fd=norm(jac_fd)
        rel_err = maximum(abs.(jac_fwd .- jac_fd)) / maximum(abs, jac_fwd)
        @test rel_err < 1e-3
    end

    @testset "NONLIN+ForwardDiff is rejected" begin
        solver_nl = Solver(body_aero; solver_type=NONLIN)
        @test_throws ErrorException VortexStepMethod.linearize(
            solver_nl, body_aero, y0;
            theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6,
            aero_coeffs=true, backend=AutoForwardDiff())
    end

    @testset "AutoForwardDiff matches AutoFiniteDiff (LOOP, POLAR_MATRICES)" begin
        data_dir = joinpath(dirname(dirname(@__DIR__)), "data", "ram_air_kite")
        body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
        foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
        cp(joinpath(data_dir, "ram_air_kite_body.obj"), body_path; force=true)
        cp(joinpath(data_dir, "ram_air_kite_foil.dat"), foil_path; force=true)
        for kind in ("cl", "cd", "cm")
            name = "ram_air_kite_foil_$(kind)_polar.csv"
            cp(joinpath(data_dir, name), joinpath(tempdir(), name); force=true)
        end

        ram_wing = ram_air_matrix_wing(; n_panels=8, n_sections=4,
            alpha_range=deg2rad.(-5:5:15),
            delta_range=deg2rad.(-3:3:3),
        )
        ram_body = BodyAerodynamics([ram_wing])
        ram_solver = Solver(ram_body;
            aerodynamic_model_type=VSM,
            is_with_artificial_damping=false,
            rtol=1e-7,
            solver_type=LOOP,
            use_gamma_prev=false,
        )

        v_a = 15.0
        aoa_rad = deg2rad(7.5)  # off-grid (grid is every 1°) to avoid piecewise-linear node discontinuities
        y_op = [zeros(4);
                [cos(aoa_rad), 0.0, sin(aoa_rad)] * v_a;
                zeros(3)]

        jac_fwd, _, conv_fwd = VortexStepMethod.linearize(
            ram_solver, ram_body, y_op;
            theta_idxs=1:4, va_idxs=5:7, omega_idxs=8:10,
            aero_coeffs=true, backend=AutoForwardDiff())
        @test conv_fwd

        jac_fd, _, conv_fd = VortexStepMethod.linearize(
            ram_solver, ram_body, y_op;
            theta_idxs=1:4, va_idxs=5:7, omega_idxs=8:10,
            aero_coeffs=true,
            backend=AutoFiniteDiff(absstep=1e-5, relstep=1e-5))
        @test conv_fd

        @info "POLAR_MATRICES linearize jacobian norms" norm_fwd=norm(jac_fwd) norm_fd=norm(jac_fd)
        rel_err = maximum(abs.(jac_fwd .- jac_fd)) / maximum(abs, jac_fwd)
        @test rel_err < 1e-3
    end
end

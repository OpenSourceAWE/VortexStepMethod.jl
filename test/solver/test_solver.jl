using VortexStepMethod
using LinearAlgebra
using Test
if !@isdefined(test_data_path)
    include("../test_data_utils.jl")
end

@testset "Solver Constructor Tests" begin
    @testset "Solver Constructor with VSMSettings" begin
        # Use module-specific test data files
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)

        try
            # Test Solver constructor with VSMSettings
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            refine!(wing)
            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            # Verify solver properties match settings
            @test solver.aerodynamic_model_type == VSM
            @test solver.density == 1.225

            # Test that the solver can solve
            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)
            sol = solve!(solver, body_aero)
            @test sol isa VSMSolution

        finally
            # Cleanup
            rm(settings_file; force=true)
        end
    end
end

@testset "NONLIN solve! re-runs across calls" begin
    settings_file = create_temp_wing_settings(
        "solver", "solver_test_wing.yaml";
        alpha=5.0, beta=0.0, wind_speed=10.0,
    )
    try
        settings = VSMSettings(settings_file)
        wing = Wing(settings)
        refine!(wing)

        body_aero = BodyAerodynamics([wing])
        solver = Solver(
            body_aero;
            solver_type=NONLIN,
            aerodynamic_model_type=VSM,
            type_initial_gamma_distribution=ELLIPTIC,
        )

        set_va!(body_aero, [10.0, 0.0, 0.0])
        solve!(solver, body_aero)
        gamma_low = copy(solver.sol.gamma_distribution)

        set_va!(body_aero, [10.0, 0.0, 5.0])
        solve!(solver, body_aero)
        gamma_high = copy(solver.sol.gamma_distribution)

        @test !isapprox(gamma_low, gamma_high; atol=1e-6, rtol=1e-4)
        @test norm(gamma_high .- gamma_low) >
              1e-3 * max(norm(gamma_low), norm(gamma_high))
    finally
        rm(settings_file; force=true)
    end
end

calc_forces_allocs(solver, body_aero) =
    (calc_forces!(solver, body_aero); @allocated calc_forces!(solver, body_aero))

@testset "calc_forces! is zero-alloc" begin
    settings_file = create_temp_wing_settings(
        "solver", "solver_test_wing.yaml";
        alpha=5.0, beta=0.0, wind_speed=10.0,
    )
    try
        settings = VSMSettings(settings_file)
        wing = Wing(settings)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])
        solver = Solver(body_aero, settings)
        set_va!(body_aero, [10.0, 0.0, 0.0])
        solve!(solver, body_aero)

        @test calc_forces_allocs(solver, body_aero) == 0
    finally
        rm(settings_file; force=true)
    end
end

@testset "Spanwise Laplacian tip closures" begin
    # Interior three-point stencil plus the Eq. 15 tip closures.
    n = 5
    laplacian = zeros(n, n)
    VortexStepMethod.build_spanwise_laplacian!(laplacian, n)

    @test laplacian[3, :] == [0.0, 1.0, -2.0, 1.0, 0.0]
    @test laplacian[1, :] == [-4.0, 4.0/3.0, 0.0, 0.0, 0.0]
    @test laplacian[end, :] == [0.0, 0.0, 0.0, 4.0/3.0, -4.0]
    @test sum(laplacian[3, :]) == 0.0

    # Degenerate spans stay all-zero (no regularization possible).
    small = zeros(2, 2)
    VortexStepMethod.build_spanwise_laplacian!(small, 2)
    @test all(small .== 0.0)
end

@testset "Artificial viscosity stabilizes post-stall" begin
    # Paper case b: AR=10, Cl'=-pi; plain iteration diverges, implicit converges.
    AR, clp, N, V, c = 10.0, -pi, 50, 1.0, 2.0
    b = AR * c
    dz = b / N
    alpha_g = deg2rad(5.0)
    idx = collect(0:N-1)
    induction = [1.0 / (pi * dz) / (4.0 * (j - i)^2 - 1.0) for i in idx, j in idx]

    residual(gamma) =
        0.5 .* V .* c .* (clp .* ((alpha_g .+ (induction * gamma) ./ V) .- alpha_g) .+ 1.0)

    laplacian = zeros(N, N)
    VortexStepMethod.build_spanwise_laplacian!(laplacian, N)
    mu = -0.035 * N^2 / AR * clp                  # Eq. (16), > 0 since clp < 0
    system = Matrix(1.0I, N, N) .- mu .* laplacian

    function iterate_loop(use_viscosity, omega; max_iter=20000)
        gamma = zeros(N)
        for _ in 1:max_iter
            prev = gamma
            nxt = use_viscosity ? (system \ residual(prev)) : residual(prev)
            gamma = (1 - omega) .* prev .+ omega .* nxt
            ref = max(maximum(abs.(gamma)), 1e-4)
            maximum(abs.(gamma .- prev)) / ref < 1e-6 && return true, gamma
        end
        return false, gamma
    end

    base_converged, _ = iterate_loop(false, 0.1)
    av_converged, gamma_av = iterate_loop(true, 0.5)

    @test !base_converged    # base fixed point diverges (sawtooth)
    @test av_converged
    peak = maximum(abs.(gamma_av))
    sawtooth = sum(abs.(gamma_av[1:end-2] .- 2 .* gamma_av[2:end-1] .+
        gamma_av[3:end])) / (length(gamma_av) - 2) / peak
    @test sawtooth < 0.02
end

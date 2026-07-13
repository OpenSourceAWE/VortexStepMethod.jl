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

# A deep-alpha polar is required: the viscosity only fires where the local lift
# slope is negative, and `cl` clamps flat past the tabulated alpha range. The
# default 15-deg range never stalls, so we extend it to 40 deg.
poststall_wing = ram_air_matrix_wing(; n_panels=20, n_sections=4,
    alpha_range=deg2rad.(-5:2:40), delta_range=deg2rad.(-3:3:3))
refine!(poststall_wing)

roughness(v) = sum(abs, @views v[1:end-2] .- 2 .* v[2:end-1] .+ v[3:end])

@testset "apply_artificial_viscosity! smooths post-stall, no-op attached" begin
    body_aero = BodyAerodynamics([poststall_wing])
    panels = body_aero.panels
    n = length(panels)

    laplacian = zeros(n, n)
    VortexStepMethod.build_spanwise_laplacian!(laplacian, n)
    viscosity_matrix = zeros(n, n)
    lift_slope = zeros(n)
    mu_array = zeros(n)
    gamma_target = zeros(n)
    planform_area = sum(p.width * p.chord for p in panels)

    spiky() = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    attached = fill(deg2rad(2.0), n)
    post_stall = fill(deg2rad(22.0), n)

    # Attached flow: the local slope is positive everywhere, so mu stays zero,
    # the solve is skipped, and gamma is returned untouched.
    gamma_attached = spiky()
    fired_attached = VortexStepMethod.apply_artificial_viscosity!(gamma_attached,
        panels, attached, laplacian, viscosity_matrix, lift_slope, mu_array,
        gamma_target, planform_area, 0.035)
    @test !fired_attached
    @test gamma_attached == spiky()

    # Post-stall: the solve fires and smooths the sawtooth.
    gamma_stalled = spiky()
    rough_before = roughness(gamma_stalled)
    fired_stalled = VortexStepMethod.apply_artificial_viscosity!(gamma_stalled,
        panels, post_stall, laplacian, viscosity_matrix, lift_slope, mu_array,
        gamma_target, planform_area, 0.035)
    @test fired_stalled
    @test roughness(gamma_stalled) < rough_before

    # The attached (hot) path must not allocate.
    gamma_alloc = spiky()
    VortexStepMethod.apply_artificial_viscosity!(gamma_alloc, panels, attached,
        laplacian, viscosity_matrix, lift_slope, mu_array, gamma_target,
        planform_area, 0.035)
    allocs = @allocated VortexStepMethod.apply_artificial_viscosity!(gamma_alloc,
        panels, attached, laplacian, viscosity_matrix, lift_slope, mu_array,
        gamma_target, planform_area, 0.035)
    @test allocs == 0
end

@testset "solve! artificial viscosity: attached no-op, post-stall finite" begin
    body_aero = BodyAerodynamics([poststall_wing])
    solver_off = Solver(body_aero; solver_type=LOOP, aerodynamic_model_type=VSM,
        is_with_artificial_viscosity=false)
    solver_on = Solver(body_aero; solver_type=LOOP, aerodynamic_model_type=VSM,
        is_with_artificial_viscosity=true)

    # Attached flow: viscosity never fires, so results are bit-identical.
    set_va!(body_aero, [10.0, 0.0, 0.0])
    gamma_off = copy(solve!(solver_off, body_aero).gamma_distribution)
    set_va!(body_aero, [10.0, 0.0, 0.0])
    gamma_on = copy(solve!(solver_on, body_aero).gamma_distribution)
    @test gamma_on == gamma_off

    # High angle of attack: the viscosity path runs and stays finite.
    set_va!(body_aero, [10.0 * cosd(25), 0.0, 10.0 * sind(25)])
    sol = solve!(solver_on, body_aero)
    @test all(isfinite, sol.gamma_distribution)
end

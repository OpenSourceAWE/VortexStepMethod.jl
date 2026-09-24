using VortexStepMethod
using VortexStepMethod.AirfoilAero: lei_poly_coeffs
using ForwardDiff
using LinearAlgebra
using Test
if !@isdefined(test_data_path)
    include("../test_data_utils.jl")
end

@testset "Solver Constructor Tests" begin
    @testset "Solver Constructor with VSMSettings" begin
        # Use module-specific test data files
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml";
            alpha=5.0, beta=0.0, va=10.0)

        try
            # Test Solver constructor with VSMSettings
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            refine!(wing)
            body_aero = BodyAerodynamics([wing])
            solver = Solver(settings)

            # Verify solver properties match settings
            @test solver.aerodynamic_model_type == VSM
            @test solver.density == 1.225
            @test Solver(settings; density=1.0).density == 1.0

            # Test that the solver can solve
            va_vec = [10.0, 0.0, 0.0]
            set_va!(body_aero, va_vec)
            sol = solve!(solver, body_aero)
            @test sol isa VSMSolution

            @testset "body_aero constructors warn and match Solver(settings)" begin
                with_settings = r"`Solver\(body_aero, settings\)` is deprecated"
                solver_from_body = @test_logs((:warn, with_settings),
                    Solver(body_aero, settings))
                @test solver_from_body isa typeof(solver)
                @test solver_from_body.density == solver.density
                @test solve!(solver_from_body, body_aero).force ≈ sol.force

                with_kwargs = r"`Solver\(body_aero; kwargs...\)` is deprecated"
                solver_from_body = @test_logs((:warn, with_kwargs),
                    Solver(body_aero; density=1.0))
                @test solver_from_body isa typeof(solver)
                @test solver_from_body.density == 1.0
            end

            @testset "Solver from panel and section counts" begin
                n_sections = wing.n_unrefined_sections
                solver_from_counts = Solver(wing.n_panels, n_sections)
                @test solver_from_counts isa typeof(solver)
                @test solve!(solver_from_counts, body_aero).force ≈ sol.force
            end

            @testset "solve refuses a body_aero sized for another solver" begin
                n_sections = wing.n_unrefined_sections
                for other in (Solver(wing.n_panels + 1, n_sections),
                              Solver(wing.n_panels, n_sections + 1))
                    @test_throws DimensionMismatch solve!(other, body_aero)
                    @test_throws "Solver built for" solve!(other, body_aero)
                    @test_throws DimensionMismatch solve(other, body_aero)
                end
            end
        finally
            # Cleanup
            rm(settings_file; force=true)
        end
    end
end

@testset "NONLIN solve! re-runs across calls" begin
    settings_file = create_temp_wing_settings(
        "solver", "solver_test_wing.yaml";
        alpha=5.0, beta=0.0, va=10.0,
    )
    try
        settings = VSMSettings(settings_file)
        wing = Wing(settings)
        refine!(wing)

        body_aero = BodyAerodynamics([wing])
        solver = Solver(
            wing.n_panels, wing.n_unrefined_sections;
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

@testset "NONLIN converges past stall, where LOOP already did" begin
    settings_file = create_temp_wing_settings(
        "solver", "solver_test_wing.yaml";
        alpha=5.0, beta=0.0, va=10.0,
    )
    try
        settings = VSMSettings(settings_file)
        wing = Wing(settings)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])
        va_vec = [10.0, 0.0, 5.0]   # 26.6 deg angle of attack, past stall
        nonlin = Solver(wing.n_panels, wing.n_unrefined_sections; solver_type=NONLIN,
            aerodynamic_model_type=VSM, type_initial_gamma_distribution=ELLIPTIC)
        loop = Solver(wing.n_panels, wing.n_unrefined_sections; solver_type=LOOP,
            aerodynamic_model_type=VSM, type_initial_gamma_distribution=ELLIPTIC)

        set_va!(body_aero, va_vec)
        sol_nonlin = solve!(nonlin, body_aero)
        gamma_nonlin = copy(sol_nonlin.gamma_distribution)
        set_va!(body_aero, va_vec)
        sol_loop = solve!(loop, body_aero)

        @test sol_nonlin.solver_status == FEASIBLE
        @test sol_loop.solver_status == FEASIBLE
        @test isapprox(gamma_nonlin, sol_loop.gamma_distribution; rtol=1e-3)
    finally
        rm(settings_file; force=true)
    end
end

"""
    unrelaxed_step(body_aero, gamma)

One unrelaxed fixed-point step `F(gamma)` of the LOOP iteration, so that
`F(gamma) - gamma` is the residual of `gamma`.
"""
function unrelaxed_step(body_aero, gamma)
    # An infinite rtol accepts the single step, so solve_base! skips its retry.
    wing = only(body_aero.wings)
    probe = Solver(wing.n_panels, wing.n_unrefined_sections; solver_type=LOOP,
        aerodynamic_model_type=VSM, relaxation_factor=1.0, max_iterations=1, rtol=Inf)
    VortexStepMethod.solve_base!(probe, body_aero, gamma)
    return copy(probe.lr.gamma_new)
end

@testset "LOOP converges on the residual, not on the relaxed step" begin
    settings_file = create_temp_wing_settings(
        "solver", "solver_test_wing.yaml";
        alpha=5.0, beta=0.0, va=10.0,
    )
    try
        settings = VSMSettings(settings_file)
        wing = Wing(settings)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])
        solver = Solver(wing.n_panels, wing.n_unrefined_sections; solver_type=LOOP,
            aerodynamic_model_type=VSM, type_initial_gamma_distribution=ELLIPTIC)

        # 0 deg, and 26.6 deg past stall
        for va_vec in ([10.0, 0.0, 0.0], [10.0, 0.0, 5.0])
            set_va!(body_aero, va_vec)
            gamma = copy(solve!(solver, body_aero).gamma_distribution)
            @test solver.lr.converged
            residual = maximum(abs, unrelaxed_step(body_aero, gamma) .- gamma)
            @test residual < solver.rtol * maximum(abs, gamma)
        end
    finally
        rm(settings_file; force=true)
    end
end

calc_forces_allocs(solver, body_aero) =
    (calc_forces!(solver, body_aero); @allocated calc_forces!(solver, body_aero))

@testset "calc_forces! is zero-alloc" begin
    settings_file = create_temp_wing_settings(
        "solver", "solver_test_wing.yaml";
        alpha=5.0, beta=0.0, va=10.0,
    )
    try
        settings = VSMSettings(settings_file)
        wing = Wing(settings)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])
        solver = Solver(settings)
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

"""
    flat_plate_wing(; n_panels=20, span=20.0, chord=1.0)

A rectangular wing with a Breukels flat-plate (zero-camber) LEI polar that genuinely
stalls: the local lift slope is positive in attached flow and turns negative past
stall, which is what the artificial-viscosity path keys on. The planar rectangular
planform keeps the stall behaviour predictable across the span, unlike a drooped kite.
"""
function flat_plate_wing(; n_panels=20, span=20.0, chord=1.0)
    wing = Wing(n_panels; spanwise_distribution=LINEAR)
    coeffs = lei_poly_coeffs(0.1, 0.0)
    add_section!(wing, [0.0,  span/2, 0.0], [chord,  span/2, 0.0],
                 LEI_AIRFOIL_BREUKELS, coeffs)
    add_section!(wing, [0.0, -span/2, 0.0], [chord, -span/2, 0.0],
                 LEI_AIRFOIL_BREUKELS, coeffs)
    refine!(wing)
    return wing
end
poststall_wing = flat_plate_wing()

roughness(v) = sum(abs, @views v[1:end-2] .- 2 .* v[2:end-1] .+ v[3:end])

@testset "apply_artificial_viscosity! smooths post-stall, no-op attached" begin
    body_aero = BodyAerodynamics([poststall_wing])
    panels = body_aero.panels
    n = length(panels)

    laplacian = zeros(n, n)
    VortexStepMethod.build_spanwise_laplacian!(laplacian, n)
    viscosity_matrix = zeros(n, n)
    lift_slope = zeros(n)
    mu_dist = zeros(n)
    gamma_target = zeros(n)
    planform_area = sum(p.width * p.chord for p in panels)

    spiky() = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    attached = fill(deg2rad(2.0), n)
    post_stall = fill(deg2rad(16.0), n)

    # Attached flow: the local slope is positive everywhere, so mu stays zero,
    # the solve is skipped, and gamma is returned untouched.
    gamma_attached = spiky()
    fired_attached = VortexStepMethod.apply_artificial_viscosity!(gamma_attached,
        panels, attached, laplacian, viscosity_matrix, lift_slope, mu_dist,
        gamma_target, planform_area, 0.035)
    @test !fired_attached
    @test gamma_attached == spiky()

    # Post-stall: the solve fires and smooths the sawtooth.
    gamma_stalled = spiky()
    rough_before = roughness(gamma_stalled)
    fired_stalled = VortexStepMethod.apply_artificial_viscosity!(gamma_stalled,
        panels, post_stall, laplacian, viscosity_matrix, lift_slope, mu_dist,
        gamma_target, planform_area, 0.035)
    @test fired_stalled
    @test roughness(gamma_stalled) < rough_before

    # The attached (hot) path must not allocate.
    gamma_alloc = spiky()
    VortexStepMethod.apply_artificial_viscosity!(gamma_alloc, panels, attached,
        laplacian, viscosity_matrix, lift_slope, mu_dist, gamma_target,
        planform_area, 0.035)
    allocs = @allocated VortexStepMethod.apply_artificial_viscosity!(gamma_alloc,
        panels, attached, laplacian, viscosity_matrix, lift_slope, mu_dist,
        gamma_target, planform_area, 0.035)
    @test allocs == 0
end

@testset "solve! artificial viscosity: attached no-op, post-stall finite" begin
    body_aero = BodyAerodynamics([poststall_wing])
    solver_off = Solver(poststall_wing.n_panels, poststall_wing.n_unrefined_sections;
        solver_type=LOOP, aerodynamic_model_type=VSM, is_with_artificial_viscosity=false)
    solver_on = Solver(poststall_wing.n_panels, poststall_wing.n_unrefined_sections;
        solver_type=LOOP, aerodynamic_model_type=VSM, is_with_artificial_viscosity=true)

    # Attached flow: viscosity never fires, so results are bit-identical.
    set_va!(body_aero, [10.0, 0.0, 0.0])
    gamma_off = copy(solve!(solver_off, body_aero).gamma_distribution)
    set_va!(body_aero, [10.0, 0.0, 0.0])
    gamma_on = copy(solve!(solver_on, body_aero).gamma_distribution)
    @test gamma_on == gamma_off

    # High angle of attack: the viscosity path runs and stays finite.
    set_va!(body_aero, [10.0 * cosd(20), 0.0, 10.0 * sind(20)])
    sol = solve!(solver_on, body_aero)
    @test all(isfinite, sol.gamma_distribution)
end

@testset "solve! reports a solve that missed the tolerances" begin
    body_aero = BodyAerodynamics([poststall_wing])
    solver = Solver(poststall_wing.n_panels, poststall_wing.n_unrefined_sections;
        solver_type=LOOP, aerodynamic_model_type=VSM, max_iterations=1)
    set_va!(body_aero, [10.0, 0.0, 0.0])

    sol = solve!(solver, body_aero)
    @test !solver.lr.converged
    @test sol.solver_status == FAILURE

    @test_throws SolveFailure solve!(solver, body_aero; throw_on_fail=true)
    @test_throws "did not converge in 1 iterations" solve!(solver, body_aero;
        throw_on_fail=true)

    converged = Solver(poststall_wing.n_panels, poststall_wing.n_unrefined_sections;
        solver_type=LOOP, aerodynamic_model_type=VSM)
    @test solve!(converged, body_aero; throw_on_fail=true) isa VSMSolution
end

@testset "finite_full sees a Dual's partials, not just its value" begin
    @test VortexStepMethod.finite_full(1.0)
    @test !VortexStepMethod.finite_full(NaN)
    @test VortexStepMethod.finite_full(ForwardDiff.Dual(1.0, 2.0))
    @test !VortexStepMethod.finite_full(ForwardDiff.Dual(1.0, Inf))
end

@testset "alpha is the crossflow angle on a swept, tapered wing" begin
    polar = (deg2rad.([-20.0, 0.0, 20.0]), [-1.8, 0.2, 2.2], fill(0.02, 3), fill(-0.05, 3))
    wing = Wing(8; spanwise_distribution=LINEAR)
    add_section!(wing, [1.5, 4.0, 0.3], [2.2, 4.0, 0.1], POLAR_VECTORS, polar)
    add_section!(wing, [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], POLAR_VECTORS, polar)
    add_section!(wing, [1.5, -4.0, 0.3], [2.2, -4.0, 0.1], POLAR_VECTORS, polar)
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    va = [15.0, 1.0, 1.5]
    set_va!(body_aero, va)

    function crossflow_alpha(panel, velocity)
        chord_in_plane = panel.x_airf .- dot(panel.x_airf, panel.y_airf) .* panel.y_airf
        return atan(dot(velocity, panel.z_airf), dot(velocity, normalize(chord_in_plane)))
    end
    induced(AIC, gamma, i) = [dot(AIC[i, :, k], gamma) for k in 1:3]

    for correct_aoa in (false, true)
        solver = Solver(length(body_aero.panels), 3; correct_aoa)
        sol = solve!(solver, body_aero)
        @test sol.solver_status == FEASIBLE
        gamma = solver.lr.gamma_new
        for (i, panel) in enumerate(body_aero.panels)
            @test solver.lr.alpha_dist[i] ≈
                  crossflow_alpha(panel, va .+ induced(body_aero.AIC, gamma, i)) atol = 1e-6
            @test sol.alpha_geometric_dist[i] ≈ crossflow_alpha(panel, va) atol = 1e-10
            correct_aoa || continue
            @test sol.alpha_dist[i] ≈ crossflow_alpha(panel,
                va .+ induced(body_aero.AIC_aero_center, gamma, i)) atol = 1e-8
        end
    end
end

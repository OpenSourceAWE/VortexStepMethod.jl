using VortexStepMethod
using VortexStepMethod: attached_trailed_loads, calc_forces!, solver_kwargs
using DifferentiationInterface
using LinearAlgebra
using Test

"""
    curved_elliptic_wing(; sweep, anhedral, aspect_ratio=4.0, n_sections=61)

Elliptic wing of unit semi-span whose quarter-chord line runs back by `sweep` and down
by `anhedral` at the tips, parabolic in span, with each section perpendicular to that
line. Returns the wing and its aspect ratio.
"""
function curved_elliptic_wing(; sweep, anhedral, aspect_ratio=4.0, n_sections=61)
    max_chord = 8 / (π * aspect_ratio)
    wing = Wing(n_sections - 1; spanwise_distribution=UNCHANGED)
    for theta in range(0, π, length=n_sections)
        eta = cos(theta) * (1 - 1e-5)
        quarter_chord = [sweep * eta^2, eta, -anhedral * eta^2]
        tangent = normalize([2sweep * eta, 1.0, -2anhedral * eta])
        chordwise = normalize([1.0, 0.0, 0.0] .- tangent[1] .* tangent)
        chord = max_chord * sqrt(1 - eta^2)
        add_section!(wing, quarter_chord .- 0.25chord .* chordwise,
            quarter_chord .+ 0.75chord .* chordwise, INVISCID)
    end
    refine!(wing)
    return wing, aspect_ratio
end

"""
    lift_and_oswald(wing, aspect_ratio, model, is_with_attached_trailed_force)

Wing `CL` and blade-force Oswald factor `CL²/(π AR CD)` at α = 10°.
"""
function lift_and_oswald(wing, aspect_ratio, model, is_with_attached_trailed_force)
    body_aero = BodyAerodynamics([wing])
    solver = Solver(wing.n_panels, wing.n_unrefined_sections; use_gamma_prev=false,
        aerodynamic_model_type=model, is_with_attached_trailed_force)
    set_va!(body_aero, 20.0 .* [cosd(10), 0.0, sind(10)])
    results = solve(solver, body_aero)
    return results["cl"], results["cl"]^2 / (π * aspect_ratio * results["cd"])
end

@testset "Attached trailed vortex force" begin
    @testset "defaults to off and follows the solver settings" begin
        @test VortexStepMethod.SolverSettings().is_with_attached_trailed_force == false
        @test Solver(4, 2).is_with_attached_trailed_force == false
        settings = VortexStepMethod.SolverSettings(is_with_attached_trailed_force=true)
        @test solver_kwargs(settings).is_with_attached_trailed_force
    end

    # Gaunaa, Li & Pirrung, TORQUE 2026, §5.1-5.2: omitting the force overpredicts CL
    # by about 2 % on a swept, anhedral wing and barely changes a straight one.
    @testset "lowers lift by about 2 % on a swept, anhedral wing ($model)" for
            model in (VSM, LLT)
        straight, aspect_ratio = curved_elliptic_wing(; sweep=0.0, anhedral=0.0)
        CL_off, _ = lift_and_oswald(straight, aspect_ratio, model, false)
        CL_on, _ = lift_and_oswald(straight, aspect_ratio, model, true)
        @test abs(CL_on / CL_off - 1) < 0.005

        curved, aspect_ratio = curved_elliptic_wing(; sweep=0.1, anhedral=0.25)
        CL_off, oswald_off = lift_and_oswald(curved, aspect_ratio, model, false)
        CL_on, oswald_on = lift_and_oswald(curved, aspect_ratio, model, true)
        @test -0.03 < CL_on / CL_off - 1 < -0.015
        @test oswald_on > oswald_off
    end

    wing, _ = curved_elliptic_wing(; sweep=0.1, anhedral=0.25, n_sections=21)
    va_vec = 20.0 .* [cosd(10), 0.0, sind(10)]

    @testset "adds attached_trailed_loads to every panel in solve! and solve" begin
        body_aero = BodyAerodynamics([wing])
        set_va!(body_aero, va_vec)
        reference_point = [0.1, 0.0, -0.2]
        solver_off = Solver(wing.n_panels, wing.n_unrefined_sections; use_gamma_prev=false,
            reference_point)
        solver_on = Solver(wing.n_panels, wing.n_unrefined_sections; use_gamma_prev=false,
            reference_point, is_with_attached_trailed_force=true)
        force_off = copy(solve!(solver_off, body_aero).f_body_3D)
        moment_off = copy(solver_off.sol.m_body_3D)
        solve!(solver_on, body_aero)
        for i in eachindex(body_aero.panels)
            attached = attached_trailed_loads(body_aero, i, solver_on.lr.gamma_new,
                solver_on.density, solver_on.core_radius_fraction, reference_point)
            @test solver_on.sol.f_body_3D[:, i] ≈ force_off[:, i] .+ attached.force
            @test solver_on.sol.m_body_3D[:, i] ≈ moment_off[:, i] .+ attached.moment
        end
        results = solve(solver_on, body_aero)
        @test results["F_distribution"] ≈ solver_on.sol.f_body_3D
        @test results["M_distribution"] ≈ solver_on.sol.m_body_3D
        @test (@allocated calc_forces!(solver_on, body_aero)) == 0
    end

    @testset "ForwardDiff matches FiniteDiff with the force on" begin
        body_aero = BodyAerodynamics([wing])
        set_va!(body_aero, va_vec)
        solver = Solver(wing.n_panels, wing.n_unrefined_sections; use_gamma_prev=false,
            type_initial_gamma_distribution=ELLIPTIC, is_with_attached_trailed_force=true)
        jac_fwd, _ = VortexStepMethod.linearize(solver, body_aero, va_vec;
            theta_idxs=nothing, va_vec_idxs=1:3, aero_coeffs=true,
            backend=AutoForwardDiff())
        jac_fd, _ = VortexStepMethod.linearize(solver, body_aero, va_vec;
            theta_idxs=nothing, va_vec_idxs=1:3, aero_coeffs=true,
            backend=AutoFiniteDiff(absstep=1e-5, relstep=1e-5))
        @test maximum(abs.(jac_fwd .- jac_fd)) / maximum(abs, jac_fwd) < 1e-3
    end
end

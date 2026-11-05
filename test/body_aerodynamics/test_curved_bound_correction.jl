using LinearAlgebra
using Test
using VortexStepMethod
using VortexStepMethod: curved_bound_correction!, quarter_chord_AIC!, reinit!,
    thin_airfoil_strips

"""Constant-chord wing of `n_panels` panels from the section points of `le_points`."""
function wing_through(le_points, n_panels)
    wing = Wing(n_panels; spanwise_distribution=LINEAR)
    for le in le_points
        add_section!(wing, le, le .+ [1.0, 0.0, 0.0], INVISCID)
    end
    refine!(wing)
    return wing
end

arc_points(radius=4.0, half_angle=55.0, n_sections=21) = [
    [0.0, radius * sind(angle), radius * cosd(angle)]
    for angle in range(half_angle, -half_angle; length=n_sections)]

@testset "Curved bound vortex correction" begin
    @testset "the strips hold a thin-airfoil load centred on the quarter chord" begin
        stations, weights = thin_airfoil_strips()
        @test length(stations) == 24
        @test sum(weights) ≈ 1
        @test sum(stations .* weights) ≈ 0.25 atol = 2e-3
        @test all(0 .< stations .< 1)
    end

    n_panels = 20
    va_dist = fill(15.0, n_panels)

    @testset "it vanishes on a straight wing" begin
        straight = wing_through([[0.0, 4.0, 0.0], [0.0, -4.0, 0.0]], n_panels)
        body_aero = BodyAerodynamics([straight])
        correction = zeros(n_panels, n_panels, 3)
        curved_bound_correction!(correction, body_aero, 0.05, va_dist)
        @test maximum(abs, correction) < 1e-10
    end

    @testset "on an arc it is kept per shape and rebuilt after a shape change" begin
        wing = wing_through(arc_points(), n_panels)
        body_aero = BodyAerodynamics([wing])
        solver = Solver(n_panels, wing.n_unrefined_sections)
        set_va!(body_aero, [15.0, 0.0, 2.0])
        solve!(solver, body_aero)
        @test maximum(abs, body_aero.AIC_curved_correction) > 1e-3
        @test body_aero.curved_correction_core_radius == solver.core_radius_fraction

        body_aero.AIC_curved_correction .= 0.0
        solve!(solver, body_aero)
        @test iszero(body_aero.AIC_curved_correction)

        reinit!(body_aero)
        @test isnan(body_aero.curved_correction_core_radius)
        set_va!(body_aero, [15.0, 0.0, 2.0])
        solve!(solver, body_aero)
        @test maximum(abs, body_aero.AIC_curved_correction) > 1e-3
    end

    @testset "a new inflow on the same shape reuses it, as a fresh body agrees" begin
        wing = wing_through(arc_points(), n_panels)
        body_aero = BodyAerodynamics([wing])
        solver = Solver(n_panels, wing.n_unrefined_sections)
        set_va!(body_aero, [10.0, 0.0, 1.0])
        solve!(solver, body_aero, nothing)
        kept = copy(body_aero.AIC_curved_correction)

        second_inflow = [20.0, 2.0, 2.0]
        body_aero.AIC_curved_correction .= 2 .* kept
        set_va!(body_aero, second_inflow)
        solve!(solver, body_aero, nothing)
        @test body_aero.AIC_curved_correction == 2 .* kept

        body_aero.AIC_curved_correction .= kept
        sol = solve!(solver, body_aero, nothing)
        fresh_wing = wing_through(arc_points(), n_panels)
        fresh_body = BodyAerodynamics([fresh_wing])
        set_va!(fresh_body, second_inflow)
        fresh_solver = Solver(n_panels, fresh_wing.n_unrefined_sections)
        fresh = solve!(fresh_solver, fresh_body, nothing)
        @test sol.cl ≈ fresh.cl rtol = 1e-12
        @test sol.cd ≈ fresh.cd rtol = 1e-12
    end

    @testset "rebuilding it allocates a small setup, not a vector per filament" begin
        arc = wing_through(arc_points(), n_panels)
        body_aero = BodyAerodynamics([arc])
        correction = zeros(n_panels, n_panels, 3)
        rebuild() = @allocated curved_bound_correction!(correction, body_aero, 0.05,
            va_dist)
        rebuild()
        @test rebuild() < 2_000
    end

    @testset "a second wing does not enter another wing's correction" begin
        low = wing_through([point .+ [0.0, 0.0, -6.0] for point in arc_points()], 10)
        high = wing_through(arc_points(), 10)
        body_aero = BodyAerodynamics([low, high])
        correction = zeros(20, 20, 3)
        curved_bound_correction!(correction, body_aero, 0.05, fill(15.0, 20))
        @test iszero(correction[1:10, 11:20, :])
        @test iszero(correction[11:20, 1:10, :])
        @test !iszero(correction[1:10, 1:10, :])
    end
end

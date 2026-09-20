using VortexStepMethod
using VortexStepMethod: coeffs_at_angles
using Test

"""
    trimmable_wing(cm; kwargs...)

A rectangular wing with lift slope 2π, `cl = 0.1` at zero angle of attack and a constant
section `cm`, as a `BodyAerodynamics` and a `Solver` sized for it. Keyword arguments go to
the solver.
"""
function trimmable_wing(cm; kwargs...)
    chord, span = 1.0, 6.0
    alpha_range = deg2rad.(-10.0:5.0:20.0)
    polar = (alpha_range, 2π .* alpha_range .+ 0.1, fill(0.02, length(alpha_range)),
             fill(cm, length(alpha_range)))
    wing = Wing(10)
    add_section!(wing, [0.0, span / 2, 0.0], [chord, span / 2, 0.0], POLAR_VECTORS, polar)
    add_section!(wing, [0.0, -span / 2, 0.0], [chord, -span / 2, 0.0], POLAR_VECTORS,
                 polar)
    refine!(wing)
    return BodyAerodynamics([wing]),
           Solver(wing.n_panels, wing.n_unrefined_sections; kwargs...)
end

@testset "stability_derivatives match central differences of solve!" begin
    body_aero, solver = trimmable_wing(0.05; reference_point=[0.25, 0.5, 0.1],
                                       use_gamma_prev=false)
    alpha, beta, wind_speed, step = deg2rad(4.0), deg2rad(3.0), 20.0, 1e-4

    derivatives = stability_derivatives(solver, body_aero, alpha, beta, wind_speed)
    @test derivatives.converged
    @test derivatives.coeffs ≈
          coeffs_at_angles(solver, body_aero, alpha, beta, wind_speed)

    central_difference(coeffs_plus, coeffs_minus) = (coeffs_plus - coeffs_minus) / 2step
    dalpha = central_difference(
        coeffs_at_angles(solver, body_aero, alpha + step, beta, wind_speed),
        coeffs_at_angles(solver, body_aero, alpha - step, beta, wind_speed))
    dbeta = central_difference(
        coeffs_at_angles(solver, body_aero, alpha, beta + step, wind_speed),
        coeffs_at_angles(solver, body_aero, alpha, beta - step, wind_speed))
    @test !iszero(dbeta)
    @test derivatives.dalpha ≈ dalpha rtol = 1e-4 atol = 1e-6
    @test derivatives.dbeta ≈ dbeta rtol = 1e-4 atol = 1e-6
end

@testset "trim_angle finds where CMy changes sign" begin
    beta, wind_speed = 0.0, 20.0

    @testset "moments about the leading edge: stable trim" begin
        body_aero, solver = trimmable_wing(0.05)
        trims = trim_angle(solver, body_aero, beta, wind_speed)
        @test length(trims) == 1
        trim = only(trims)
        trim_coeffs = coeffs_at_angles(solver, body_aero, trim.alpha, beta, wind_speed)
        @test abs(trim_coeffs[5]) < 1e-5
        @test trim.dCMy_dalpha < 0
        derivatives = stability_derivatives(solver, body_aero, trim.alpha, beta, wind_speed)
        @test trim.dCMy_dalpha ≈ derivatives.dalpha[5]
    end

    @testset "moments about the trailing edge: unstable trim" begin
        body_aero, solver = trimmable_wing(-0.05; reference_point=[1.0, 0.0, 0.0])
        trim = only(trim_angle(solver, body_aero, beta, wind_speed))
        trim_coeffs = coeffs_at_angles(solver, body_aero, trim.alpha, beta, wind_speed)
        @test abs(trim_coeffs[5]) < 1e-5
        @test trim.dCMy_dalpha > 0
    end

    @testset "a NONLIN solver with backend=nothing finds the same trim" begin
        loop_body, loop_solver = trimmable_wing(0.05)
        body_aero, solver = trimmable_wing(0.05; solver_type=NONLIN)
        trim_loop = only(trim_angle(loop_solver, loop_body, beta, wind_speed))
        trim = only(trim_angle(solver, body_aero, beta, wind_speed; backend=nothing))
        @test trim.alpha ≈ trim_loop.alpha atol = 1e-4
        @test trim.dCMy_dalpha ≈ trim_loop.dCMy_dalpha rtol = 1e-4
    end

    @testset "a solve that misses the tolerances throws" begin
        body_aero, solver = trimmable_wing(0.05; max_iterations=2)
        @test_throws SolveFailure trim_angle(solver, body_aero, beta, wind_speed)
    end

    @testset "no sign change in alpha_range: no trim" begin
        body_aero, solver = trimmable_wing(0.05)
        @test isempty(trim_angle(solver, body_aero, beta, wind_speed;
                                 alpha_range=deg2rad.(4:2:12)))
    end
end

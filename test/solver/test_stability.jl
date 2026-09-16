using VortexStepMethod
using VortexStepMethod: apparent_wind
using Test

"""
    trimmable_wing_aero(cm)

A rectangular wing with lift slope 2π, `cl = 0.1` at zero angle of attack and a constant
section `cm`.
"""
function trimmable_wing_aero(cm)
    chord, span = 1.0, 6.0
    alpha_range = deg2rad.(-10.0:5.0:20.0)
    polar = (alpha_range, 2π .* alpha_range .+ 0.1, fill(0.02, length(alpha_range)),
             fill(cm, length(alpha_range)))
    wing = Wing(10)
    add_section!(wing, [0.0, span / 2, 0.0], [chord, span / 2, 0.0], POLAR_VECTORS, polar)
    add_section!(wing, [0.0, -span / 2, 0.0], [chord, -span / 2, 0.0], POLAR_VECTORS,
                 polar)
    refine!(wing)
    return BodyAerodynamics([wing])
end

function coeffs_at(solver, body_aero, alpha, beta, wind_speed)
    set_va!(body_aero, apparent_wind(alpha, beta, wind_speed))
    sol = solve!(solver, body_aero)
    return [sol.force_coeffs; sol.moment_coeffs]
end

@testset "stability_derivatives match central differences of solve!" begin
    body_aero = trimmable_wing_aero(0.05)
    solver = Solver(body_aero; reference_point=[0.25, 0.5, 0.1], use_gamma_prev=false)
    alpha, beta, wind_speed, step = deg2rad(4.0), deg2rad(3.0), 20.0, 1e-4

    derivatives = stability_derivatives(solver, body_aero, alpha, beta, wind_speed)
    @test derivatives.converged
    @test derivatives.coeffs ≈ coeffs_at(solver, body_aero, alpha, beta, wind_speed)

    central_difference(coeffs_plus, coeffs_minus) = (coeffs_plus - coeffs_minus) / 2step
    dalpha = central_difference(
        coeffs_at(solver, body_aero, alpha + step, beta, wind_speed),
        coeffs_at(solver, body_aero, alpha - step, beta, wind_speed))
    dbeta = central_difference(
        coeffs_at(solver, body_aero, alpha, beta + step, wind_speed),
        coeffs_at(solver, body_aero, alpha, beta - step, wind_speed))
    @test !iszero(dbeta)
    @test derivatives.dalpha ≈ dalpha rtol = 1e-4 atol = 1e-6
    @test derivatives.dbeta ≈ dbeta rtol = 1e-4 atol = 1e-6
end

@testset "rate derivatives match central differences of solve! about reference_point" begin
    body_aero = trimmable_wing_aero(0.05)
    reference_point = [0.25, 0.5, 0.1]
    solver = Solver(body_aero; reference_point, use_gamma_prev=false)
    alpha, beta, wind_speed, step = deg2rad(4.0), deg2rad(3.0), 20.0, 1e-4
    va = apparent_wind(alpha, beta, wind_speed)
    omega = [0.1, -0.05, 0.08]
    set_va!(body_aero, va, omega)

    derivatives = stability_derivatives(solver, body_aero, alpha, beta, wind_speed)
    @test derivatives.converged
    @test body_aero.reference_point == reference_point

    function coeffs_at_rate(rate)
        set_va!(body_aero, va, rate; reference_point)
        sol = solve!(solver, body_aero)
        return [sol.force_coeffs; sol.moment_coeffs]
    end
    rate_scales = 2wind_speed ./ [body_aero.wings[1].span, body_aero.c_ref,
                                  body_aero.wings[1].span]
    rate_derivatives = map(1:3) do axis
        rate_step = step .* (1:3 .== axis)
        (coeffs_at_rate(omega + rate_step) - coeffs_at_rate(omega - rate_step)) /
            2step * rate_scales[axis]
    end
    @test all(!iszero, rate_derivatives)
    @test derivatives.dp ≈ rate_derivatives[1] rtol = 1e-4 atol = 1e-6
    @test derivatives.dq ≈ rate_derivatives[2] rtol = 1e-4 atol = 1e-6
    @test derivatives.dr ≈ rate_derivatives[3] rtol = 1e-4 atol = 1e-6
end

@testset "trim_angle finds where CMy changes sign" begin
    beta, wind_speed = 0.0, 20.0

    @testset "moments about the leading edge: stable trim" begin
        body_aero = trimmable_wing_aero(0.05)
        solver = Solver(body_aero)
        trims = trim_angle(solver, body_aero, beta, wind_speed)
        @test length(trims) == 1
        trim = only(trims)
        @test abs(coeffs_at(solver, body_aero, trim.alpha, beta, wind_speed)[5]) < 1e-5
        @test trim.dCMy_dalpha < 0
        derivatives = stability_derivatives(solver, body_aero, trim.alpha, beta, wind_speed)
        @test trim.dCMy_dalpha ≈ derivatives.dalpha[5]
    end

    @testset "moments about the trailing edge: unstable trim" begin
        body_aero = trimmable_wing_aero(-0.05)
        solver = Solver(body_aero; reference_point=[1.0, 0.0, 0.0])
        trim = only(trim_angle(solver, body_aero, beta, wind_speed))
        @test abs(coeffs_at(solver, body_aero, trim.alpha, beta, wind_speed)[5]) < 1e-5
        @test trim.dCMy_dalpha > 0
    end

    @testset "pitching body: trim turning about the reference point" begin
        body_aero = trimmable_wing_aero(-0.05)
        reference_point = [1.0, 0.0, 0.0]
        solver = Solver(body_aero; reference_point)
        omega = [0.0, 0.5, 0.0]
        set_va!(body_aero, apparent_wind(0.0, beta, wind_speed), omega)
        trim = only(trim_angle(solver, body_aero, beta, wind_speed))
        set_va!(body_aero, apparent_wind(trim.alpha, beta, wind_speed), omega;
                reference_point)
        cmy = solve!(solver, body_aero).moment_coeffs[2]
        @test abs(cmy) < 1e-5 * abs(trim.dCMy_dalpha)
    end

    @testset "no sign change in alpha_range: no trim" begin
        body_aero = trimmable_wing_aero(0.05)
        solver = Solver(body_aero)
        @test isempty(trim_angle(solver, body_aero, beta, wind_speed;
                                 alpha_range=deg2rad.(4:2:12)))
    end
end

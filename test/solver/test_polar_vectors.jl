using VortexStepMethod
using LinearAlgebra
using Test

"""
    linear_polar(slope, cd, cm)

POLAR_VECTORS data with `cl = slope * alpha` and constant `cd`, `cm`, for alpha
from -20° to 30°.
"""
function linear_polar(slope, cd, cm)
    alphas = collect(deg2rad.(-20.0:1.0:30.0))
    return (alphas, slope .* alphas, fill(cd, length(alphas)), fill(cm, length(alphas)))
end

"""
    solve_three_panel_wing(aero_model, polar_pos_y=nothing, polar_neg_y=polar_pos_y)

Solve a flat rectangular wing of three 1 m × 1 m panels at 5° incidence, whose two
sections on the +y side carry `polar_pos_y` and whose two on the -y side carry
`polar_neg_y`.
"""
function solve_three_panel_wing(aero_model, polar_pos_y=nothing, polar_neg_y=polar_pos_y)
    wing = Wing(3)
    for (y, polar) in ((1.5, polar_pos_y), (0.5, polar_pos_y),
                       (-0.5, polar_neg_y), (-1.5, polar_neg_y))
        add_section!(wing, [0.0, y, 0.0], [1.0, y, 0.0], aero_model, polar)
    end
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, 15.0 .* [cosd(5), 0.0, sind(5)])
    solver = Solver(wing.n_panels, wing.n_unrefined_sections; rtol=1e-10)
    return solve!(solver, body_aero)
end

@testset "POLAR_VECTORS wing" begin
    @testset "each panel solves on the average of its two sections' polars" begin
        sol = solve_three_panel_wing(POLAR_VECTORS,
            linear_polar(4π, 0.1, -0.02), linear_polar(2π, 0.05, -0.01))
        @test sol.solver_status == FEASIBLE
        @test sol.cl_dist ≈ [4π, 3π, 2π] .* sol.alpha_dist
        @test sol.cd_dist ≈ [0.1, 0.075, 0.05]
        @test sol.cm_dist ≈ [-0.02, -0.015, -0.01]
        @testset "the stronger +y side rolls the wing positive about x" begin
            @test sol.moment[1] > 0
        end
    end

    @testset "a 2π·alpha polar matches INVISCID plus its profile drag" begin
        polar = linear_polar(2π, 0.05, 0.0)
        sol_polar = solve_three_panel_wing(POLAR_VECTORS, polar)
        sol_inviscid = solve_three_panel_wing(INVISCID)
        @test sol_polar.gamma_distribution ≈ sol_inviscid.gamma_distribution rtol=1e-12
        @test sol_polar.lift_dist ≈ sol_inviscid.lift_dist rtol=1e-12
        @test sol_polar.cd_dist ≈ fill(0.05, 3)
        profile_drag = eachcol(sol_polar.f_body_3D .- sol_inviscid.f_body_3D)
        @test norm.(profile_drag) ≈ sol_polar.drag_dist .* sol_polar.width_dist
        @test abs(sol_polar.moment[1]) < 1e-10
    end
end

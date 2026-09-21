using VortexStepMethod
using LinearAlgebra
using Test

"""
    rectangular_wing(rotation, offset)

A rectangular wing spanning `y`, with every point turned by `rotation` and moved by
`offset`, and its spanwise direction turned with it.
"""
function rectangular_wing(rotation, offset)
    chord, span = 1.5, 6.0
    alpha_range = deg2rad.([-10.0, 0.0, 10.0])
    polar = (alpha_range, [-0.6, 0.4, 1.4], fill(0.02, 3), fill(-0.08, 3))
    wing = Wing(10; spanwise_direction=rotation * [0.0, 1.0, 0.0])
    for y in (span / 2, -span / 2)
        add_section!(wing, rotation * [0.0, y, 0.0] + offset,
            rotation * [chord, y, 0.0] + offset, POLAR_VECTORS, polar)
    end
    refine!(wing)
    return wing
end

@testset "each wing's panels take their own wing's spanwise direction" begin
    # 90° about x leaves the inflow along x unchanged; the second wing is far enough
    # away that the two wings do not induce on each other
    rotation = [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0]
    va = [20.0, 0.0, 0.0]
    solo_aero = BodyAerodynamics([rectangular_wing(I(3), zeros(3))])
    pair_aero = BodyAerodynamics([rectangular_wing(I(3), zeros(3)),
        rectangular_wing(rotation, [0.0, 1000.0, 0.0])])
    set_va!(solo_aero, va)
    set_va!(pair_aero, va)
    n_panels = length(solo_aero.panels)
    rotated = n_panels .+ (1:n_panels)

    @testset "solve!" begin
        solo, pair = map((solo_aero, pair_aero)) do body_aero
            solver = Solver(length(body_aero.panels), 2length(body_aero.wings))
            solve!(solver, body_aero)
        end
        @test pair.solver_status == FEASIBLE
        @test pair.f_body_3D[:, 1:n_panels] ≈ solo.f_body_3D rtol = 1e-5
        @test pair.f_body_3D[:, rotated] ≈ rotation * solo.f_body_3D rtol = 1e-5
        @test pair.lift_dist[rotated] ≈ solo.lift_dist rtol = 1e-5
        @test pair.drag_dist[rotated] ≈ solo.drag_dist rtol = 1e-5
    end

    @testset "solve" begin
        solo, pair = map((solo_aero, pair_aero)) do body_aero
            solver = Solver(length(body_aero.panels), 2length(body_aero.wings))
            solve(solver, body_aero)
        end
        @test pair["F_distribution"][:, rotated] ≈ rotation * solo["F_distribution"] rtol = 1e-5
        for key in ("cl_distribution", "cd_distribution", "cs_distribution")
            @test pair[key][rotated] ≈ solo[key] rtol = 1e-5 atol = 1e-8
        end
    end
end

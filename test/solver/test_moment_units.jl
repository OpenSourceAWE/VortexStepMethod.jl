using VortexStepMethod
using LinearAlgebra
using Test

"""
    scaled_wing_aero(scale)

A rectangular wing with a constant nonzero `cm`, every length multiplied by `scale`,
at an inflow that does not depend on it.
"""
function scaled_wing_aero(scale)
    chord, span = 1.5scale, 6.0scale
    alpha_range = deg2rad.([-10.0, 0.0, 10.0])
    polar = (alpha_range, [-0.6, 0.4, 1.4], fill(0.02, 3), fill(-0.08, 3))
    wing = Wing(10)
    add_section!(wing, [0.0, span / 2, 0.0], [chord, span / 2, 0.0], POLAR_VECTORS, polar)
    add_section!(wing, [0.0, -span / 2, 0.0], [chord, -span / 2, 0.0], POLAR_VECTORS,
                 polar)
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, [20.0, 0.0, 2.0])
    return body_aero
end

@testset "moments carry N·m: lengths scaled by k scale them by k³" begin
    # at fixed inflow and coefficients, N scales as k², N·m as k³, and N·m per
    # unit span as k²; a moment short of one chord factor scales as k² instead
    k = 2.0
    reference_point(scale) = scale .* [-0.4, 0.3, 0.2]

    @testset "solve!" begin
        small, large = map((1.0, k)) do scale
            body_aero = scaled_wing_aero(scale)
            wing = only(body_aero.wings)
            solver = Solver(wing.n_panels, wing.n_unrefined_sections;
                reference_point=reference_point(scale))
            solve!(solver, body_aero)
        end
        @test small.solver_status == large.solver_status == FEASIBLE
        @test all(!iszero, small.panel_moment_dist)
        @test large.force ≈ k^2 .* small.force rtol = 1e-6
        @test large.moment ≈ k^3 .* small.moment rtol = 1e-6
        @test large.m_body_3D ≈ k^3 .* small.m_body_3D rtol = 1e-6
        @test large.moment_dist ≈ k^3 .* small.moment_dist rtol = 1e-6
        @test large.moment_unrefined_dist ≈ k^3 .* small.moment_unrefined_dist rtol = 1e-6
        @test large.panel_moment_dist ≈ k^2 .* small.panel_moment_dist rtol = 1e-6
        @test large.moment_coeffs ≈ small.moment_coeffs rtol = 1e-6
        @test large.moment_coeff_dist ≈ small.moment_coeff_dist rtol = 1e-6
    end
end

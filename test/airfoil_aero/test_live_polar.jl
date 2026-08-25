using Test
using LinearAlgebra
using Statistics
using VortexStepMethod
using VortexStepMethod.AirfoilAero
using VortexStepMethod: Panel, calculate_cl, calculate_cd, calculate_cm,
                        set_taylor_polar!

@testset "TAYLOR panel polar" begin
    panel = Panel{Float64}()
    set_taylor_polar!(panel, deg2rad(5.0), [0.6, 5.0, -2.0], [0.02, 0.1, 1.0],
                      [-0.1, 0.2, 0.0])
    @test panel.aero_model == TAYLOR
    alpha = deg2rad(6.0)
    d = alpha - deg2rad(5.0)
    @test calculate_cl(panel, alpha) ≈ 0.6 + 5.0d - 2.0d^2
    @test calculate_cd(panel, alpha) ≈ 0.02 + 0.1d + 1.0d^2
    @test calculate_cm(panel, alpha) ≈ -0.1 + 0.2d

    # A local expansion has no flap axis: delta must not change the answer.
    @test calculate_cl(panel, alpha, deg2rad(9.0)) == calculate_cl(panel, alpha)

    coeffs = panel.cl_coeffs
    set_taylor_polar!(panel, 0.0, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    @test panel.cl_coeffs === coeffs  # same order refits in place
end

@testset "Kulfan deformation" begin
    basis = KulfanBasis()
    base = KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    camber = @. 0.05 * basis.x^2 * (1 - basis.x)
    deformed = deform_kulfan(basis, base, camber)

    shape = AirfoilAero.class_function(basis.x) .*
            AirfoilAero.bernstein_basis(basis.x, 7)
    for (weights, base_weights) in ((deformed.upper_weights, base.upper_weights),
                                    (deformed.lower_weights, base.lower_weights))
        residual = shape * (weights .- base_weights) .- camber
        @test maximum(abs, residual) < 0.01 * maximum(abs, camber)
    end
    @test deformed.leading_edge_weight == base.leading_edge_weight
    @test deformed.TE_thickness == base.TE_thickness

    # Deforming the surfaces apart changes thickness, not just camber.
    apart = deform_kulfan(basis, base, camber, -camber)
    @test apart.upper_weights ≉ apart.lower_weights .+
          (base.upper_weights .- base.lower_weights)

    @test_throws ArgumentError deform_kulfan(basis, base, camber[1:end-1])
end

@testset "control point deflection" begin
    basis = KulfanBasis()
    sampled = control_point_deflection(basis, [0.0, 0.7, 0.8, 1.0],
                                       [0.0, 0.03, 0.025, 0.0])
    @test sampled[1] ≈ 0.0 atol = 1e-12
    @test sampled[end] ≈ 0.0 atol = 1e-12
    @test maximum(sampled) ≈ 0.03 atol = 1e-3

    # Unsorted input is the same curve as sorted input.
    shuffled = control_point_deflection(basis, [0.8, 0.0, 1.0, 0.7],
                                        [0.025, 0.0, 0.0, 0.03])
    @test shuffled ≈ sampled
    @test_throws ArgumentError control_point_deflection(basis, [0.0, 0.0], [1.0, 2.0])
end

@testset "live polars" begin
    base = KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    n_panels = 3
    panels = [Panel{Float64}() for _ in 1:n_panels]
    live = LivePolars(fill(base, n_panels))
    confidence = refresh_live_polars!(live, panels, deg2rad(6.0), 3e6)
    @test 0.0 < confidence <= 1.0
    @test all(p -> p.aero_model == TAYLOR, panels)

    # The fit tracks a direct NeuralFoil sweep across its own window.
    errors = [abs(calculate_cl(panels[1], deg2rad(6.0 + d)) -
                  neuralfoil_aero(base, 6.0 + d, 3e6).CL[1])
              for d in range(-4, 4, 17)]
    @test maximum(errors) < 0.02

    @test polar_drift(live, fill(deg2rad(6.0), n_panels)) ≈ 0.0 atol = 1e-12
    @test polar_drift(live, fill(deg2rad(10.0), n_panels)) ≈ 1.0

    # Deforming the camber up raises lift at the same angle of attack.
    flat = calculate_cl(panels[1], deg2rad(6.0))
    camber = @. 0.02 * live.basis.x * (1 - live.basis.x)
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6;
                         deflection=fill(camber, n_panels))
    @test calculate_cl(panels[1], deg2rad(6.0)) > flat

    @test_throws ArgumentError refresh_live_polars!(live, panels[1:2], 0.0, 3e6)
    @test_throws ArgumentError LivePolars(fill(base, 2);
        settings=LivePolarSettings(; order=4, n_samples=3))
end

@testset "TAYLOR spanwise blend" begin
    left = (deg2rad(4.0), [0.5, 5.0, 0.0], [0.02, 0.0, 0.0], [-0.1, 0.0, 0.0])
    right = (deg2rad(8.0), [0.9, 5.0, 0.0], [0.06, 0.0, 0.0], [-0.3, 0.0, 0.0])
    blended = VortexStepMethod.calculate_new_aero_data(
        (TAYLOR, TAYLOR), (left, right), 1, 0.25, 0.75)
    @test blended[1] ≈ deg2rad(7.0)
    @test blended[2] ≈ [0.8, 5.0, 0.0]

    section_left = Section([0.0, 1.0, 0.0], [1.0, 1.0, 0.0], TAYLOR, left)
    section_right = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], TAYLOR, right)
    panel = Panel{Float64}()
    VortexStepMethod.init_aero!(panel, section_left, section_right)
    @test panel.alpha_ref ≈ deg2rad(6.0)
    @test panel.cl_coeffs ≈ [0.7, 5.0, 0.0]
end

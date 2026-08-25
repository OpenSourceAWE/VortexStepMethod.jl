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

@testset "TAYLOR window is continued linearly" begin
    window = deg2rad(4.0)
    panel = Panel{Float64}()
    set_taylor_polar!(panel, 0.0, [0.6, 5.0, -20.0], [0.02, 0.0, 0.0],
                      [-0.1, 0.0, 0.0]; window)
    # Inside the window the polynomial is untouched.
    @test calculate_cl(panel, 0.5window) ≈ 0.6 + 5.0 * 0.5window - 20.0 * (0.5window)^2
    # At the edge value and slope match, so the continuation is smooth.
    edge = calculate_cl(panel, window)
    slope = (calculate_cl(panel, window + 1e-7) - edge) / 1e-7
    inner = (edge - calculate_cl(panel, window - 1e-7)) / 1e-7
    @test slope ≈ inner rtol = 1e-4
    # Drag is the exception: its fit has a minimum inside the window, so continuing
    # the lower edge's slope would take it through zero. A negative drag coefficient
    # is an energy source, and on the SK100 it drove the apparent wind from 12 to
    # 21.8 m/s inside one step before the solve died.
    drag = Panel{Float64}()
    set_taylor_polar!(drag, 0.0, [0.5, 0.0, 0.0], [0.02, 0.4, 6.0],
                      [-0.05, 0.0, 0.0]; window)
    @test calculate_cd(drag, -0.5window) < calculate_cd(drag, 0.0)  # still falling
    @test calculate_cd(drag, -window) >= 0.0
    @test calculate_cd(drag, -4window) >= 0.0
    @test calculate_cd(drag, -20window) >= 0.0
    # Above the window the continuation is untouched, drag grows.
    @test calculate_cd(drag, 4window) > calculate_cd(drag, window) > 0.0

    # And it stays linear rather than falling off with the parabola's arm, on both
    # sides: the unbounded quadratic would be 0.6 + 5·d − 20·d² far out.
    @test calculate_cl(panel, 4window) ≈ edge + slope * 3window rtol = 1e-5
    lower_edge = calculate_cl(panel, -window)
    lower_slope = 5.0 - 2 * 20.0 * (-window)
    @test calculate_cl(panel, -4window) ≈ lower_edge - lower_slope * 3window rtol = 1e-5
    @test calculate_cl(panel, 4window) > 0.6 + 5.0 * 4window - 20.0 * (4window)^2
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

@testset "an unrepresentable chord line is removed, not projected" begin
    basis = KulfanBasis()
    base = KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    camber = @. 0.02 * basis.x * (1 - basis.x)

    # A deflection whose ends are off zero is a chord rotation and translation. The
    # basis cannot hold either, and projecting one anyway used to answer a 2.5% chord
    # deflection with weights four times the ones it was correcting.
    tilted = camber .+ 0.03 .* basis.x .+ 0.01
    offset, slope = chord_line(basis, tilted)
    @test offset ≈ 0.01 atol = 1e-12
    @test slope ≈ 0.03 atol = 1e-12
    residual = chord_residual(basis, tilted)
    @test residual ≈ camber atol = 1e-12
    @test abs(residual[1]) < 1e-12 && abs(residual[end]) < 1e-12

    # So the tilted deflection has to give the same airfoil as the pure camber, and
    # weights of the order of the camber rather than of the tilt.
    pure = deform_kulfan(basis, base, camber)
    tilt = deform_kulfan(basis, base, tilted)
    @test tilt.upper_weights ≈ pure.upper_weights
    @test maximum(abs, tilt.upper_weights .- base.upper_weights) <
          10 * maximum(abs, camber)

    # A pure chord rotation is no shape change at all.
    rotation = deform_kulfan(basis, base, 0.05 .* basis.x)
    @test rotation.upper_weights ≈ base.upper_weights atol = 1e-12
    @test rotation.lower_weights ≈ base.lower_weights atol = 1e-12
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

    # The panel keeps the object that was sampled, not a copy of it: a plot of
    # live_shape has to be a picture of what the solve flew.
    @test panels[1].live_shape === live.deformed[1]
    @test panels[1].live_shape.upper_weights !=
          live.base[1].upper_weights           # the deformation reached it
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6)
    @test panels[1].live_shape === live.deformed[1]
    @test panels[1].live_shape.upper_weights ≈ live.base[1].upper_weights

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

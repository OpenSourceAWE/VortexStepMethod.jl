using Test
using LinearAlgebra
using Statistics
using VortexStepMethod
using VortexStepMethod.AirfoilAero
using VortexStepMethod: Panel, Section, calculate_cl, calculate_cd, calculate_cm,
                        set_polar!, panel_interp_types, init_aero!

"""Panels carrying the interpolations a `POLAR_VECTORS` table is rewritten into."""
function polar_panels(n_panels; n_samples=5)
    alphas = collect(range(deg2rad(-6.0), deg2rad(6.0), n_samples))
    data = (alphas, zeros(n_samples), fill(0.01, n_samples), zeros(n_samples))
    section = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], POLAR_VECTORS, data)
    CL, CD, CM, CP = panel_interp_types(section, true)
    return map(1:n_panels) do _
        panel = Panel{Float64, CL, CD, CM, CP}()
        init_aero!(panel, section, section)
        panel
    end
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

@testset "rewritten panel polar" begin
    panel = only(polar_panels(1))
    alphas = deg2rad.([-6.0, -3.0, 0.0, 3.0, 6.0])
    set_polar!(panel, alphas, [0.0, 0.3, 0.6, 0.9, 0.8],
                       [0.03, 0.02, 0.02, 0.03, 0.06], [-0.1, -0.1, -0.1, -0.1, 0.0])
    @test panel.aero_model == POLAR_VECTORS
    @test panel.alpha_ref ≈ 0.0
    @test panel.alpha_window ≈ deg2rad(6.0)

    # A sample is the polar at its own angle, and between two it interpolates.
    @test calculate_cl(panel, deg2rad(3.0)) ≈ 0.9
    @test calculate_cl(panel, deg2rad(1.5)) ≈ 0.75
    @test calculate_cd(panel, deg2rad(-4.5)) ≈ 0.025
    @test calculate_cm(panel, deg2rad(4.5)) ≈ -0.05

    # Which is what a fit cannot do: the peak is held rather than averaged away.
    @test calculate_cl(panel, deg2rad(6.0)) < calculate_cl(panel, deg2rad(3.0))

    # Past either end the last sample is held, so a solve that steps out reads a
    # value the network gave rather than an extrapolation that runs away.
    @test calculate_cl(panel, deg2rad(40.0)) ≈ 0.8
    @test calculate_cl(panel, deg2rad(-40.0)) ≈ 0.0
    @test calculate_cd(panel, deg2rad(-40.0)) ≈ 0.03

    # A sampled polar has no flap axis: delta must not change the answer.
    @test calculate_cl(panel, deg2rad(1.5), deg2rad(9.0)) ==
          calculate_cl(panel, deg2rad(1.5))

    knots, coeffs = panel.alpha_knots, panel.cl_coeffs
    set_polar!(panel, alphas .+ deg2rad(2.0), [0.2, 0.5, 0.8, 1.1, 1.0],
                       [0.03, 0.02, 0.02, 0.03, 0.06], zeros(5))
    @test panel.alpha_knots === knots && panel.cl_coeffs === coeffs

    @test_throws ArgumentError set_polar!(panel, alphas, [0.0], [0.0], [0.0])
    @test_throws ArgumentError set_polar!(panel, reverse(alphas), zeros(5),
                                          zeros(5), zeros(5))
    @test_throws ArgumentError set_polar!(panel, alphas[1:1], [0.5], [0.02], [-0.1])

    # A panel built without interpolations says so rather than failing on a field type.
    @test_throws ArgumentError set_polar!(Panel{Float64}(), alphas, zeros(5), zeros(5),
                                          zeros(5))
end

@testset "a rewritten polar lands on a matrix panel too" begin
    # A wing whose offline polars carry a delta axis still has to take a live polar:
    # the deflection is in the shape it was generated from, so the table says the same
    # at every delta rather than having no answer for one.
    alpha = collect(deg2rad.(-15.0:1.0:15.0))
    delta = collect(deg2rad.(-40.0:10.0:40.0))
    na, nd = length(alpha), length(delta)
    grid = [0.1 * i for i in 1:na, _ in 1:nd]
    data = (alpha, delta, grid, fill(0.02, na, nd), fill(-0.05, na, nd))
    section = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], POLAR_MATRICES, data)
    CL, CD, CM, CP = panel_interp_types(section, true)
    panel = Panel{Float64, CL, CD, CM, CP}()
    init_aero!(panel, section, section)

    knots = collect(deg2rad.(-12.0:3.0:12.0))
    cl = collect(range(0.0, 1.2, length(knots)))
    set_polar!(panel, knots, cl, fill(0.03, 9), fill(-0.1, 9))

    @test panel.aero_model == POLAR_MATRICES
    @test calculate_cl(panel, knots[3]) ≈ cl[3]
    # The flap axis is gone from the answer: every delta reads the same polar.
    @test calculate_cl(panel, knots[3], deg2rad(-40.0)) ≈
          calculate_cl(panel, knots[3], deg2rad(40.0)) ≈ cl[3]
    # And it is held past the window it was sampled over, as on a vector panel.
    @test calculate_cl(panel, deg2rad(40.0)) ≈ cl[end]
end

@testset "live polars" begin
    base = KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    n_panels = 3
    panels = polar_panels(n_panels)
    live = LivePolars(fill(base, n_panels))
    confidence = refresh_live_polars!(live, panels, deg2rad(6.0), 3e6)
    @test 0.0 < confidence <= 1.0
    @test all(p -> p.aero_model == POLAR_VECTORS, panels)

    # Every sample is the network's own answer, not a fit through it.
    for (k, offset) in enumerate(live.settings.offsets)
        @test calculate_cl(panels[1], deg2rad(6.0) + offset) ≈
              neuralfoil_aero(base, 6.0 + rad2deg(offset), 3e6).CL[1] rtol = 1e-6
    end
    # And between them it tracks a direct sweep. Two bounds, because they are two
    # different claims: away from the knee the polar is nearly straight and the
    # samples pin it down, while at the knee the error is set by how sharp the
    # corner is rather than by the spacing — halving the spacing there buys almost
    # nothing, and the value stays bounded either way, which is what the solve needs.
    sweep(range_deg) = maximum(abs(calculate_cl(panels[1], deg2rad(6.0 + d)) -
                                   neuralfoil_aero(base, 6.0 + d, 3e6).CL[1])
                               for d in range_deg)
    @test sweep(range(-12, 6, 37)) < 0.02
    @test sweep(range(-12, 12, 49)) < 0.08

    # The confidence is kept per panel, so the comparison can pick the panel worth
    # asking about. It reads the polar the panel flies rather than re-running the
    # network, and on an undeformed section XFoil has to agree with it.
    @test length(live.confidence) == n_panels
    @test all(0.0 .< live.confidence .<= 1.0)
    check = compare_live_polar(live, panels[1], 1, deg2rad(6.0), 3e6)
    @test check.confidence == live.confidence[1]
    @test check.cl[1] ≈ calculate_cl(panels[1], deg2rad(6.0))
    @test isapprox(check.cl[1], check.cl[2]; atol=0.05)

    @test polar_drift(live, fill(deg2rad(6.0), n_panels)) ≈ 0.0 atol = 1e-12
    @test polar_drift(live, fill(deg2rad(18.0), n_panels)) ≈ 1.0

    # Deforming the camber up raises lift at the same angle of attack.
    flat = calculate_cl(panels[1], deg2rad(6.0))
    camber = @. 0.02 * live.basis.x * (1 - live.basis.x)
    knots = panels[1].alpha_knots
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6;
                         deflection=fill(camber, n_panels))
    @test calculate_cl(panels[1], deg2rad(6.0)) > flat
    @test panels[1].alpha_knots === knots   # a refresh writes in place

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
        settings=LivePolarSettings(; offsets=deg2rad.([1.0, 2.0, 3.0])))
    @test_throws ArgumentError LivePolars(fill(base, 2);
        settings=LivePolarSettings(; offsets=deg2rad.([4.0, -4.0])))
end

@testset "live surface pressure" begin
    base = KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    n_panels = 2
    panels = polar_panels(n_panels)
    settings = LivePolarSettings(; model_size="large")
    live = LivePolars(fill(base, n_panels); settings)
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6)

    # The contour the traction pattern is spread over: Selig order, nose at argmin.
    contour_x, contour_y = kulfan_to_coordinates(base)
    contour_x = collect(float.(contour_x))
    contour_y = collect(float.(contour_y))
    leading_edge = argmin(contour_x)
    cp = [zeros(length(contour_x)) for _ in 1:n_panels]
    same = cp[1]
    refresh_live_pressure!(cp, live, fill(contour_x, n_panels),
                           fill(contour_y, n_panels),
                           fill(leading_edge, n_panels), deg2rad(6.0), 3e6)
    @test cp[1] === same                    # written in place
    @test all(isfinite, cp[1])
    @test cp[1] ≈ cp[2]                     # same shape, same angle

    # This is the same pressure the offline tables were generated with, so the live
    # path and the table generator must not be two different answers for one shape.
    solver = NeuralFoilSolver(; model_size="large", n_crit=settings.n_crit)
    section = analyze_section(solver, DeformedSection(base, contour_x, contour_y),
        deg2rad(6.0), 3e6)
    @test cp[1] ≈ section.cp rtol = 1e-8

    # Suction over the upper surface, pressure under the lower one.
    @test minimum(cp[1][1:leading_edge]) < -0.5
    loading(p) = sum(view(p, (leading_edge + 1):length(p))) /
                 (length(p) - leading_edge) - sum(view(p, 1:leading_edge)) /
                 leading_edge
    @test loading(cp[1]) > 0.0

    # Bending the camber up loads the section harder, which is the whole claim: the
    # pattern moves with the deformation. Read as net loading, not as peak suction —
    # a cambered nose meets the flow better, so the leading-edge peak weakens while
    # the section as a whole carries more.
    camber = @. 0.02 * live.basis.x * (1 - live.basis.x)
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6;
                         deflection=fill(camber, n_panels))
    cambered = deepcopy(cp)
    refresh_live_pressure!(cambered, live, fill(contour_x, n_panels),
                           fill(contour_y, n_panels),
                           fill(leading_edge, n_panels), deg2rad(6.0), 3e6)
    @test loading(cambered[1]) > loading(cp[1])

    @test_throws ArgumentError refresh_live_pressure!(cp[1:1], live,
        fill(contour_x, n_panels), fill(contour_y, n_panels),
        fill(leading_edge, n_panels), 0.0, 3e6)
end

@testset "the contour follows the deformation" begin
    base = KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    panels = polar_panels(2)
    live = LivePolars(fill(base, 2))
    shape = [contour_shape_matrix(live.basis.x, 8) for _ in 1:2]
    offset = [zeros(length(live.basis.x)) for _ in 1:2]
    same = offset[1]

    # An undeformed panel's contour is the reference one, exactly.
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6)
    live_shape_offset!(offset, live, shape)
    @test offset[1] === same
    @test all(iszero, offset[1])

    # A camber increment comes back as itself: this is the same offset the airfoil
    # the network saw was built from, so the contour and the pressures agree.
    camber = @. 0.02 * live.basis.x * (1 - live.basis.x)
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6;
                         deflection=fill(camber, 2))
    live_shape_offset!(offset, live, shape)
    @test maximum(abs, offset[1] .- camber) < 0.02 * maximum(abs, camber)
    @test offset[1] ≈ offset[2]

    # A chord rotation is not a shape change, so it moves no contour node.
    refresh_live_polars!(live, panels, deg2rad(6.0), 3e6;
                         deflection=fill(0.05 .* live.basis.x, 2))
    live_shape_offset!(offset, live, shape)
    @test maximum(abs, offset[1]) < 1e-12

    @test_throws ArgumentError live_shape_offset!(offset[1:1], live, shape)
end

@testset "live skin friction" begin
    contour_x = [1.0, 0.5, 0.0, 0.5, 1.0]
    cf = [zeros(5), zeros(5)]
    same = cf[1]
    live_surface_friction!(cf, [contour_x, contour_x], [3e6, 1e6])
    @test cf[1] === same
    @test cf[1] ≈ [AirfoilAero.flat_plate_cf(x, 3e6) for x in contour_x]
    # Lower Reynolds is more friction, and it is highest at the nose.
    @test all(cf[2] .> cf[1])
    @test argmax(cf[1]) == 3
end

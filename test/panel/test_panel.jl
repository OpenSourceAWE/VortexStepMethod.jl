using VortexStepMethod: Panel, Section, calculate_relative_alpha_and_relative_velocity, calculate_cl, calculate_cd_cm, reinit!, INVISCID, POLAR_VECTORS, MVec3
using VortexStepMethod: panel_axes, panel_span_vector
using Interpolations: linear_interpolation, Line
using LinearAlgebra
using Test

function create_panel(section1::Section, section2::Section)
    # Calculate panel geometry
    section = Dict(
        "p1" => section1.LE_point,
        "p2" => section2.LE_point,
        "p3" => section2.TE_point,
        "p4" => section1.TE_point
    )
    
    bound_1 = section["p1"] .* 3/4 .+ section["p4"] .* 1/4
    bound_2 = section["p2"] .* 3/4 .+ section["p3"] .* 1/4

    mid_LE_point = section2.LE_point .+ 0.5 .* (section1.LE_point .- section2.LE_point)
    mid_TE_point = section2.TE_point .+ 0.5 .* (section1.TE_point .- section2.TE_point)
    mid_LE_vector = mid_TE_point .- mid_LE_point
    aero_center = bound_1 .+ 0.5 .* (bound_2 .- bound_1)
    control_point = aero_center .+ 0.5 .* mid_LE_vector

    LLpoint = aero_center
    VSMpoint = control_point
    z_airf = cross(VSMpoint .- LLpoint, section["p2"] .- section["p1"])
    z_airf = z_airf ./ norm(z_airf)

    # TANGENTIAL y_airf defined parallel to the chord-line
    x_airf = VSMpoint .- LLpoint
    x_airf = x_airf ./ norm(x_airf)

    # SPAN y_airf along the LE
    y_airf = bound_2 .- bound_1
    y_airf = y_airf ./ norm(y_airf)

    CL, CD, CM, CP = VortexStepMethod.panel_interp_types(section1, true)
    panel = Panel{Float64, CL, CD, CM, CP}()
    reinit!(
        panel,
        section1,
        section2,
        aero_center,
        control_point,
        bound_1,
        bound_2,
        x_airf,
        y_airf,
        z_airf,
        0.0,
        zeros(MVec3)
    )
    return panel
end

@testset "Panel Tests" begin
    @testset "Basic Panel Properties" begin
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], INVISCID)
        panel = create_panel(section1, section2)    

        # Test panel initialization
        @test panel isa Panel

        # Test TE/LE points
        @test isapprox(panel.TE_point_1, [1.0, 0.0, 0.0])
        @test isapprox(panel.TE_point_2, [1.0, 10.0, 0.0])
        @test isapprox(panel.LE_point_1, [0.0, 0.0, 0.0])
        @test isapprox(panel.LE_point_2, [0.0, 10.0, 0.0])

        # Test corner points
        expected_corners = [
        #   LE1  TE1   TE2   TE1
            0.0  1.0   1.0   0.0;  # x coordinates
            0.0  0.0  10.0  10.0;  # y coordinates
            0.0  0.0   0.0   0.0   # z coordinates
        ]
        @test all(isapprox.(panel.corner_points, expected_corners))

        # Test chord calculation
        rib_1 = panel.TE_point_1 .- panel.LE_point_1  # Vector from LE to TE for first section
        rib_2 = panel.TE_point_2 .- panel.LE_point_2  # Vector from LE to TE for second section
        expected_chord = (norm(rib_1) + norm(rib_2)) / 2
        @test isapprox(panel.chord, expected_chord)
    end

    @testset "Panel Reference Frame" begin
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], INVISCID)
        panel = create_panel(section1, section2)

        # Test reference frame vectors
        @test isapprox(panel.z_airf, [0.0, 0.0, 1.0])
        @test isapprox(panel.x_airf, [1.0, 0.0, 0.0])
        @test isapprox(panel.y_airf, [0.0, 1.0, 0.0])
    end

    @testset "Velocity Calculations" begin
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], INVISCID)
        panel = create_panel(section1, section2)

        # Test relative velocity calculations
        panel.va = [10.0, 0.0, 0.0]
        induced_velocity = [1.0, 1.0, 1.0]
        alpha, rel_vel = calculate_relative_alpha_and_relative_velocity(panel, induced_velocity)
        
        # Verify calculations
        norm_airf = panel.z_airf
        tan_airf = panel.x_airf
        relative_velocity = panel.va .+ induced_velocity
        vn = dot(norm_airf, relative_velocity)
        vtan = dot(tan_airf, relative_velocity)
        expected_alpha = atan(vn / vtan)

        @test isapprox(alpha, expected_alpha)
        @test isapprox(rel_vel, relative_velocity)
    end

    # Generate mock polar data
    alphas = deg2rad.(collect(-10:1:25))
    n_points = length(alphas)
    polar_data = (zeros(n_points), zeros(n_points), zeros(n_points), zeros(n_points))
    big_polar_data = (zeros(n_points), zeros(n_points), zeros(n_points), zeros(n_points))
    
    # Fill polar data with realistic values
    for (i, alpha) in enumerate(alphas)
        cd = 0.015 + 0.015 * abs(alpha/10)^1.5  # Drag increases with angle
        cl = 0.1 * alpha + 0.01 * alpha^2 * exp(-alpha/20)  # Lift with stall behavior
        cm = -0.02 * alpha  # Linear pitching moment
        polar_data[1][i] = alpha
        polar_data[2][i] = cl
        polar_data[3][i] = cd
        polar_data[4][i] = cm
        big_polar_data[1][i] = alpha
        big_polar_data[2][i] = 1.1cl
        big_polar_data[3][i] = 1.1cd
        big_polar_data[4][i] = 1.1cm
    end
    
    # Create two sections with slightly different polar data
    section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], POLAR_VECTORS, polar_data)
    section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], POLAR_VECTORS, big_polar_data)
    
    # Create panel
    panel = create_panel(section1, section2)
    
    @testset "Coefficient Interpolation" begin
        test_alphas = deg2rad.([-5.0, 0.0, 5.0, 10.0])
        for alpha in test_alphas
            # Calculate coefficients using panel methods
            cl = calculate_cl(panel, alpha)
            cd, cm = calculate_cd_cm(panel, alpha)
            
            # Calculate expected values through interpolation
            expected_cl = linear_interpolation(polar_data[1], polar_data[2], extrapolation_bc=Line())(alpha)
            expected_cd = linear_interpolation(polar_data[1], polar_data[3], extrapolation_bc=Line())(alpha)
            expected_cm = linear_interpolation(polar_data[1], polar_data[4], extrapolation_bc=Line())(alpha)
            
            # Average with slightly different data (1.1 factor)
            expected_cl = (expected_cl + 1.1 * expected_cl) / 2
            expected_cd = (expected_cd + 1.1 * expected_cd) / 2
            expected_cm = (expected_cm + 1.1 * expected_cm) / 2
            
            @test isapprox(cl, expected_cl, rtol=1e-5)
            @test isapprox(cd, expected_cd, rtol=1e-5)
            @test isapprox(cm, expected_cm, rtol=1e-5)
        end
    end
end

@testset "Panel frame on a swept, twisted panel" begin
    twist = deg2rad(35)
    le_1 = [0.0, 0.0, 0.0]
    te_1 = [2.0, 0.0, 0.0]
    le_2 = [1.2, 1.0, 0.0]
    te_2 = le_2 .+ 1.4 .* [cos(twist), 0.0, -sin(twist)]
    axes = panel_axes(le_1, te_1, le_2, te_2)
    span_vec = panel_span_vector(le_1, te_1, le_2, te_2)

    leading_edge_normal = cross(axes.x_airf, le_1 .- le_2)
    leading_edge_normal ./= norm(leading_edge_normal)
    separation = rad2deg(acos(clamp(abs(dot(leading_edge_normal, axes.z_airf)), -1, 1)))

    @testset "the geometry discriminates the two definitions" begin
        @test separation > 5
    end

    @testset "the normal is square to the bound vortex and the chord" begin
        @test abs(dot(axes.z_airf, span_vec)) < 1e-12
        @test abs(dot(axes.z_airf, axes.y_airf)) < 1e-12
        @test abs(dot(axes.z_airf, axes.x_airf)) < 1e-12
    end

    @testset "the frame closes as z = x cross y" begin
        closure = cross(axes.x_airf, axes.y_airf)
        @test isapprox(closure ./ norm(closure), axes.z_airf; atol=1e-12)
    end

    @testset "orient flips y and z together" begin
        flipped = panel_axes(le_1, te_1, le_2, te_2, 0.5, -1)
        @test isapprox(flipped.y_airf, -axes.y_airf; atol=1e-12)
        @test isapprox(flipped.z_airf, -axes.z_airf; atol=1e-12)
        @test isapprox(flipped.x_airf, axes.x_airf; atol=1e-12)
    end
end

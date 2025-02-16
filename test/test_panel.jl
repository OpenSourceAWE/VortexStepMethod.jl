using VortexStepMethod: Panel, Section, calculate_relative_alpha_and_relative_velocity
using LinearAlgebra
using Test
using BenchmarkTools

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
    aerodynamic_center = bound_1 .+ 0.5 .* (bound_2 .- bound_1)
    control_point = aerodynamic_center .+ 0.5 .* mid_LE_vector

    LLpoint = aerodynamic_center
    VSMpoint = control_point
    x_airf = cross(VSMpoint .- LLpoint, section["p2"] .- section["p1"])
    x_airf = x_airf ./ norm(x_airf)

    # TANGENTIAL y_airf defined parallel to the chord-line
    y_airf = VSMpoint .- LLpoint
    y_airf = y_airf ./ norm(y_airf)

    # SPAN z_airf along the LE
    z_airf = bound_2 .- bound_1
    z_airf = z_airf ./ norm(z_airf)

    Panel(
        section1,
        section2,
        aerodynamic_center,
        control_point,
        bound_1,
        bound_2,
        x_airf,
        y_airf,
        z_airf
    )
end

@testset "Panel Tests" begin
    @testset "Basic Panel Properties" begin
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "inviscid")
        section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], "inviscid")
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
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "inviscid")
        section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], "inviscid")
        panel = create_panel(section1, section2)

        # Test reference frame vectors
        @test isapprox(panel.x_airf, [0.0, 0.0, 1.0])
        @test isapprox(panel.y_airf, [1.0, 0.0, 0.0])
        @test isapprox(panel.z_airf, [0.0, 1.0, 0.0])
    end

    @testset "Velocity Calculations" begin
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "inviscid")
        section2 = Section([0.0, 10.0, 0.0], [1.0, 10.0, 0.0], "inviscid")
        panel = create_panel(section1, section2)

        # Test relative velocity calculations
        panel.va = [10.0, 0.0, 0.0]
        induced_velocity = [1.0, 1.0, 1.0]
        alpha, rel_vel = calculate_relative_alpha_and_relative_velocity(panel, induced_velocity)
        
        # Verify calculations
        norm_airf = panel.x_airf
        tan_airf = panel.y_airf
        relative_velocity = panel.va .+ induced_velocity
        vn = dot(norm_airf, relative_velocity)
        vtan = dot(tan_airf, relative_velocity)
        expected_alpha = atan(vn / vtan)

        @test isapprox(alpha, expected_alpha)
        @test isapprox(rel_vel, relative_velocity)
    end
end
using Test
using LinearAlgebra
using VortexStepMethod: Wing, Section, add_section!, refine_mesh_by_splitting_provided_sections, refine_aerodynamic_mesh
import Base: ==

"""
    ==(a::Section, b::Section)

Define equality for Section type. Two sections are equal if their leading edge points,
trailing edge points, and aerodynamic inputs are equal.

Points are compared with approximate equality to handle floating point differences.
"""
function ==(a::Section, b::Section)
    return (isapprox(a.LE_point, b.LE_point; rtol=1e-5, atol=1e-5) &&
            isapprox(a.TE_point, b.TE_point; rtol=1e-5, atol=1e-5) &&
            a.aero_model === b.aero_model &&
            all(a.aero_data .== b.aero_data))
end

@testset "Wing Geometry Tests" begin
    @testset "Wing initialization" begin
        example_wing = Wing(10; spanwise_panel_distribution=LINEAR)
        @test example_wing.n_panels == 10
        @test example_wing.spanwise_panel_distribution == LINEAR
        @test example_wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
        @test length(example_wing.sections) == 0
    end

    @testset "Add section" begin
        example_wing = Wing(10)
        add_section!(example_wing, [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], INVISCID)
        @test length(example_wing.sections) == 1
        
        section = example_wing.sections[1]
        @test section.LE_point ≈ [0.0, 0.0, 0.0]
        @test section.TE_point ≈ [-1.0, 0.0, 0.0]
        @test section.aero_model === INVISCID
    end

    @testset "Robustness left to right" begin
        example_wing = Wing(10)
        # Test correct order
        add_section!(example_wing, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], INVISCID)
        add_section!(example_wing, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], INVISCID)
        add_section!(example_wing, [0.0, -1.5, 0.0], [0.0, -1.5, 0.0], INVISCID)
        sections = refine_aerodynamic_mesh(example_wing)

        # Test right to left order
        example_wing_1 = Wing(10)
        add_section!(example_wing_1, [0.0, -1.5, 0.0], [0.0, -1.5, 0.0], INVISCID)
        add_section!(example_wing_1, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], INVISCID)
        add_section!(example_wing_1, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], INVISCID)
        sections_1 = refine_aerodynamic_mesh(example_wing_1)

        # Test random order
        example_wing_2 = Wing(10)
        add_section!(example_wing_2, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], INVISCID)
        add_section!(example_wing_2, [0.0, -1.5, 0.0], [0.0, -1.5, 0.0], INVISCID)
        add_section!(example_wing_2, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], INVISCID)
        sections_2 = refine_aerodynamic_mesh(example_wing_2)

        for i in eachindex(sections)
            @test sections[i].LE_point ≈ sections_1[i].LE_point
            @test sections[i].TE_point ≈ sections_1[i].TE_point
            @test sections[i].LE_point ≈ sections_2[i].LE_point
            @test sections[i].TE_point ≈ sections_2[i].TE_point
        end
    end

    @testset "Refine aerodynamic mesh" begin
        n_panels = 4
        span = 20.0

        # Test linear distribution
        wing = Wing(n_panels; spanwise_panel_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)
        sections = refine_aerodynamic_mesh(wing)

        @test length(sections) == wing.n_panels + 1

        for i in eachindex(sections)
            expected_LE = [0.0, span/2 - (i-1)*span/n_panels, 0.0]
            expected_TE = [-1.0, span/2 - (i-1)*span/n_panels, 0.0]
            @test isapprox(sections[i].LE_point, expected_LE; rtol=1e-5)
            @test isapprox(sections[i].TE_point, expected_TE; rtol=1e-4)
        end

        # Test cosine distribution
        wing = Wing(n_panels; spanwise_panel_distribution=COSINE)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)
        sections = refine_aerodynamic_mesh(wing)
        
        @test length(sections) == wing.n_panels + 1

        theta = range(0, π; length=n_panels+1)
        expected_LE_y = span/2 .* cos.(theta)
        
        for i in eachindex(sections)
            expected_LE = [0.0, expected_LE_y[i], 0.0]
            expected_TE = [-1.0, expected_LE_y[i], 0.0]
            @test isapprox(sections[i].LE_point, expected_LE; atol=1e-8)
            @test isapprox(sections[i].TE_point, expected_TE; atol=1e-8)
        end
    end

    @testset "Single panel" begin
        n_panels = 1
        span = 20.0

        wing = Wing(n_panels; spanwise_panel_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)

        sections = refine_aerodynamic_mesh(wing)
        @test length(sections) == wing.n_panels + 1
        @test sections[1].LE_point ≈ [0.0, span/2, 0.0]
        @test sections[1].TE_point ≈ [-1.0, span/2, 0.0]
    end

    @testset "Two panels" begin
        n_panels = 2
        span = 20.0

        wing = Wing(n_panels; spanwise_panel_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)

        sections = refine_aerodynamic_mesh(wing)
        @test length(sections) == wing.n_panels + 1
        @test sections[1].LE_point ≈ [0.0, span/2, 0.0]
        @test sections[2].LE_point ≈ [0.0, 0.0, 0.0]
        @test sections[3].LE_point ≈ [0.0, -span/2, 0.0]
    end

    @testset "More sections than panels" begin
        n_panels = 2
        span = 20.0

        wing = Wing(n_panels; spanwise_panel_distribution=LINEAR)
        y_coords = [span/2, span/4, 0.0, -span/4, -span/3, -span/2]
        for y in y_coords
            add_section!(wing, [0.0, y, 0.0], [-1.0, y, 0.0], INVISCID)
        end

        sections = refine_aerodynamic_mesh(wing)
        @test length(sections) == wing.n_panels + 1

        for i in eachindex(sections)
            expected_LE = [0.0, span/2 - (i-1)*span/n_panels, 0.0]
            expected_TE = [-1.0, span/2 - (i-1)*span/n_panels, 0.0]
            @test isapprox(sections[i].LE_point, expected_LE; rtol=1e-5)
            @test isapprox(sections[i].TE_point, expected_TE; rtol=1e-4)
        end
    end

    @testset "Symmetrical wing" begin
        n_panels = 2
        span = 10.0  # Total span from -5 to 5

        wing = Wing(n_panels; spanwise_panel_distribution=LINEAR)
        add_section!(wing, [0.0, 5.0, 0.0], [-1.0, 5.0, 0.0], INVISCID)
        add_section!(wing, [0.0, -5.0, 0.0], [-1.0, -5.0, 0.0], INVISCID)

        sections = refine_aerodynamic_mesh(wing)

        # Calculate expected quarter-chord points
        qc_start = [-0.25, 5.0, 0.0]
        qc_end = [-0.25, -5.0, 0.0]
        expected_qc_y = range(qc_start[2], qc_end[2]; length=n_panels+1)

        for (i, section) in enumerate(sections)
            # Calculate expected quarter-chord point
            expected_qc = [-0.25, expected_qc_y[i], 0.0]

            # Calculate expected chord vector
            chord_start = [-1.0, 5.0, 0.0] - [0.0, 5.0, 0.0]
            chord_end = [-1.0, -5.0, 0.0] - [0.0, -5.0, 0.0]
            t = (expected_qc_y[i] - qc_start[2]) / (qc_end[2] - qc_start[2])

            # Normalize chord vectors
            chord_start_norm = chord_start ./ norm(chord_start)
            chord_end_norm = chord_end ./ norm(chord_end)

            # Interpolate direction
            avg_direction = (1-t) .* chord_start_norm .+ t .* chord_end_norm
            avg_direction = avg_direction ./ norm(avg_direction)

            # Interpolate length
            chord_start_length = norm(chord_start)
            chord_end_length = norm(chord_end)
            avg_length = (1-t) * chord_start_length + t * chord_end_length

            expected_chord = avg_direction .* avg_length

            # Calculate expected LE and TE points
            expected_LE = expected_qc .- 0.25 .* expected_chord
            expected_TE = expected_qc .+ 0.75 .* expected_chord

            @debug "Section $i:" LE=section.LE_point expected_LE=expected_LE TE=section.TE_point expected_TE=expected_TE

            @test isapprox(section.LE_point, expected_LE; rtol=1e-5, atol=1e-5)
            @test isapprox(section.TE_point, expected_TE; rtol=1e-5, atol=1e-5)
        end

        @test length(sections) == n_panels + 1
        @test isapprox(sections[1].LE_point[2], 5.0; atol=1e-5)
        @test isapprox(sections[end].LE_point[2], -5.0; atol=1e-5)
        @test isapprox(sections[1].TE_point[2], 5.0; atol=1e-5)
        @test isapprox(sections[end].TE_point[2], -5.0; atol=1e-5)
    end

    @testset "LEI airfoil interpolation" begin
        n_panels = 4
        span = 20.0

        wing = Wing(n_panels; spanwise_panel_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], LEI_AIRFOIL_BREUKELS, (0.0, 0.0))
        add_section!(wing, [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], LEI_AIRFOIL_BREUKELS, (2.0, 0.5))
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], LEI_AIRFOIL_BREUKELS, (4.0, 1.0))

        sections = refine_aerodynamic_mesh(wing)
        @test length(sections) == wing.n_panels + 1

        expected_tube_diameter = range(0, 4; length=n_panels+1)
        expected_chamber_height = range(0, 1; length=n_panels+1)

        for (i, section) in enumerate(sections)
            expected_LE = [0.0, span/2 - (i-1)*span/n_panels, 0.0]
            expected_TE = [-1.0, span/2 - (i-1)*span/n_panels, 0.0]

            @test isapprox(section.LE_point, expected_LE; rtol=1e-5)
            @test isapprox(section.TE_point, expected_TE; rtol=1e-4)

            aero_model = section.aero_model
            aero_data = section.aero_data
            @test aero_model === LEI_AIRFOIL_BREUKELS
            @test isapprox(aero_data[1], expected_tube_diameter[i])
            @test isapprox(aero_data[2], expected_chamber_height[i])
        end
    end

    @testset "Split provided sections" begin
        section1 = Section([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        section2 = Section([0.0, 1.0, 0.0], [1.0, 1.0, 0.0], INVISCID)
        section3 = Section([0.0, 2.0, 0.0], [1.0, 2.0, 0.0], INVISCID)

        wing = Wing(6; spanwise_panel_distribution=SPLIT_PROVIDED)
        add_section!(wing, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        add_section!(wing, [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], INVISCID)
        add_section!(wing, [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], INVISCID)

        new_sections = refine_mesh_by_splitting_provided_sections(wing)

        @test length(new_sections) - 1 == 6
        @test new_sections[1] == section1
        @test new_sections[4] == section2
        @test new_sections[end] == section3

        @test 0.0 < new_sections[2].LE_point[2] < 1.0
        @test 0.0 < new_sections[3].LE_point[2] < 1.0
        @test 1.0 < new_sections[5].LE_point[2] < 2.0

        for section in new_sections
            @test section.aero_model === INVISCID
        end
    end
end
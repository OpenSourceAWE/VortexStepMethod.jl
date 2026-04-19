using Test
using LinearAlgebra
using VortexStepMethod
using VortexStepMethod: Wing, Section, add_section!, refine_mesh_by_splitting_provided_sections!, refine!
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
            a.aero_model == b.aero_model &&
            all(a.aero_data .== b.aero_data))
end

@testset "Wing Geometry Tests" begin
    @testset "Wing initialization" begin
        example_wing = Wing(10; spanwise_distribution=LINEAR)
        @test example_wing.n_panels == 10
        @test example_wing.spanwise_distribution == LINEAR
        @test example_wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
        @test length(example_wing.unrefined_sections) == 0
    end

    @testset "Add section" begin
        example_wing = Wing(10)
        add_section!(example_wing, [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], INVISCID)
        @test length(example_wing.unrefined_sections) == 1
        
        section = example_wing.unrefined_sections[1]
        @test section.LE_point ≈ [0.0, 0.0, 0.0]
        @test section.TE_point ≈ [-1.0, 0.0, 0.0]
        @test section.aero_model === INVISCID
    end

    @testset "Vector polar data" begin
        wing = Wing(2; remove_nan=true)
        
        # Create vector data with NaNs
        alpha = range(-10, 10, 21)
        cl = cos.(alpha)
        cd = sin.(alpha)
        cm = zeros(21)
        
        # Insert NaNs
        cl[5] = NaN
        cd[10] = NaN
        cm[15] = NaN
        
        aero_data = (collect(alpha), cl, cd, cm)
        add_section!(wing, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], POLAR_VECTORS, aero_data)
        
        # Check if NaNs were removed consistently
        cleaned_data = wing.unrefined_sections[1].aero_data
        @test length(cleaned_data[1]) == 18  # 21 - 3 NaN positions
        @test !any(isnan, cleaned_data[2])   # cl
        @test !any(isnan, cleaned_data[3])   # cd
        @test !any(isnan, cleaned_data[4])   # cm
        @test all(diff(cleaned_data[1]) .> 0) # alpha still monotonic
    end

    @testset "Matrix polar data" begin
        wing = Wing(2)
        
        # Create matrix data with NaNs
        alpha = collect(range(-10, 10, 21))
        delta = collect(range(-5, 5, 11))
        cl = [cos(a) * cos(b) for a in alpha, b in delta]
        cd = [sin(a) * sin(b) for a in alpha, b in delta]
        cm = zeros(21, 11)
        
        # Insert NaNs at various positions
        cl[5,3] = NaN
        cd[10,5] = NaN
        cm[15,7] = NaN
        
        aero_data = (alpha, delta, cl, cd, cm)
        add_section!(wing, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], POLAR_MATRICES, aero_data)
        
        # Check if NaNs were removed consistently
        cleaned_data = wing.unrefined_sections[1].aero_data
        @test !any(isnan, cleaned_data[3])   # cl
        @test !any(isnan, cleaned_data[4])   # cd
        @test !any(isnan, cleaned_data[5])   # cm
        @test all(diff(cleaned_data[1]) .> 0) # alpha still monotonic
        @test all(diff(cleaned_data[2]) .> 0) # delta still monotonic
    end

    @testset "Matrix polar interpolation" begin
        alpha = [-0.2, 0.0, 0.2]
        delta = [-0.1, 0.1]

        cl_left = [1.0 2.0; 3.0 4.0; 5.0 6.0]
        cd_left = [0.1 0.2; 0.3 0.4; 0.5 0.6]
        cm_left = [-0.1 -0.2; -0.3 -0.4; -0.5 -0.6]

        cl_right = cl_left .+ 6.0
        cd_right = cd_left .+ 0.6
        cm_right = cm_left .- 0.6

        polar_left = (alpha, delta, cl_left, cd_left, cm_left)
        polar_right = (copy(alpha), copy(delta), cl_right, cd_right, cm_right)

        result = VortexStepMethod.calculate_new_aero_data(
            (POLAR_MATRICES, POLAR_MATRICES),
            (polar_left, polar_right),
            1,
            0.25,
            0.75,
        )

        @test result[1] === alpha
        @test result[2] === delta
        @test result[3] ≈ cl_left .* 0.25 .+ cl_right .* 0.75
        @test result[4] ≈ cd_left .* 0.25 .+ cd_right .* 0.75
        @test result[5] ≈ cm_left .* 0.25 .+ cm_right .* 0.75

        bad_alpha_right = [-0.2, 0.05, 0.2]
        @test_throws(
            ArgumentError("Alpha steps must be identical."),
            VortexStepMethod.calculate_new_aero_data(
                (POLAR_MATRICES, POLAR_MATRICES),
                (polar_left, (bad_alpha_right, delta, cl_right, cd_right, cm_right)),
                1,
                0.5,
                0.5,
            )
        )

        bad_delta_right = [-0.1, 0.15]
        @test_throws(
            ArgumentError("Delta steps must be identical."),
            VortexStepMethod.calculate_new_aero_data(
                (POLAR_MATRICES, POLAR_MATRICES),
                (polar_left, (copy(alpha), bad_delta_right, cl_right, cd_right, cm_right)),
                1,
                0.5,
                0.5,
            )
        )
    end

    @testset "Robustness left to right" begin
        example_wing = Wing(10)
        # Test correct order
        add_section!(example_wing, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], INVISCID)
        add_section!(example_wing, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], INVISCID)
        add_section!(example_wing, [0.0, -1.5, 0.0], [0.0, -1.5, 0.0], INVISCID)
        refine!(example_wing)
        sections = example_wing.refined_sections

        # Test right to left order
        example_wing_1 = Wing(10)
        add_section!(example_wing_1, [0.0, -1.5, 0.0], [0.0, -1.5, 0.0], INVISCID)
        add_section!(example_wing_1, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], INVISCID)
        add_section!(example_wing_1, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], INVISCID)
        refine!(example_wing_1)
        sections_1 = example_wing_1.refined_sections

        # Test random order
        example_wing_2 = Wing(10)
        add_section!(example_wing_2, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], INVISCID)
        add_section!(example_wing_2, [0.0, -1.5, 0.0], [0.0, -1.5, 0.0], INVISCID)
        add_section!(example_wing_2, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], INVISCID)
        refine!(example_wing_2)
        sections_2 = example_wing_2.refined_sections

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
        wing = Wing(n_panels; spanwise_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)
        refine!(wing)
        sections = wing.refined_sections

        @test length(sections) == wing.n_panels + 1

        for i in eachindex(sections)
            expected_LE = [0.0, span/2 - (i-1)*span/n_panels, 0.0]
            expected_TE = [-1.0, span/2 - (i-1)*span/n_panels, 0.0]
            @test isapprox(sections[i].LE_point, expected_LE; rtol=1e-5)
            @test isapprox(sections[i].TE_point, expected_TE; rtol=1e-4)
        end

        # Test cosine distribution
        wing = Wing(n_panels; spanwise_distribution=COSINE)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)
        refine!(wing)
        sections = wing.refined_sections
        
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

        wing = Wing(n_panels; spanwise_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)

        refine!(wing)
        sections = wing.unrefined_sections
        @test length(sections) == wing.n_panels + 1
        @test sections[1].LE_point ≈ [0.0, span/2, 0.0]
        @test sections[1].TE_point ≈ [-1.0, span/2, 0.0]
    end

    @testset "Two panels" begin
        n_panels = 2
        span = 20.0

        wing = Wing(n_panels; spanwise_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], INVISCID)
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], INVISCID)

        refine!(wing)
        sections = wing.refined_sections
        @test length(sections) == wing.n_panels + 1
        @test sections[1].LE_point ≈ [0.0, span/2, 0.0]
        @test sections[2].LE_point ≈ [0.0, 0.0, 0.0]
        @test sections[3].LE_point ≈ [0.0, -span/2, 0.0]
    end

    @testset "More sections than panels" begin
        n_panels = 2
        span = 20.0

        wing = Wing(n_panels; spanwise_distribution=LINEAR)
        y_coords = [span/2, span/4, 0.0, -span/4, -span/3, -span/2]
        for y in y_coords
            add_section!(wing, [0.0, y, 0.0], [-1.0, y, 0.0], INVISCID)
        end

        refine!(wing)
        sections = wing.refined_sections
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

        wing = Wing(n_panels; spanwise_distribution=LINEAR)
        add_section!(wing, [0.0, 5.0, 0.0], [-1.0, 5.0, 0.0], INVISCID)
        add_section!(wing, [0.0, -5.0, 0.0], [-1.0, -5.0, 0.0], INVISCID)

        refine!(wing)
        sections = wing.refined_sections

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

        wing = Wing(n_panels; spanwise_distribution=LINEAR)
        add_section!(wing, [0.0, span/2, 0.0], [-1.0, span/2, 0.0], LEI_AIRFOIL_BREUKELS, (0.0, 0.0))
        add_section!(wing, [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], LEI_AIRFOIL_BREUKELS, (2.0, 0.5))
        add_section!(wing, [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0], LEI_AIRFOIL_BREUKELS, (4.0, 1.0))

        refine!(wing)
        sections = wing.refined_sections
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
        section2 = Section([0.0, -1.0, 0.0], [1.0, -1.0, 0.0], INVISCID)
        section3 = Section([0.0, -2.0, 0.0], [1.0, -2.0, 0.0], INVISCID)

        wing = Wing(6; spanwise_distribution=SPLIT_PROVIDED)
        add_section!(wing, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        add_section!(wing, [0.0, -1.0, 0.0], [1.0, -1.0, 0.0], INVISCID)
        add_section!(wing, [0.0, -2.0, 0.0], [1.0, -2.0, 0.0], INVISCID)

        refine!(wing)
        new_sections = wing.refined_sections

        @test length(new_sections) - 1 == 6
        @test new_sections[1] == section1
        @test new_sections[4] == section2
        @test new_sections[end] == section3

        @test -1.0 < new_sections[2].LE_point[2] < 0.0
        @test -1.0 < new_sections[3].LE_point[2] < 0.0
        @test -2.0 < new_sections[5].LE_point[2] < -1.0

        for section in new_sections
            @test section.aero_model == INVISCID
        end
    end

    @testset "use_prior_polar reuses refined aero data" begin
        alpha = deg2rad.([-5.0, 0.0, 5.0])
        left_aero = (collect(alpha), [1.0, 1.2, 1.4], [0.1, 0.11, 0.12], [0.01, 0.02, 0.03])
        right_aero = (collect(alpha), [3.0, 3.2, 3.4], [0.3, 0.31, 0.32], [0.03, 0.04, 0.05])

        wing = Wing(4; spanwise_distribution=LINEAR, use_prior_polar=true)
        add_section!(wing, [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], POLAR_VECTORS, left_aero)
        add_section!(wing, [0.0, -1.0, 0.0], [1.0, -1.0, 0.0], POLAR_VECTORS, right_aero)
        refine!(wing)

        baseline_section_aero = deepcopy(wing.refined_sections[3].aero_data)
        baseline_le = copy(wing.refined_sections[1].LE_point)

        wing.unrefined_sections[1].aero_data = (collect(alpha), fill(99.0, 3), fill(9.0, 3), fill(0.9, 3))
        wing.unrefined_sections[2].aero_data = (collect(alpha), fill(-99.0, 3), fill(8.0, 3), fill(-0.8, 3))
        wing.unrefined_sections[1].LE_point .+= [0.2, 0.0, 0.0]
        wing.unrefined_sections[1].TE_point .+= [0.2, 0.0, 0.0]
        refine!(wing)

        @test wing.refined_sections[3].aero_data == baseline_section_aero
        @test !isapprox(wing.refined_sections[1].LE_point[1], baseline_le[1])

        wing_no_reuse = Wing(4; spanwise_distribution=LINEAR, use_prior_polar=false)
        add_section!(wing_no_reuse, [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], POLAR_VECTORS, left_aero)
        add_section!(wing_no_reuse, [0.0, -1.0, 0.0], [1.0, -1.0, 0.0], POLAR_VECTORS, right_aero)
        refine!(wing_no_reuse)
        no_reuse_baseline = deepcopy(wing_no_reuse.refined_sections[3].aero_data)

        wing_no_reuse.unrefined_sections[1].aero_data = (collect(alpha), fill(99.0, 3), fill(9.0, 3), fill(0.9, 3))
        wing_no_reuse.unrefined_sections[2].aero_data = (collect(alpha), fill(-99.0, 3), fill(8.0, 3), fill(-0.8, 3))
        refine!(wing_no_reuse)

        @test wing_no_reuse.refined_sections[3].aero_data != no_reuse_baseline
    end

    @testset "Refined panel mapping" begin
        # Test that refined panel mapping actually maps each panel to its closest unrefined panel

        @testset "LINEAR distribution" begin
            n_panels = 20
            span = 10.0

            wing = Wing(n_panels; spanwise_distribution=LINEAR)
            # 3 sections
            add_section!(wing, [0.0, span/2, 0.0], [1.0, span/2, 0.0], INVISCID)
            add_section!(wing, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
            add_section!(wing, [0.0, -span/2, 0.0], [1.0, -span/2, 0.0], INVISCID)

            refine!(wing)

            @test length(wing.refined_panel_mapping) == n_panels

            # Manually verify each refined panel is mapped to its closest unrefined section
            n_unrefined_sections = length(wing.unrefined_sections)
            for refined_panel_idx in 1:n_panels
                # Calculate refined panel center
                le_mid = (wing.refined_sections[refined_panel_idx].LE_point +
                         wing.refined_sections[refined_panel_idx+1].LE_point) / 2
                te_mid = (wing.refined_sections[refined_panel_idx].TE_point +
                         wing.refined_sections[refined_panel_idx+1].TE_point) / 2
                refined_center = (le_mid + te_mid) / 2

                # Find closest unrefined section manually
                min_dist = Inf
                closest_idx = 1
                for unrefined_section_idx in 1:n_unrefined_sections
                    le_point = wing.unrefined_sections[unrefined_section_idx].LE_point
                    te_point = wing.unrefined_sections[unrefined_section_idx].TE_point
                    unrefined_center = (le_point + te_point) / 2

                    dist = norm(refined_center - unrefined_center)
                    if dist < min_dist
                        min_dist = dist
                        closest_idx = unrefined_section_idx
                    end
                end

                # Verify mapping is correct
                @test wing.refined_panel_mapping[refined_panel_idx] == closest_idx
            end
        end

        @testset "COSINE distribution" begin
            n_panels = 30
            span = 20.0

            wing = Wing(n_panels; spanwise_distribution=COSINE)
            # 4 sections
            add_section!(wing, [0.0, span/2, 0.0], [1.0, span/2, 0.0], INVISCID)
            add_section!(wing, [0.0, span/6, 0.0], [1.0, span/6, 0.0], INVISCID)
            add_section!(wing, [0.0, -span/6, 0.0], [1.0, -span/6, 0.0], INVISCID)
            add_section!(wing, [0.0, -span/2, 0.0], [1.0, -span/2, 0.0], INVISCID)

            refine!(wing)

            @test length(wing.refined_panel_mapping) == n_panels

            # Verify each panel is mapped to its closest unrefined section
            n_unrefined_sections = length(wing.unrefined_sections)
            for refined_panel_idx in 1:n_panels
                # Calculate refined panel center
                le_mid = (wing.refined_sections[refined_panel_idx].LE_point +
                         wing.refined_sections[refined_panel_idx+1].LE_point) / 2
                te_mid = (wing.refined_sections[refined_panel_idx].TE_point +
                         wing.refined_sections[refined_panel_idx+1].TE_point) / 2
                refined_center = (le_mid + te_mid) / 2

                # Find closest unrefined section manually
                min_dist = Inf
                closest_idx = 1
                for unrefined_section_idx in 1:n_unrefined_sections
                    le_point = wing.unrefined_sections[unrefined_section_idx].LE_point
                    te_point = wing.unrefined_sections[unrefined_section_idx].TE_point
                    unrefined_center = (le_point + te_point) / 2

                    dist = norm(refined_center - unrefined_center)
                    if dist < min_dist
                        min_dist = dist
                        closest_idx = unrefined_section_idx
                    end
                end

                # Verify mapping is correct
                @test wing.refined_panel_mapping[refined_panel_idx] == closest_idx
            end
        end

        @testset "SPLIT_PROVIDED distribution" begin
            n_panels = 12

            wing = Wing(n_panels; spanwise_distribution=SPLIT_PROVIDED)
            # 4 sections
            add_section!(wing, [0.0, 6.0, 0.0], [1.0, 6.0, 0.0], INVISCID)
            add_section!(wing, [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], INVISCID)
            add_section!(wing, [0.0, -2.0, 0.0], [1.0, -2.0, 0.0], INVISCID)
            add_section!(wing, [0.0, -6.0, 0.0], [1.0, -6.0, 0.0], INVISCID)

            refine!(wing)

            @test length(wing.refined_panel_mapping) == n_panels

            # Verify each panel is mapped to its closest unrefined section
            n_unrefined_sections = length(wing.unrefined_sections)
            for refined_panel_idx in 1:n_panels
                # Calculate refined panel center
                le_mid = (wing.refined_sections[refined_panel_idx].LE_point +
                         wing.refined_sections[refined_panel_idx+1].LE_point) / 2
                te_mid = (wing.refined_sections[refined_panel_idx].TE_point +
                         wing.refined_sections[refined_panel_idx+1].TE_point) / 2
                refined_center = (le_mid + te_mid) / 2

                # Find closest unrefined section manually
                min_dist = Inf
                closest_idx = 1
                for unrefined_section_idx in 1:n_unrefined_sections
                    le_point = wing.unrefined_sections[unrefined_section_idx].LE_point
                    te_point = wing.unrefined_sections[unrefined_section_idx].TE_point
                    unrefined_center = (le_point + te_point) / 2

                    dist = norm(refined_center - unrefined_center)
                    if dist < min_dist
                        min_dist = dist
                        closest_idx = unrefined_section_idx
                    end
                end

                # Verify mapping is correct
                @test wing.refined_panel_mapping[refined_panel_idx] == closest_idx
            end
        end

    end

    @testset "UNCHANGED preserves sections exactly" begin
        n_panels = 4
        span = 10.0
        wing = Wing(n_panels; spanwise_distribution=UNCHANGED)
        ys = range(span/2, -span/2; length=n_panels+1)
        for y in ys
            add_section!(wing,
                [0.0, y, 0.0], [-1.0, y, 0.0], INVISCID)
        end
        refine!(wing)
        @test length(wing.refined_sections) == n_panels + 1
        for (i, sec) in enumerate(wing.refined_sections)
            @test isapprox(sec.LE_point,
                wing.unrefined_sections[i].LE_point;
                atol=1e-14)
            @test isapprox(sec.TE_point,
                wing.unrefined_sections[i].TE_point;
                atol=1e-14)
            chord = norm(sec.TE_point - sec.LE_point)
            @test isapprox(chord, 1.0; atol=1e-14)
        end
        # Re-refine preserves same result
        refine!(wing)
        for (i, sec) in enumerate(wing.refined_sections)
            @test isapprox(sec.LE_point,
                wing.unrefined_sections[i].LE_point;
                atol=1e-14)
        end
    end

    @testset "COSINE concentrates panels at tips" begin
        n_panels = 20
        span = 10.0
        wing = Wing(n_panels; spanwise_distribution=COSINE)
        add_section!(wing,
            [0.0, span/2, 0.0], [-1.0, span/2, 0.0],
            INVISCID)
        add_section!(wing,
            [0.0, -span/2, 0.0], [-1.0, -span/2, 0.0],
            INVISCID)
        refine!(wing)
        sections = wing.refined_sections
        @test length(sections) == n_panels + 1
        # Endpoints match
        @test isapprox(sections[1].LE_point[2],
            span/2; atol=1e-10)
        @test isapprox(sections[end].LE_point[2],
            -span/2; atol=1e-10)
        # Tip spacing smaller than root spacing
        tip_gap = abs(sections[1].LE_point[2] -
            sections[2].LE_point[2])
        mid = div(n_panels, 2) + 1
        root_gap = abs(sections[mid].LE_point[2] -
            sections[mid+1].LE_point[2])
        @test tip_gap < root_gap
        # All chords preserved
        for sec in sections
            @test isapprox(
                norm(sec.TE_point - sec.LE_point),
                1.0; atol=1e-10)
        end
    end

    @testset "SPLIT_PROVIDED chord preservation" begin
        n_panels = 12
        n_ribs = 4
        span = 9.0
        wing = Wing(n_panels;
            spanwise_distribution=SPLIT_PROVIDED)
        for i in 1:n_ribs
            y = span/2 - (i-1) * span / (n_ribs - 1)
            add_section!(wing,
                [0.0, y, 0.0], [-1.0, y, 0.0], INVISCID)
        end
        refine!(wing)
        @test length(wing.refined_sections) == n_panels + 1
        for sec in wing.refined_sections
            @test isapprox(
                norm(sec.TE_point - sec.LE_point),
                1.0; atol=1e-10)
        end
        # Rib endpoints preserved
        for i in 1:n_ribs
            found = false
            for sec in wing.refined_sections
                if isapprox(sec.LE_point,
                        wing.unrefined_sections[i].LE_point;
                        atol=1e-10)
                    @test isapprox(sec.TE_point,
                        wing.unrefined_sections[i].TE_point;
                        atol=1e-10)
                    found = true
                    break
                end
            end
            @test found
        end
    end
end

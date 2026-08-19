using VortexStepMethod
using LinearAlgebra
using Test

@testset "YAML Wing Deformation Tests" begin
    @testset "Simple Wing Deformation" begin
        # Load existing simple_wing.yaml
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        # Store original TE point for comparison
        i = length(body_aero.panels) ÷ 2
        original_te_point = copy(body_aero.panels[i].TE_point_1)
        original_le_point = copy(body_aero.panels[i].LE_point_1)

        # Apply deformation with non-zero angles (panel-level)
        theta_dist = fill(deg2rad(30.0), wing.n_panels)  # 30 degrees twist per panel
        delta_dist = fill(deg2rad(5.0), wing.n_panels)   # 5 degrees TE deflection per panel

        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero)

        # Check if TE point changed after deformation
        deformed_te_point = copy(body_aero.panels[i].TE_point_1)
        deformed_le_point = copy(body_aero.panels[i].LE_point_1)

        # TE point should change significantly due to twist and deflection
        @test !isapprox(original_te_point, deformed_te_point, atol=1e-2)
        @test deformed_te_point[3] < original_te_point[3]  # TE should move down with positive twist

        # LE point stays fixed (deform! rotates TE around LE)
        @test isapprox(original_le_point, deformed_le_point, atol=1e-4)

        # Check delta is set correctly
        @test body_aero.panels[i].delta ≈ deg2rad(5.0)

        # Reset deformation with zero angles
        zero_theta_dist = zeros(wing.n_panels)
        zero_delta_dist = zeros(wing.n_panels)

        VortexStepMethod.deform!(wing, zero_theta_dist, zero_delta_dist)
        VortexStepMethod.reinit!(body_aero)

        # Check if TE point returned to original position
        reset_te_point = copy(body_aero.panels[i].TE_point_1)
        reset_le_point = copy(body_aero.panels[i].LE_point_1)
        @test original_te_point ≈ reset_te_point atol=1e-4
        @test original_le_point ≈ reset_le_point atol=1e-4
        @test body_aero.panels[i].delta ≈ 0.0 atol=1e-4
    end

    @testset "Complex Wing Deformation" begin
        # Load existing complex_wing.yaml with multiple sections
        complex_wing_file = test_data_path("yaml_geometry", "complex_wing.yaml")
        wing = Wing(complex_wing_file; n_panels=12)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        # Store original points for multiple panels
        original_points = []
        test_indices = [1, length(body_aero.panels) ÷ 2, length(body_aero.panels)]
        for i in test_indices
            push!(original_points, (
                LE=copy(body_aero.panels[i].LE_point_1),
                TE=copy(body_aero.panels[i].TE_point_1)
            ))
        end

        # Apply spanwise-varying deformation (panel-level)
        n = wing.n_panels
        theta_dist = [deg2rad(10.0 * i / n) for i in 1:n]  # Linear twist distribution
        delta_dist = [deg2rad(-5.0 + 10.0 * i / n) for i in 1:n]  # Varying deflection

        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero)

        # Check that different panels have different deformations
        for (idx, i) in enumerate(test_indices)
            deformed_te = body_aero.panels[i].TE_point_1
            deformed_le = body_aero.panels[i].LE_point_1

            # TE points should have changed due to deformation
            @test !isapprox(original_points[idx].TE, deformed_te, atol=1e-2)
            # LE points stay fixed (deform! rotates TE around LE)
            @test isapprox(original_points[idx].LE, deformed_le, atol=1e-4)
        end

        # Check that the deformation is applied correctly
        # First panel should have smaller theta, last panel should have larger theta
        @test body_aero.panels[1].delta < body_aero.panels[end].delta

        # Reset and verify
        VortexStepMethod.deform!(wing, zeros(wing.n_panels), zeros(wing.n_panels))
        VortexStepMethod.reinit!(body_aero)

        for (idx, i) in enumerate(test_indices)
            reset_te = body_aero.panels[i].TE_point_1
            reset_le = body_aero.panels[i].LE_point_1
            @test original_points[idx].TE ≈ reset_te atol=1e-4
            @test original_points[idx].LE ≈ reset_le atol=1e-4
            @test body_aero.panels[i].delta ≈ 0.0 atol=1e-4
        end
    end

    @testset "Multiple Reinit Calls with NTuple aero_data" begin
        # This test specifically checks the NTuple handling fix
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)
        refine!(wing)

        # Verify that sections have NTuple aero_data (for wings with simple polars)
        # or other valid AeroData types
        @test wing.unrefined_sections[1].aero_data !== nothing

        # Perform multiple reinit! calls to ensure NTuple handling works
        for _ in 1:5
            VortexStepMethod.reinit!(wing)
        end

        # Wing should still be valid after multiple reinits
        @test wing.unrefined_sections[1].aero_data !== nothing
        # Verify refined_sections and non_deformed_sections are properly populated
        @test length(wing.refined_sections) == wing.n_panels + 1
        @test length(wing.non_deformed_sections) == wing.n_panels + 1
    end

    @testset "Deformation with BodyAerodynamics Reinit" begin
        # Test that reinit! on BodyAerodynamics properly handles deformed wings
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        # Apply deformation
        theta_dist = fill(deg2rad(15.0), wing.n_panels)
        delta_dist = fill(deg2rad(3.0), wing.n_panels)
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero)

        # Store state after deformation
        i = length(body_aero.panels) ÷ 2

        # Multiple reinit calls should work without errors
        for _ in 1:3
            VortexStepMethod.reinit!(body_aero;
                va=zeros(3),
                omega=zeros(3),
                init_aero=true,
                
            )
        end

        # Panel should maintain deformation
        @test body_aero.panels[i].delta ≈ deg2rad(3.0) atol=1e-6
    end

    @testset "Edge Cases" begin
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=2)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        # Test zero deformation
        VortexStepMethod.deform!(wing, zeros(wing.n_panels), zeros(wing.n_panels))
        VortexStepMethod.reinit!(body_aero)
        @test all(p.delta ≈ 0.0 for p in body_aero.panels)

        # Test large deformation angles
        theta_dist = fill(deg2rad(60.0), wing.n_panels)
        delta_dist = fill(deg2rad(30.0), wing.n_panels)

        # Should not error even with large angles
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero)
        @test all(p.delta ≈ deg2rad(30.0) for p in body_aero.panels)

        # Test negative angles
        theta_dist = fill(deg2rad(-20.0), wing.n_panels)
        delta_dist = fill(deg2rad(-10.0), wing.n_panels)
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero)
        @test all(p.delta ≈ deg2rad(-10.0) for p in body_aero.panels)
    end

    @testset "Panel to Section Angle Mapping" begin
        # Test that panel angles are correctly averaged to section angles
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        # Create varying panel angles
        theta_panel = [10.0, 20.0, 30.0, 40.0]  # degrees per panel
        delta_panel = [5.0, 10.0, 15.0, 20.0]   # degrees per panel

        # Apply deformation
        VortexStepMethod.deform!(wing, deg2rad.(theta_panel), deg2rad.(delta_panel))

        # Verify section angles by checking the geometry
        # Section 1: should use panel 1 angle = 10°
        # Section 2: should avg panels 1,2 = (10+20)/2 = 15°
        # Section 3: should avg panels 2,3 = (20+30)/2 = 25°
        # Section 4: should avg panels 3,4 = (30+40)/2 = 35°
        # Section 5: should use panel 4 angle = 40°

        # We can verify this by checking that delta values are correct
        VortexStepMethod.reinit!(body_aero)

        # Each panel gets its delta directly
        @test body_aero.panels[1].delta ≈ deg2rad(5.0) atol=1e-6
        @test body_aero.panels[2].delta ≈ deg2rad(10.0) atol=1e-6
        @test body_aero.panels[3].delta ≈ deg2rad(15.0) atol=1e-6
        @test body_aero.panels[4].delta ≈ deg2rad(20.0) atol=1e-6
    end

    @testset "unrefined_deform! Linear Interpolation" begin
        # unrefined_deform! linearly interpolates angles from unrefined sections to
        # refined sections; panel-level dist values are the average of adjacent
        # refined-section values.
        complex_wing_file = test_data_path("yaml_geometry", "complex_wing.yaml")
        wing = Wing(complex_wing_file; n_panels=12)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        @test wing.n_unrefined_sections == 7

        theta_unrefined = deg2rad.([10.0, 15.0, 20.0, 25.0, 20.0, 15.0, 10.0])
        delta_unrefined = deg2rad.([5.0, 7.5, 10.0, 12.5, 10.0, 7.5, 5.0])

        VortexStepMethod.unrefined_deform!(wing, theta_unrefined, delta_unrefined)
        VortexStepMethod.reinit!(body_aero)

        # Endpoint refined sections take exact unrefined endpoint values, so the
        # first/last panel value is the average of an endpoint and its neighbor.
        n_sec = wing.n_panels + 1
        @test wing.refined_section_weight[1] ≈ 1.0
        @test wing.refined_section_weight[n_sec] ≈ 0.0
        @test wing.refined_section_left_idx[1] == 1
        @test wing.refined_section_left_idx[n_sec] == wing.n_unrefined_sections - 1

        # All panel deltas lie within the input range (no overshoot from interp).
        delta_min = minimum(delta_unrefined)
        delta_max = maximum(delta_unrefined)
        for p in body_aero.panels
            @test delta_min - 1e-9 ≤ p.delta ≤ delta_max + 1e-9
        end
    end

    @testset "unrefined_deform! Endpoint and Linearity" begin
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=40)
        refine!(wing)
        @test wing.n_unrefined_sections == 2

        delta_input = deg2rad.([0.0, 10.0])
        VortexStepMethod.unrefined_deform!(wing, nothing, delta_input)
        delta_panels = copy(wing.delta_dist)

        # With two unrefined sections, the refined-section interpolation is exactly
        # linear in arc-length, so panel-level values (averages of adjacent refined
        # sections) are linear and monotonic — no discontinuities.
        for i in 1:(wing.n_panels - 1)
            @test delta_panels[i + 1] ≥ delta_panels[i] - 1e-12
        end
        max_gradient = maximum(abs.(diff(delta_panels)))
        @test max_gradient < (delta_input[2] - delta_input[1]) / wing.n_panels * 1.5

        # Refined-section endpoints match unrefined endpoints exactly.
        @test wing.refined_section_weight[1] ≈ 1.0
        @test wing.refined_section_weight[end] ≈ 0.0
    end

    @testset "Twist rotates a swept chord rigidly" begin
        # A swept or dihedral section's chord has a component along the twist
        # axis, which only a full rotation carries through unchanged.
        wing = Wing(8; n_unrefined_sections=3)
        add_section!(wing, [0.0, 5.0, 1.0], [1.0, 4.6, 1.0], INVISCID)
        add_section!(wing, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], INVISCID)
        add_section!(wing, [0.0, -5.0, 1.0], [1.0, -4.6, 1.0], INVISCID)
        refine!(wing)
        chords = [norm(s.TE_point - s.LE_point) for s in wing.refined_sections]
        for theta in (0.05, 0.2, 0.5)
            VortexStepMethod.unrefined_deform!(wing, fill(theta, 3), nothing)
            twisted = [norm(s.TE_point - s.LE_point) for s in wing.refined_sections]
            @test maximum(abs.(twisted .- chords) ./ chords) < 1e-12
        end
    end
end

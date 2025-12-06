using VortexStepMethod
using LinearAlgebra
using Test

@testset "YAML Wing Deformation Tests" begin
    @testset "Simple Wing Deformation" begin
        # Load existing simple_wing.yaml
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)
        body_aero = BodyAerodynamics([wing])

        # Store original TE point for comparison
        i = length(body_aero.panels) ÷ 2
        original_te_point = copy(body_aero.panels[i].TE_point_1)
        original_le_point = copy(body_aero.panels[i].LE_point_1)

        # Apply deformation with non-zero angles (panel-level)
        theta_dist = fill(deg2rad(30.0), wing.n_panels)  # 30 degrees twist per panel
        delta_dist = fill(deg2rad(5.0), wing.n_panels)   # 5 degrees TE deflection per panel

        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

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
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

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
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

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
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

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
        body_aero = BodyAerodynamics([wing])

        # Apply deformation
        theta_dist = fill(deg2rad(15.0), wing.n_panels)
        delta_dist = fill(deg2rad(3.0), wing.n_panels)
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

        # Store state after deformation
        i = length(body_aero.panels) ÷ 2

        # Multiple reinit calls should work without errors
        for _ in 1:3
            VortexStepMethod.reinit!(body_aero;
                va=zeros(3),
                omega=zeros(3),
                init_aero=true,
                refine_mesh=false
            )
        end

        # Panel should maintain deformation
        @test body_aero.panels[i].delta ≈ deg2rad(3.0) atol=1e-6
    end

    @testset "Edge Cases" begin
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=2)
        body_aero = BodyAerodynamics([wing])

        # Test zero deformation
        VortexStepMethod.deform!(wing, zeros(wing.n_panels), zeros(wing.n_panels))
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)
        @test all(p.delta ≈ 0.0 for p in body_aero.panels)

        # Test large deformation angles
        theta_dist = fill(deg2rad(60.0), wing.n_panels)
        delta_dist = fill(deg2rad(30.0), wing.n_panels)

        # Should not error even with large angles
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)
        @test all(p.delta ≈ deg2rad(30.0) for p in body_aero.panels)

        # Test negative angles
        theta_dist = fill(deg2rad(-20.0), wing.n_panels)
        delta_dist = fill(deg2rad(-10.0), wing.n_panels)
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)
        @test all(p.delta ≈ deg2rad(-10.0) for p in body_aero.panels)
    end

    @testset "Panel to Section Angle Mapping" begin
        # Test that panel angles are correctly averaged to section angles
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)
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
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

        # Each panel gets its delta directly
        @test body_aero.panels[1].delta ≈ deg2rad(5.0) atol=1e-6
        @test body_aero.panels[2].delta ≈ deg2rad(10.0) atol=1e-6
        @test body_aero.panels[3].delta ≈ deg2rad(15.0) atol=1e-6
        @test body_aero.panels[4].delta ≈ deg2rad(20.0) atol=1e-6
    end

    @testset "unrefined_deform! Maps to Panels" begin
        # Test that unrefined_deform! correctly maps unrefined sections to panels
        # Use complex_wing which has 7 unrefined sections
        complex_wing_file = test_data_path("yaml_geometry", "complex_wing.yaml")
        wing = Wing(complex_wing_file; n_panels=12)
        body_aero = BodyAerodynamics([wing])

        # Verify we have 7 unrefined sections
        @test wing.n_unrefined_sections == 7

        # Create unrefined section angles (7 sections)
        # These will be mapped to panels via refined_panel_mapping
        theta_unrefined = deg2rad.([10.0, 15.0, 20.0, 25.0, 20.0, 15.0, 10.0])
        delta_unrefined = deg2rad.([5.0, 7.5, 10.0, 12.5, 10.0, 7.5, 5.0])

        # Apply using unrefined_deform!
        VortexStepMethod.unrefined_deform!(wing, theta_unrefined, delta_unrefined)
        VortexStepMethod.reinit!(body_aero; refine_mesh=false)

        # Each panel should have the delta from its mapped unrefined section
        for i in 1:wing.n_panels
            unrefined_idx = wing.refined_panel_mapping[i]
            expected_delta = delta_unrefined[unrefined_idx]
            @test body_aero.panels[i].delta ≈ expected_delta atol=1e-6
        end
    end

    @testset "Smooth vs Non-Smooth Deformation" begin
        # Create test wing with 2 unrefined sections, refined to 40 panels
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=40)
        @test wing.n_unrefined_sections == 2

        # Define varying input angles at unrefined section level
        delta_input = deg2rad.([0.0, 10.0])

        # Test 1: Non-smooth deformation has step-wise discontinuities
        VortexStepMethod.unrefined_deform!(wing, nothing, delta_input; smooth=false)
        delta_nonsmooth = copy(wing.delta_dist)

        # Verify step-wise pattern: panels in same unrefined section have identical angles
        for i in 1:wing.n_panels
            unrefined_idx = wing.refined_panel_mapping[i]
            @test delta_nonsmooth[i] ≈ delta_input[unrefined_idx] atol=1e-10
        end

        # Verify discontinuities exist at unrefined section boundaries
        # Find boundary indices (where panel mapping changes)
        max_gradient_nonsmooth = 0.0
        for i in 1:(wing.n_panels-1)
            if wing.refined_panel_mapping[i] != wing.refined_panel_mapping[i+1]
                gradient = abs(delta_nonsmooth[i+1] - delta_nonsmooth[i])
                max_gradient_nonsmooth = max(max_gradient_nonsmooth, gradient)
            end
        end
        @test max_gradient_nonsmooth > deg2rad(5.0)  # Should have large jumps

        # Test 2: Smooth deformation is continuous
        VortexStepMethod.unrefined_deform!(wing, nothing, delta_input; smooth=true)
        delta_smooth = copy(wing.delta_dist)

        # Verify gradients between adjacent panels are small
        max_gradient_smooth = maximum(abs.(diff(delta_smooth)))
        @test max_gradient_smooth < deg2rad(3.0)  # Should be smooth

        # Verify no sharp discontinuities
        for i in 1:(wing.n_panels-1)
            @test abs(delta_smooth[i+1] - delta_smooth[i]) < deg2rad(3.0)
        end

        # Test 3: Smoothing reduces maximum gradient
        @test max_gradient_smooth < max_gradient_nonsmooth

        # Test 4: Angles match input at unrefined section centers (both modes)
        # For non-smooth: extract angle at center panel of each unrefined section
        VortexStepMethod.unrefined_deform!(wing, nothing, delta_input; smooth=false)
        for i in 1:wing.n_unrefined_sections
            # Find panels belonging to this unrefined section
            panel_indices = findall(==(i), wing.refined_panel_mapping)
            center_panel_idx = panel_indices[div(length(panel_indices), 2) + 1]
            @test wing.delta_dist[center_panel_idx] ≈ delta_input[i] atol=1e-10
        end

        # For smooth: angles at center should be close to input (tolerance larger due to smoothing)
        VortexStepMethod.unrefined_deform!(wing, nothing, delta_input; smooth=true)
        for i in 1:wing.n_unrefined_sections
            panel_indices = findall(==(i), wing.refined_panel_mapping)
            center_panel_idx = panel_indices[div(length(panel_indices), 2) + 1]
            # Smoothing may shift values slightly, use absolute tolerance for small angles
            @test wing.delta_dist[center_panel_idx] ≈ delta_input[i] atol=deg2rad(2.0)
        end
    end
end

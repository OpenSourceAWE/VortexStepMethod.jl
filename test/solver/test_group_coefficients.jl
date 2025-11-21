using VortexStepMethod
using LinearAlgebra
using Test

@testset "Unrefined Coefficient Arrays Tests" begin
    @testset "Unrefined coefficients aggregation" begin
        # Create a simple wing with unrefined sections
        n_panels = 20
        n_unrefined_sections = 5  # This gives 4 unrefined panels

        # Create a test wing settings file
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)

        try
            # Modify settings to use specific panel configuration
            settings = VSMSettings(settings_file)
            settings.wings[1].n_panels = n_panels
            settings.solver_settings.n_panels = n_panels

            # Create wing and solver
            wing = Wing(settings)
            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            # Set conditions and solve
            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)
            sol = solve!(solver, body_aero)

            n_unrefined_panels = wing.n_unrefined_sections - 1

            # Test 1: Unrefined arrays exist and have correct size
            @test length(sol.cl_unrefined_array) == n_unrefined_panels
            @test length(sol.cd_unrefined_array) == n_unrefined_panels
            @test length(sol.cm_unrefined_array) == n_unrefined_panels

            # Test 2: Unrefined arrays are not all zeros (solver computed them)
            @test !all(sol.cl_unrefined_array .== 0.0)
            @test !all(sol.cd_unrefined_array .== 0.0)

            # Test 3: Verify unrefined coefficients are aggregated from refined panels
            # using refined_panel_mapping
            for unrefined_idx in 1:n_unrefined_panels
                # Find all refined panels that map to this unrefined panel
                refined_panel_indices = findall(x -> x == unrefined_idx, wing.refined_panel_mapping)

                if !isempty(refined_panel_indices)
                    # Calculate expected average from refined panel coefficients
                    expected_cl = sum(sol.cl_array[refined_panel_indices]) / length(refined_panel_indices)
                    expected_cd = sum(sol.cd_array[refined_panel_indices]) / length(refined_panel_indices)
                    expected_cm = sum(sol.cm_array[refined_panel_indices]) / length(refined_panel_indices)

                    # Check if unrefined coefficients match expected averages
                    # Handle NaN values that can occur in INVISCID models
                    if isnan(expected_cl)
                        @test isnan(sol.cl_unrefined_array[unrefined_idx])
                    else
                        @test isapprox(sol.cl_unrefined_array[unrefined_idx], expected_cl, rtol=1e-10)
                    end
                    if isnan(expected_cd)
                        @test isnan(sol.cd_unrefined_array[unrefined_idx])
                    else
                        @test isapprox(sol.cd_unrefined_array[unrefined_idx], expected_cd, rtol=1e-10)
                    end
                    if isnan(expected_cm)
                        @test isnan(sol.cm_unrefined_array[unrefined_idx])
                    else
                        @test isapprox(sol.cm_unrefined_array[unrefined_idx], expected_cm, rtol=1e-10)
                    end
                end
            end

            # Test 4: Verify physical consistency (lift coefficients should be positive at positive AoA)
            # Skip test if values are NaN
            if !any(isnan.(sol.cl_unrefined_array))
                @test all(sol.cl_unrefined_array .> 0.0)
            end

        finally
            rm(settings_file; force=true)
        end
    end

    @testset "Unrefined coefficients with different panel counts" begin
        # Test with various panel/section combinations
        test_cases = [
            (n_panels=40, n_unrefined_expected=21),  # From YAML file sections
            (n_panels=30, n_unrefined_expected=21),
            (n_panels=24, n_unrefined_expected=21),
        ]

        for (n_panels, n_unrefined_expected) in test_cases
            settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)

            try
                settings = VSMSettings(settings_file)
                settings.wings[1].n_panels = n_panels
                settings.solver_settings.n_panels = n_panels

                wing = Wing(settings)
                body_aero = BodyAerodynamics([wing])
                solver = Solver(body_aero, settings)

                va = [10.0, 0.0, 0.0]
                set_va!(body_aero, va)
                sol = solve!(solver, body_aero)

                n_unrefined_panels = wing.n_unrefined_sections - 1

                # Verify arrays have correct size
                @test length(sol.cl_unrefined_array) == n_unrefined_panels
                @test length(sol.cd_unrefined_array) == n_unrefined_panels
                @test length(sol.cm_unrefined_array) == n_unrefined_panels

                # Verify unrefined coefficients are computed correctly using mapping
                for unrefined_idx in 1:n_unrefined_panels
                    refined_panel_indices = findall(x -> x == unrefined_idx, wing.refined_panel_mapping)

                    if !isempty(refined_panel_indices)
                        expected_cl = sum(sol.cl_array[refined_panel_indices]) / length(refined_panel_indices)
                        expected_cd = sum(sol.cd_array[refined_panel_indices]) / length(refined_panel_indices)
                        expected_cm = sum(sol.cm_array[refined_panel_indices]) / length(refined_panel_indices)

                        # Handle NaN for all coefficients
                        if isnan(expected_cl)
                            @test isnan(sol.cl_unrefined_array[unrefined_idx])
                        else
                            @test isapprox(sol.cl_unrefined_array[unrefined_idx], expected_cl, rtol=1e-10)
                        end
                        if isnan(expected_cd)
                            @test isnan(sol.cd_unrefined_array[unrefined_idx])
                        else
                            @test isapprox(sol.cd_unrefined_array[unrefined_idx], expected_cd, rtol=1e-10)
                        end
                        if isnan(expected_cm)
                            @test isnan(sol.cm_unrefined_array[unrefined_idx])
                        else
                            @test isapprox(sol.cm_unrefined_array[unrefined_idx], expected_cm, rtol=1e-10)
                        end
                    end
                end

            finally
                rm(settings_file; force=true)
            end
        end
    end
end

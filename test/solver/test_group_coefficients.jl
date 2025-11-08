using VortexStepMethod
using LinearAlgebra
using Test

@testset "Group Coefficient Arrays Tests" begin
    @testset "Group coefficients with EQUAL_SIZE method" begin
        # Create a simple wing with groups
        n_panels = 20
        n_groups = 4

        # Create a test wing settings file
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)

        try
            # Modify settings to use specific panel/group configuration
            settings = VSMSettings(settings_file)
            settings.wings[1].n_panels = n_panels
            settings.wings[1].n_groups = n_groups
            settings.wings[1].grouping_method = EQUAL_SIZE
            settings.solver_settings.n_panels = n_panels
            settings.solver_settings.n_groups = n_groups

            # Create wing and solver
            wing = Wing(settings)
            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            # Set conditions and solve
            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)
            sol = solve!(solver, body_aero)

            # Test 1: Group arrays exist and have correct size
            @test length(sol.cl_group_array) == n_groups
            @test length(sol.cd_group_array) == n_groups
            @test length(sol.cm_group_array) == n_groups

            # Test 2: Group arrays are not all zeros (solver computed them)
            @test !all(sol.cl_group_array .== 0.0)
            @test !all(sol.cd_group_array .== 0.0)

            # Test 3: Verify group coefficients are averages of panel coefficients
            panels_per_group = n_panels ÷ n_groups
            for group_idx in 1:n_groups
                panel_start = (group_idx - 1) * panels_per_group + 1
                panel_end = group_idx * panels_per_group

                # Calculate expected average from panel coefficients
                expected_cl = sum(sol.cl_array[panel_start:panel_end]) / panels_per_group
                expected_cd = sum(sol.cd_array[panel_start:panel_end]) / panels_per_group
                expected_cm = sum(sol.cm_array[panel_start:panel_end]) / panels_per_group

                # Check if group coefficients match expected averages
                # Handle NaN values that can occur in INVISCID models
                if isnan(expected_cl)
                    @test isnan(sol.cl_group_array[group_idx])
                else
                    @test isapprox(sol.cl_group_array[group_idx], expected_cl, rtol=1e-10)
                end
                if isnan(expected_cd)
                    @test isnan(sol.cd_group_array[group_idx])
                else
                    @test isapprox(sol.cd_group_array[group_idx], expected_cd, rtol=1e-10)
                end
                if isnan(expected_cm)
                    @test isnan(sol.cm_group_array[group_idx])
                else
                    @test isapprox(sol.cm_group_array[group_idx], expected_cm, rtol=1e-10)
                end
            end

            # Test 4: Verify physical consistency (lift coefficients should be positive at positive AoA)
            # Skip test if values are NaN
            if !any(isnan.(sol.cl_group_array))
                @test all(sol.cl_group_array .> 0.0)
            end

        finally
            rm(settings_file; force=true)
        end
    end

    @testset "Group coefficients with n_groups=0 (no grouping)" begin
        # Create a wing with no groups
        n_panels = 20
        n_groups = 0

        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)

        try
            settings = VSMSettings(settings_file)
            settings.wings[1].n_panels = n_panels
            settings.wings[1].n_groups = n_groups
            settings.solver_settings.n_panels = n_panels
            settings.solver_settings.n_groups = n_groups

            wing = Wing(settings)
            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)
            sol = solve!(solver, body_aero)

            # Test: Group arrays should be empty when n_groups=0
            @test length(sol.cl_group_array) == 0
            @test length(sol.cd_group_array) == 0
            @test length(sol.cm_group_array) == 0

        finally
            rm(settings_file; force=true)
        end
    end

    @testset "Group coefficients with different group sizes" begin
        # Test with various panel/group combinations
        test_cases = [
            (n_panels=40, n_groups=8),
            (n_panels=30, n_groups=5),
            (n_panels=24, n_groups=6),
        ]

        for (n_panels, n_groups) in test_cases
            settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)

            try
                settings = VSMSettings(settings_file)
                settings.wings[1].n_panels = n_panels
                settings.wings[1].n_groups = n_groups
                settings.wings[1].grouping_method = EQUAL_SIZE
                settings.solver_settings.n_panels = n_panels
                settings.solver_settings.n_groups = n_groups

                wing = Wing(settings)
                body_aero = BodyAerodynamics([wing])
                solver = Solver(body_aero, settings)

                va = [10.0, 0.0, 0.0]
                set_va!(body_aero, va)
                sol = solve!(solver, body_aero)

                # Verify arrays have correct size
                @test length(sol.cl_group_array) == n_groups
                @test length(sol.cd_group_array) == n_groups
                @test length(sol.cm_group_array) == n_groups

                # Verify group coefficients are computed correctly
                panels_per_group = n_panels ÷ n_groups
                for group_idx in 1:n_groups
                    panel_start = (group_idx - 1) * panels_per_group + 1
                    panel_end = group_idx * panels_per_group

                    expected_cl = sum(sol.cl_array[panel_start:panel_end]) / panels_per_group
                    expected_cd = sum(sol.cd_array[panel_start:panel_end]) / panels_per_group
                    expected_cm = sum(sol.cm_array[panel_start:panel_end]) / panels_per_group

                    # Handle NaN for all coefficients
                    if isnan(expected_cl)
                        @test isnan(sol.cl_group_array[group_idx])
                    else
                        @test isapprox(sol.cl_group_array[group_idx], expected_cl, rtol=1e-10)
                    end
                    if isnan(expected_cd)
                        @test isnan(sol.cd_group_array[group_idx])
                    else
                        @test isapprox(sol.cd_group_array[group_idx], expected_cd, rtol=1e-10)
                    end
                    if isnan(expected_cm)
                        @test isnan(sol.cm_group_array[group_idx])
                    else
                        @test isapprox(sol.cm_group_array[group_idx], expected_cm, rtol=1e-10)
                    end
                end

            finally
                rm(settings_file; force=true)
            end
        end
    end
end

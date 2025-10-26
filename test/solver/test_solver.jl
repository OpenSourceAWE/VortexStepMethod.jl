using VortexStepMethod
using LinearAlgebra
using Test

@testset "Solver Constructor Tests" begin
    @testset "Solver Constructor with VSMSettings" begin
        # Use module-specific test data files
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml"; alpha=5.0, beta=0.0, wind_speed=10.0)
        
        try
            # Test Solver constructor with VSMSettings
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)
            
            # Verify solver properties match settings
            @test solver.aerodynamic_model_type == VSM
            @test solver.density == 1.225
            
            # Test that the solver can solve 
            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)
            sol = solve!(solver, body_aero)
            @test sol isa VSMSolution
            
        finally
            # Cleanup
            rm(settings_file; force=true)
        end
    end

    @testset "Solver with n_groups=0" begin
        # Test that solver works correctly when n_groups=0 (no group functionality)
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml";
            alpha=5.0, beta=0.0, wind_speed=10.0, n_groups=0)

        try
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            @test wing.n_groups == 0

            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            # Verify solver has zero groups
            @test length(solver.sol.group_moment_dist) == 0
            @test length(solver.sol.group_moment_coeff_dist) == 0

            # Test that solve! works without errors
            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)
            sol = solve!(solver, body_aero)

            @test sol isa VSMSolution
            @test sol.solver_status == FEASIBLE

            # Verify that group moments are empty
            @test length(sol.group_moment_dist) == 0
            @test length(sol.group_moment_coeff_dist) == 0

            # But force and moment should still be computed
            @test !all(sol.force .== 0.0)
            @test norm(sol.force) > 0

        finally
            rm(settings_file; force=true)
        end
    end

    @testset "Linearize with n_groups=0" begin
        # Test that linearize works correctly when n_groups=0
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml";
            alpha=5.0, beta=0.0, wind_speed=10.0, n_groups=0, n_panels=4)

        try
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            @test wing.n_groups == 0

            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            # Set velocity
            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)

            # Test linearize with only velocity (no theta or delta since n_groups=0)
            y = va  # Only velocity in input vector
            jac, results = VortexStepMethod.linearize(
                solver,
                body_aero,
                y;
                theta_idxs=nothing,
                delta_idxs=nothing,
                va_idxs=1:3,
                omega_idxs=nothing
            )

            # Results should only have 6 elements (force + moment, no group moments)
            @test length(results) == 6
            @test size(jac) == (6, 3)  # 6 outputs, 3 inputs (vx, vy, vz)

            # Verify forces are non-zero
            @test norm(results[1:3]) > 0

            # Test that using theta_idxs with n_groups=0 throws an error
            @test_throws ArgumentError VortexStepMethod.linearize(
                solver,
                body_aero,
                [0.0, 10.0, 0.0, 0.0];  # Invalid: trying to use theta
                theta_idxs=1:1,
                va_idxs=2:4
            )

            # Test that using delta_idxs with n_groups=0 throws an error
            @test_throws ArgumentError VortexStepMethod.linearize(
                solver,
                body_aero,
                [0.0, 10.0, 0.0, 0.0];  # Invalid: trying to use delta
                delta_idxs=1:1,
                va_idxs=2:4
            )

        finally
            rm(settings_file; force=true)
        end
    end

    @testset "Linearize theta_idxs validation" begin
        # Test that theta_idxs length must match n_groups
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml";
            alpha=5.0, beta=0.0, wind_speed=10.0, n_groups=2, n_panels=4)

        try
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            @test wing.n_groups == 2

            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)

            # Test with correct number of theta angles (2)
            y_correct = [0.0, 0.0, 10.0, 0.0, 0.0]  # 2 theta + 3 va
            jac, results = VortexStepMethod.linearize(
                solver,
                body_aero,
                y_correct;
                theta_idxs=1:2,
                va_idxs=3:5
            )
            @test size(jac, 1) == 8  # 6 + 2 group moments

            # Test with wrong number of theta angles (should throw error)
            y_wrong = [0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0]  # 4 theta + 3 va
            @test_throws ArgumentError VortexStepMethod.linearize(
                solver,
                body_aero,
                y_wrong;
                theta_idxs=1:4,  # Wrong: 4 angles but only 2 groups
                va_idxs=5:7
            )

            # Test with wrong number of delta angles
            @test_throws ArgumentError VortexStepMethod.linearize(
                solver,
                body_aero,
                y_wrong;
                delta_idxs=1:4,  # Wrong: 4 angles but only 2 groups
                va_idxs=5:7
            )

        finally
            rm(settings_file; force=true)
        end
    end

    @testset "Wing type cannot deform" begin
        # Test that regular Wing type (YAML-based) cannot use group_deform!
        settings_file = create_temp_wing_settings("solver", "solver_test_wing.yaml";
            alpha=5.0, beta=0.0, wind_speed=10.0, n_groups=2, n_panels=4)

        try
            settings = VSMSettings(settings_file)
            wing = Wing(settings)
            @test wing isa Wing
            @test wing.n_groups == 2

            body_aero = BodyAerodynamics([wing])
            solver = Solver(body_aero, settings)

            va = [10.0, 0.0, 0.0]
            set_va!(body_aero, va)

            # Test that trying to use theta angles with Wing throws an error
            y = [0.0, 0.0, 10.0, 0.0, 0.0]  # 2 theta + 3 va
            @test_throws ArgumentError VortexStepMethod.linearize(
                solver,
                body_aero,
                y;
                theta_idxs=1:2,
                va_idxs=3:5
            )

            # But linearize should work fine with only velocity (no deformation)
            y_velocity_only = [10.0, 0.0, 0.0]
            jac, results = VortexStepMethod.linearize(
                solver,
                body_aero,
                y_velocity_only;
                theta_idxs=nothing,
                va_idxs=1:3
            )
            @test size(jac, 1) == 8  # 6 + 2 group moments

        finally
            rm(settings_file; force=true)
        end
    end
end

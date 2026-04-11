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
            refine!(wing)
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
end

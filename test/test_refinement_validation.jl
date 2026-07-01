using VortexStepMethod
using Test

@testset "Refinement Validation" begin
    @testset "Error when refinement forgotten" begin
        wing = Wing(20; spanwise_distribution=LINEAR)
        add_section!(wing, [0,10,0], [1,10,0], INVISCID)
        add_section!(wing, [0,-10,0], [1,-10,0], INVISCID)

        # Should error without refinement
        @test_throws ArgumentError BodyAerodynamics([wing])

        # Should work after refinement
        refine!(wing)
        body_aero = BodyAerodynamics([wing])
        @test length(body_aero.panels) == 20
    end

    @testset "Multiple refine! calls work correctly" begin
        wing = Wing(20; spanwise_distribution=LINEAR)
        add_section!(wing, [0,10,0], [1,10,0], INVISCID)
        add_section!(wing, [0,-10,0], [1,-10,0], INVISCID)

        # Multiple calls should work (not idempotent - re-refines each time)
        refine!(wing)
        n1 = length(wing.refined_sections)

        refine!(wing)
        n2 = length(wing.refined_sections)

        @test n1 == n2 == 21
    end

    @testset "YAML wing deformation support after refinement" begin
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        wing = Wing(simple_wing_file; n_panels=4)

        # After refinement, should have non_deformed_sections
        refine!(wing)
        @test !isempty(wing.non_deformed_sections)
        @test length(wing.non_deformed_sections) == 5

        # unrefined_deform! should work
        @test_nowarn VortexStepMethod.unrefined_deform!(wing, [0.1, 0.2], nothing)
    end

    @testset "Converted-wing deformation support" begin
        wing = ram_air_matrix_wing(; n_panels=20, n_sections=2)

        # ObjWing creates non_deformed_sections in constructor (no refine! needed)
        @test !isempty(wing.non_deformed_sections)
        @test length(wing.non_deformed_sections) == 21

        # unrefined_deform! should work
        @test_nowarn VortexStepMethod.unrefined_deform!(wing, deg2rad.([5.0, -5.0]), nothing)
    end

    @testset "Refinement populates all required fields" begin
        wing = Wing(10; spanwise_distribution=COSINE)
        add_section!(wing, [0,5,0], [1,5,0], INVISCID)
        add_section!(wing, [0,-5,0], [1,-5,0], INVISCID)

        refine!(wing)

        # Check refined_sections populated
        @test length(wing.refined_sections) == 11

        # Check non_deformed_sections populated
        @test length(wing.non_deformed_sections) == 11

        # Check theta_dist and delta_dist resized
        @test length(wing.theta_dist) == 10
        @test length(wing.delta_dist) == 10

        # Check all zeros initially
        @test all(wing.theta_dist .== 0.0)
        @test all(wing.delta_dist .== 0.0)
    end
end

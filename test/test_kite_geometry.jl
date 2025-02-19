
using Test
using VortexStepMethod
using VortexStepMethod: create_interpolations, find_circle_center_and_radius, calculate_inertia_tensor, calculate_com, read_faces
using LinearAlgebra

@testset "KiteGeometry Tests" begin
    # Test data
    test_obj_path = joinpath(@__DIR__, "data", "test.obj")
    
    @testset "OBJ File Reading" begin
        # Create minimal test OBJ file
        test_vertices = """
        v 0.0 0.0 0.0
        v 1.0 0.0 0.0
        v 0.0 1.0 0.0
        f 1 2 3
        """
        write(test_obj_path, test_vertices)
        
        vertices, faces = read_faces(test_obj_path)
        
        @test length(vertices) == 3
        @test length(faces) == 1
        @test vertices[1] ≈ [0.0, 0.0, 0.0]
        @test vertices[2] ≈ [1.0, 0.0, 0.0]
        @test vertices[3] ≈ [0.0, 1.0, 0.0]
        @test faces[1] == [1, 2, 3]
    end
    
    @testset "Center of Mass Calculation" begin
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        faces = [[1, 2, 3]]
        
        com = calculate_com(vertices, faces)
        expected_com = [1/3, 1/3, 0.0]
        
        @test isapprox(com, expected_com, rtol=1e-5)
    end
    
    @testset "Inertia Tensor Calculation" begin
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        faces = [[1, 2, 3]]
        mass = 1.0
        com = [1/3, 1/3, 0.0]
        
        I = calculate_inertia_tensor(vertices, faces, mass, com)
        
        # Test properties of inertia tensor
        @test size(I) == (3,3)
        @test isapprox(I, I', rtol=1e-10)  # Symmetric
        @test all(diag(I) .≥ 0)  # Non-negative diagonal
    end
    
    @testset "Circle Fitting" begin
        # Create simple curved wing vertices
        r = 5.0
        z_center = 2.0
        vertices = []
        for θ in range(-π/4, π/4, length=100)
            push!(vertices, [0.0, r*sin(θ), z_center + r*cos(θ)])
        end
        
        z, radius, gamma_tip = find_circle_center_and_radius(vertices)
        
        @test isapprox(z, z_center, rtol=1e-2)
        @test isapprox(radius, r, rtol=1e-2)
        @test gamma_tip ≈ -π/4 rtol=1e-2
    end
    
    @testset "Interpolation Creation" begin
        vertices = []
        r = 5.0
        z_center = 2.0
        for θ in range(-π/4, π/4, length=10)
            push!(vertices, [0.0, r*sin(θ), z_center + r*cos(θ)])
            push!(vertices, [1.0, r*sin(θ), z_center + r*cos(θ)])
        end
        
        le_interp, te_interp, area_interp, gammas, max_xs, min_xs = 
            create_interpolations(vertices, z_center, r, π/4)
            
        # Test interpolation at middle point
        @test le_interp(0.0) ≈ 0.0 rtol=1e-2
        @test te_interp(0.0) ≈ 1.0 rtol=1e-2
    end
    
    @testset "KiteWing Construction" begin
        # Create minimal test wing
        wing = KiteWing(test_obj_path)
        
        @test wing.n_panels == 54  # Default value
        @test wing.spanwise_panel_distribution == "linear"
        @test wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
        @test isempty(wing.sections)
        @test wing.mass ≈ 1.0
        @test length(wing.center_of_mass) == 3
        @test typeof(wing.le_interp) <: Function
        @test typeof(wing.te_interp) <: Function
        @test typeof(wing.area_interp) <: Function
    end
    
    @testset "Section Addition" begin
        wing = KiteWing(test_obj_path)
        gamma = 0.0
        aero_input = "inviscid"
        
        add_section!(wing, gamma, aero_input)
        
        @test length(wing.sections) == 1
        @test wing.sections[1].aero_input == aero_input
    end

    rm(test_obj_path)
end
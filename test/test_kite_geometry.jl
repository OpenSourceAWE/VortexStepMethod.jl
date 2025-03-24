
using Test
using VortexStepMethod
using VortexStepMethod: create_interpolations, find_circle_center_and_radius, calculate_inertia_tensor, center_to_com!, read_faces, calc_inertia_y_rotation
using LinearAlgebra
using Interpolations
using Serialization

@testset "Kite Geometry Tests" begin
    # Test data
    test_obj_path = joinpath(tempdir(), "test.obj")
    test_dat_path = joinpath(tempdir(), "test.dat")
    
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
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        faces = [[1, 2, 3]]
        
        com = center_to_com!(vertices, faces)
        expected_com = [1/3, 0.0, 1/3]
        
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
        @test gamma_tip ≈ π/4 rtol=1e-2
    end
    
    r = 5.0
    @testset "Interpolation Creation" begin
        vertices = []
        z_center = 2.0
        Δθ = π/2 / 1000
        for θ in range(-π/4, π/4-Δθ, 1000)
            frac = 1.0 - abs(4θ/π)  # Goes from 0 to 1 to 0
            x_le = 0.1 - 0.1 * frac      # Goes from 0.1 to 0.0 to 0.1
            x_te = 0.9 + 0.1 * frac # Goes from 0.9 to 1.0 to 0.9
        
            push!(vertices, [x_le, r*sin(θ), z_center + r*cos(θ)])
            push!(vertices, [x_te, r*sin(θ+0.5Δθ), z_center + r*cos(θ+0.5Δθ)])
            push!(vertices, [x_le, r*sin(θ+Δθ), z_center + r*cos(θ+Δθ)])
        end
        
        # Create test airfoil data file
        test_dat_path = joinpath(tempdir(), "test.dat")
        write(test_dat_path, "1.0 0.0\n0.0 0.0\n-1.0 0.0\n")
        
        # Create polar data
        alphas = -1.0:1.0:1.0
        d_trailing_edge_angles = -1.0:1.0:1.0
        cl_matrix = zeros(length(alphas), length(d_trailing_edge_angles))
        cd_matrix = zeros(length(alphas), length(d_trailing_edge_angles))
        cm_matrix = zeros(length(alphas), length(d_trailing_edge_angles))
        
        # Fill matrices with sample data
        for i in eachindex(alphas)
            for j in eachindex(d_trailing_edge_angles)
                cl_matrix[i,j] = sin(deg2rad(alphas[i]))
                cd_matrix[i,j] = 0.01 + 0.1*sin(deg2rad(alphas[i]))^2
                cm_matrix[i,j] = -0.1*sin(deg2rad(alphas[i]))
            end
        end
        cl_matrix[end] = NaN
        cd_matrix[end] = NaN
        cm_matrix[end] = NaN
        
        # Serialize polar data
        polar_path = test_dat_path[1:end-4] * "_polar.bin"
        serialize(polar_path, (alphas, d_trailing_edge_angles, cl_matrix, cd_matrix, cm_matrix))
        
        # Create and serialize obj file
        faces = [[i, i+1, i+2] for i in 1:3:length(vertices)-2]
        open(test_obj_path, "w") do io
            for v in vertices
                println(io, "v $(v[1]) $(v[2]) $(v[3])")
            end
            for f in faces
                println(io, "f $(f[1]) $(f[2]) $(f[3])")
            end
        end
        
        # Create info file
        info_path = test_obj_path[1:end-4] * "_info.bin"
        le_interp, te_interp, area_interp = create_interpolations(vertices, z_center, r, π/4, I(3))
        center_of_mass = center_to_com!(vertices, faces)
        inertia_tensor = calculate_inertia_tensor(vertices, faces, 1.0, zeros(3))
        
        serialize(info_path, (inertia_tensor, center_of_mass, r, π/4, 
            le_interp, te_interp, area_interp))
        
        # Test interpolation at middle point
        @test isapprox([le_interp[i](0.0) for i in 1:3], [0.0, 0.0, r+z_center], atol=0.03)
        @test isapprox([te_interp[i](0.0) for i in 1:3], [1.0, 0.0, r+z_center], atol=0.03)
    end

    @testset "Alignment to principal frame" begin
        vertices, faces = read_faces(test_obj_path)
        center_of_mass = center_to_com!(vertices, faces)
        inertia_tensor_b = calculate_inertia_tensor(vertices, faces, 1.0, zeros(3))
        inertia_tensor_p, R_b_p = calc_inertia_y_rotation(inertia_tensor_b)
        for v in vertices
            v .= R_b_p * v
        end
        inertia_tensor_b2 = calculate_inertia_tensor(vertices, faces, 1.0, zeros(3))
        inertia_tensor_p2, R_b_p2 = calc_inertia_y_rotation(inertia_tensor_b2)
        @test inertia_tensor_p ≈ inertia_tensor_p2
        @test R_b_p2 ≈ I(3)
    end
    
    @testset "RamAirWing Construction" begin
        wing = RamAirWing(test_obj_path, test_dat_path; remove_nan=true)
        
        @test wing.n_panels == 54  # Default value
        @test wing.spanwise_distribution == UNCHANGED
        @test wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
        @test length(wing.sections) > 0  # Should have sections now
        @test wing.mass ≈ 1.0
        @test wing.radius ≈ r rtol=1e-2
        @test wing.gamma_tip ≈ π/4 rtol=1e-2
        @test !isnan(wing.sections[1].aero_data[3][end])
        @test !isnan(wing.sections[1].aero_data[4][end])
        @test !isnan(wing.sections[1].aero_data[5][end])

        wing = RamAirWing(test_obj_path, test_dat_path; remove_nan=false)
        @test isnan(wing.sections[1].aero_data[3][end])
        @test isnan(wing.sections[1].aero_data[4][end])
        @test isnan(wing.sections[1].aero_data[5][end])
    end

    @testset "Wing Deformation" begin
        # Create a RamAirWing for testing
        wing = RamAirWing(test_obj_path, test_dat_path; remove_nan=true)
        body_aero = BodyAerodynamics([wing])
        
        # Store original TE point for comparison
        i = length(body_aero.panels) ÷ 2
        original_te_point = copy(body_aero.panels[i].TE_point_1)
        
        # Apply deformation with non-zero angles
        theta_dist = fill(deg2rad(30.0), wing.n_panels)  # 10 degrees twist
        delta_dist = fill(deg2rad(5.0), wing.n_panels)   # 5 degrees trailing edge deflection
        
        VortexStepMethod.deform!(wing, theta_dist, delta_dist)
        VortexStepMethod.init!(body_aero)
        
        # Check if TE point changed after deformation
        deformed_te_point = copy(body_aero.panels[i].TE_point_1)
        @test !isapprox(original_te_point, deformed_te_point, atol=1e-2)
        @test deformed_te_point[3] < original_te_point[3] # right hand rule
        @test deformed_te_point[2] ≈ original_te_point[2] atol=1e-2 # right hand rule
        @test deformed_te_point[1] < original_te_point[1] # right hand rule
        @test body_aero.panels[i].delta ≈ deg2rad(5.0)
        
        # Reset deformation with zero angles
        zero_theta_dist = zeros(wing.n_panels)
        zero_delta_dist = zeros(wing.n_panels)
        
        VortexStepMethod.deform!(wing, zero_theta_dist, zero_delta_dist)
        VortexStepMethod.init!(body_aero)
        
        # Check if TE point returned to original position
        reset_te_point = copy(body_aero.panels[i].TE_point_1)
        @test original_te_point ≈ reset_te_point atol=1e-4
    end
    
    rm(test_obj_path)
    rm(test_dat_path)
end
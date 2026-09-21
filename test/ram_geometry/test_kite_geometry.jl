
using Test
using VortexStepMethod
using VortexStepMethod: read_aero_matrix
using VortexStepMethod.AirfoilAero: write_aero_matrix
using VortexStepMethod.ObjAdapter: create_interpolations, find_circle_center_and_radius,
    read_faces
using LinearAlgebra

@testset "Kite Geometry Tests" begin
    work_dir = mktempdir()
    test_obj_path = joinpath(work_dir, "test.obj")
    test_dat_path = joinpath(work_dir, "test.dat")
    
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
        @test vertices isa Vector{Vector{Float64}}
        @test faces isa Vector{Vector{Int64}}
    end
    
    @testset "Circle Fitting" begin
        # Create simple curved wing vertices
        r = 5.0
        z_center = 2.0
        vertices = [[0.0, r*sin(θ), z_center + r*cos(θ)]
                    for θ in range(-π/4, π/4, length=100)]
        
        z, radius, gamma_tip = find_circle_center_and_radius(vertices)
        
        @test isapprox(z, z_center, rtol=1e-2)
        @test isapprox(radius, r, rtol=1e-2)
        @test gamma_tip ≈ π/4 rtol=1e-2
    end
    
    r = 5.0
    @testset "Interpolation Creation" begin
        vertices = Vector{Float64}[]
        z_center = 2.0
        Δθ = π/2 / 1000
        for θ in range(-π/4, π/4-Δθ, 1000)
            frac = 1.0 - abs(4θ/π)  # Goes from 0 to 1 to 0
            x_le = 0.1 - 0.1 * frac # Goes from 0.1 to 0.0 to 0.1
            x_te = 0.9 + 0.1 * frac # Goes from 0.9 to 1.0 to 0.9
        
            push!(vertices, [x_le, r*sin(θ), z_center + r*cos(θ)])
            push!(vertices, [x_te, r*sin(θ+0.5Δθ), z_center + r*cos(θ+0.5Δθ)])
            push!(vertices, [x_le, r*sin(θ+Δθ), z_center + r*cos(θ+Δθ)])
        end
        
        # Create test airfoil data file
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
        
        cl_polar_path = test_dat_path[1:end-4] * "_cl_polar.csv"
        cd_polar_path = test_dat_path[1:end-4] * "_cd_polar.csv"
        cm_polar_path = test_dat_path[1:end-4] * "_cm_polar.csv"
        
        # Write matrices to CSV
        write_aero_matrix(cl_polar_path, cl_matrix, deg2rad.(alphas), deg2rad.(d_trailing_edge_angles), "C_l")
        write_aero_matrix(cd_polar_path, cd_matrix, deg2rad.(alphas), deg2rad.(d_trailing_edge_angles), "C_d")
        write_aero_matrix(cm_polar_path, cm_matrix, deg2rad.(alphas), deg2rad.(d_trailing_edge_angles), "C_m")
        
        # Test reading back the matrices
        cl_read, alphas_read, deltas_read = read_aero_matrix(cl_polar_path)
        # write_aero_matrix stores coefficients rounded to 4 decimals
        @test maximum(abs.(cl_read[1:end-1,:] .- cl_matrix[1:end-1,:])) <= 5e-5
        @test isnan(cl_read[end,end])
        @test alphas_read ≈ deg2rad.(alphas)
        @test deltas_read ≈ deg2rad.(d_trailing_edge_angles)
        
        le_interp, te_interp, area_interp = create_interpolations(vertices, z_center, r, π/4, I(3))

        # Test interpolation at middle point
        @test isapprox([le_interp[i](0.0) for i in 1:3], [0.0, 0.0, r+z_center], atol=0.03)
        @test isapprox([te_interp[i](0.0) for i in 1:3], [1.0, 0.0, r+z_center], atol=0.03)
    end

    @testset "Converted-wing construction and deformation" begin
        # TODO: redesign. These previously tested ObjWing internals (radius,
        # gamma_tip, UNCHANGED distribution, obj deform\!) that were dropped when
        # ObjWing was replaced by convert-then-load (obj_to_matrix_yaml -> Wing).
        # Rebuild against ram_air_matrix_wing() geometry once its numerics are set.
        @test_skip false
    end
end

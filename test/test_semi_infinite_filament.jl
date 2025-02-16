using VortexStepMethod: SemiInfiniteFilament, velocity_3D_trailing_vortex_semiinfinite!
using LinearAlgebra
using Test
using BenchmarkTools

function create_test_filament()
    x1 = [0.0, 0.0, 0.0]
    direction = [1.0, 0.0, 0.0]
    filament_direction = 1
    vel_mag = 1.0
    SemiInfiniteFilament(x1, direction, vel_mag, filament_direction)
end

function analytical_solution(control_point, gamma, x1, direction, filament_direction, vel_mag)
    gamma = -gamma  # Sign convention difference
    r1 = control_point - x1
    r1_cross_direction = cross(r1, direction)
    r_perp = dot(r1, direction) * direction
    
    alpha0 = 1.25643
    nu = 1.48e-5
    epsilon = sqrt(4 * alpha0 * nu * norm(r_perp) / vel_mag)

    if norm(r1_cross_direction) > epsilon
        K = (gamma / (4π * norm(r1_cross_direction)^2)) * 
            (1 + dot(r1, direction) / norm(r1))
        return K * r1_cross_direction * filament_direction
    else
        r1_proj = dot(r1, direction) * direction + epsilon * 
                  (r1/norm(r1) - direction) / norm(r1/norm(r1) - direction)
        r1_cross_direction_proj = cross(r1_proj, direction)
        K_proj = (gamma / (4π * norm(r1_cross_direction_proj)^2)) * 
                 (1 + dot(r1_proj, direction) / norm(r1_proj))
        return K_proj * r1_cross_direction_proj * filament_direction
    end
end

@testset "SemiInfiniteFilament Tests" begin
    gamma = 1.0
    core_radius_fraction = 0.01
    work_vectors = ntuple(_ -> Vector{Float64}(undef, 3), 10)

    @testset "Allocation Tests" begin
        filament = create_test_filament()
        control_point = [0.5, 0.5, 2.0]
        induced_velocity = zeros(3)
        
        b = @benchmarkable velocity_3D_trailing_vortex_semiinfinite!(
            $induced_velocity,
            $filament,
            $filament.direction,
            $control_point,
            $gamma,
            $filament.vel_mag,
            $work_vectors
        )
        result = run(b)
        @test result.allocs == 0
        @test result.memory == 0
    end

    @testset "Calculate Induced Velocity" begin
        filament = create_test_filament()
        control_point = [0.5, 0.5, 2.0]
        induced_velocity = zeros(3)
        
        velocity_3D_trailing_vortex_semiinfinite!(
            induced_velocity,
            filament,
            filament.direction,
            control_point,
            gamma,
            filament.vel_mag,
            work_vectors
        )
        
        analytical = analytical_solution(
            control_point, gamma, filament.x1, filament.direction,
            filament.filament_direction, filament.vel_mag  
        )
        
        @test isapprox(induced_velocity, analytical, rtol=1e-6)
    end

    @testset "Point on Filament" begin
        filament = create_test_filament()
        test_points = [
            [0.0, 0.0, 0.0],  # Start point
            [0.5, 0.0, 0.0],  # Along filament
            [5.0, 0.0, 0.0],  # Further along
        ]
        induced_velocity = zeros(3)
        
        for point in test_points
            velocity_3D_trailing_vortex_semiinfinite!(
                induced_velocity,
                filament,
                filament.direction,
                point,
                gamma,
                filament.vel_mag,
                work_vectors
            )
            @test all(isapprox.(induced_velocity, zeros(3), atol=1e-6))
        end
    end

    @testset "Different Gamma Values" begin
        filament = create_test_filament()
        control_point = [0.5, 1.0, 0.0]
        v1 = zeros(3)
        v2 = zeros(3)
        v4 = zeros(3)
        
        velocity_3D_trailing_vortex_semiinfinite!(v1, filament, filament.direction, control_point, 1.0, filament.vel_mag, work_vectors)
        velocity_3D_trailing_vortex_semiinfinite!(v2, filament, filament.direction, control_point, 2.0, filament.vel_mag, work_vectors)
        velocity_3D_trailing_vortex_semiinfinite!(v4, filament, filament.direction, control_point, 4.0, filament.vel_mag, work_vectors)
        
        @test isapprox(v4, 2 * v2)
        @test isapprox(v4, 4 * v1)
    end

    @testset "Symmetry" begin
        filament = create_test_filament()
        vel_pos = zeros(3)
        vel_neg = zeros(3)
        
        velocity_3D_trailing_vortex_semiinfinite!(
            vel_pos,
            filament,
            filament.direction,
            [0.0, 1.0, 0.0],
            gamma,
            filament.vel_mag,
            work_vectors
        )
        velocity_3D_trailing_vortex_semiinfinite!(
            vel_neg,
            filament,
            filament.direction,
            [0.0, -1.0, 0.0],
            gamma,
            filament.vel_mag,
            work_vectors
        )
        
        @test isapprox(vel_pos, -vel_neg)
    end
end
using VortexStepMethod: BoundFilament, velocity_3D_bound_vortex!
using LinearAlgebra
using BenchmarkTools: @benchmarkable
using Test

# Test helper functions
function create_test_filament()
    BoundFilament([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
end

function analytical_solution(control_point, gamma)
    A = [0.0, 0.0, 0.0]
    B = [1.0, 0.0, 0.0]
    P = control_point
    
    r0 = B - A
    r1 = P - A
    r2 = P - B
    
    r1Xr2 = cross(r1, r2)
    norm_r1Xr2 = norm(r1Xr2)
    
    return (gamma / (4π)) * (r1Xr2 / (norm_r1Xr2^2)) * 
           dot(r0, r1/norm(r1) - r2/norm(r2))
end

@testset "BoundFilament Tests" begin
    gamma = 1.0
    core_radius_fraction = 0.01
    induced_velocity = zeros(3)
    work_vectors = ntuple(_ -> Vector{Float64}(undef, 3), 10)

    @testset "Non-allocating" begin
        filament = create_test_filament()
        control_point = [0.5, 1.0, 0.0]
        induced_velocity = zeros(3)
        
        b = @benchmarkable velocity_3D_bound_vortex!(
            $induced_velocity,
            $filament,
            $control_point,
            $gamma,
            $core_radius_fraction,
            $work_vectors
        )
        result = run(b)
        @test result.allocs == 0
        @test result.memory == 0
    end

    @testset "Calculate Induced Velocity" begin
        filament = create_test_filament()
        control_point = [0.5, 1.0, 0.0]
        
        velocity_3D_bound_vortex!(
            induced_velocity,
            filament,
            control_point,
            gamma,
            core_radius_fraction,
            work_vectors
        )
        analytical_velocity = analytical_solution(control_point, gamma)
        
        @test isapprox(induced_velocity, analytical_velocity, rtol=1e-6)
    end

    @testset "Point Exactly on Filament" begin
        filament = create_test_filament()
        test_points = [
            [0.0, 0.0, 0.0],  # start point
            [1.0, 0.0, 0.0],  # end point
            [0.5, 0.0, 0.0],  # middle point
        ]
        
        for point in test_points
            velocity_3D_bound_vortex!(
                induced_velocity,
                filament,
                point,
                gamma,
                core_radius_fraction,
                work_vectors
            )
            @test all(isapprox.(induced_velocity, zeros(3), atol=1e-10))
            @test !any(isnan.(induced_velocity))
        end
    end

    @testset "Long Filament" begin
        filament = BoundFilament([0.0, 0.0, 0.0], [1e6, 0.0, 0.0])
        control_point = [5e5, 1.0, 0.0]
        
        velocity_3D_bound_vortex!(
            induced_velocity,
            filament,
            control_point,
            gamma,
            core_radius_fraction,
            work_vectors
        )
        @test !any(isnan.(induced_velocity))
        @test isapprox(induced_velocity[1], 0.0, atol=1e-8)
        @test abs(induced_velocity[2]) < 1e-8
        @test isapprox(induced_velocity[3], 0.0)
    end

    @testset "Point Far from Filament" begin
        filament = create_test_filament()
        control_point = [0.5, 1e6, 0.0]
        
        velocity_3D_bound_vortex!(
            induced_velocity,
            filament,
            control_point,
            gamma,
            core_radius_fraction,
            work_vectors
        )
        
        @test !any(isnan.(induced_velocity))
        @test all(isapprox.(induced_velocity, zeros(3), atol=1e-12))
    end

    v1, v2, v4 = zeros(3), zeros(3), zeros(3)

    @testset "Different Gamma Values" begin
        filament = create_test_filament()
        control_point = [0.5, 1.0, 0.0]
        
        velocity_3D_bound_vortex!(v1, filament, control_point, 1.0, core_radius_fraction, work_vectors)
        velocity_3D_bound_vortex!(v2, filament, control_point, 2.0, core_radius_fraction, work_vectors)
        velocity_3D_bound_vortex!(v4, filament, control_point, 4.0, core_radius_fraction, work_vectors)
        
        @test isapprox(v4, 2 * v2)
        @test isapprox(v4, 4 * v1)
    end

    @testset "Symmetry" begin
        filament = BoundFilament([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0])

        velocity_3D_bound_vortex!(v1, filament, [0.0, 1.0, 0.0], gamma, core_radius_fraction, work_vectors)
        velocity_3D_bound_vortex!(v2, filament, [0.0, -1.0, 0.0], gamma, core_radius_fraction, work_vectors)
        
        @test isapprox(v1, -v2)
    end

    @testset "Around Core Radius" begin
        filament = create_test_filament()
        delta = 1e-5
        
        points = [
            [0.5, core_radius_fraction - delta, 0.0],
            [0.5, core_radius_fraction, 0.0],
            [0.5, core_radius_fraction + delta, 0.0]
        ]
        
        velocities = [zeros(3) for p in points]
        [
            velocity_3D_bound_vortex!(velocities[i], filament, p, gamma, core_radius_fraction, work_vectors)
            for (i, p) in enumerate(points)
        ]
        
        # Check for NaN and finite values
        @test all(!any(isnan.(v)) for v in velocities)
        @test all(all(isfinite.(v)) for v in velocities)
        
        # Check magnitude is maximum at core radius
        @test norm(velocities[2]) > norm(velocities[1])
        @test norm(velocities[2]) > norm(velocities[3])
        
        # Check continuity around core radius
        @test isapprox(velocities[1], velocities[2], rtol=1e-2)
        
        # Check non-zero velocities
        @test !all(isapprox.(velocities[1], zeros(3), atol=1e-10))
        @test !all(isapprox.(velocities[2], zeros(3), atol=1e-10))
        @test !all(isapprox.(velocities[3], zeros(3), atol=1e-10))
        
        # Check symmetry
        v_neg = zeros(3)
        velocity_3D_bound_vortex!(
            v_neg,
            filament,
            [0.5, -core_radius_fraction, 0.0],
            gamma,
            core_radius_fraction,
            work_vectors
        )
        @test isapprox(velocities[2], -v_neg)
    end
end
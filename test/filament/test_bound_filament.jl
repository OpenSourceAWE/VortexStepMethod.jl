using VortexStepMethod: BoundFilament, velocity_3D_bound_vortex!, reinit!
using LinearAlgebra
using Test

# Test helper functions
function create_test_filament()
    fil = BoundFilament{Float64}()
    reinit!(fil, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
    return fil
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
        filament = BoundFilament{Float64}()
        reinit!(filament, [0.0, 0.0, 0.0], [1e6, 0.0, 0.0])
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
        @test isapprox(induced_velocity[2], 0.0, atol=1e-8)
        @test abs(induced_velocity[3]) > 1e-10
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
        filament = BoundFilament{Float64}()
        reinit!(filament, [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0])

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

    @testset "Velocity is azimuthal (perpendicular to axis and radius)" begin
        filament = create_test_filament()
        r0 = [1.0, 0.0, 0.0]

        for d in (0.25, 0.5, 0.99, 1.0, 1.01, 2.0) .* core_radius_fraction
            for phi in (0.0, π/4, π/2, π, -π/3)
                p = [0.5, d * cos(phi), d * sin(phi)]
                v = zeros(3)
                velocity_3D_bound_vortex!(
                    v, filament, p, gamma,
                    core_radius_fraction, work_vectors)

                r_radial = [0.0, p[2], p[3]]
                @test isapprox(dot(v, r0), 0.0; atol=1e-10)
                @test isapprox(dot(v, r_radial), 0.0; atol=1e-8)
            end
        end
    end

    @testset "Solid-body rotation inside core" begin
        filament = create_test_filament()

        v_half = zeros(3)
        v_quarter = zeros(3)
        velocity_3D_bound_vortex!(
            v_half, filament,
            [0.5, 0.5 * core_radius_fraction, 0.0],
            gamma, core_radius_fraction, work_vectors)
        velocity_3D_bound_vortex!(
            v_quarter, filament,
            [0.5, 0.25 * core_radius_fraction, 0.0],
            gamma, core_radius_fraction, work_vectors)

        @test isapprox(norm(v_quarter), 0.5 * norm(v_half); rtol=1e-3)
        @test isapprox(normalize(v_quarter), normalize(v_half); atol=1e-8)
    end

    @testset "Branch 1 / Branch 3 agree at core boundary" begin
        filament = create_test_filament()
        delta = 0.05 * core_radius_fraction
        v_inside = zeros(3)
        v_outside = zeros(3)
        velocity_3D_bound_vortex!(
            v_inside, filament,
            [0.5, core_radius_fraction - delta, 0.0],
            gamma, core_radius_fraction, work_vectors)
        velocity_3D_bound_vortex!(
            v_outside, filament,
            [0.5, core_radius_fraction + delta, 0.0],
            gamma, core_radius_fraction, work_vectors)

        @test isapprox(norm(v_inside), norm(v_outside); rtol=0.1)
        @test isapprox(normalize(v_inside), normalize(v_outside); atol=1e-2)
    end

    @testset "Matches Rankine vortex profile inside core" begin
        L = 1e6
        r_c = 0.01 * L
        filament_long = BoundFilament{Float64}()
        reinit!(filament_long, [0.0, 0.0, 0.0], [L, 0.0, 0.0])

        for r in (1.0, 10.0, 100.0, 1000.0)
            v = zeros(3)
            velocity_3D_bound_vortex!(
                v, filament_long, [L/2, r, 0.0],
                gamma, 0.01, work_vectors)

            expected_mag = gamma * r / (2π * r_c^2)

            @test isapprox(v[1], 0.0; atol=expected_mag * 1e-4)
            @test isapprox(v[2], 0.0; atol=expected_mag * 1e-4)
            @test isapprox(abs(v[3]), expected_mag; rtol=1e-3)
            @test v[3] > 0
        end
    end
end
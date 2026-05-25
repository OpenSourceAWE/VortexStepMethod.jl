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
        # A real vortex's induced velocity is purely azimuthal:
        # perpendicular to BOTH the filament axis r0 AND the radial
        # direction from the axis to the field point. The old Branch 3
        # projected the field point in the azimuthal direction instead
        # of the radial direction, producing a velocity with a nonzero
        # radial component. This test catches that.
        filament = create_test_filament()
        r0 = [1.0, 0.0, 0.0]

        # Probe at multiple distances spanning both branches (in-core
        # and outside-core) and azimuthal positions.
        for d in (0.25, 0.5, 0.99, 1.0, 1.01, 2.0) .* core_radius_fraction
            for phi in (0.0, π/4, π/2, π, -π/3)
                p = [0.5, d * cos(phi), d * sin(phi)]
                v = zeros(3)
                velocity_3D_bound_vortex!(
                    v, filament, p, gamma,
                    core_radius_fraction, work_vectors)

                # Radial direction at p (perpendicular to axis)
                r_radial = [0.0, p[2], p[3]]

                # v must be perpendicular to both axis and radial
                @test isapprox(dot(v, r0), 0.0; atol=1e-10)
                @test isapprox(dot(v, r_radial), 0.0; atol=1e-8)
            end
        end
    end

    @testset "Solid-body rotation inside core" begin
        # Inside the core, velocity magnitude scales linearly with
        # d_perp (Branch 3 = linear ramp from 0 on axis to peak at
        # boundary). Direction stays constant along a radial line.
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

        # Magnitudes scale linearly with d_perp
        @test isapprox(norm(v_quarter), 0.5 * norm(v_half); rtol=1e-3)
        # Directions are parallel
        @test isapprox(normalize(v_quarter), normalize(v_half); atol=1e-8)
    end

    @testset "Branch 1 / Branch 3 agree at core boundary" begin
        # The two branches are designed to be continuous at d_perp =
        # epsilon. Probe just inside and just outside the boundary with
        # enough delta to clearly land in different branches; results
        # must agree in both magnitude AND direction. Old code returned
        # orthogonal directions here.
        filament = create_test_filament()
        delta = 0.05 * core_radius_fraction  # 5% in/out
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

        # Magnitudes within a few % of each other (peak is at boundary)
        @test isapprox(norm(v_inside), norm(v_outside); rtol=0.1)
        # Directions match — would fail with old azimuthal-projection bug
        @test isapprox(normalize(v_inside), normalize(v_outside); atol=1e-2)
    end
end
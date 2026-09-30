using VortexStepMethod: SemiInfiniteFilament, velocity_3D_trailing_vortex_semiinfinite!,
    reinit!, ALPHA0, NU
using ForwardDiff
using LinearAlgebra
using Test

function create_test_filament2()
    x1 = [0.0, 0.0, 0.0]
    direction = [1.0, 0.0, 0.0]
    filament_direction = 1
    va = 1.0
    filament = SemiInfiniteFilament{Float64}()
    reinit!(filament, x1, direction, va, filament_direction)
    return filament
end

function analytical_solution(control_point, gamma, x1, direction, filament_direction)
    gamma = -gamma  # Sign convention difference
    r1 = control_point - x1
    r1_cross_direction = cross(r1, direction)
    K = (gamma / (4π * norm(r1_cross_direction)^2)) * (1 + dot(r1, direction) / norm(r1))
    return K * r1_cross_direction * filament_direction
end

"""
    core_radius(axial_distance, va)

Lamb–Oseen core radius [m] of a trailing vortex `axial_distance` [m] downstream of its
start.
"""
core_radius(axial_distance, va) = sqrt(4 * ALPHA0 * NU * axial_distance / va)

"""
    off_axis_velocity(offset, gamma)

z velocity [m/s] induced by the unit-speed trailing filament from the origin along x at
the point `offset` [m] off its axis, half a metre downstream.
"""
function off_axis_velocity(offset, gamma)
    T = typeof(offset)
    filament = SemiInfiniteFilament{T}()
    reinit!(filament, zeros(T, 3), T[1, 0, 0], one(T), 1)
    velocity = zeros(T, 3)
    velocity_3D_trailing_vortex_semiinfinite!(velocity, filament, filament.direction,
        [0.5, offset, 0.0], gamma, filament.va, ntuple(_ -> zeros(T, 3), 10))
    return velocity[3]
end

@testset "SemiInfiniteFilament Tests" begin
    gamma = 1.0
    core_radius_fraction = 0.01
    work_vectors = ntuple(_ -> Vector{Float64}(undef, 3), 10)
    @testset "Calculate Induced Velocity" begin
        filament = create_test_filament2()
        control_point = [0.5, 0.5, 2.0]
        induced_velocity = zeros(3)
        
        velocity_3D_trailing_vortex_semiinfinite!(
            induced_velocity,
            filament,
            filament.direction,
            control_point,
            gamma,
            filament.va,
            work_vectors
        )
        
        analytical = analytical_solution(
            control_point, gamma, filament.x1, filament.direction,
            filament.filament_direction
        )
        
        @test isapprox(induced_velocity, analytical, rtol=1e-6)
    end

    @testset "Point on Filament" begin
        filament = create_test_filament2()
        start_point = [0.0, 0.0, 0.0]
        test_points = [
            [0.5, 0.0, 0.0],  # Along filament
            [5.0, 0.0, 0.0],  # Further along
        ]
        induced_velocity = zeros(3)

        velocity_3D_trailing_vortex_semiinfinite!(
            induced_velocity,
            filament,
            filament.direction,
            start_point,
            gamma,
            filament.va,
            work_vectors
        )
        @test induced_velocity == zeros(3)

        for point in test_points
            velocity_3D_trailing_vortex_semiinfinite!(
                induced_velocity,
                filament,
                filament.direction,
                point,
                gamma,
                filament.va,
                work_vectors
            )
            @test all(isapprox.(induced_velocity, zeros(3), atol=1e-5))
        end
    end

    @testset "Different Gamma Values" begin
        filament = create_test_filament2()
        control_point = [0.5, 1.0, 0.0]
        v1 = zeros(3)
        v2 = zeros(3)
        v4 = zeros(3)
        
        velocity_3D_trailing_vortex_semiinfinite!(v1, filament, filament.direction,
            control_point, 1.0, filament.va, work_vectors)
        velocity_3D_trailing_vortex_semiinfinite!(v2, filament, filament.direction,
            control_point, 2.0, filament.va, work_vectors)
        velocity_3D_trailing_vortex_semiinfinite!(v4, filament, filament.direction,
            control_point, 4.0, filament.va, work_vectors)
        
        @test isapprox(v4, 2 * v2)
        @test isapprox(v4, 4 * v1)
    end

    @testset "Symmetry" begin
        filament = create_test_filament2()
        vel_pos = zeros(3)
        vel_neg = zeros(3)

        velocity_3D_trailing_vortex_semiinfinite!(
            vel_pos,
            filament,
            filament.direction,
            [0.0, 1.0, 0.0],
            gamma,
            filament.va,
            work_vectors
        )
        velocity_3D_trailing_vortex_semiinfinite!(
            vel_neg,
            filament,
            filament.direction,
            [0.0, -1.0, 0.0],
            gamma,
            filament.va,
            work_vectors
        )

        @test isapprox(vel_pos, -vel_neg)
    end

    @testset "Velocity is azimuthal (perpendicular to axis and radius)" begin
        filament = create_test_filament2()
        direction = filament.direction

        for d in (1e-4, 1e-3, 5e-3, 1e-2, 1e-1)
            for phi in (0.0, π/4, π/2, π, -π/3)
                p = [0.5, d * cos(phi), d * sin(phi)]
                v = zeros(3)
                velocity_3D_trailing_vortex_semiinfinite!(
                    v, filament, direction, p, gamma,
                    filament.va, work_vectors)

                r_radial = [0.0, p[2], p[3]]
                @test isapprox(dot(v, direction), 0.0; atol=1e-10)
                @test isapprox(dot(v, r_radial) / norm(v), 0.0; atol=1e-6)
            end
        end
    end

    @testset "Constant azimuthal direction inside core" begin
        filament = create_test_filament2()
        va = filament.va

        d_inside = 1e-4
        v1 = zeros(3); v2 = zeros(3)
        velocity_3D_trailing_vortex_semiinfinite!(
            v1, filament, filament.direction,
            [0.5, d_inside, 0.0], gamma, va, work_vectors)
        velocity_3D_trailing_vortex_semiinfinite!(
            v2, filament, filament.direction,
            [0.5, 2 * d_inside, 0.0], gamma, va, work_vectors)

        @test isapprox(normalize(v2), normalize(v1); atol=1e-8)
    end

    @testset "Velocity scales linearly with distance inside core" begin
        epsilon = core_radius(0.5, 1.0)
        v_half = off_axis_velocity(0.5 * epsilon, gamma)
        v_quarter = off_axis_velocity(0.25 * epsilon, gamma)

        @test v_quarter ≈ 0.5 * v_half rtol = 1e-12
        @test off_axis_velocity(-0.25 * epsilon, gamma) ≈ -v_quarter rtol = 1e-12
    end

    @testset "Velocity is continuous at the core boundary" begin
        epsilon = core_radius(0.5, 1.0)
        v_inside = off_axis_velocity(epsilon * (1 - 1e-9), gamma)
        v_outside = off_axis_velocity(epsilon * (1 + 1e-9), gamma)

        @test v_inside ≈ v_outside rtol = 1e-6
    end

    @testset "ForwardDiff sees the core's slope on the axis" begin
        velocity_at(offset) = off_axis_velocity(offset, gamma)
        slope = ForwardDiff.derivative(velocity_at, 0.0)

        @test slope != 0
        @test slope ≈ velocity_at(1e-5) / 1e-5 rtol = 1e-9
    end
end

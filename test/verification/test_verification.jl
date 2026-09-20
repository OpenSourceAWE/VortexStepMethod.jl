using DelimitedFiles
using Interpolations
using LinearAlgebra
using Test
using VortexStepMethod
using VortexStepMethod: BoundFilament, SemiInfiniteFilament, reinit!,
    velocity_3D_bound_vortex!, velocity_3D_trailing_vortex!,
    velocity_3D_trailing_vortex_semiinfinite!
using VortexStepMethod.AirfoilAero: lei_poly_coeffs

include("../utils.jl")

"""
    read_columns(file_name)

Columns of the comma-separated table `file_name` in `test/verification/data`,
header row skipped.
"""
function read_columns(file_name)
    table = readdlm(joinpath(@__DIR__, "data", file_name), ',', Float64; skipstart=1)
    return eachcol(table)
end

"""
    polar_vectors(alpha, cl, cd=zero(cl), cm=zero(cl))

`POLAR_VECTORS` aero data from a section polar with `alpha` in degrees.
"""
polar_vectors(alpha, cl, cd=zero(cl), cm=zero(cl)) =
    (deg2rad.(alpha), collect(cl), collect(cd), collect(cm))

"""
    wing_from_coordinates(coordinates, aero_model, aero_data=nothing;
                          n_panels, spanwise_distribution=UNCHANGED)

Refined wing with one section per leading/trailing-edge row pair of the 2N×3
`coordinates`, which run from -y to +y.
"""
function wing_from_coordinates(coordinates, aero_model, aero_data=nothing;
                               n_panels=size(coordinates, 1) ÷ 2 - 1,
                               spanwise_distribution=UNCHANGED)
    ordered = flip_created_coord_in_pairs(coordinates)
    wing = Wing(n_panels; spanwise_distribution)
    for i in 1:2:size(ordered, 1)
        add_section!(wing, ordered[i, :], ordered[i + 1, :], aero_model, aero_data)
    end
    refine!(wing)
    return wing
end

"""
    lift_drag_polar(wing, model, alphas; wind_speed, relaxation_factor)

Wing `CL` and `CD` at each angle of attack in `alphas` [deg], solved with the
settings of the Python verification cases.
"""
function lift_drag_polar(wing, model, alphas; wind_speed, relaxation_factor)
    body_aero = BodyAerodynamics([wing])
    solver = Solver(wing.n_panels, wing.n_unrefined_sections; aerodynamic_model_type=model,
                    relaxation_factor, core_radius_fraction=1e-20)
    CL = zeros(length(alphas))
    CD = zeros(length(alphas))
    for (i, alpha) in enumerate(alphas)
        set_va!(body_aero, wind_speed .* [cosd(alpha), 0.0, sind(alpha)])
        results = solve(solver, body_aero)
        CL[i] = results["cl"]
        CD[i] = results["cd"]
    end
    return CL, CD
end

"""
    bound_filament(x1, x2)

Straight vortex filament from `x1` to `x2`.
"""
function bound_filament(x1, x2)
    filament = BoundFilament{Float64}()
    reinit!(filament, x1, x2)
    return filament
end

max_error(actual, expected) = maximum(abs.(actual .- expected))

@testset "Python verification cases" begin
    @testset "elliptic wing matches lifting-line theory" begin
        span = 15.709
        max_chord = 1.0
        aspect_ratio = span^2 / (π * span * max_chord / 4)
        alphas = [3.0, 9.0]
        coordinates = generate_coordinates_el_wing(max_chord, span, 40, "cos")
        wing = wing_from_coordinates(coordinates, INVISCID)
        CL_llt, CD_llt = lift_drag_polar(wing, LLT, alphas; wind_speed=20.0,
                                         relaxation_factor=0.05)
        CL_vsm, CD_vsm = lift_drag_polar(wing, VSM, alphas; wind_speed=20.0,
                                         relaxation_factor=0.05)
        CL_theory = 2π .* deg2rad.(alphas) ./ (1 + 2 / aspect_ratio)
        CD_theory = CL_theory .^ 2 ./ (π * aspect_ratio)

        @test max_error(CL_llt, CL_theory) < 1e-2
        @test max_error(CD_llt, CD_theory) < 1e-3
        @test max_error(CL_vsm, CL_theory) < 1e-1
        @test max_error(CD_vsm, CD_theory) < 2e-3
        @test max_error(CL_llt, CL_vsm) < 1e-1
        @test max_error(CD_llt, CD_vsm) < 1e-2
    end

    @testset "curved Clark Y wing matches RANS" begin
        alphas = [9.0]
        coordinates = generate_coordinates_curved_wing(2.18, 6.969, π / 4, 4.673, 60,
                                                       "lin")
        polar = polar_vectors(read_columns("clarky_polar.csv")...)
        wing = wing_from_coordinates(coordinates, POLAR_VECTORS, polar)
        CL_llt, CD_llt = lift_drag_polar(wing, LLT, alphas; wind_speed=20.0,
                                         relaxation_factor=0.03)
        CL_vsm, CD_vsm = lift_drag_polar(wing, VSM, alphas; wind_speed=20.0,
                                         relaxation_factor=0.03)
        alpha_rans, CL_rans, CD_rans, _ = read_columns("curved_wing_rans.csv")

        @test max_error(CL_vsm, linear_interpolation(alpha_rans, CL_rans).(alphas)) < 1e-1
        @test max_error(CD_vsm, linear_interpolation(alpha_rans, CD_rans).(alphas)) < 1e-1
        @test max_error(CL_llt, CL_vsm) < 2e-1
        @test max_error(CD_llt, CD_vsm) < 4e-2
    end

    @testset "rectangular AR 12 NACA 4415 wing matches CFD" begin
        n_sections = 60
        alphas = [3.0, 6.0, 9.0]
        coordinates = generate_coordinates_rect_wing(ones(n_sections), 12.0,
            zeros(n_sections), zeros(n_sections), n_sections, "lin")
        polar = polar_vectors(read_columns("naca4415_cfd_polar.csv")...)
        wing = wing_from_coordinates(coordinates, POLAR_VECTORS, polar)
        CL_llt, CD_llt = lift_drag_polar(wing, LLT, alphas; wind_speed=20.0,
                                         relaxation_factor=0.03)
        CL_vsm, CD_vsm = lift_drag_polar(wing, VSM, alphas; wind_speed=20.0,
                                         relaxation_factor=0.03)
        alpha_cfd, CL_cfd = read_columns("rectangular_wing_ar12_cfd.csv")

        @test max_error(CL_vsm, linear_interpolation(alpha_cfd, CL_cfd).(alphas)) < 1e-1
        @test max_error(CL_llt, CL_vsm) < 2e-1
        @test max_error(CD_llt, CD_vsm) < 4e-2
    end

    @testset "V3 kite with Breukels sections matches RANS" begin
        leading_edges_mm = [
            859.580 -4139.660 5654.227
            -17.623 -3967.978 6471.622
            -237.683 -3134.335 7476.759
            -383.733 -1959.729 8078.914
            -456.562 -664.252 8339.101
            -456.562 664.252 8339.101
            -383.733 1959.729 8078.914
            -237.683 3134.335 7476.759
            -17.623 3967.978 6471.622
            859.580 4139.660 5654.227
        ]
        trailing_edges_mm = [
            1538.773 -4113.307 5530.496
            1703.467 -3955.506 6467.819
            2002.516 -3116.753 7454.254
            2110.145 -1946.587 8041.739
            2158.559 -660.721 8294.064
            2158.559 660.721 8294.064
            2110.145 1946.587 8041.739
            2002.516 3116.753 7454.254
            1703.467 3955.506 6467.819
            1538.773 4113.307 5530.496
        ]
        coordinates = zeros(20, 3)
        coordinates[1:2:end, :] .= leading_edges_mm ./ 1000
        coordinates[2:2:end, :] .= trailing_edges_mm ./ 1000
        alphas = [3.0, 6.0, 9.0]
        wing = wing_from_coordinates(coordinates, POLY, lei_poly_coeffs(0.1, 0.095);
                                     n_panels=36, spanwise_distribution=SPLIT_PROVIDED)
        CL_vsm, CD_vsm = lift_drag_polar(wing, VSM, alphas; wind_speed=22.0,
                                         relaxation_factor=0.03)
        alpha_cl, CL_rans = read_columns("v3_kite_rans_cl.csv")
        alpha_cd, CD_rans = read_columns("v3_kite_rans_cd.csv")

        @test length(wing.refined_sections) == 37
        @test max_error(CL_vsm, linear_interpolation(alpha_cl, CL_rans).(alphas)) < 1e-1
        @test max_error(CD_vsm, linear_interpolation(alpha_cd, CD_rans).(alphas)) < 1e-1
    end

    @testset "three horseshoe vortices match Biot-Savart" begin
        work_vectors = ntuple(_ -> zeros(3), 10)
        flow_direction = [1.0, 0.0, 0.0]
        wind_speed = 1.0
        trailing_length = 100.0
        evaluation_point = zeros(3)
        horseshoes = (
            (gamma=2.0, left=[0.0, -1.0, 0.0], right=[0.0, -3.0, 0.0],
             velocity=[0.0, 0.0, 0.1061]),
            (gamma=10.0, left=[0.0, 1.0, 0.0], right=[0.0, -1.0, 0.0],
             velocity=[0.0, 0.0, -1.5915]),
            (gamma=5.0, left=[0.0, 3.0, 0.0], right=[0.0, 1.0, 0.0],
             velocity=[0.0, 0.0, 0.2653]),
        )
        for (gamma, left, right, velocity) in horseshoes
            bound_velocity = zeros(3)
            velocity_3D_bound_vortex!(bound_velocity, bound_filament(right, left),
                evaluation_point, gamma, 1e-5, work_vectors)

            left_velocity = zeros(3)
            right_velocity = zeros(3)
            left_leg = SemiInfiniteFilament{Float64}()
            right_leg = SemiInfiniteFilament{Float64}()
            reinit!(left_leg, left, flow_direction, wind_speed, 1)
            reinit!(right_leg, right, flow_direction, wind_speed, -1)
            velocity_3D_trailing_vortex_semiinfinite!(left_velocity, left_leg,
                flow_direction, evaluation_point, gamma, wind_speed, work_vectors)
            velocity_3D_trailing_vortex_semiinfinite!(right_velocity, right_leg,
                flow_direction, evaluation_point, gamma, wind_speed, work_vectors)
            @test bound_velocity + left_velocity + right_velocity ≈ velocity atol=1e-4

            wake_offset = trailing_length * flow_direction
            velocity_3D_trailing_vortex!(left_velocity,
                bound_filament(left, left + wake_offset), evaluation_point, gamma,
                wind_speed, work_vectors)
            velocity_3D_trailing_vortex!(right_velocity,
                bound_filament(right + wake_offset, right), evaluation_point, gamma,
                wind_speed, work_vectors)
            @test bound_velocity + left_velocity + right_velocity ≈ velocity atol=1e-4
        end
    end
end

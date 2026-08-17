using VortexStepMethod
using VortexStepMethod: flow_curvature_cm
using LinearAlgebra
using Test

@testset "section_pitch_rate reduces to the rigid-body rate" begin
    # a rigid section rotating at q about y_airf must give back exactly q, which
    # is what lets one expression serve both rigid and deforming wings
    x_airf, y_airf, z_airf = [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]
    chord, q = 1.7, 0.8
    omega = q * y_airf
    leading = [0.0, 0.0, 0.0]
    trailing = chord * x_airf
    @test section_pitch_rate(cross(omega, leading), cross(omega, trailing),
                             z_airf, chord) ≈ q
    @test section_pitch_rate([0.0, 0, 0], [0.0, 0, 0], z_airf, chord) == 0.0
    @test section_pitch_rate([0.0, 0, 0], [0.0, 0, 1.0], z_airf, 0.0) == 0.0
end

@testset "Flow curvature pitch-rate moment" begin
    chord, span, V, q = 1.0, 8.0, 20.0, 0.5

    wing = Wing(12)
    for y in range(span / 2, -span / 2, length=5)
        add_section!(wing, [0.0, y, 0.0], [chord, y, 0.0], INVISCID)
    end
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, [V, 0.0, 0.0], [0.0, q, 0.0])

    panel = body_aero.panels[6]
    q_local = dot(body_aero.omega, panel.y_airf)

    @testset "increment matches thin airfoil theory" begin
        @test body_aero.pitch_rate_dist[6] ≈ q_local
        @test flow_curvature_cm(q_local, chord, V) ≈
              -0.25π * q_local * chord / (2V)
        @test flow_curvature_cm(0.0, chord, V) == 0.0
        @test flow_curvature_cm(q_local, chord, 0.0) == 0.0
    end

    @testset "positive rate about y_airf raises aft incidence" begin
        # nose-up rotation loads the aft chord more, which is what makes the
        # negative increment a damping rather than a driving moment
        aft = panel.aero_center + 0.1chord * panel.x_airf
        dv = cross(body_aero.omega, panel.aero_center) -
             cross(body_aero.omega, aft)
        @test sign(dot(dv, panel.z_airf)) == sign(q_local)
    end

    solver_off = Solver(body_aero; flow_curvature=false)
    solver_on = Solver(body_aero; flow_curvature=true)

    function moment_at(solver, omega)
        set_va!(body_aero, [V, 0.0, 0.0], omega)
        solve!(solver, body_aero)
        return copy(solver.sol.moment)
    end

    @testset "opposes the rotation" begin
        # closed form of the summed increment, pi*rho*V*S*c^2/16 per rad/s,
        # exact in the limit where the induced velocity is small against V
        expected = π * 1.225 * V * span * chord * chord^2 / 16 * q
        for omega in ([0.0, q, 0.0], [0.0, -q, 0.0])
            dM = moment_at(solver_on, omega) .- moment_at(solver_off, omega)
            @test dot(dM, omega) < 0
            @test norm(dM) ≈ expected rtol = 0.05
        end
    end

    @testset "no effect without rotation" begin
        dM = moment_at(solver_on, zeros(3)) .- moment_at(solver_off, zeros(3))
        @test all(iszero, dM)
    end

    @testset "defaults to off" begin
        @test Solver(body_aero).flow_curvature == false
        @test VortexStepMethod.SolverSettings().flow_curvature == false
    end

    @testset "distributed rates drive a deformation mode" begin
        n = length(body_aero.panels)
        va_dist = repeat([V 0.0 0.0], n)

        set_va!(body_aero, va_dist)
        @test all(iszero, body_aero.pitch_rate_dist)
        solve!(solver_on, body_aero)
        base = copy(solver_on.sol.cm_dist)

        # antisymmetric twist rate: no rigid-body omega can express this
        rates = [panel.aero_center[2] > 0 ? 1.0 : -1.0 for panel in body_aero.panels]
        set_va!(body_aero, va_dist; pitch_rate_dist=rates)
        @test body_aero.pitch_rate_dist ≈ rates
        solve!(solver_on, body_aero)

        for i in 1:n
            @test solver_on.sol.cm_dist[i] - base[i] ≈
                  flow_curvature_cm(rates[i], solver_on.sol._chord_dist[i],
                                    solver_on.lr.v_a_dist[i])
        end
        # the two half-wings must be driven in opposite senses
        @test sign(solver_on.sol.cm_dist[1] - base[1]) ==
              -sign(solver_on.sol.cm_dist[n] - base[n])

        @test_throws ArgumentError set_va!(body_aero, va_dist;
                                           pitch_rate_dist=rates[1:end-1])
    end
end

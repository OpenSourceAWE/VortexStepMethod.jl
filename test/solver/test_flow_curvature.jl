using VortexStepMethod
using VortexStepMethod: flow_curvature_cm
using LinearAlgebra
using Test

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
        @test flow_curvature_cm(body_aero.omega, panel, chord, V) ≈
              -0.25π * q_local * chord / (2V)
        @test flow_curvature_cm(zeros(3), panel, chord, V) == 0.0
        @test flow_curvature_cm(body_aero.omega, panel, chord, 0.0) == 0.0
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
end

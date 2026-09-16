using VortexStepMethod
using VortexStepMethod: spanwise_flow_drag, calc_forces!
using LinearAlgebra
using Test

@testset "spanwise_flow_drag is the Gaunaa et al. 2024 model" begin
    # Gaunaa et al. 2024, J. Phys.: Conf. Ser. 2767 022068, with Re and both
    # coefficients built on the speed normal to the span:
    # f0 = 0.062 Re^(-1/7), ΔCd = f0 (cos(β)^(-5/7) - 1), C_par = f0 tan(β) cos(β)^(-5/7)
    density, mu, chord, speed = 1.225, 1.81e-5, 0.8, 25.0
    for beta in deg2rad.((-40.0, -5.0, 0.0, 10.0, 35.0))
        v_normal = speed * cos(beta)
        f0 = 0.062 * (density * v_normal * chord / mu)^(-1 / 7)
        correction = spanwise_flow_drag(v_normal, speed * sin(beta), chord, density, mu)
        @test correction.delta_cd ≈ f0 * (cos(beta)^(-5 / 7) - 1) atol = 1e-14
        @test correction.c_span ≈ f0 * tan(beta) * cos(beta)^(-5 / 7) atol = 1e-14
    end
    @test spanwise_flow_drag(speed, 0.0, chord, density, mu) == (delta_cd=0.0, c_span=0.0)
end

@testset "Spanwise-flow viscous drag correction in the solver" begin
    chord, span, V = 1.0, 8.0, 20.0
    alpha, beta = deg2rad(6.0), deg2rad(12.0)

    wing = Wing(12)
    for y in range(span / 2, -span / 2, length=5)
        add_section!(wing, [0.0, y, 0.0], [chord, y, 0.0], INVISCID)
    end
    refine!(wing)
    body_aero = BodyAerodynamics([wing])

    solver_off = Solver(body_aero; use_gamma_prev=false)
    solver_on = Solver(body_aero; use_gamma_prev=false,
                       is_with_viscous_drag_correction=true)
    va_sideslip = V .* [cos(alpha) * cos(beta), sin(beta), sin(alpha) * cos(beta)]
    va_straight = V .* [cos(alpha), 0.0, sin(alpha)]

    function force_dist_at(solver, va)
        set_va!(body_aero, va)
        solve!(solver, body_aero)
        return copy(solver.sol.f_body_3D)
    end

    @testset "adds the model's drag and spanwise force to each panel" begin
        density, mu = solver_on.density, solver_on.mu
        for va in (va_sideslip, va_straight)
            delta_force = force_dist_at(solver_on, va) .- force_dist_at(solver_off, va)
            for (i, panel) in enumerate(body_aero.panels)
                v_normal = solver_on.lr.v_a_dist[i]
                v_span = solver_on.lr.v_span_dist[i]
                cos_beta = v_normal / hypot(v_normal, v_span)
                f0 = 0.062 * (density * v_normal * panel.chord / mu)^(-1 / 7)
                q_chord_width = 0.5 * density * v_normal^2 * panel.chord * panel.width
                span_force = dot(delta_force[:, i], panel.y_airf)
                drag_force = norm(delta_force[:, i] .- span_force .* panel.y_airf)
                @test span_force ≈ q_chord_width * f0 * (v_span / v_normal) *
                                   cos_beta^(-5 / 7) atol = 1e-12
                @test drag_force ≈ q_chord_width * f0 * (cos_beta^(-5 / 7) - 1) atol = 1e-12
            end
        end
    end

    @testset "sideslip drives the spanwise flow and raises the drag" begin
        delta_force = force_dist_at(solver_on, va_sideslip) .-
                      force_dist_at(solver_off, va_sideslip)
        for (i, panel) in enumerate(body_aero.panels)
            @test solver_on.lr.v_span_dist[i] ≈ dot(va_sideslip, panel.y_airf) rtol = 0.05
            @test dot(delta_force[:, i], va_sideslip) > 0
        end
    end

    @testset "solve reports the corrected forces" begin
        set_va!(body_aero, va_sideslip)
        solve!(solver_on, body_aero)
        results = solve(solver_on, body_aero)
        @test [results["Fx"], results["Fy"], results["Fz"]] ≈ solver_on.sol.force
        @test results["F_distribution"] ≈ solver_on.sol.f_body_3D
    end

    @testset "linearize reports the corrected forces" begin
        y = [va_sideslip; zeros(3)]
        results_for(solver) = VortexStepMethod.linearize(solver, body_aero, y;
            theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6)[2]
        results_on, results_off = results_for(solver_on), results_for(solver_off)
        @test results_on[1:3] ≈ vec(sum(force_dist_at(solver_on, va_sideslip); dims=2))
        @test !(results_on[1:3] ≈ results_off[1:3])
    end

    @testset "calc_forces! stays zero-alloc" begin
        set_va!(body_aero, va_sideslip)
        solve!(solver_on, body_aero)
        calc_forces!(solver_on, body_aero)
        @test (@allocated calc_forces!(solver_on, body_aero)) == 0
    end

    @testset "defaults to off" begin
        @test Solver(body_aero).is_with_viscous_drag_correction == false
        @test VortexStepMethod.SolverSettings().is_with_viscous_drag_correction == false
    end
end

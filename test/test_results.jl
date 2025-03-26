using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!
using LinearAlgebra
using Test
using Logging

include("utils.jl")

if !@isdefined ram_wing
    cp("data/ram_air_kite_body.obj", "/tmp/ram_air_kite_body.obj"; force=true)
    cp("data/ram_air_kite_foil.dat", "/tmp/ram_air_kite_foil.dat"; force=true)
    ram_wing = RamAirWing("/tmp/ram_air_kite_body.obj", "/tmp/ram_air_kite_foil.dat"; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1))
end

@testset "Calculate results against output results" begin
    # Setup
    density = 1.225
    N = 40
    max_chord = 1.0
    span = 15.709  # AR = 20
    v_a = 20.0
    AR = span^2 / (π * span * max_chord / 4)
    aoa = deg2rad(5)
    v_a = [cos(aoa), 0.0, sin(aoa)] .* v_a
    model = VSM

    # Setup wing geometry
    dist = "cos"
    core_radius_fraction = 1e-20
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N; spanwise_distribution=UNCHANGED)
    for idx in 1:2:length(coord_left_to_right[:, 1])
        @debug "coord_left_to_right[$idx] = $(coord_left_to_right[idx,:])"
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            INVISCID
        )
    end
    
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, v_a)

    # Run analysis
    P = length(body_aero.panels)
    loop_solver = Solver{P}(
        aerodynamic_model_type=model, 
        core_radius_fraction=core_radius_fraction,
        solver_type=LOOP,
        atol=1e-8,
        rtol=1e-8
    )
    nonlin_solver = Solver{P}(
        aerodynamic_model_type=model, 
        core_radius_fraction=core_radius_fraction,
        solver_type=NONLIN,
        atol=1e-8,
        rtol=1e-8
    )
    results_NEW = solve(loop_solver, body_aero; reference_point=[0,1,0])
    # println(results_NEW)

    @test results_NEW isa Dict

    @testset "Loop and nonlin solve!" begin
        loop_sol = solve!(loop_solver, body_aero; reference_point=[0,1,0])
        nonlin_sol = solve!(nonlin_solver, body_aero; reference_point=[0,1,0])

        @test all(isapprox.(nonlin_sol.gamma_distribution, loop_sol.gamma_distribution; atol=1e-4))

        @test loop_sol.force.x ≈ -117.96518414816444 atol=1e-4
        @test loop_sol.force.y ≈ 0.0 atol=1e-10
        @test loop_sol.force.z ≈ 1481.996390329679 atol=1e-4 rtol= 1e-4

        @test loop_sol.moment.x ≈ -1481.996390329678 atol=1e-4 rtol= 1e-4
        @test loop_sol.moment.y ≈ 0.0 atol=1e-10
        @test loop_sol.moment.z ≈ -117.9651841481644 atol=1e-4

        @test loop_sol.force_coefficients[1] ≈ -0.039050322560956294 atol=1e-4 # CFx
        @test loop_sol.force_coefficients[2] ≈ 0.0                   atol=1e-4 # CFy
        @test loop_sol.force_coefficients[3] ≈ 0.49055973654418716   atol=1e-4 # CFz
        @test loop_sol.force_coefficients[3] / loop_sol.force_coefficients[1] ≈ loop_sol.force[3] / loop_sol.force[1]
        @test loop_sol.moment_distribution[1] ≈ -0.0006683569356186426 atol=1e-8
        @test loop_sol.moment_coeff_dist[1] ≈ -2.212405554436003e-7 atol=1e-10
        @test loop_sol.moment_distribution[1] / loop_sol.moment_distribution[2] ≈ loop_sol.moment_coeff_dist[1] / loop_sol.moment_coeff_dist[2]

        @test loop_sol.solver_status == FEASIBLE
    end

    # Calculate forces using uncorrected alpha
    alpha = results_NEW["alpha_uncorrected"]
    dyn_visc = 0.5 * density * norm(v_a)^2
    n_panels = length(body_aero.panels)
    lift = zeros(n_panels)
    drag = zeros(n_panels)
    moment = zeros(n_panels)
    
    for (i, panel) in enumerate(body_aero.panels)
        lift[i] = dyn_visc * calculate_cl(panel, alpha[i]) * panel.chord
        cd_cm = calculate_cd_cm(panel, alpha[i])
        drag[i] = dyn_visc * cd_cm[1] * panel.chord
        moment[i] = dyn_visc * cd_cm[2] * panel.chord^2
        # @info "lift: $lift, drag: $drag, moment: $moment"
    end
    Fmag = hcat(lift, drag, moment)

    # Calculate coefficients using corrected alpha
    alpha = results_NEW["alpha_at_ac"]
    aero_coeffs = hcat(
        [alpha[i] for (i, panel) in enumerate(body_aero.panels)],
        [calculate_cl(panel, alpha[i]) for (i, panel) in enumerate(body_aero.panels)],
        [calculate_cd_cm(panel, alpha[i])[1] for (i, panel) in enumerate(body_aero.panels)],
        [calculate_cd_cm(panel, alpha[i])[2] for (i, panel) in enumerate(body_aero.panels)]
    )
    
    ringvec = [Dict("r0" => panel.width * panel.y_airf) for panel in body_aero.panels]
    controlpoints = [Dict("tangential" => panel.x_airf, "normal" => panel.z_airf) 
                    for panel in body_aero.panels]
    Atot = calculate_projected_area(wing)

    F_rel_ref, F_gl_ref, Ltot_ref, Dtot_ref, CL_ref, CD_ref, CS_ref = 
        output_results(Fmag, aero_coeffs, ringvec, v_a, controlpoints, Atot)

    # Compare results
    @info "Comparing results"
    @info "cl_calculated: $(results_NEW["cl"]), CL_ref: $CL_ref"
    @info "cd_calculated: $(results_NEW["cd"]), CD_ref: $CD_ref"
    @info "cs_calculated: $(results_NEW["cs"]), CS_ref: $CS_ref"
    @info "L_calculated: $(results_NEW["lift"]), Ltot_ref: $Ltot_ref"
    @info "D_calculated: $(results_NEW["drag"]), Dtot_ref: $Dtot_ref"

    # Assert results
    @test isapprox(results_NEW["cl"], CL_ref, rtol=1e-4)
    @test isapprox(results_NEW["cd"], CD_ref, rtol=1e-4)
    @test isapprox(results_NEW["cs"], CS_ref, rtol=1e-4)
    @test isapprox(results_NEW["lift"], Ltot_ref, rtol=1e-4)
    @test isapprox(results_NEW["drag"], Dtot_ref, rtol=1e-4)
    @test isapprox(results_NEW["Fx"], results_NEW["Mz"], rtol=1e-4) # 1 meter arm
    @test isapprox(results_NEW["My"], 0.0, atol=1e-3)
    @test isapprox(results_NEW["Fz"], -results_NEW["Mx"], rtol=1e-4) # 1 meter arm

    # Check array shapes
    @test length(results_NEW["cl_distribution"]) == length(body_aero.panels)
    @test length(results_NEW["cd_distribution"]) == length(body_aero.panels)
end

@testset "Linearized model vs nonlinear model" begin
    body_aero = BodyAerodynamics([ram_wing])
    
    
end

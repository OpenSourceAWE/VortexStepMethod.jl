using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!
using LinearAlgebra
using Test
using Logging

include("../utils.jl")
if !@isdefined(create_temp_wing_settings)
    include("../test_data_utils.jl")
end

"""
    inviscid_wing(section_y; n_panels=length(section_y) - 1)

A refined flat `INVISCID` wing of unit chord [m] in the plane z = 0, with one section at
each spanwise position in `section_y` [m].
"""
function inviscid_wing(section_y; n_panels=length(section_y) - 1)
    wing = Wing(n_panels)
    for y in section_y
        add_section!(wing, [0.0, y, 0.0], [1.0, y, 0.0], INVISCID)
    end
    refine!(wing)
    return wing
end

@testset "Induction Matrix Creation" begin
    # Setup
    n_panels = 3
    N = n_panels + 1  # number of SECTIONS
    max_chord = 1.0
    span = 2.36
    AR = span^2 / (π * span * max_chord / 4)
    dist = "cos"
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)

    Atot = max_chord / 2 * span / 2 * π
    @debug "N: $N"
    @debug "size(coord): $(size(coord))"

    va = 20.0
    aoa = 5.7106 * π / 180
    va_vec = [cos(aoa), 0.0, sin(aoa)] .* va

    # Create wing geometry
    core_radius_fraction = 1e-20
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(n_panels; spanwise_distribution=UNCHANGED)
    for idx in 1:2:size(coord_left_to_right, 1)
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            INVISCID
        )
    end
    
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, va_vec)

    # Calculate reference matrices using thesis functions
    controlpoints, rings, bladepanels, ringvec, coord_L =
        create_geometry_general(coord, va_vec, N, "5fil", LLT)
    
    # Test LLT matrices
    @testset "LLT Matrices" begin
        # Calculate reference matrices
        MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
            deepcopy(ringvec),
            deepcopy(controlpoints),
            deepcopy(rings),
            deepcopy(va_vec),
            zeros(N-1),
            nothing,  # data_airf not needed
            nothing,  # conv_crit not needed
            LLT
        )

        # Calculate new matrices
        va_dist = fill(norm(va_vec), length(body_aero.panels))
        va_unit_dist = repeat(reshape(va_vec ./ norm(va_vec), 1, 3),
                              length(body_aero.panels))
        calculate_AIC_matrices!(
            body_aero,
            LLT,
            core_radius_fraction,
            va_dist,
            va_unit_dist
        )
        AIC_x, AIC_y, AIC_z = @views body_aero.AIC[:, :, 1], body_aero.AIC[:, :, 2], body_aero.AIC[:, :, 3]

        # Compare matrices
        @test isapprox(MatrixU, AIC_x, atol=1e-5)
        @test isapprox(MatrixV, -AIC_y, atol=1e-5)
        @test isapprox(MatrixW, AIC_z, atol=1e-5)
    end

    # Test VSM matrices
    @testset "VSM Matrices" begin
        # Calculate reference matrices for VSM
        controlpoints, rings, bladepanels, ringvec, coord_L = 
            create_geometry_general(coord, va_vec, N, "5fil", VSM)
        
        MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
            deepcopy(ringvec),
            deepcopy(controlpoints),
            deepcopy(rings),
            deepcopy(va_vec),
            zeros(N-1),
            nothing,
            nothing,
            VSM
        )

        # Calculate new matrices
        va_dist = fill(norm(va_vec), length(body_aero.panels))
        va_unit_dist = repeat(reshape(va_vec ./ norm(va_vec), 1, 3),
                              length(body_aero.panels))
        calculate_AIC_matrices!(
            body_aero,
            VSM,
            core_radius_fraction,
            va_dist,
            va_unit_dist
        )
        AIC_x, AIC_y, AIC_z = body_aero.AIC[:, :, 1], body_aero.AIC[:, :, 2], body_aero.AIC[:, :, 3]

        # Compare matrices with higher precision for VSM
        @test isapprox(MatrixU, AIC_x, atol=1e-8)
        @test isapprox(MatrixV, -AIC_y, atol=1e-8)
        @test isapprox(MatrixW, AIC_z, atol=1e-8)
    end
end


@testset "Wing Geometry Creation" begin
    @testset "Origin Translation" begin
        wing = inviscid_wing([0.0, 1.0, 2.0])

        # Test non-zero origin translation
        origin = MVec3(1.0, 2.0, 3.0)
        body_aero = BodyAerodynamics([wing]; kite_body_origin=origin)
        
        # Check if sections are correctly translated
        @test wing.unrefined_sections[3].LE_point ≈ [-1.0, -2.0, -3.0]
        @test wing.unrefined_sections[3].TE_point ≈ [0.0, -2.0, -3.0]
        @test wing.unrefined_sections[2].LE_point ≈ [-1.0, -1.0, -3.0]
        @test wing.unrefined_sections[2].TE_point ≈ [0.0, -1.0, -3.0]
        @test wing.unrefined_sections[1].LE_point ≈ [-1.0, 0.0, -3.0]
        @test wing.unrefined_sections[1].TE_point ≈ [0.0, 0.0, -3.0]
    end

    function create_geometry(; model=VSM, wing_type=:rectangular, plotting=false, N=40)
        max_chord = 1.0
        span = 17.0
        AR = span^2 / (π * span * max_chord / 4)
        @debug "AR: $AR"
        va = 20.0
        aoa = 5.7106 * π / 180
        va_vec = [cos(aoa), 0.0, sin(aoa)] .* va
    
        coord = if wing_type === :rectangular
            theta = range(-0.5, 0.5, length=N)
            beta = range(-2, 2, length=N)
            generate_coordinates_rect_wing(
                fill(max_chord, N),
                span,
                theta,
                beta,
                N,
                "lin"
            )
        elseif wing_type === :curved
            generate_coordinates_curved_wing(
                max_chord, span, π/4, 5, N, "cos"
            )
        elseif wing_type === :elliptical
            generate_coordinates_el_wing(max_chord, span, N, "cos")
        else
            error("Invalid wing type")
        end
    
        coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
        wing = Wing(N-1; spanwise_distribution=UNCHANGED)
        for i in 1:2:size(coord_left_to_right, 1)
            add_section!(
                wing,
                coord_left_to_right[i,:],
                coord_left_to_right[i+1,:],
                INVISCID
            )
        end
        refine!(wing)
        body_aero = BodyAerodynamics([wing])
        set_va!(body_aero, va_vec)

        return body_aero, coord, va_vec, model
    end

    for model in [VSM, LLT]
        @debug "model: $model"
        for wing_type in [:rectangular, :curved, :elliptical]
            @debug "wing_type: $wing_type"
            body_aero, coord, va_vec, model = create_geometry(
                model=model, wing_type=wing_type
            )
            
            # Generate geometry
            expected_controlpoints, expected_rings, expected_bladepanels, 
                expected_ringvec, expected_coord_L = create_geometry_general(
                coord, va_vec, div(size(coord,1), 2), "5fil", model
            )

            for i in 1:length(body_aero.panels)
                @debug "i: $i"
                # Handle control points
                index_reversed = length(body_aero.panels) - i + 1
                panel = body_aero.panels[index_reversed]
                
                evaluation_point = if model === VSM
                    panel.control_point
                else  # LLT
                    panel.aero_center
                end

                @test isapprox(evaluation_point, expected_controlpoints[i]["coordinates"], atol=1e-4)
                @test isapprox(panel.chord, expected_controlpoints[i]["chord"], atol=1e-4)
                @test isapprox(panel.z_airf, expected_controlpoints[i]["normal"], atol=1e-4)
                @test isapprox(panel.x_airf, expected_controlpoints[i]["tangential"], atol=1e-4)
                @test isapprox(
                    hcat(panel.z_airf, panel.x_airf, panel.y_airf),
                    expected_controlpoints[i]["airf_coord"],
                    atol=1e-4
                )
                
                if model === VSM
                    @test isapprox(
                        panel.aero_center,
                        expected_controlpoints[i]["coordinates_aoa"],
                        atol=1e-4
                    )
                end

                # Handle rings
                expected_ring_i = expected_rings[i]
                expected_ring_i_list = [
                    expected_ring_i[1],
                    expected_ring_i[2],
                    expected_ring_i[3],
                    expected_ring_i[4],
                    expected_ring_i[5]
                ]

                filaments = panel.filaments
                filament_list = [
                    filaments[1],
                    filaments[3],
                    filaments[5],
                    filaments[2],
                    filaments[4]
                ]

                for (j, fil) in enumerate(filament_list)
                    if j == 1  # bound filaments
                        @test isapprox(fil.x1, expected_ring_i_list[j]["x1"], atol=1e-4)
                        @test isapprox(fil.x2, expected_ring_i_list[j]["x2"], atol=1e-4)
                    elseif j ∈ (2, 4)  # trailing filaments
                        @test isapprox(fil.x1, expected_ring_i_list[j]["x1"], atol=1e-4)
                        @test isapprox(fil.x2, expected_ring_i_list[j]["x2"], atol=1e-4)
                    else  # semi-infinite filaments
                        @test isapprox(fil.x1, expected_ring_i_list[j]["x1"], atol=1e-4)
                    end
                end

                # Handle bladepanels
                exp_bladepanels = expected_bladepanels[i]
                @test isapprox(panel.LE_point_2, exp_bladepanels["p1"], atol=1e-4)
                @test isapprox(panel.LE_point_1, exp_bladepanels["p2"], atol=1e-4)
                @test isapprox(panel.TE_point_1, exp_bladepanels["p3"], atol=1e-4)
                @test isapprox(panel.TE_point_2, exp_bladepanels["p4"], atol=1e-4)

                # Handle ringvec
                exp_ringvec = expected_ringvec[i]
                r0 = panel.bound_point_1 - panel.bound_point_2
                r3 = evaluation_point - (panel.bound_point_1 + panel.bound_point_2) / 2
                @test isapprox(r0, exp_ringvec["r0"], atol=1e-4)
                @test isapprox(r3, exp_ringvec["r3"], atol=1e-4)

                # Handle coord_L
                @test all(isapprox.(panel.aero_center, expected_coord_L[:, i]))
            end
        end
    end
end

@testset "Calculate results against output results" begin
    # Setup
    density = 1.225
    N = 40
    max_chord = 1.0
    span = 15.709  # AR = 20
    va = 20.0
    AR = span^2 / (π * span * max_chord / 4)
    aoa = deg2rad(5)
    va_vec = [cos(aoa), 0.0, sin(aoa)] .* va
    model = VSM

    # Setup wing geometry
    dist = "cos"
    core_radius_fraction = 1e-20
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N-1; spanwise_distribution=UNCHANGED)
    for idx in 1:2:length(coord_left_to_right[:, 1])
        @debug "coord_left_to_right[$idx] = $(coord_left_to_right[idx,:])"
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            INVISCID
        )
    end
    
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, va_vec)

    # Run analysis
    loop_solver = Solver(body_aero;
        aerodynamic_model_type=model,
        core_radius_fraction=core_radius_fraction,
        solver_type=LOOP,
        correct_aoa=true,
        atol=1e-8,
        rtol=1e-8
    )
    nonlin_solver = Solver(body_aero;
        aerodynamic_model_type=model,
        core_radius_fraction=core_radius_fraction,
        solver_type=NONLIN,
        correct_aoa=true,
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

        @test loop_sol.force_coeffs[1] ≈ -0.039050322560956294 atol=1e-4 # CFx
        @test loop_sol.force_coeffs[2] ≈ 0.0                   atol=1e-4 # CFy
        @test loop_sol.force_coeffs[3] ≈ 0.49055973654418716   atol=3e-4 # CFz
        @test loop_sol.force_coeffs[3] / loop_sol.force_coeffs[1] ≈ loop_sol.force[3] / loop_sol.force[1]
        @test loop_sol.moment_dist[1] ≈ -0.0006683569356186426 atol=1e-8
        @test loop_sol.moment_coeff_dist[1] ≈ -2.212405554436003e-7 atol=1e-9
        @test loop_sol.moment_dist[1] / loop_sol.moment_dist[2] ≈ loop_sol.moment_coeff_dist[1] / loop_sol.moment_coeff_dist[2]

        @test loop_sol.solver_status == FEASIBLE

    end

    # Calculate forces using uncorrected alpha
    alpha = results_NEW["alpha_uncorrected"]
    dyn_visc = 0.5 * density * norm(va_vec)^2
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
        output_results(Fmag, aero_coeffs, ringvec, va_vec, controlpoints, Atot)

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

@testset "set_va! with VSMSettings" begin
    settings_file = create_temp_wing_settings("body_aerodynamics", "test_wing.yaml";
                                              alpha=10.0, beta=5.0, wind_speed=15.0)
    try
        settings   = VSMSettings(settings_file)
        wing       = Wing(settings)
        refine!(wing)
        body_aero = BodyAerodynamics([wing])

        set_va!(body_aero, settings)

        α, β, wind_speed = deg2rad(10.0), deg2rad(5.0), 15.0
        expected_va_vec = wind_speed .* [cos(α)*cos(β), sin(β), sin(α)*cos(β)]

        for p in body_aero.panels
            @test p.va ≈ expected_va_vec atol=1e-10
        end
        @test body_aero._va ≈ expected_va_vec atol=1e-10
    finally
        isfile(settings_file) && rm(settings_file; force=true)
    end
end

@testset "set_va! with distributed inflow blocks body_aero.va access" begin
    body_aero = BodyAerodynamics([inviscid_wing([0.0, 1.0, 2.0])])

    va_vec_dist = [
        10.0 0.0 0.0
        9.0 0.0 1.0
    ]
    set_va!(body_aero, va_vec_dist)

    @test body_aero.has_distributed_va
    try
        body_aero.va
        @test false
    catch err
        @test err isa ArgumentError
        @test occursin("distributed inflow", sprint(showerror, err))
    end

    set_va!(body_aero, [11.0, 0.0, 0.0])
    @test !body_aero.has_distributed_va
    @test body_aero.va ≈ [11.0, 0.0, 0.0]
end

@testset "set_va! with omega on multi-wing body" begin
    body_aero = BodyAerodynamics([inviscid_wing([0.0, 1.0, 2.0]),
                                  inviscid_wing([10.0, 11.0, 12.0])])

    va_vec = [10.0, 0.0, 0.0]
    omega = [0.0, 0.0, 1.0]
    set_va!(body_aero, va_vec, omega)

    for panel in body_aero.panels
        expected_va_vec = va_vec .+ (-omega × panel.control_point)
        @test panel.va ≈ expected_va_vec atol=1e-12
    end
    @test body_aero.omega ≈ omega
    @test !body_aero.has_distributed_va
    @test body_aero.va ≈ va_vec

    new_omega = [0.0, 0.0, 2.0]
    @test body_aero._va ≈ va_vec
    body_aero.omega = new_omega

    for panel in body_aero.panels
        expected_va_vec = va_vec .+ (-new_omega × panel.control_point)
        @test panel.va ≈ expected_va_vec atol=1e-12
    end
    @test body_aero.omega ≈ new_omega
end

"""
    solve_wings(wings)

The `BodyAerodynamics` built from `wings` in a 10 m/s inflow and its `solve!` solution.
"""
function solve_wings(wings)
    body_aero = BodyAerodynamics(wings; va=[10.0, 0.0, 1.0])
    return body_aero, solve!(Solver(body_aero), body_aero)
end

"""
    wing_pair(section_y, n_panels, offset)

Two `inviscid_wing`s of `n_panels` panels, the second shifted by `-offset` [m] in y.
"""
function wing_pair(section_y, n_panels, offset)
    return [inviscid_wing(section_y; n_panels),
            inviscid_wing(section_y .- offset; n_panels)]
end

"""
    linearize_body(body_aero; kwargs...)

`linearize` of `body_aero` at zero twist and deflection over the twist and the deflection of
every unrefined section, then the inflow and the angular rate.
"""
function linearize_body(body_aero; kwargs...)
    n_sections = sum(wing -> wing.n_unrefined_sections, body_aero.wings)
    solver = Solver(body_aero; use_gamma_prev=false, rtol=1e-10)
    y0 = [zeros(2n_sections); body_aero.va; zeros(3)]
    return VortexStepMethod.linearize(solver, body_aero, y0;
        theta_idxs=1:n_sections, delta_idxs=n_sections+1:2n_sections,
        va_idxs=2n_sections+1:2n_sections+3, omega_idxs=2n_sections+4:2n_sections+6,
        kwargs...)
end

@testset "solve! on a two-wing body" begin
    n_panels = 6
    section_y = [2.0, 0.0, -2.0]
    span = section_y[1] - section_y[end]
    single_aero, single = solve_wings([inviscid_wing(section_y; n_panels)])
    @test single.solver_status == FEASIBLE

    @testset "wings far apart each act as the isolated wing" begin
        body_aero, sol = solve_wings(wing_pair(section_y, n_panels, 1e4))

        @test length(body_aero.panels) == 2n_panels
        @test sol.solver_status == FEASIBLE
        @test body_aero.projected_area ≈ 2single_aero.projected_area
        @test sol.gamma_distribution ≈ repeat(single.gamma_distribution, 2) rtol=1e-6
        @test sol.cl_unrefined_dist ≈ repeat(single.cl_unrefined_dist, 2) rtol=1e-6
        @test sol.force ≈ 2single.force rtol=1e-6
        @test sol.force_coeffs ≈ single.force_coeffs rtol=1e-6
    end

    @testset "wings one chord apart induce on each other" begin
        _, sol = solve_wings(wing_pair(section_y, n_panels, span + 1.0))
        gamma = sol.gamma_distribution

        @test sol.solver_status == FEASIBLE
        @test gamma ≈ reverse(gamma) rtol=1e-6
        @test gamma[n_panels] > 1.05single.gamma_distribution[n_panels]
        @test sol.force[3] > 2single.force[3]
    end

    n_wing_sections = length(section_y)
    first_sections = 1:n_wing_sections
    second_sections = n_wing_sections+1:2n_wing_sections
    @testset "unrefined_deform! hands each wing its own run of angles" begin
        body_aero = BodyAerodynamics(wing_pair(section_y, n_panels, span + 1.0))
        isolated = wing_pair(section_y, n_panels, span + 1.0)
        theta = deg2rad.(1.0:2n_wing_sections)
        delta = -2theta
        VortexStepMethod.unrefined_deform!(body_aero, theta, delta)

        for (wing_idx, section_range) in enumerate((first_sections, second_sections))
            wing = isolated[wing_idx]
            VortexStepMethod.unrefined_deform!(wing, theta[section_range],
                delta[section_range])
            @test VortexStepMethod.unrefined_section_range(body_aero, wing_idx) ==
                  section_range
            @test body_aero.wings[wing_idx].theta_dist == wing.theta_dist
            @test body_aero.wings[wing_idx].delta_dist == wing.delta_dist
        end
    end

    @testset "linearize: theta of each wing moves that wing's sections" begin
        body_aero, _ = solve_wings(wing_pair(section_y, n_panels, 1e4))
        jac, _, converged = linearize_body(body_aero)
        own_first = jac[6 .+ first_sections, first_sections]

        @test converged
        @test norm(own_first) > 0
        @test jac[6 .+ second_sections, second_sections] ≈ own_first rtol=1e-9
        @test norm(jac[6 .+ first_sections, second_sections]) < 1e-7norm(own_first)
        @test norm(jac[6 .+ second_sections, first_sections]) < 1e-7norm(own_first)
    end

    @testset "linearize: AutoForwardDiff matches AutoFiniteDiff" begin
        body_aero, _ = solve_wings(wing_pair(section_y, n_panels, span + 1.0))
        jac_fwd, _, fwd_converged = linearize_body(body_aero)
        jac_fd, _, fd_converged = linearize_body(body_aero; backend=nothing,
            fd_absstep=1e-6, fd_relstep=1e-6)

        @test fwd_converged
        @test fd_converged
        @test norm(jac_fwd[:, 1:2n_wing_sections]) > 0
        @test jac_fwd ≈ jac_fd rtol=1e-4
    end
end

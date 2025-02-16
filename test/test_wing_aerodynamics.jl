# using VortexStepMethod: Wing, WingAerodynamics, BoundFilament, SemiInfiniteFilament, add_section!, set_va!, solve, calculate_cl
using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices
using LinearAlgebra
using Test
using Logging

# ENV["JULIA_DEBUG"] = "all"
# global_logger(ConsoleLogger(stderr, Logging.Debug))

include("utils.jl")

@testset "Calculate results against output results" begin
    # Setup
    density = 1.225
    N = 40
    max_chord = 1.0
    span = 15.709  # AR = 20
    Umag = 20.0
    AR = span^2 / (π * span * max_chord / 4)
    aoa = deg2rad(5)
    Uinf = [cos(aoa), 0.0, sin(aoa)] .* Umag
    model = "VSM"

    # Setup wing geometry
    dist = "cos"
    core_radius_fraction = 1e-20
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N; spanwise_panel_distribution="unchanged")
    for idx in 1:2:length(coord_left_to_right[:, 1])
        @debug "coord_left_to_right[$idx] = $(coord_left_to_right[idx,:])"
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            "inviscid"
        )
    end
    
    wing_aero = WingAerodynamics([wing])
    set_va!(wing_aero, (Uinf, 0.0))

    # Run analysis
    solver_object = Solver(
        aerodynamic_model_type=model, 
        core_radius_fraction=core_radius_fraction
    )
    results_NEW = solve(solver_object, wing_aero)

    @test results_NEW isa Dict

    # Calculate forces using uncorrected alpha
    alpha = results_NEW["alpha_uncorrected"]
    dyn_visc = 0.5 * density * norm(Uinf)^2
    n_panels = length(wing_aero.panels)
    lift = zeros(n_panels)
    drag = zeros(n_panels)
    moment = zeros(n_panels)
    
    for (i, panel) in enumerate(wing_aero.panels)
        lift[i] = dyn_visc * calculate_cl(panel, alpha[i]) * panel.chord
        cd_cm = calculate_cd_cm(panel, alpha[i])
        drag[i] = dyn_visc * cd_cm[1] * panel.chord
        moment[i] = dyn_visc * cd_cm[2] * panel.chord^2
        @info "lift: $lift, drag: $drag, moment: $moment"
    end
    Fmag = hcat(lift, drag, moment)

    # Calculate coefficients using corrected alpha
    alpha = results_NEW["alpha_at_ac"]
    aero_coeffs = hcat(
        [alpha[i] for (i, panel) in enumerate(wing_aero.panels)],
        [calculate_cl(panel, alpha[i]) for (i, panel) in enumerate(wing_aero.panels)],
        [calculate_cd_cm(panel, alpha[i])[1] for (i, panel) in enumerate(wing_aero.panels)],
        [calculate_cd_cm(panel, alpha[i])[2] for (i, panel) in enumerate(wing_aero.panels)]
    )
    
    ringvec = [Dict("r0" => panel.width * panel.z_airf) for panel in wing_aero.panels]
    controlpoints = [Dict("tangential" => panel.y_airf, "normal" => panel.x_airf) 
                    for panel in wing_aero.panels]
    Atot = calculate_projected_area(wing)

    F_rel_ref, F_gl_ref, Ltot_ref, Dtot_ref, CL_ref, CD_ref, CS_ref = 
        output_results(Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot)

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

    # Check array shapes
    @test length(results_NEW["cl_distribution"]) == length(wing_aero.panels)
    @test length(results_NEW["cd_distribution"]) == length(wing_aero.panels)
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

    Umag = 20.0
    aoa = 5.7106 * π / 180
    Uinf = [cos(aoa), 0.0, sin(aoa)] .* Umag

    # Create wing geometry
    core_radius_fraction = 1e-20
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(n_panels; spanwise_panel_distribution="unchanged")
    for idx in 1:2:size(coord_left_to_right, 1)
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            "inviscid"
        )
    end
    
    wing_aero = WingAerodynamics([wing])
    set_va!(wing_aero, (Uinf, 0.0))

    # Calculate reference matrices using thesis functions
    controlpoints, rings, bladepanels, ringvec, coord_L = 
        create_geometry_general(coord, Uinf, N, "5fil", "LLT")
    
    # Test LLT matrices
    @testset "LLT Matrices" begin
        # Calculate reference matrices
        MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
            deepcopy(ringvec),
            deepcopy(controlpoints),
            deepcopy(rings),
            deepcopy(Uinf),
            zeros(N-1),
            nothing,  # data_airf not needed
            nothing,  # conv_crit not needed
            "LLT"
        )

        # Calculate new matrices
        va_norm_array = fill(norm(Uinf), length(coord))
        va_unit_array = repeat(reshape(Uinf ./ norm(Uinf), 1, 3), length(coord))
        AIC_x, AIC_y, AIC_z = calculate_AIC_matrices(
            wing_aero,
            "LLT",
            core_radius_fraction,
            va_norm_array,
            va_unit_array
        )

        # Compare matrices
        @test isapprox(MatrixU, AIC_x, atol=1e-5)
        @test isapprox(MatrixV, -AIC_y, atol=1e-5)
        @test isapprox(MatrixW, AIC_z, atol=1e-5)
    end

    # Test VSM matrices
    @testset "VSM Matrices" begin
        # Calculate reference matrices for VSM
        controlpoints, rings, bladepanels, ringvec, coord_L = 
            create_geometry_general(coord, Uinf, N, "5fil", "VSM")
        
        MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
            deepcopy(ringvec),
            deepcopy(controlpoints),
            deepcopy(rings),
            deepcopy(Uinf),
            zeros(N-1),
            nothing,
            nothing,
            "VSM"
        )

        # Calculate new matrices
        va_norm_array = fill(norm(Uinf), length(coord))
        va_unit_array = repeat(reshape(Uinf ./ norm(Uinf), 1, 3), length(coord))
        AIC_x, AIC_y, AIC_z = calculate_AIC_matrices(
            wing_aero,
            "VSM",
            core_radius_fraction,
            va_norm_array,
            va_unit_array
        )

        # Compare matrices with higher precision for VSM
        @test isapprox(MatrixU, AIC_x, atol=1e-8)
        @test isapprox(MatrixV, -AIC_y, atol=1e-8)
        @test isapprox(MatrixW, AIC_z, atol=1e-8)
        @show MatrixU MatrixV MatrixW
    end
end
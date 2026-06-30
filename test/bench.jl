using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using BenchmarkTools
using StaticArrays
using VortexStepMethod
using VortexStepMethod: calculate_AIC_matrices!, gamma_loop!, calculate_results,
                       update_effective_angle_of_attack!, calculate_projected_area,
                       calculate_cl, calculate_cd_cm,
                       calculate_velocity_induced_single_ring_semiinfinite!,
                       calculate_velocity_induced_bound_2D!,
                       velocity_3D_bound_vortex!,
                       velocity_3D_trailing_vortex!,
                       velocity_3D_trailing_vortex_semiinfinite!,
                       Panel,
                       reinit!,
                       Wing, add_section!, refine!,
                       apply_billowing_to_pair!,
                       billowing_arc_length, MVec3
using Test
using LinearAlgebra

@testset "Function Allocation Tests" begin
    # Define wing parameters
    n_panels = 20          # Number of panels
    span = 20.0            # Wing span [m]
    chord = 1.0            # Chord length [m]
    v_a = 20.0             # Magnitude of inflow velocity [m/s]
    density = 1.225        # Air density [kg/m³]
    alpha_deg = 30.0       # Angle of attack [degrees]
    alpha = deg2rad(alpha_deg)
    
    wing = Wing(n_panels, spanwise_distribution=LINEAR)
    add_section!(wing, 
        [0.0, span/2, 0.0],    # Left tip LE 
        [chord, span/2, 0.0],  # Left tip TE
        INVISCID)
    add_section!(wing, 
        [0.0, -span/2, 0.0],   # Right tip LE
        [chord, -span/2, 0.0], # Right tip TE
        INVISCID)

    unchanged_wing = Wing(2, spanwise_distribution=UNCHANGED)
    add_section!(unchanged_wing, 
        [0.0, span/2, 0.0],    # Left tip LE 
        [chord, span/2, 0.0],  # Left tip TE
        INVISCID)
    add_section!(unchanged_wing, 
        [0.0, 0.0, 0.0],    # Left tip LE 
        [chord, 0.0, 0.0],  # Left tip TE
        INVISCID)
    add_section!(unchanged_wing, 
        [0.0, -span/2, 0.0],   # Right tip LE
        [chord, -span/2, 0.0], # Right tip TE
        INVISCID)
    
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    refine!(unchanged_wing)
    unchanged_body_aero = BodyAerodynamics([unchanged_wing])
    reinit!(unchanged_body_aero)

    @testset "Re-initialization" begin
        result = @benchmark reinit!($unchanged_body_aero; init_aero=false) samples=1 evals=1
        @info "Re-initializing Allocations: $(result.allocs) \t Memory: $(result.memory)"
        @test result.allocs ≤ 50
    end

    vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
    set_va!(body_aero, vel_app)

    # Initialize solvers for both LLT and VSM methods
    solver = Solver(body_aero)
    nonlin_solver = Solver(body_aero; solver_type=NONLIN)

    # Pre-allocate arrays
    gamma = rand(n_panels)
    gamma_new = similar(gamma)
    AIC_x = rand(n_panels, n_panels)
    AIC_y = similar(AIC_x)
    AIC_z = similar(AIC_x)
    v_ind = zeros(3)
    point = rand(3)
    va_norm_array = ones(n_panels)
    va_unit_array = ones(n_panels, 3)
    
    models = [VSM, LLT]
    core_radius_fractions = [0.001, 10.0]

    @testset "AIC Matrix Calculation" begin
        @info "AIC Matrix Calculation"
        for model in models
            for frac in core_radius_fractions
                @testset "Model $model Core Radius Fraction $frac" begin
                    result = @benchmark calculate_AIC_matrices!($body_aero, $model, $frac, $va_norm_array, $va_unit_array) samples=1 evals=1
                    @test result.allocs ≤ 30
                    @info "Model: $(model) \t Core radius fraction: $(frac) \t Allocations: $(result.allocs) \t Memory: $(result.memory)"
                end
            end
        end
    end
    
    @testset "Gamma Loop" begin
        @info "Gamma Loop"
        # Pre-allocate arrays
        gamma_new = zeros(n_panels)
        va_array = zeros(n_panels, 3)
        chord_array = zeros(n_panels)
        x_airf_array = zeros(n_panels, 3)
        y_airf_array = zeros(n_panels, 3)
        z_airf_array = zeros(n_panels, 3)
        
        # Fill arrays with data
        for (i, panel) in enumerate(body_aero.panels)
            va_array[i, :] .= panel.va
            chord_array[i] = panel.chord
            x_airf_array[i, :] .= panel.x_airf
            y_airf_array[i, :] .= panel.y_airf
            z_airf_array[i, :] .= panel.z_airf
        end

        alphas = collect(-20.0:30.0)
        cls = [2π * α for α in alphas]
        cds = [0.01 + 0.05 * α^2 for α in alphas]
        cms = [-0.1 * α for α in alphas]

        for model in models
            for (aero_model, aero_data) in [(INVISCID, nothing), (POLAR_VECTORS, (alphas, cls, cds, cms))]
                wing = Wing(n_panels, spanwise_distribution=LINEAR)
                add_section!(wing, 
                    [0.0, span/2, 0.0],    # Left tip LE 
                    [chord, span/2, 0.0],  # Left tip TE
                    aero_model,
                    aero_data)
                add_section!(wing, 
                    [0.0, -span/2, 0.0],   # Right tip LE
                    [chord, -span/2, 0.0], # Right tip TE
                    aero_model,
                    aero_data)
                refine!(wing)
                body_aero = BodyAerodynamics([wing])
                
                solver = Solver(body_aero;
                    aerodynamic_model_type=model
                )
                solver.sol._va_dist .= va_array
                solver.sol._chord_dist .= chord_array
                solver.sol._x_airf_dist .= x_airf_array
                solver.sol._y_airf_dist .= y_airf_array
                solver.sol._z_airf_dist .= z_airf_array
                result = @benchmark gamma_loop!(
                    $solver,
                    $body_aero,
                    $body_aero.panels,
                    0.5;
                    log = false
                ) samples=1 evals=1
                @test result.allocs ≤ 10
                @info "Model: $model \t Aero_model: $aero_model \t Allocations: $(result.allocs) Memory: $(result.memory)"
            end
        end
    end
    
    @testset "Results Calculation" begin
        # Pre-allocate arrays
        alpha_array = zeros(n_panels)
        v_a_array = zeros(n_panels)
        chord_array = zeros(n_panels)
        x_airf_array = zeros(n_panels, 3)
        y_airf_array = zeros(n_panels, 3)
        z_airf_array = zeros(n_panels, 3)
        va_array = zeros(n_panels, 3)
        va_norm_array = zeros(n_panels)
        va_unit_array = zeros(n_panels, 3)
        reference_point = zeros(3)
        

        set_va!(body_aero, vel_app)
        # Fill arrays with panel data to satisfy calculate_results preconditions.
        for (i, panel) in enumerate(body_aero.panels)
            chord_array[i] = panel.chord
            x_airf_array[i, :] .= panel.x_airf
            y_airf_array[i, :] .= panel.y_airf
            z_airf_array[i, :] .= panel.z_airf
            va_array[i, :] .= panel.va
            va_norm_array[i] = norm(panel.va)
            va_unit_array[i, :] .= va_norm_array[i] > 0.0 ? panel.va ./ va_norm_array[i] : [1.0, 0.0, 0.0]
            v_a_array[i] = va_norm_array[i]
        end
        results = @MVector zeros(3)
        
        result = @benchmark calculate_results(
            $body_aero,
            $gamma,
            $reference_point,
            $density,
            1e-20,
            0.0,
            $alpha_array,
            $v_a_array,
            $chord_array,
            $x_airf_array,
            $z_airf_array,
            $va_array,
            $va_norm_array,
            $va_unit_array,
            $body_aero.panels,
            false
        ) samples=1 evals=1
        @info "Calculate Results Allocations: $(result.allocs) Memory: $(result.memory)"
        @test result.allocs ≤ 700
    end

    @testset "Allocation Tests for solve() and solve!()" begin
        result = @benchmark  solve_base!($solver, $body_aero, nothing) samples=1 evals=1
        @test result.allocs <= 55
        # time Python: 32.0 ms  Ryzen 7950x
        # time Julia:   0.45 ms Ryzen 7950x
        result = @benchmark  sol = solve!($solver, $body_aero, nothing) samples=1 evals=1
        @test result.allocs <= 110

        # Step 5: Solve using both methods
        result = @benchmark  solve_base!($nonlin_solver, $body_aero, nothing) samples=1 evals=1
        @test result.allocs <= 55
        result = @benchmark  sol = solve!($nonlin_solver, $body_aero, nothing) samples=1 evals=1
        @test result.allocs <= 110
    end

    @testset "Refinement Allocation Tests" begin
        n = 20
        s = 10.0
        c = 1.0

        function make_wing(dist, n_panels, n_ribs;
                           billowing_pct=0.0)
            w = Wing(n_panels;
                spanwise_distribution=dist,
                billowing_percentage=billowing_pct)
            for i in 1:n_ribs
                y = s * (i - 1) / (n_ribs - 1)
                add_section!(w,
                    [0.0, y, 0.0],
                    [c, y, 0.0],
                    INVISCID)
            end
            return w
        end

        @testset "LINEAR" begin
            w = make_wing(LINEAR, n, 2)
            refine!(w)
            result = @benchmark refine!($w) samples=1 evals=1
            @info "LINEAR refine! allocs: $(result.allocs)" *
                  " memory: $(result.memory)"
            @test result.allocs ≤ 70
        end

        @testset "COSINE" begin
            w = make_wing(COSINE, n, 2)
            refine!(w)
            result = @benchmark refine!($w) samples=1 evals=1
            @info "COSINE refine! allocs: $(result.allocs)" *
                  " memory: $(result.memory)"
            @test result.allocs ≤ 70
        end

        @testset "SPLIT_PROVIDED" begin
            w = make_wing(SPLIT_PROVIDED, n, 5)
            refine!(w)
            result = @benchmark refine!($w) samples=1 evals=1
            @info "SPLIT_PROVIDED refine! allocs: " *
                  "$(result.allocs) memory: $(result.memory)"
            @test result.allocs ≤ 100
        end

        @testset "UNCHANGED" begin
            w = make_wing(UNCHANGED, 4, 5)
            refine!(w)
            result = @benchmark refine!($w) samples=1 evals=1
            @info "UNCHANGED refine! allocs: $(result.allocs)" *
                  " memory: $(result.memory)"
            @test result.allocs ≤ 5
        end

        @testset "BILLOWING" begin
            w = make_wing(BILLOWING, 4 * 8, 5;
                          billowing_pct=10.0)
            refine!(w)
            result = @benchmark refine!($w) samples=1 evals=1
            @info "BILLOWING refine! allocs: $(result.allocs)" *
                  " memory: $(result.memory)"
            @test result.allocs ≤ 140
        end

        @testset "billowing_arc_length" begin
            w = make_wing(BILLOWING, 4 * 8, 5;
                          billowing_pct=10.0)
            refine!(w)
            y_hat = MVec3(0.0, -1.0, 0.0)
            le_ref = MVec3(w.unrefined_sections[2].LE_point)
            te_l = MVec3(w.unrefined_sections[1].TE_point)
            te_r = MVec3(w.unrefined_sections[2].TE_point)
            span_len = norm(
                w.unrefined_sections[1].LE_point -
                w.unrefined_sections[2].LE_point)
            # warmup
            billowing_arc_length(
                w.refined_sections, 2, 9,
                y_hat, span_len, le_ref,
                te_l, te_r, 0.3)
            result = @benchmark billowing_arc_length(
                $(w.refined_sections), $2, $9,
                $y_hat, $span_len, $le_ref,
                $te_l, $te_r, $0.3) samples=1 evals=1
            @info "billowing_arc_length allocs: " *
                  "$(result.allocs) memory: $(result.memory)"
            @test result.allocs == 0
        end

        @testset "apply_billowing_to_pair!" begin
            w = make_wing(BILLOWING, 4 * 8, 5;
                          billowing_pct=10.0)
            refine!(w)
            # Reset TE to linear interpolation for re-application
            w2 = make_wing(SPLIT_PROVIDED, 4 * 8, 5)
            refine!(w2)
            for (s1, s2) in zip(
                    w.refined_sections, w2.refined_sections)
                s1.TE_point .= s2.TE_point
            end
            y_hat = MVec3(0.0, -1.0, 0.0)
            le_ref = MVec3(w.unrefined_sections[2].LE_point)
            te_l = MVec3(w.unrefined_sections[1].TE_point)
            te_r = MVec3(w.unrefined_sections[2].TE_point)
            span_len = norm(
                w.unrefined_sections[1].LE_point -
                w.unrefined_sections[2].LE_point)
            # warmup
            apply_billowing_to_pair!(
                w.refined_sections, 2, 9,
                y_hat, span_len, le_ref,
                te_l, te_r, 10.0)
            # Reset again
            for (s1, s2) in zip(
                    w.refined_sections, w2.refined_sections)
                s1.TE_point .= s2.TE_point
            end
            result = @benchmark apply_billowing_to_pair!(
                $(w.refined_sections), $2, $9,
                $y_hat, $span_len, $le_ref,
                $te_l, $te_r, $10.0) samples=1 evals=1
            @info "apply_billowing_to_pair! allocs: " *
                  "$(result.allocs) memory: $(result.memory)"
            @test result.allocs == 0
        end
    end
end

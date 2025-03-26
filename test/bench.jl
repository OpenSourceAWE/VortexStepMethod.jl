using Pkg
if !("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
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
                       init!
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
    
    unchanged_wing = Wing(n_panels, spanwise_distribution=UNCHANGED)
    wing = Wing(n_panels, spanwise_distribution=LINEAR)
    add_section!(wing, 
        [0.0, span/2, 0.0],    # Left tip LE 
        [chord, span/2, 0.0],  # Left tip TE
        INVISCID)
    add_section!(wing, 
        [0.0, -span/2, 0.0],   # Right tip LE
        [chord, -span/2, 0.0], # Right tip TE
        INVISCID)
    
    body_aero = BodyAerodynamics([wing])
    unchanged_body_aero = BodyAerodynamics([unchanged_wing])
    init!(unchanged_body_aero)
    
    @testset "Re-initialization" begin
        result = @benchmark init!($unchanged_body_aero; init_aero=false) samples=1 evals=1
        @info "Re-initializing Allocations: $(result.allocs) \t Memory: $(result.memory)"
        @test result.allocs < 100
    end

    vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
    set_va!(body_aero, vel_app)

    # Initialize solvers for both LLT and VSM methods
    P = length(body_aero.panels)
    solver = Solver{P}()

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
                    @test result.allocs ≤ 100
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

        n_angles = 5
        alphas = collect(range(-deg2rad(10), deg2rad(10), n_angles))
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
                body_aero = BodyAerodynamics([wing])
                
                P = length(body_aero.panels)
                solver = Solver{P}(
                    aerodynamic_model_type=model
                )
                solver.sol._va_array .= va_array
                solver.sol._chord_array .= chord_array
                solver.sol._x_airf_array .= x_airf_array
                solver.sol._y_airf_array .= y_airf_array
                solver.sol._z_airf_array .= z_airf_array
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
        

        # # Fill arrays with data
        # for (i, panel) in enumerate(body_aero.panels)
        #     chord_array[i] = panel.chord
        #     x_airf_array[i, :] .= panel.x_airf
        #     y_airf_array[i, :] .= panel.y_airf
        #     z_airf_array[i, :] .= panel.z_airf
        #     va_array[i, :] .= panel.va
        # end
        set_va!(body_aero, vel_app)
        results = @MVector zeros(3)
        
        result = @benchmark calculate_results(
            $body_aero,
            $gamma,
            $reference_point,
            $density,
            VSM,
            1e-20,
            0.0,
            $alpha_array,
            $v_a_array,
            $chord_array,
            $x_airf_array,
            $y_airf_array,
            $z_airf_array,
            $va_array,
            $va_norm_array,
            $va_unit_array,
            $body_aero.panels,
            false
        ) samples=1 evals=1
        @test result.allocs ≤ 300
    end
end


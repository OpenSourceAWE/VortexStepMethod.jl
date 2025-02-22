using BenchmarkTools
using VortexStepMethod
using VortexStepMethod: calculate_AIC_matrices!, gamma_loop, calculate_results,
                       update_effective_angle_of_attack_if_VSM, calculate_projected_area,
                       calculate_cl, calculate_cd_cm,
                       calculate_velocity_induced_single_ring_semiinfinite,
                       calculate_velocity_induced_bound_2D,
                       velocity_3D_bound_vortex!,
                       velocity_3D_trailing_vortex!,
                       velocity_3D_trailing_vortex_semiinfinite!,
                       Panel
using Test
using LinearAlgebra

@testset "Function Allocation Tests" begin
    # Define wing parameters
    n_panels = 20          # Number of panels
    span = 20.0            # Wing span [m]
    chord = 1.0            # Chord length [m]
    v_a = 20.0            # Magnitude of inflow velocity [m/s]
    density = 1.225        # Air density [kg/m³]
    alpha_deg = 30.0       # Angle of attack [degrees]
    alpha = deg2rad(alpha_deg)
    
    # Create test panels
    panels = []
    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    add_section!(wing, 
        [0.0, span/2, 0.0],   # Left tip LE 
        [chord, span/2, 0.0],  # Left tip TE
        "inviscid")
    add_section!(wing, 
        [0.0, -span/2, 0.0],  # Right tip LE
        [chord, -span/2, 0.0], # Right tip TE
        "inviscid")
    
    wing_aero = WingAerodynamics([wing])

    vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
    set_va!(wing_aero, (vel_app, 0.0))  # Second parameter is yaw rate

    # Initialize solvers for both LLT and VSM methods
    solver = Solver()

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
    
    models = ["VSM", "LLT"]
    core_radius_fractions = [0.001, 10.0]

    @testset "AIC Matrix Calculation" begin
        for model in models
            for frac in core_radius_fractions
                @testset "Model $model Core Radius Fraction $frac" begin
                    result = @benchmark calculate_AIC_matrices!($wing_aero, $model, $frac, $va_norm_array, $va_unit_array)
                    @test result.allocs ≤ 100  # Allow some allocations for matrix setup
                end
            end
        end
    end
    
    # @testset "Gamma Loop" begin
    #     result = @benchmark gamma_loop($solver, $wing, $gamma_new, $AIC_x, $AIC_y, $AIC_z)
    #     @test result.allocs ≤ 50  # Main iteration should be mostly allocation-free
    # end
    
    # @testset "Results Calculation" begin
    #     result = @benchmark calculate_results($wing, $gamma)
    #     @test result.allocs ≤ 20  # Allow minimal allocations for results
    # end
    
    # @testset "Angle of Attack Update" begin
    #     result = @benchmark update_effective_angle_of_attack_if_VSM($wing, $gamma)
    #     @test result.allocs == 0  # Should be allocation-free
    # end
    
    # @testset "Area Calculations" begin
    #     result = @benchmark calculate_projected_area($wing)
    #     @test result.allocs ≤ 10  # Geometric calculations may need some allocations
    # end
    
    # @testset "Aerodynamic Coefficients" begin
    #     panel = panels[1]
    #     alpha = 0.1
        
    #     @test (@ballocated calculate_cl($panel, $alpha)) == 0
    #     @test (@ballocated calculate_cd_cm($panel, $alpha)) == 0
    # end
    
    # @testset "Induced Velocity Calculations" begin
    #     # Test single ring velocity calculation
    #     @test (@ballocated calculate_velocity_induced_single_ring_semiinfinite(
    #         $point, $panels[1], $gamma[1])) == 0
            
    #     # Test 2D bound vortex
    #     @test (@ballocated calculate_velocity_induced_bound_2D(
    #         $point, $panels[1], $gamma[1])) == 0
            
    #     # Test 3D velocity components
    #     @test (@ballocated velocity_3D_bound_vortex!(
    #         $v_ind, $point, $panels[1], $gamma[1])) == 0
            
    #     @test (@ballocated velocity_3D_trailing_vortex!(
    #         $v_ind, $point, $panels[1], $gamma[1])) == 0
            
    #     @test (@ballocated velocity_3D_trailing_vortex_semiinfinite!(
    #         $v_ind, $point, $panels[1], $gamma[1])) == 0
    # end
end
using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!, init!
using LinearAlgebra
using Test
using Logging

if !@isdefined ram_wing
    cp("data/ram_air_kite_body.obj", "/tmp/ram_air_kite_body.obj"; force=true)
    cp("data/ram_air_kite_foil.dat", "/tmp/ram_air_kite_foil.dat"; force=true)
    ram_wing = RamAirWing("/tmp/ram_air_kite_body.obj", "/tmp/ram_air_kite_foil.dat"; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1))
end

# @testset "Nonlinear vs Linear" begin
#     va = [15.0, 0.0, 0.0]
#     theta = zeros(4)
#     delta = zeros(4)
#     omega = zeros(3)

#     body_aero = BodyAerodynamics([ram_wing]; va)
#     solver = Solver(body_aero;
#         aerodynamic_model_type=VSM,
#         is_with_artificial_damping=false,
#         solver_type=NONLIN
#     )

#     jac, lin_res = VortexStepMethod.linearize(
#         solver, 
#         body_aero, 
#         ram_wing, 
#         [theta; va; omega; delta]; 
#         theta_idxs=1:4, 
#         va_idxs=5:7, 
#         omega_idxs=8:10,
#         delta_idxs=11:14,
#         moment_frac=0.1
#     )
#     nonlin_res = VortexStepMethod.solve!(solver, body_aero; log=true)
#     nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]

#     @test nonlin_res ≈ lin_res

#     dva = [0.01, 0.01, 0.01]
#     init!(body_aero; init_aero=false, va=va+dva)
#     nonlin_res = VortexStepMethod.solve!(solver, body_aero; log=true)
#     nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
#     @test lin_res ≉ nonlin_res atol=1e-2
#     @test lin_res + jac * [zeros(4); dva; zeros(3); zeros(4)] ≈ nonlin_res atol=1e-2
#     # Test accuracy
#     @test norm(lin_res + jac * [zeros(4); dva; zeros(3); zeros(4)] - nonlin_res) /
#         norm(lin_res - nonlin_res) < 2e-3
# end

@testset "Nonlinear vs Linear - Comprehensive Input Testing" begin
    # Initialize with base parameters
    va = [15.0, 0.0, 0.0]
    theta = zeros(4)
    delta = zeros(4)
    omega = zeros(3)
    
    # Define perturbation magnitudes
    dva_magnitudes = [0.01, 0.01, 0.01]  # Velocity perturbations (m/s)
    dtheta_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0), deg2rad(2.0)] # Angle perturbations (rad)
    ddelta_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0), deg2rad(2.0)] # Trailing edge perturbations (rad)
    domega_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0)]  # Angular rate perturbations (rad/s)
    
    # Create body aerodynamics and solver
    VortexStepMethod.group_deform!(ram_wing, theta, delta; smooth=false)
    body_aero = BodyAerodynamics([ram_wing]; va, omega)
    solver = Solver(body_aero;
        aerodynamic_model_type=VSM,
        is_with_artificial_damping=false,
        atol=1e-8,
        rtol=1e-8,
        solver_type=NONLIN
    )
    
    # Get initial Jacobian and linearization result
    base_inputs = [theta; va; omega; delta]
    jac, lin_res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        base_inputs; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        delta_idxs=11:14,
        moment_frac=0.1
    )
    
    # Verify that linearization results match nonlinear results at operating point
    baseline_res = VortexStepMethod.solve!(solver, body_aero; log=false)
    baseline_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
    @test baseline_res ≈ lin_res
    
    # Define test cases
    test_cases = [
        ("va", 5:7, dva_magnitudes),
        ("theta", 1:4, dtheta_magnitudes),
        ("omega", 8:10, domega_magnitudes),
        ("delta", 11:14, ddelta_magnitudes)
    ]
    
    # For each test case
    for (input_name, indices, magnitudes) in test_cases
        @testset "Input: $input_name" begin
            max_error_ratio = 0.0
            
            # Select one index at a time for focused testing
            for idx in indices
                # Test each magnitude
                for mag in magnitudes
                    # Prepare perturbation vector
                    perturbation = zeros(length(base_inputs))
                    
                    if length(indices) > 1 && input_name != "theta" && input_name != "delta"
                        # For vector quantities (va, omega), perturb all components together
                        perturbation[indices] .= mag
                    else
                        # For control angles, perturb one at a time
                        perturbation[idx] = mag
                    end
                    
                    # Reset to baseline
                    reset_va = copy(va)
                    reset_theta = copy(theta)
                    reset_delta = copy(delta)
                    reset_omega = zeros(3)
                    
                    # Apply the perturbation to the nonlinear model
                    if input_name == "va"
                        reset_va = va + perturbation[5:7]
                    elseif input_name == "theta"
                        reset_theta = theta + perturbation[1:4]
                    elseif input_name == "delta"
                        reset_delta = delta + perturbation[11:14]
                    elseif input_name == "omega"
                        reset_omega = perturbation[8:10]
                    else
                        throw(ArgumentError())
                    end
                    VortexStepMethod.group_deform!(ram_wing, reset_theta, reset_delta; smooth=false)
                    init!(body_aero; init_aero=false, va=reset_va, omega=reset_omega)
                    
                    # Get nonlinear solution
                    nonlin_res = VortexStepMethod.solve!(solver, body_aero, nothing; log=false)
                    nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
                    @test nonlin_res ≉ baseline_res
                    
                    # Compute linearized prediction
                    lin_prediction = lin_res + jac * perturbation
                    
                    # Calculate error ratio for this case
                    prediction_error = norm(lin_prediction - nonlin_res)
                    baseline_difference = norm(baseline_res - nonlin_res)
                    error_ratio = prediction_error / baseline_difference

                    # Update max error ratio
                    max_error_ratio = max(max_error_ratio, error_ratio)
                    
                    # For small perturbations, test that error ratio is small
                    if idx == first(indices) && mag == first(magnitudes)
                        @test error_ratio < 2e-3
                    end
                end
            end
            
            @info "Max error ratio for $input_name: $max_error_ratio"
        end
    end
    
    # Test combinations of input variations
    @testset "Combined Input Variations" begin
        # Sample some key combinations
        combination_tests = [
            ("Theta + VA", [dtheta_magnitudes[1] * ones(4); dva_magnitudes[1] * ones(3); zeros(3); zeros(4)]),
            ("Theta + Omega", [dtheta_magnitudes[1] * ones(4); zeros(3); domega_magnitudes[1] * ones(3); zeros(4)]),
            ("VA + Omega", [zeros(4); dva_magnitudes[1] * ones(3); domega_magnitudes[1] * ones(3); zeros(4)]),
            ("Delta + VA", [zeros(4); dva_magnitudes[1] * ones(3); zeros(3); ddelta_magnitudes[1] * ones(4)]),
            ("All Inputs", [dtheta_magnitudes[1] * ones(4); dva_magnitudes[1] * ones(3); 
                           domega_magnitudes[1] * ones(3); ddelta_magnitudes[1] * ones(4)])
        ]
        
        for (combo_name, perturbation) in combination_tests
            # Reset to baseline
            VortexStepMethod.group_deform!(ram_wing, theta, delta; smooth=false)
            init!(body_aero; init_aero=false, va=va, omega=zeros(3))
            
            # Extract components
            perturbed_theta = theta + perturbation[1:4]
            perturbed_va = va + perturbation[5:7]
            perturbed_omega = perturbation[8:10]
            perturbed_delta = delta + perturbation[11:14]
            
            # Apply to nonlinear model
            VortexStepMethod.group_deform!(ram_wing, perturbed_theta, perturbed_delta; smooth=false)
            init!(body_aero; init_aero=false, va=perturbed_va, omega=perturbed_omega)
            
            # Get nonlinear solution
            nonlin_res = VortexStepMethod.solve!(solver, body_aero; log=false)
            nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
            
            # Compute linearized prediction
            lin_prediction = lin_res + jac * perturbation
            
            # Calculate error ratio
            prediction_error = norm(lin_prediction - nonlin_res)
            baseline_difference = norm(lin_res - nonlin_res)
            error_ratio = prediction_error / baseline_difference
            
            @info "$combo_name error ratio: $error_ratio"
            
            # Test combinations
            @test lin_res ≉ nonlin_res atol=1e-2
            @test lin_prediction ≈ nonlin_res rtol=0.2 atol=0.05
            @test error_ratio < 2e-3
        end
    end
end
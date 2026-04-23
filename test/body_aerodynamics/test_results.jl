using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!, reinit!
using LinearAlgebra
using Test
using Logging

# Helper to find the data directory for test files
function _find_ram_data_dir()
    data_dir = joinpath(dirname(dirname(@__DIR__)), "data", "ram_air_kite")
    if isdir(data_dir)
        return data_dir
    else
        # Fallback for case where test is run from different working directory
        return joinpath(@__DIR__, "..", "..", "data", "ram_air_kite")
    end
end

if !@isdefined ram_wing_results
    data_dir = _find_ram_data_dir()
    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
    
    body_src = joinpath(data_dir, "ram_air_kite_body.obj")
    foil_src = joinpath(data_dir, "ram_air_kite_foil.dat")
    
    # Check if source files exist before copying
    if isfile(body_src) && isfile(foil_src)
        cp(body_src, body_path; force=true)
        cp(foil_src, foil_path; force=true)
    else
        error("Required data files not found: $body_src or $foil_src")
    end
    
    ram_wing = ObjWing(body_path, foil_path; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1), n_unrefined_sections=4)
end

@testset "Nonlinear vs Linear - Comprehensive Input Testing" begin
    # Initialize with base parameters
    va = [15.0, 1.0, 0.5]
    theta = deg2rad.([2.0, 1.0, -1.0, -2.0])
    delta = deg2rad.([1.0, 0.5, -0.5, -1.0])
    omega = [0.0, 0.1, 0.0]
    
    # Define perturbation magnitudes
    dva_magnitudes = [0.01, 0.01, 0.01]  # Velocity perturbations (m/s)
    dtheta_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0), deg2rad(2.0)] # Angle perturbations (rad)
    ddelta_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0), deg2rad(2.0)] # Trailing edge perturbations (rad)
    domega_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0)]  # Angular rate perturbations (rad/s)
    
    # Create body aerodynamics and solver
    VortexStepMethod.unrefined_deform!(ram_wing, theta, delta; smooth=false)
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

    coeff_jac, coeff_lin_res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        base_inputs; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        delta_idxs=11:14,
        moment_frac=0.1,
        aero_coeffs=true
    )
    
    # Verify that linearization results match nonlinear results at operating point
    baseline_res = VortexStepMethod.solve!(solver, body_aero; log=false)
    baseline_res = [solver.sol.force; solver.sol.moment; solver.sol.moment_unrefined_dist]
    coeff_baseline_res = [solver.sol.force_coeffs; solver.sol.moment_coeffs; solver.sol.cm_unrefined_dist]
    @test baseline_res ≈ lin_res
    @test coeff_baseline_res ≈ coeff_lin_res
    
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
            coeff_max_error_ratio = 0.0
            
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
                    reset_omega = copy(omega)
                    
                    # Apply the perturbation to the nonlinear model
                    if input_name == "va"
                        reset_va = va + perturbation[5:7]
                    elseif input_name == "theta"
                        reset_theta = theta + perturbation[1:4]
                    elseif input_name == "delta"
                        reset_delta = delta + perturbation[11:14]
                    elseif input_name == "omega"
                        reset_omega = omega + perturbation[8:10]
                    else
                        throw(ArgumentError())
                    end
                    VortexStepMethod.unrefined_deform!(ram_wing, reset_theta, reset_delta; smooth=false)
                    reinit!(body_aero; init_aero=false, va=reset_va, omega=reset_omega)
                    
                    # Get nonlinear solution
                    nonlin_res = VortexStepMethod.solve!(solver, body_aero, nothing; log=false)
                    nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.moment_unrefined_dist]
                    coeff_nonlin_res = [solver.sol.force_coeffs; solver.sol.moment_coeffs; solver.sol.cm_unrefined_dist]
                    @test nonlin_res ≉ baseline_res
                    @test coeff_nonlin_res ≉ baseline_res
                    
                    # Compute linearized prediction
                    lin_prediction = lin_res + jac * perturbation
                    coeff_lin_prediction = coeff_lin_res + coeff_jac * perturbation
                    
                    # Calculate error ratio for this case
                    prediction_error = norm(lin_prediction - nonlin_res)
                    baseline_difference = norm(baseline_res - nonlin_res)
                    error_ratio = prediction_error / baseline_difference
                    coeff_prediction_error = norm(coeff_lin_prediction - coeff_nonlin_res)
                    coeff_baseline_difference = norm(coeff_baseline_res - coeff_nonlin_res)
                    coeff_error_ratio = coeff_prediction_error / coeff_baseline_difference

                    # Update max error ratio
                    max_error_ratio = max(max_error_ratio, error_ratio)
                    coeff_max_error_ratio = max(coeff_max_error_ratio, coeff_error_ratio)
                    
                    # For small perturbations, test that error ratio is small
                    if idx == first(indices) && mag == first(magnitudes)
                        @test error_ratio < 2e-3
                    end
                end
            end
            
            @info "Max error ratio for $input_name: $max_error_ratio"
            @info "Max coeff error ratio for $input_name: $coeff_max_error_ratio"
        end
    end
    
    # Test combinations of input variations
    @testset "Combined Input Variations" begin
        # Smaller perturbations for combination tests to stay within linear regime.
        # The new deform! averages theta at section boundaries, introducing coupling
        # that makes combined perturbations less linear than individual ones.
        combo_scale = 0.25
        combo_dtheta = dtheta_magnitudes .* combo_scale
        combo_dva = dva_magnitudes .* combo_scale
        combo_domega = domega_magnitudes .* combo_scale
        combo_ddelta = ddelta_magnitudes .* combo_scale

        # For combination testing, we'll create targeted test cases
        # that use only the specific indices we want to test together
        combination_tests = [
            # Name, active indices, perturbation vector, indices mappings
            (
                "Theta + VA",
                # Use only theta and va indices
                [1:4; 5:7],
                # Perturbation values for these indices
                [combo_dtheta; combo_dva],
                # Mappings for linearize function
                (theta_idxs=1:4, va_idxs=5:7, omega_idxs=nothing, delta_idxs=nothing)
            ),
            (
                "Theta + Omega",
                [1:4; 8:10],
                [combo_dtheta; combo_domega],
                (theta_idxs=1:4, va_idxs=nothing, omega_idxs=5:7, delta_idxs=nothing)
            ),
            (
                "VA + Omega",
                [5:7; 8:10],
                [combo_dva; combo_domega],
                (theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6, delta_idxs=nothing)
            ),
            (
                "Delta + VA",
                [5:7; 11:14],
                [combo_dva; combo_ddelta],
                (theta_idxs=nothing, va_idxs=1:3, omega_idxs=nothing, delta_idxs=4:7)
            ),
            (
                "All Inputs",
                [1:4; 5:7; 8:10; 11:14],
                [combo_dtheta; combo_dva;
                combo_domega; combo_ddelta],
                (theta_idxs=1:4, va_idxs=5:7, omega_idxs=8:10, delta_idxs=11:14)
            )
        ]
        
        for (combo_name, active_indices, perturbation, idx_mappings) in combination_tests
            @testset "$combo_name" begin
                # Start with a fresh model for each combination test
                VortexStepMethod.unrefined_deform!(ram_wing, theta, delta; smooth=false)
                reinit!(body_aero; init_aero=false, va, omega)
                
                # Create the appropriate input vector for this combination
                input_vec = Vector{Float64}(undef, length(active_indices))
                
                # Fill with base values
                if !isnothing(idx_mappings.theta_idxs)
                    input_vec[idx_mappings.theta_idxs] .= theta
                end
                if !isnothing(idx_mappings.va_idxs)
                    input_vec[idx_mappings.va_idxs] .= va
                end
                if !isnothing(idx_mappings.omega_idxs)
                    input_vec[idx_mappings.omega_idxs] .= omega
                end
                if !isnothing(idx_mappings.delta_idxs)
                    input_vec[idx_mappings.delta_idxs] .= delta
                end
                
                # Get the Jacobian and linearization result for this specific combination
                jac_combo, lin_res_combo = VortexStepMethod.linearize(
                    solver, 
                    body_aero, 
                    input_vec; 
                    idx_mappings...
                )
                
                # Get baseline results
                baseline_res = VortexStepMethod.solve!(solver, body_aero; log=false)
                baseline_res = [solver.sol.force; solver.sol.moment; solver.sol.moment_unrefined_dist]
                
                # Should match the linearization result
                @test baseline_res ≈ lin_res_combo
                
                # Apply perturbation using the appropriate indices
                perturbed_input = copy(input_vec) + perturbation
                
                # Extract components based on the combination being tested
                perturbed_theta = !isnothing(idx_mappings.theta_idxs) ? 
                    perturbed_input[idx_mappings.theta_idxs] : theta
                
                perturbed_va = !isnothing(idx_mappings.va_idxs) ? 
                    perturbed_input[idx_mappings.va_idxs] : va
                
                perturbed_omega = !isnothing(idx_mappings.omega_idxs) ? 
                    perturbed_input[idx_mappings.omega_idxs] : omega
                
                perturbed_delta = !isnothing(idx_mappings.delta_idxs) ? 
                    perturbed_input[idx_mappings.delta_idxs] : delta
                
                # Apply to nonlinear model
                VortexStepMethod.unrefined_deform!(ram_wing, perturbed_theta, perturbed_delta; smooth=false)
                reinit!(body_aero; init_aero=false, va=perturbed_va, omega=perturbed_omega)
                
                # Get nonlinear solution with perturbation
                nonlin_res = VortexStepMethod.solve!(solver, body_aero; log=false)
                nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.moment_unrefined_dist]
                
                # Compute linearized prediction using our specialized Jacobian
                lin_prediction = lin_res_combo + jac_combo * perturbation
                
                # Ensure perturbation had an effect
                @test baseline_res ≉ nonlin_res atol=1e-3
                
                # Calculate error ratio
                prediction_error = norm(lin_prediction - nonlin_res)
                baseline_difference = norm(baseline_res - nonlin_res)
                error_ratio = prediction_error / baseline_difference
                
                @info "$combo_name error metrics" prediction_error baseline_difference error_ratio
                
                # Validate the prediction
                @test lin_prediction ≈ nonlin_res rtol=0.05 atol=1e-3
                @test error_ratio < 0.05
            end
        end
    end
end

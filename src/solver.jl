
"""
    Solver

Main solver structure for the Vortex Step Method.
"""
struct Solver
    # General settings
    aerodynamic_model_type::String
    density::Float64
    max_iterations::Int
    allowed_error::Float64
    tol_reference_error::Float64
    relaxation_factor::Float64
    
    # Damping settings
    is_with_artificial_damping::Bool
    artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}}
    
    # Additional settings
    type_initial_gamma_distribution::String
    core_radius_fraction::Float64
    mu::Float64
    is_only_f_and_gamma_output::Bool

    function Solver(;
        aerodynamic_model_type::String="VSM",
        density::Float64=1.225,
        max_iterations::Int=1500,
        allowed_error::Float64=1e-5,
        tol_reference_error::Float64=0.001,
        relaxation_factor::Float64=0.03,
        is_with_artificial_damping::Bool=false,
        artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}}=(k2=0.1, k4=0.0),
        type_initial_gamma_distribution::String="elliptic",
        core_radius_fraction::Float64=1e-20,
        mu::Float64=1.81e-5,
        is_only_f_and_gamma_output::Bool=false
    )
        new(
            aerodynamic_model_type,
            density,
            max_iterations,
            allowed_error,
            tol_reference_error,
            relaxation_factor,
            is_with_artificial_damping,
            artificial_damping,
            type_initial_gamma_distribution,
            core_radius_fraction,
            mu,
            is_only_f_and_gamma_output
        )
    end
end

"""
    solve(solver::Solver, wing_aero::WingAerodynamics, gamma_distribution=nothing)

Main solving routine for the aerodynamic model.
"""
function solve(solver::Solver, wing_aero::WingAerodynamics, gamma_distribution=nothing)
    isnothing(wing_aero.va) && throw(ArgumentError("Inflow conditions are not set"))

    # Initialize variables
    panels = wing_aero.panels
    n_panels = length(panels)
    alpha_array = zeros(n_panels)
    relaxation_factor = solver.relaxation_factor

    # Preallocate arrays
    x_airf_array = zeros(n_panels, 3)
    y_airf_array = zeros(n_panels, 3)
    z_airf_array = zeros(n_panels, 3)
    va_array = zeros(n_panels, 3)
    chord_array = zeros(n_panels)

    # Fill arrays from panels
    for (i, panel) in enumerate(panels)
        x_airf_array[i, :] .= panel.x_airf
        y_airf_array[i, :] .= panel.y_airf
        z_airf_array[i, :] .= panel.z_airf
        va_array[i, :] .= panel.va
        chord_array[i] = panel.chord
    end

    # Calculate unit vectors
    va_norm_array = norm.(eachrow(va_array))
    va_unit_array = va_array ./ va_norm_array

    # Calculate AIC matrices
    AIC_x, AIC_y, AIC_z = calculate_AIC_matrices(
        wing_aero,
        solver.aerodynamic_model_type,
        solver.core_radius_fraction,
        va_norm_array,
        va_unit_array
    )

    # Initialize gamma distribution
    gamma_initial = if isnothing(gamma_distribution)
        if solver.type_initial_gamma_distribution == "elliptic"
            calculate_circulation_distribution_elliptical_wing(wing_aero)
        else
            zeros(n_panels)
        end
    else
        length(gamma_distribution) == n_panels || 
            throw(ArgumentError("gamma_distribution length must match number of panels"))
        gamma_distribution
    end

    @debug "Initial gamma_new: $gamma_initial"

    # Run main iteration loop
    converged, gamma_new, alpha_array, Umag_array = gamma_loop(
        solver,
        gamma_initial,
        AIC_x,
        AIC_y,
        AIC_z,
        va_array,
        chord_array,
        x_airf_array,
        y_airf_array,
        z_airf_array,
        panels,
        relaxation_factor
    )

    # Try again with reduced relaxation factor if not converged
    if !converged && relaxation_factor > 1e-3
        @warn "Running again with half the relaxation_factor = $(relaxation_factor/2)"
        converged, gamma_new, alpha_array, Umag_array = gamma_loop(
            solver,
            gamma_initial,
            AIC_x,
            AIC_y,
            AIC_z,
            va_array,
            chord_array,
            x_airf_array,
            y_airf_array,
            z_airf_array,
            panels,
            relaxation_factor/2
        )
    end

    # Calculate final results
    results = calculate_results(
        wing_aero,
        gamma_new,
        solver.density,
        solver.aerodynamic_model_type,
        solver.core_radius_fraction,
        solver.mu,
        alpha_array,
        Umag_array,
        chord_array,
        x_airf_array,
        y_airf_array,
        z_airf_array,
        va_array,
        va_norm_array,
        va_unit_array,
        panels,
        solver.is_only_f_and_gamma_output
    )
    return results
end

"""
    gamma_loop(solver::Solver, gamma_new::Vector{Float64}, AIC_x::Matrix{Float64}, 
              AIC_y::Matrix{Float64}, AIC_z::Matrix{Float64}, va_array::Matrix{Float64}, 
              chord_array::Vector{Float64}, x_airf_array::Matrix{Float64}, 
              y_airf_array::Matrix{Float64}, z_airf_array::Matrix{Float64}, 
              panels::Vector{Panel}, relaxation_factor::Float64)

Main iteration loop for calculating circulation distribution.
"""
function gamma_loop(
    solver::Solver,
    gamma_new::Vector{Float64},
    AIC_x::Matrix{Float64},
    AIC_y::Matrix{Float64},
    AIC_z::Matrix{Float64},
    va_array::Matrix{Float64},
    chord_array::Vector{Float64},
    x_airf_array::Matrix{Float64},
    y_airf_array::Matrix{Float64},
    z_airf_array::Matrix{Float64},
    panels::Vector{Panel},
    relaxation_factor::Float64
)
    converged = false
    n_panels = length(panels)

    alpha_array = zeros(n_panels)
    Umag_array = zeros(n_panels)

    iters = 0
    for i in 1:solver.max_iterations
        iters += 1
        gamma = copy(gamma_new)
        
        # Calculate induced velocities
        induced_velocity_all = hcat(
            AIC_x * gamma,
            AIC_y * gamma,
            AIC_z * gamma
        )
        
        relative_velocity_array = va_array .+ induced_velocity_all
        relative_velocity_crossz = cross.(eachrow(relative_velocity_array), eachrow(z_airf_array))
        Uinfcrossz_array = cross.(eachrow(va_array), eachrow(z_airf_array))

        v_normal_array = vec(sum(x_airf_array .* relative_velocity_array, dims=2))
        v_tangential_array = vec(sum(y_airf_array .* relative_velocity_array, dims=2))
        alpha_array .= atan.(v_normal_array, v_tangential_array)
        
        Umag_array = norm.(relative_velocity_crossz)
        Umagw_array = norm.(Uinfcrossz_array)
        
        cl_array = [calculate_cl(panel, alpha) for (panel, alpha) in zip(panels, alpha_array)]
        gamma_new = 0.5 .* Umag_array.^2 ./ Umagw_array .* cl_array .* chord_array

        # Apply damping if needed
        if solver.is_with_artificial_damping
            damp, is_damping_applied = smooth_circulation(gamma, 0.1, 0.5)
            @debug "damp: $damp"
        else
            damp = zeros(length(gamma))
            is_damping_applied = false
        end

        # Update gamma with relaxation and damping
        gamma_new = (1 - relaxation_factor) .* gamma + 
                    relaxation_factor .* gamma_new .+ damp

        # Check convergence
        reference_error = maximum(abs.(gamma_new))
        reference_error = max(reference_error, solver.tol_reference_error)
        error = maximum(abs.(gamma_new - gamma))
        normalized_error = error / reference_error

        @debug "Iteration: $i, normalized_error: $normalized_error, is_damping_applied: $is_damping_applied"

        if normalized_error < solver.allowed_error
            converged = true
            break
        end
    end

    if converged
        @info "Converged after $iters iterations"
    else
        @warn "NO convergence after $(solver.max_iterations) iterations"
    end

    return converged, gamma_new, alpha_array, Umag_array
end

"""
    calculate_artificial_damping(solver::Solver, gamma::Vector{Float64}, 
                               alpha::Vector{Float64}, stall_angle_list::Vector{Float64})

Calculate artificial damping for numerical stability.
"""
function calculate_artificial_damping(
    solver::Solver,
    gamma::Vector{Float64},
    alpha::Vector{Float64},
    stall_angle_list::Vector{Float64}
)
    # Check for stall condition
    is_stalled = false
    if solver.aerodynamic_model_type == "LLT" || 
       (solver.artificial_damping.k2 == 0 && solver.artificial_damping.k4 == 0)
        return zeros(length(gamma)), is_stalled
    end

    for (ia, alpha_i) in enumerate(alpha)
        if alpha_i > stall_angle_list[ia]
            is_stalled = true
            break
        end
    end

    !is_stalled && return zeros(length(gamma)), is_stalled

    # Calculate damping
    n_gamma = length(gamma)
    damp = zeros(n_gamma)

    for ig in 1:n_gamma
        # Handle boundary cases
        gim2, gim1, gi, gip1, gip2 = if ig == 1
            gamma[1], gamma[1], gamma[1], gamma[2], gamma[3]
        elseif ig == 2
            gamma[1], gamma[1], gamma[2], gamma[3], gamma[4]
        elseif ig == n_gamma - 1
            gamma[n_gamma-3], gamma[n_gamma-2], gamma[n_gamma-1], gamma[n_gamma], gamma[n_gamma]
        elseif ig == n_gamma
            gamma[n_gamma-2], gamma[n_gamma-1], gamma[n_gamma], gamma[n_gamma], gamma[n_gamma]
        else
            gamma[ig-2], gamma[ig-1], gamma[ig], gamma[ig+1], gamma[ig+2]
        end

        dif2 = (gip1 - gi) - (gi - gim1)
        dif4 = (gip2 - 3.0 * gip1 + 3.0 * gi - gim1) - 
               (gip1 - 3.0 * gi + 3.0 * gim1 - gim2)
        
        damp[ig] = solver.artificial_damping.k2 * dif2 - 
                   solver.artificial_damping.k4 * dif4
    end

    return damp, is_stalled
end

"""
    smooth_circulation(circulation::Vector{Float64}, 
                      smoothness_factor::Float64, 
                      damping_factor::Float64)

Smooth circulation distribution if needed.

Returns:
- Tuple of smoothed circulation and boolean indicating if smoothing was applied
"""
function smooth_circulation(
    circulation::Vector{Float64},
    smoothness_factor::Float64,
    damping_factor::Float64
)
    # Calculate mean circulation excluding endpoints
    circulation_mean = mean(circulation[2:end-1])
    smoothness_threshold = smoothness_factor * circulation_mean

    # Calculate differences between adjacent points
    differences = diff(circulation[2:end-1])
    @debug "circulation_mean: $circulation_mean, diff: $differences"

    # Check smoothness
    if isempty(differences)
        return zeros(length(circulation)), false
    end

    if maximum(abs.(differences)) <= smoothness_threshold
        return zeros(length(circulation)), false
    end

    # Apply smoothing
    smoothed = copy(circulation)
    for i in 2:length(circulation)-1
        left = circulation[i-1]
        center = circulation[i]
        right = circulation[i+1]
        avg = (left + right) / 2
        smoothed[i] = center + damping_factor * (avg - center)
    end

    # Preserve total circulation
    total_original = sum(circulation)
    total_smoothed = sum(smoothed)
    smoothed .*= total_original / total_smoothed

    damp = smoothed - circulation
    return damp, true
end
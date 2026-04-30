
"""
    VSMSolution

Struct for storing the solution of the [solve!](@ref) function. Must contain all info needed by `KiteModels.jl`.

# Naming Convention
- Variables ending in `_dist`: Per-panel distributions (length P, one value per panel)
- Variables ending in `_unrefined_dist`: Per-unrefined-section distributions (length U, averaged values per unrefined section)

# Attributes
- `width_dist`::Vector{Float64}: Width of the panels [m]
- `alpha_dist`::Vector{Float64}: Angle of attack of each panel relative to the apparent wind [rad]
- cl_dist::Vector{Float64}: Lift coefficients of the panels [-]
- cd_dist::Vector{Float64}: Drag coefficients of the panels [-]
- cm_dist::Vector{Float64}: Pitching moment coefficients of the panels [-]
- lift_dist::Vector{Float64}: Lift force of the panels [N]
- drag_dist::Vector{Float64}: Drag force of the panels [N]
- panel_moment_dist::Vector{Float64}: Pitching moment around the spanwise vector of the panels [Nm]
- `f_body_3D`::Matrix{Float64}: Matrix of the aerodynamic forces (x, y, z vectors) [N]
- `m_body_3D`::Matrix{Float64}: Matrix of the aerodynamic moments [Nm]
- `gamma_distribution`::Union{Nothing, Vector{Float64}}: Vector containing the panel circulations.
- force::MVec3: Aerodynamic force vector in KB reference frame [N]
- moment::MVec3: Aerodynamic moments [Mx, My, Mz] around the reference point [Nm]
- force_coeffs::MVec3: Aerodynamic force coefficients [CFx, CFy, CFz] [-]
- `moment_coeffs`::MVec3: Aerodynamic moment coefficients [CMx, CMy, CMz] [-]
- `moment_dist`::Vector{Float64}: Pitching moments around the spanwise vector of each panel. [Nm]
- `moment_coeff_dist`::Vector{Float64}: Pitching moment coefficient around the spanwise vector of each panel. [-]
- `moment_unrefined_dist`::MVector{U, Float64}: Averaged moments for unrefined sections [Nm]
- `cl_unrefined_dist`::MVector{U, Float64}: Averaged lift coefficients for unrefined sections [-]
- `cd_unrefined_dist`::MVector{U, Float64}: Averaged drag coefficients for unrefined sections [-]
- `cm_unrefined_dist`::MVector{U, Float64}: Averaged moment coefficients for unrefined sections [-]
- `alpha_unrefined_dist`::MVector{U, Float64}: Averaged angles of attack for unrefined sections [rad]
- `solver_status`::SolverStatus: enum, see [SolverStatus](@ref)
"""
@with_kw mutable struct VSMSolution{P,U}
    ### private vectors of solve_base!
    _x_airf_dist::Matrix{Float64} = zeros(P, 3)
    _y_airf_dist::Matrix{Float64} = zeros(P, 3)
    _z_airf_dist::Matrix{Float64} = zeros(P, 3)
    _va_dist::Matrix{Float64} = zeros(P, 3)
    _chord_dist::Vector{Float64} = zeros(P)
    ### end of private vectors
    width_dist::Vector{Float64} = zeros(P)
    alpha_dist::Vector{Float64} = zeros(P)
    alpha_geometric_dist::Vector{Float64} = zeros(P)
    cl_dist::Vector{Float64} = zeros(P)
    cd_dist::Vector{Float64} = zeros(P)
    cm_dist::Vector{Float64} = zeros(P)
    lift_dist::Vector{Float64} = zeros(P)
    drag_dist::Vector{Float64} = zeros(P)
    panel_moment_dist::Vector{Float64} = zeros(P)
    f_body_3D::Matrix{Float64} = zeros(3, P)
    m_body_3D::Matrix{Float64} = zeros(3, P)
    gamma_distribution::Union{Nothing, Vector{Float64}} = nothing
    force::MVec3 = zeros(MVec3)          
    moment::MVec3 = zeros(MVec3)       
    force_coeffs::MVec3 = zeros(MVec3)  
    moment_coeffs::MVec3 = zeros(MVec3)  
    center_of_pressure::Union{Nothing, MVec3} = nothing
    panel_cp_locations::Vector{MVec3} = MVec3[]
    moment_dist::MVector{P, Float64} = zeros(P)
    moment_coeff_dist::MVector{P, Float64} = zeros(P)
    moment_unrefined_dist::MVector{U, Float64} = zeros(U)
    cl_unrefined_dist::MVector{U, Float64} = zeros(U)
    cd_unrefined_dist::MVector{U, Float64} = zeros(U)
    cm_unrefined_dist::MVector{U, Float64} = zeros(U)
    alpha_unrefined_dist::MVector{U, Float64} = zeros(U)
    x_airf_unrefined_dist::Vector{MVec3} = [MVec3(0,0,0) for _ in 1:U]
    y_airf_unrefined_dist::Vector{MVec3} = [MVec3(0,0,0) for _ in 1:U]
    z_airf_unrefined_dist::Vector{MVec3} = [MVec3(0,0,0) for _ in 1:U]
    va_unrefined_dist::Vector{MVec3} = [MVec3(0,0,0) for _ in 1:U]
    chord_unrefined_dist::MVector{U, Float64} = zeros(U)
    width_unrefined_dist::MVector{U, Float64} = zeros(U)
    solver_status::SolverStatus = FAILURE
end

# Output of the function gamma_loop!
@with_kw mutable struct LoopResult{P}
    converged::Bool              = false
    gamma_new::MVector{P, Float64}   = zeros(P)
    alpha_dist::MVector{P, Float64} = zeros(P) # TODO: Is this different from BodyAerodynamics.alpha_dist ?
    v_a_dist::MVector{P, Float64}   = zeros(P)
end

@with_kw struct BaseResult{P}
    va_norm_dist::MVector{P, Float64} = zeros(P)
    va_unit_dist::Matrix{Float64} = zeros(P, 3)
end

@inline function check_reference_point(reference_point)
    msg = "reference_point must be a list/array with 3 numbers."
    reference_point isa AbstractVector || throw(ArgumentError(msg))
    length(reference_point) == 3 || throw(ArgumentError(msg))
    all(x -> x isa Number, reference_point) || throw(ArgumentError(msg))
    try
        return MVec3(Float64(reference_point[1]), Float64(reference_point[2]), Float64(reference_point[3]))
    catch
        throw(ArgumentError(msg))
    end
end

"""
    Solver

Main solver structure for the Vortex Step Method.See also: [solve](@ref)

# Attributes

## General settings
- `aerodynamic_model_type`::Model = VSM: The model type, see: [Model](@ref)
- density::Float64 = 1.225: Air density [kg/m³] 
- `max_iterations`::Int64 = 1500
- `rtol`::Float64 = 1e-5: relative error
- `tol_reference_error`::Float64 = 0.001
- `relaxation_factor`::Float64 = 0.03: Relaxation factor for convergence 

## Damping settings
- `is_with_artificial_damping`::Bool = false: Whether to apply artificial damping
- `artificial_damping`::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}} = (k2=0.1, k4=0.0): Artificial damping parameters

## Additional settings
- `type_initial_gamma_distribution`::InitialGammaDistribution = ELLIPTIC: see: [InitialGammaDistribution](@ref)
- `use_gamma_prev`::Bool = true: reuse provided previous gamma as initial guess when available
- `core_radius_fraction`::Float64 = 1e-20: 
- mu::Float64 = 1.81e-5: Dynamic viscosity [N·s/m²]
- `is_only_f_and_gamma_output`::Bool = false: Whether to only output f and gamma
- `reference_point`::MVec3 = [0.0, 0.0, 0.0]: Moment reference point in body frame

## Solution
sol::VSMSolution = VSMSolution(): The result of calling [solve!](@ref) 
"""
@with_kw mutable struct Solver{P,U}
    # General settings
    solver_type::SolverType = LOOP
    aerodynamic_model_type::Model = VSM
    density::Float64 = 1.225
    max_iterations::Int64 = 1500
    rtol::Float64 = 1e-5
    tol_reference_error::Float64 = 0.001
    relaxation_factor::Float64 = 0.03

    # Nonlin solver fields
    atol::Float64 = 1e-5
    nonlin_jac::Matrix{Float64} = zeros(P, P)
    nonlin_residual::MVector{P, Float64} = zeros(P)
    nonlin_residual_perturbed::MVector{P, Float64} = zeros(P)
    nonlin_gamma_perturbed::MVector{P, Float64} = zeros(P)
    
    # Damping settings
    is_with_artificial_damping::Bool = false
    artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}} =(k2=0.1, k4=0.0)
    
    # Additional settings
    type_initial_gamma_distribution::InitialGammaDistribution = ZEROS
    use_gamma_prev::Bool = true
    core_radius_fraction::Float64 = 0.05
    mu::Float64 = 1.81e-5
    is_only_f_and_gamma_output::Bool = false
    correct_aoa::Bool = false
    reference_point::MVec3 = zeros(MVec3)

    # Intermediate results
    lr::LoopResult{P} = LoopResult{P}()
    br::BaseResult{P} = BaseResult{P}()
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:11]
    cache_base::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}  = [LazyBufferCache()]
    cache_lin::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:4]
    
    # Solution
    sol::VSMSolution{P,U} = VSMSolution{P,U}()
end

function Solver(body_aero; reference_point=[0.0, 0.0, 0.0], kwargs...)
    P = length(body_aero.panels)
    U = sum([wing.n_unrefined_sections for wing in body_aero.wings])
    reference_point_checked = check_reference_point(reference_point)
    return Solver{P,U}(; reference_point=reference_point_checked, kwargs...)
end

function Solver(body_aero, settings::VSMSettings)
    ss = settings.solver_settings
    solver_type = ss.solver_type == "NONLIN" ? NONLIN : LOOP
    reference_point = hasproperty(ss, :reference_point) ? ss.reference_point : [0.0, 0.0, 0.0]
    Solver(body_aero;
        solver_type,
        aerodynamic_model_type=ss.aerodynamic_model_type,
        density=ss.density,
        max_iterations=ss.max_iterations,
        rtol=ss.rtol,
        tol_reference_error=ss.tol_reference_error,
        relaxation_factor=ss.relaxation_factor,
        is_with_artificial_damping=ss.artificial_damping,
        artificial_damping=(k2=ss.k2, k4=ss.k4),
        type_initial_gamma_distribution=ss.type_initial_gamma_distribution,
        use_gamma_prev=ss.use_gamma_prev,
        core_radius_fraction=ss.core_radius_fraction,
        mu=ss.mu,
        is_only_f_and_gamma_output=ss.calc_only_f_and_gamma,
        correct_aoa=ss.correct_aoa,
        reference_point=reference_point,
    )
end

"""
    solve!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=solver.sol.gamma_distribution; 
          log=false, reference_point=solver.reference_point, moment_frac=0.1)

Main solving routine for the aerodynamic model. Reference point is in the kite body (KB) frame.
This version is modifying the `solver.sol` struct and is faster than the `solve` function which returns
a dictionary.

# Arguments:
- solver::Solver: The solver to use, could be a VSM or LLT solver. See: [Solver](@ref)
- body_aero::BodyAerodynamics: The aerodynamic body. See: [BodyAerodynamics](@ref)
- gamma_distribution: Initial circulation vector or nothing; Length: Number of segments. [m²/s]

# Keyword Arguments:
- log=false: If true, print the number of iterations and other info.
- reference_point=solver.reference_point
- moment_frac=0.1: X-coordinate of normalized panel around which the moment distribution should be calculated.

# Returns
The solution of type [VSMSolution](@ref)
"""
function solve!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=solver.sol.gamma_distribution; 
        log=false, reference_point=solver.reference_point, moment_frac=0.1)

    # calculate intermediate result
    solve_base!(solver, body_aero, gamma_distribution; log)
    gamma_new = solver.lr.gamma_new
    if !isnothing(solver.sol.gamma_distribution)
        solver.sol.gamma_distribution .= gamma_new
    else
        solver.sol.gamma_distribution = gamma_new
    end

    # Initialize arrays
    cl_dist = solver.sol.cl_dist
    cd_dist = solver.sol.cd_dist
    cm_dist = solver.sol.cm_dist
    converged = solver.lr.converged
    alpha_dist = solver.lr.alpha_dist
    alpha_corrected = solver.sol.alpha_dist
    alpha_geometric_dist = solver.sol.alpha_geometric_dist
    v_a_dist = solver.lr.v_a_dist
    panels = body_aero.panels
   
    width_dist = solver.sol.width_dist
    solver.sol.moment_dist .= 0
    solver.sol.moment_coeff_dist .= 0
    moment_dist = solver.sol.moment_dist
    moment_coeff_dist = solver.sol.moment_coeff_dist
    density = solver.density
    aerodynamic_model_type = solver.aerodynamic_model_type

    # Calculate coefficients for each panel
    for (i, panel) in enumerate(panels)                                               # zero bytes
        cl_dist[i] = calculate_cl(panel, alpha_dist[i])
        cd_dist[i], cm_dist[i] = calculate_cd_cm(panel, alpha_dist[i])
        width_dist[i] = panel.width

        # Geometric AoA using panel-local axes and prescribed
        # freestream — scalar ops to avoid allocations
        begin
            va1 = solver.sol._va_dist[i,1]
            va2 = solver.sol._va_dist[i,2]
            va3 = solver.sol._va_dist[i,3]
            va_norm = sqrt(va1^2 + va2^2 + va3^2)
            x1 = solver.sol._x_airf_dist[i,1]
            x2 = solver.sol._x_airf_dist[i,2]
            x3 = solver.sol._x_airf_dist[i,3]
            x_norm = sqrt(x1^2 + x2^2 + x3^2)
            z1 = solver.sol._z_airf_dist[i,1]
            z2 = solver.sol._z_airf_dist[i,2]
            z3 = solver.sol._z_airf_dist[i,3]
            z_norm = sqrt(z1^2 + z2^2 + z3^2)
            if va_norm == 0 || x_norm == 0 || z_norm == 0
                alpha_geometric_dist[i] = NaN
            else
                inv_va = -1.0 / va_norm
                vu1 = va1 * inv_va
                vu2 = va2 * inv_va
                vu3 = va3 * inv_va
                v_tangential = (x1*vu1+x2*vu2+x3*vu3) / x_norm
                v_normal = (z1*vu1+z2*vu2+z3*vu3) / z_norm
                alpha_geometric_dist[i] = pi + atan(
                    v_normal, v_tangential)
            end
        end

    end

    # create an alias for the three vertical output vectors
    lift = solver.sol.lift_dist
    drag = solver.sol.drag_dist
    panel_moment_dist = solver.sol.panel_moment_dist

    # Compute using fused broadcasting (no intermediate allocations)
    @. lift = cl_dist * 0.5 * density * v_a_dist^2 * solver.sol._chord_dist
    @. drag = cd_dist * 0.5 * density * v_a_dist^2 * solver.sol._chord_dist
    @. panel_moment_dist = cm_dist * 0.5 * density * v_a_dist^2 * solver.sol._chord_dist^2

    # Calculate alpha corrections based on model type
    if solver.correct_aoa && aerodynamic_model_type == VSM      # 64 bytes
        update_effective_angle_of_attack!(
            alpha_corrected,
            body_aero,
            gamma_new,
            solver.core_radius_fraction,
            solver.sol._z_airf_dist,
            solver.sol._x_airf_dist,
            solver.sol._va_dist,
            solver.br.va_norm_dist,
            solver.br.va_unit_dist
        )
    else
        alpha_corrected .= alpha_dist
    end

    # Initialize result arrays
    area_all_panels = 0.0
    panel_areas = zeros(length(panels))

    # Get wing properties
    spanwise_direction = body_aero.wings[1].spanwise_direction

    # Calculate wing geometry properties
    projected_area = body_aero.projected_area
    c_ref = body_aero.c_ref
    
    wv = body_aero.work_vectors
    dir_iva = wv[1]
    dir_lift = wv[2]
    dir_drag = wv[3]
    lift_va = wv[4]
    drag_va = wv[5]
    r_vec = wv[6]
    f_tmp = wv[7]
    cross_tmp = wv[8]

    for (i, panel) in enumerate(panels)

        ### Lift and Drag ###
        panel_area = panel.chord * panel.width
        area_all_panels += panel_area
        panel_areas[i] = panel_area

        # Calculate induced velocity direction
        alpha_corrected_i = alpha_corrected[i]
        c_alpha = cos(alpha_corrected_i)
        s_alpha = sin(alpha_corrected_i)
        @inbounds for k in 1:3
            dir_iva[k] = c_alpha * panel.x_airf[k] +
                         s_alpha * panel.z_airf[k]
        end
        normalize3!(dir_iva)

        # Calculate lift and drag directions
        cross3!(dir_lift, dir_iva, panel.y_airf)
        normalize3!(dir_lift)
        cross3!(dir_drag, spanwise_direction, dir_lift)
        normalize3!(dir_drag)

        # Calculate force vectors
        li = lift[i]
        di = drag[i]
        @inbounds for k in 1:3
            lift_va[k] = li * dir_lift[k]
            drag_va[k] = di * dir_drag[k]
        end

        # Body frame forces
        width = panel.width
        @inbounds for k in 1:3
            solver.sol.f_body_3D[k, i] = (lift_va[k] +
                                          drag_va[k]) * width
        end

        # Calculate the moments
        m_scale = panel_moment_dist[i] * width
        @inbounds for k in 1:3
            r_vec[k] = panel.aero_center[k] -
                       reference_point[k]
            f_tmp[k] = solver.sol.f_body_3D[k, i]
        end
        cross3!(cross_tmp, r_vec, f_tmp)
        @inbounds for k in 1:3
            solver.sol.m_body_3D[k, i] = m_scale *
                panel.y_airf[k] + cross_tmp[k]
        end

        # Moment distribution (moment on each panel)
        arm = (moment_frac - 0.25) * panel.chord
        ftotal_dot_z = 0.0
        @inbounds for k in 1:3
            ftotal_dot_z += (lift_va[k] + drag_va[k]) *
                            panel.z_airf[k]
        end
        moment_dist[i] = (ftotal_dot_z * arm +
                          panel_moment_dist[i]) * width
    end

    # Python parity: normalize with area-weighted reference velocity for distributed inflow.
    va_ref_vector = _compute_reference_velocity_from_distribution(
        solver.sol._va_dist,
        length(panels),
        panel_areas
    )
    va_ref_mag = norm(va_ref_vector)
    va_ref_mag > 0.0 || throw(ArgumentError("Reference freestream magnitude must be positive."))
    q_ref = 0.5 * density * va_ref_mag^2
    moment_coeff_dist .= moment_dist ./ (q_ref * projected_area * c_ref)

    # Only compute unrefined arrays if there are unrefined sections
    if length(solver.sol.moment_unrefined_dist) > 0
        moment_unrefined_dist = solver.sol.moment_unrefined_dist
        cl_unrefined_dist = solver.sol.cl_unrefined_dist
        cd_unrefined_dist = solver.sol.cd_unrefined_dist
        cm_unrefined_dist = solver.sol.cm_unrefined_dist
        alpha_unrefined_dist = solver.sol.alpha_unrefined_dist
        x_airf_unrefined_dist = solver.sol.x_airf_unrefined_dist
        y_airf_unrefined_dist = solver.sol.y_airf_unrefined_dist
        z_airf_unrefined_dist = solver.sol.z_airf_unrefined_dist
        va_unrefined_dist = solver.sol.va_unrefined_dist
        chord_unrefined_dist = solver.sol.chord_unrefined_dist
        width_unrefined_dist = solver.sol.width_unrefined_dist

        # Zero all unrefined arrays
        moment_unrefined_dist .= 0.0
        cl_unrefined_dist .= 0.0
        cd_unrefined_dist .= 0.0
        cm_unrefined_dist .= 0.0
        alpha_unrefined_dist .= 0.0
        for i in eachindex(x_airf_unrefined_dist)
            x_airf_unrefined_dist[i] .= 0.0
            y_airf_unrefined_dist[i] .= 0.0
            z_airf_unrefined_dist[i] .= 0.0
            va_unrefined_dist[i] .= 0.0
        end
        chord_unrefined_dist .= 0.0
        width_unrefined_dist .= 0.0

        panel_idx = 1
        unrefined_idx = 1
        for wing in body_aero.wings
            if wing.n_unrefined_sections > 0
                # Accumulate values from refined panels to unrefined sections
                unrefined_section_counts = zeros(Int, wing.n_unrefined_sections)
                for local_panel_idx in 1:wing.n_panels
                    panel = body_aero.panels[panel_idx]
                    original_section_idx = wing.refined_panel_mapping[local_panel_idx]
                    target_unrefined_idx = unrefined_idx + original_section_idx - 1

                    # Accumulate coefficients and moments
                    moment_unrefined_dist[target_unrefined_idx] += moment_dist[panel_idx]
                    cl_unrefined_dist[target_unrefined_idx] += solver.sol.cl_dist[panel_idx]
                    cd_unrefined_dist[target_unrefined_idx] += solver.sol.cd_dist[panel_idx]
                    cm_unrefined_dist[target_unrefined_idx] += solver.sol.cm_dist[panel_idx]
                    alpha_unrefined_dist[target_unrefined_idx] += solver.sol.alpha_dist[panel_idx]

                    # Accumulate geometry
                    x_airf_unrefined_dist[target_unrefined_idx] .+= panel.x_airf
                    y_airf_unrefined_dist[target_unrefined_idx] .+= panel.y_airf
                    z_airf_unrefined_dist[target_unrefined_idx] .+= panel.z_airf
                    va_unrefined_dist[target_unrefined_idx] .+= panel.va
                    chord_unrefined_dist[target_unrefined_idx] += panel.chord
                    width_unrefined_dist[target_unrefined_idx] += panel.width

                    unrefined_section_counts[original_section_idx] += 1
                    panel_idx += 1
                end

                # Average coefficients and geometry (width stays summed)
                for i in 1:wing.n_unrefined_sections
                    target_unrefined_idx = unrefined_idx + i - 1
                    if unrefined_section_counts[i] > 0
                        count = unrefined_section_counts[i]
                        moment_unrefined_dist[target_unrefined_idx] /= count
                        cl_unrefined_dist[target_unrefined_idx] /= count
                        cd_unrefined_dist[target_unrefined_idx] /= count
                        cm_unrefined_dist[target_unrefined_idx] /= count
                        alpha_unrefined_dist[target_unrefined_idx] /= count
                        x_airf_unrefined_dist[target_unrefined_idx] ./= count
                        y_airf_unrefined_dist[target_unrefined_idx] ./= count
                        z_airf_unrefined_dist[target_unrefined_idx] ./= count
                        va_unrefined_dist[target_unrefined_idx] ./= count
                        chord_unrefined_dist[target_unrefined_idx] /= count
                        # width_unrefined_dist is NOT averaged - it is the
                        # sum of panel widths in the unrefined section
                    end
                end
                unrefined_idx += wing.n_unrefined_sections
            else
                # Skip panels for wings with no unrefined sections
                panel_idx += wing.n_panels
            end
        end
    end

    # update the result struct
    solver.sol.force .= 0.0
    solver.sol.moment .= 0.0
    @inbounds for i in 1:length(panels)
        for k in 1:3
            solver.sol.force[k] += solver.sol.f_body_3D[k, i]
            solver.sol.moment[k] += solver.sol.m_body_3D[k, i]
        end
    end
    solver.sol.force_coeffs .= solver.sol.force ./ (q_ref * projected_area)
    solver.sol.moment_coeffs .= solver.sol.moment ./ (q_ref * projected_area * c_ref)
    # Keep solve! fast: center-of-pressure is only computed in solve() dictionary path.
    solver.sol.center_of_pressure = nothing
    empty!(solver.sol.panel_cp_locations)
    if converged
        # TODO: Check if the result if feasible if converged
        solver.sol.solver_status = FEASIBLE
    else
        solver.sol.solver_status = FAILURE
    end

    return solver.sol
end

"""
    solve(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
          log=false, reference_point=solver.reference_point)

Main solving routine for the aerodynamic model. Reference point is in the kite body (KB) frame.
See also: [solve!](@ref)

# Arguments:
- solver::Solver: The solver to use, could be a VSM or LLT solver. See: [Solver](@ref)
- body_aero::BodyAerodynamics: The aerodynamic body. See: [BodyAerodynamics](@ref)
- gamma_distribution: Initial circulation vector or nothing; Length: Number of segments. [m²/s]

# Keyword Arguments:
- log=false: If true, print the number of iterations and other info.
- reference_point=solver.reference_point

# Returns
A dictionary with the results.
"""
function solve(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
    log=false, reference_point=solver.reference_point)
    reference_point_checked = check_reference_point(reference_point)
    # calculate intermediate result
    solve_base!(solver, body_aero, gamma_distribution; log)

    # Calculate final results as dictionary
    results = calculate_results(
        body_aero,
        solver.lr.gamma_new,
        reference_point_checked,
        solver.density,
        solver.core_radius_fraction,
        solver.mu,
        solver.lr.alpha_dist,
        solver.lr.v_a_dist,
        solver.sol._chord_dist,
        solver.sol._x_airf_dist,
        solver.sol._z_airf_dist,
        solver.sol._va_dist,
        solver.br.va_norm_dist,
        solver.br.va_unit_dist,
        body_aero.panels,
        solver.is_only_f_and_gamma_output;
        correct_aoa=solver.correct_aoa
    )
    # Attach geometric AoA (already computed in calculate_results) to solver.sol
    if haskey(results, "alpha_geometric")
        solver.sol.alpha_geometric_dist .= results["alpha_geometric"]
    end
    return results
end

@inline @inbounds function calc_norm_array!(va_norm_dist, va_array)
    for i in axes(va_array, 1)
        va_norm_dist[i] = sqrt(
            va_array[i,1]^2 + va_array[i,2]^2 +
            va_array[i,3]^2)
    end
end

function solve_base!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing;
               log=false)
    
    # check arguments
    isnothing(body_aero.panels[1].va) && throw(ArgumentError("Inflow conditions are not set, use set_va!(body_aero, va)"))
    
    # Initialize variables
    panels = body_aero.panels
    n_panels = length(panels)
    relaxation_factor = solver.relaxation_factor
    
    # Clear arrays
    solver.sol._x_airf_dist .= 0
    solver.sol._y_airf_dist .= 0
    solver.sol._z_airf_dist .= 0
    solver.sol._va_dist .= 0
    solver.sol._chord_dist .= 0

    # Fill arrays from panels
    for (i, panel) in enumerate(panels)
        @inbounds for k in 1:3
            solver.sol._x_airf_dist[i, k] = panel.x_airf[k]
            solver.sol._y_airf_dist[i, k] = panel.y_airf[k]
            solver.sol._z_airf_dist[i, k] = panel.z_airf[k]
            solver.sol._va_dist[i, k] = panel.va[k]
        end
        solver.sol._chord_dist[i] = panel.chord
    end

    # Calculate unit vectors
    calc_norm_array!(solver.br.va_norm_dist, solver.sol._va_dist)
    @inbounds for i in 1:n_panels
        inv_norm = 1.0 / solver.br.va_norm_dist[i]
        for k in 1:3
            solver.br.va_unit_dist[i, k] =
                solver.sol._va_dist[i, k] * inv_norm
        end
    end

    # Calculate AIC matrices
    calculate_AIC_matrices!(body_aero, solver.aerodynamic_model_type, solver.core_radius_fraction, solver.br.va_norm_dist,
                            solver.br.va_unit_dist)

    # Initialize gamma distribution
    gamma_initial = solver.cache_base[1][solver.sol._chord_dist]
    if isnothing(gamma_distribution) || !solver.use_gamma_prev
        if !isnothing(gamma_distribution) && !solver.use_gamma_prev
            @debug "Ignoring provided gamma_distribution because use_gamma_prev=false"
        end
        if solver.type_initial_gamma_distribution == ELLIPTIC
            calculate_circulation_distribution_elliptical_wing(gamma_initial, body_aero)
        else
            gamma_initial .= 0
        end
    else
        length(gamma_distribution) == n_panels || 
            throw(ArgumentError("gamma_distribution length must match number of panels"))
        gamma_initial .= gamma_distribution
    end

    @debug "Initial gamma_new: $gamma_initial"
    solver.lr.gamma_new .= gamma_initial
    # Run main iteration loop
    gamma_loop!(solver, body_aero, panels, relaxation_factor; log)
    # Try again with reduced relaxation factor if not converged
    if solver.solver_type == LOOP && !solver.lr.converged && relaxation_factor > 1e-3
        log && @warn "Running again with half the relaxation_factor = $(relaxation_factor/2)"
        solver.lr.gamma_new .= gamma_initial
        gamma_loop!(solver, body_aero, panels, relaxation_factor/2; log)
    end

    nothing
end

@inline function update_gamma_candidate!(
    gamma_out,
    gamma_in,
    solver::Solver,
    panels::Vector{Panel},
    n_panels::Int,
    AIC_x,
    AIC_y,
    AIC_z,
    velocity_view_x,
    velocity_view_y,
    velocity_view_z,
    va_array,
    induced_velocity_all,
    relative_velocity_array,
    y_airf_array,
    relative_velocity_crossz,
    v_acrossz_array,
    z_airf_array,
    x_airf_array,
    v_normal_array,
    v_tangential_array,
    va_magw_array,
    cl_dist,
    chord_array,
)
    mul!(velocity_view_x, AIC_x, gamma_in)
    mul!(velocity_view_y, AIC_y, gamma_in)
    mul!(velocity_view_z, AIC_z, gamma_in)

    relative_velocity_array .= va_array .+ induced_velocity_all
    @inbounds for i in 1:n_panels
        ax = relative_velocity_array[i,1]
        ay = relative_velocity_array[i,2]
        az = relative_velocity_array[i,3]
        bx = y_airf_array[i,1]
        by = y_airf_array[i,2]
        bz = y_airf_array[i,3]
        relative_velocity_crossz[i,1] = ay*bz - az*by
        relative_velocity_crossz[i,2] = az*bx - ax*bz
        relative_velocity_crossz[i,3] = ax*by - ay*bx
        ax = va_array[i,1]
        ay = va_array[i,2]
        az = va_array[i,3]
        v_acrossz_array[i,1] = ay*bz - az*by
        v_acrossz_array[i,2] = az*bx - ax*bz
        v_acrossz_array[i,3] = ax*by - ay*bx
    end

    @inbounds for i in 1:n_panels
        v_normal_array[i] =
            z_airf_array[i,1]*relative_velocity_array[i,1] +
            z_airf_array[i,2]*relative_velocity_array[i,2] +
            z_airf_array[i,3]*relative_velocity_array[i,3]
        v_tangential_array[i] =
            x_airf_array[i,1]*relative_velocity_array[i,1] +
            x_airf_array[i,2]*relative_velocity_array[i,2] +
            x_airf_array[i,3]*relative_velocity_array[i,3]
    end
    solver.lr.alpha_dist .= atan.(v_normal_array, v_tangential_array)

    @inbounds for i in 1:n_panels
        solver.lr.v_a_dist[i] = sqrt(
            relative_velocity_crossz[i,1]^2 +
            relative_velocity_crossz[i,2]^2 +
            relative_velocity_crossz[i,3]^2)
        va_magw_array[i] = sqrt(
            v_acrossz_array[i,1]^2 +
            v_acrossz_array[i,2]^2 +
            v_acrossz_array[i,3]^2)
    end

    for (i, (panel, alpha)) in enumerate(zip(panels, solver.lr.alpha_dist))
        cl_dist[i] = calculate_cl(panel, alpha)
    end
    gamma_out .= 0.5 .* solver.lr.v_a_dist.^2 ./ va_magw_array .* cl_dist .* chord_array
    return nothing
end

@inline function lu_solve_inplace!(
    A::AbstractMatrix{Float64},
    b::AbstractVector{Float64},
    n::Int,
)
    @inbounds for k in 1:n
        amax = abs(A[k, k])
        pivot_row = k
        for i in (k + 1):n
            v = abs(A[i, k])
            if v > amax
                amax = v
                pivot_row = i
            end
        end
        amax == 0.0 && return false
        if pivot_row != k
            for j in 1:n
                tmp = A[k, j]
                A[k, j] = A[pivot_row, j]
                A[pivot_row, j] = tmp
            end
            tmp = b[k]
            b[k] = b[pivot_row]
            b[pivot_row] = tmp
        end
        inv_pivot = 1.0 / A[k, k]
        for i in (k + 1):n
            factor = A[i, k] * inv_pivot
            A[i, k] = factor
            for j in (k + 1):n
                A[i, j] -= factor * A[k, j]
            end
            b[i] -= factor * b[k]
        end
    end
    @inbounds for k in n:-1:1
        s = b[k]
        for j in (k + 1):n
            s -= A[k, j] * b[j]
        end
        b[k] = s / A[k, k]
    end
    return true
end

"""
    gamma_loop!(solver::Solver, AIC_x::Matrix{Float64},
              AIC_y::Matrix{Float64}, AIC_z::Matrix{Float64},
              panels::Vector{Panel}, relaxation_factor::Float64; log=true)

Main iteration loop for calculating circulation distribution.
"""
function gamma_loop!(
    solver::Solver,
    body_aero::BodyAerodynamics,
    panels::Vector{Panel},
    relaxation_factor::Float64;
    log::Bool = true
)
    va_array = solver.sol._va_dist
    chord_array = solver.sol._chord_dist
    x_airf_array = solver.sol._x_airf_dist
    y_airf_array = solver.sol._y_airf_dist
    z_airf_array = solver.sol._z_airf_dist
    solver.lr.converged   = false
    n_panels    = length(body_aero.panels)
    solver.lr.alpha_dist .= body_aero.alpha_dist
    solver.lr.v_a_dist   .= body_aero.v_a_dist
    
    va_magw_array            = solver.cache[1][solver.lr.v_a_dist]
    gamma                    = solver.cache[2][solver.lr.gamma_new]
    abs_gamma_new            = solver.cache[3][solver.lr.gamma_new]
    induced_velocity_all     = solver.cache[4][va_array]
    relative_velocity_array  = solver.cache[5][va_array]
    relative_velocity_crossz = solver.cache[6][va_array]
    v_acrossz_array          = solver.cache[7][va_array]
    cl_dist                 = solver.cache[8][solver.lr.gamma_new]
    damp                     = solver.cache[9][solver.lr.gamma_new]
    v_normal_array           = solver.cache[10][solver.lr.gamma_new]
    v_tangential_array       = solver.cache[11][solver.lr.gamma_new]

    AIC_x = @view body_aero.AIC[1, :, :]
    AIC_y = @view body_aero.AIC[2, :, :]
    AIC_z = @view body_aero.AIC[3, :, :]

    velocity_view_x = @view induced_velocity_all[:, 1]
    velocity_view_y = @view induced_velocity_all[:, 2]
    velocity_view_z = @view induced_velocity_all[:, 3]

    if solver.solver_type == NONLIN
        gamma_iter = solver.lr.gamma_new
        residual = solver.nonlin_residual
        residual_perturbed = solver.nonlin_residual_perturbed
        gamma_perturbed = solver.nonlin_gamma_perturbed
        jac = solver.nonlin_jac
        relstep = sqrt(eps(Float64))
        abstep = sqrt(eps(Float64))

        update_gamma_candidate!(
            residual, gamma_iter, solver, panels, n_panels,
            AIC_x, AIC_y, AIC_z,
            velocity_view_x, velocity_view_y, velocity_view_z,
            va_array, induced_velocity_all, relative_velocity_array,
            y_airf_array, relative_velocity_crossz, v_acrossz_array,
            z_airf_array, x_airf_array,
            v_normal_array, v_tangential_array,
            va_magw_array, cl_dist, chord_array,
        )
        @inbounds for i in 1:n_panels
            residual[i] -= gamma_iter[i]
        end

        solver.lr.converged = false
        for iter in 1:solver.max_iterations
            @inbounds for j in 1:n_panels
                for i in 1:n_panels
                    gamma_perturbed[i] = gamma_iter[i]
                end
                step = max(abstep, relstep * abs(gamma_iter[j]))
                gamma_perturbed[j] += step
                update_gamma_candidate!(
                    residual_perturbed, gamma_perturbed, solver, panels, n_panels,
                    AIC_x, AIC_y, AIC_z,
                    velocity_view_x, velocity_view_y, velocity_view_z,
                    va_array, induced_velocity_all, relative_velocity_array,
                    y_airf_array, relative_velocity_crossz, v_acrossz_array,
                    z_airf_array, x_airf_array,
                    v_normal_array, v_tangential_array,
                    va_magw_array, cl_dist, chord_array,
                )
                inv_step = 1.0 / step
                for i in 1:n_panels
                    residual_perturbed[i] -= gamma_perturbed[i]
                    jac[i, j] = (residual_perturbed[i] - residual[i]) * inv_step
                end
            end

            lu_solve_inplace!(jac, residual, n_panels) || break

            @inbounds for i in 1:n_panels
                gamma_iter[i] -= residual[i]
            end

            update_gamma_candidate!(
                residual, gamma_iter, solver, panels, n_panels,
                AIC_x, AIC_y, AIC_z,
                velocity_view_x, velocity_view_y, velocity_view_z,
                va_array, induced_velocity_all, relative_velocity_array,
                y_airf_array, relative_velocity_crossz, v_acrossz_array,
                z_airf_array, x_airf_array,
                v_normal_array, v_tangential_array,
                va_magw_array, cl_dist, chord_array,
            )
            max_residual = 0.0
            @inbounds for i in 1:n_panels
                residual[i] -= gamma_iter[i]
                a = abs(residual[i])
                a > max_residual && (max_residual = a)
            end
            if max_residual < solver.atol
                solver.lr.converged = true
                break
            end
        end

        gamma .= gamma_iter
        return nothing
    end

    if solver.solver_type == LOOP
        function f_loop!(gamma_new, gamma, damp)
            gamma .= gamma_new
            update_gamma_candidate!(
                gamma_new,
                gamma,
                solver,
                panels,
                n_panels,
                AIC_x,
                AIC_y,
                AIC_z,
                velocity_view_x,
                velocity_view_y,
                velocity_view_z,
                va_array,
                induced_velocity_all,
                relative_velocity_array,
                y_airf_array,
                relative_velocity_crossz,
                v_acrossz_array,
                z_airf_array,
                x_airf_array,
                v_normal_array,
                v_tangential_array,
                va_magw_array,
                cl_dist,
                chord_array,
            )
            # Update gamma with relaxation and damping
            @. gamma_new = (1 - relaxation_factor) * gamma + 
                    relaxation_factor * gamma_new + damp
            
            # Apply damping if needed
            if solver.is_with_artificial_damping
                smooth_circulation!(damp, gamma, 0.1, 0.5)
                @debug "damp: $damp"
            else
                damp .= 0.0
            end
            return nothing
        end
        iters = 0
        for i in 1:solver.max_iterations
            iters += 1

            f_loop!(solver.lr.gamma_new, gamma, damp)

            # Check convergence
            abs_gamma_new .= abs.(solver.lr.gamma_new)
            reference_error = maximum(abs_gamma_new)
            reference_error = max(reference_error, solver.tol_reference_error)
            abs_gamma_new .= abs.(solver.lr.gamma_new .- gamma)
            error = maximum(abs_gamma_new)
            normalized_error = error / reference_error

            @debug "Iteration: $i, normalized_error: $normalized_error"

            if normalized_error < solver.rtol
                solver.lr.converged = true
                break
            end
        end
    
        if log && solver.lr.converged
            @info "Converged after $iters iterations"
        elseif log
            @warn "NO convergence after $(solver.max_iterations) iterations"
        end
        return nothing
    end
end

"""
    smooth_circulation!(damp, circulation, 
                      smoothness_factor::Float64, 
                      damping_factor::Float64)

Smooth circulation distribution if needed.

Returns:
- Tuple of smoothed circulation and boolean indicating if smoothing was applied
"""
function smooth_circulation!(
    damp,
    circulation,
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

    damp .= smoothed .- circulation
    return true
end

"""
    linearize(solver::Solver, body_aero::BodyAerodynamics, y::Vector{T};
        theta_idxs=1:4, delta_idxs=nothing, va_idxs=nothing, omega_idxs=nothing,
        aero_coeffs=false, kwargs...) where T

Compute Jacobian matrix of aerodynamic outputs with respect to control and kinematic inputs using
finite differences. Used for control system design and linear stability analysis.

The function uses automatic differentiation with finite differences to compute ∂outputs/∂inputs.
Deformations are applied to the wing's unrefined sections (the original sections before mesh
refinement), with each control angle affecting one unrefined section.

# Arguments
- `solver::Solver`: Solver instance (must be configured for the wing)
- `body_aero::BodyAerodynamics`: Body aerodynamics with exactly one wing
- `y::Vector{T}`: Input vector at operating point containing control angles and/or kinematic states

# Keyword Arguments
- `theta_idxs`: Indices in `y` for twist angles (one per unrefined section, default: 1:4)
- `delta_idxs`: Indices in `y` for trailing edge deflections (one per unrefined section, default: nothing)
- `va_idxs`: Indices in `y` for apparent wind velocity components `(vx, vy, vz)` (default: nothing)
- `omega_idxs`: Indices in `y` for angular velocity components `(ωx, ωy, ωz)` (default: nothing)
- `aero_coeffs::Bool`: Return force/moment coefficients instead of dimensional values (default: false)
- `kwargs...`: Additional arguments passed to `solve!`

# Index Validation
The lengths of `theta_idxs` and `delta_idxs` (if provided) must match `wing.n_unrefined_sections`.
Unrefined sections are the original wing sections before mesh refinement for aerodynamic analysis.

# Caching
The function caches previous deformation angles to avoid redundant `unrefined_deform!` calls during
Jacobian computation. When the same angles are encountered, geometry deformation is skipped.

# Returns
- `jac::Matrix{Float64}`: Jacobian matrix (m×n) where m = 6 + n_unrefined_sections, n = length(y)
- `results::Vector{Float64}`: Output vector at operating point
    - If `aero_coeffs=false`: `(Fx, Fy, Fz, Mx, My, Mz, moment_unrefined_dist...)`
    - If `aero_coeffs=true`: `(CFx, CFy, CFz, CMx, CMy, CMz, cm_unrefined_dist...)`

# Example
```julia
# Create deformable wing with 4 unrefined sections
wing = ObjWing("kite.obj", "airfoil.dat"; n_unrefined_sections=4)
body_aero = BodyAerodynamics([wing], va=[15.0, 0, 0])
solver = Solver(body_aero)

# Operating point: 4 twist angles + velocity + angular rates
y_op = [zeros(4);        # theta angles [rad]
        [15.0, 0.0, 0.0]; # va [m/s]
        zeros(3)]         # omega [rad/s]

# Compute Jacobian
jac, results = linearize(solver, body_aero, y_op;
    theta_idxs=1:4,
    va_idxs=5:7,
    omega_idxs=8:10,
    aero_coeffs=true
)

# jac is (10×10): [6 force/moment coeffs + 4 unrefined moment coeffs] × [4 theta + 3 va + 3 omega]
```
"""
function linearize(solver::Solver, body_aero::BodyAerodynamics, y::Vector{T};
        theta_idxs=1:4,
        delta_idxs=nothing,
        va_idxs=nothing,
        omega_idxs=nothing,
        aero_coeffs=false,
        kwargs...) where T

    !(length(body_aero.wings) == 1) && throw(ArgumentError("Linearization only works for a body_aero with one wing"))
    wing = body_aero.wings[1]

    # Validate that theta_idxs and delta_idxs match the number of unrefined sections
    if !isnothing(theta_idxs) && wing.n_unrefined_sections > 0
        length(theta_idxs) != wing.n_unrefined_sections && throw(ArgumentError(
            "Length of theta_idxs ($(length(theta_idxs))) must match number of unrefined sections ($(wing.n_unrefined_sections))"))
    end
    if !isnothing(delta_idxs) && wing.n_unrefined_sections > 0
        length(delta_idxs) != wing.n_unrefined_sections && throw(ArgumentError(
            "Length of delta_idxs ($(length(delta_idxs))) must match number of unrefined sections ($(wing.n_unrefined_sections))"))
    end
    if wing.n_unrefined_sections == 0 && (!isnothing(theta_idxs) || !isnothing(delta_idxs))
        throw(ArgumentError("Cannot use theta_idxs or delta_idxs when wing has no unrefined sections"))
    end

    init_va = body_aero.cache[1][body_aero._va]
    init_va .= body_aero._va
    init_omega = copy(body_aero.omega)
    last_theta_ref = Ref{Vector{T}}(Vector{T}(undef, 0))
    if !isnothing(theta_idxs)
        @views last_theta_ref[] = body_aero.cache[2][y[theta_idxs]]
        last_theta_ref[] .= NaN
    end
    last_delta_ref = Ref{Vector{T}}(Vector{T}(undef, 0))
    if !isnothing(delta_idxs)
        @views last_delta_ref[] = body_aero.cache[3][y[delta_idxs]]
        last_delta_ref[] .= NaN
    end

    # Function to compute forces for given control inputs
    function calc_results!(results, y)
        @views theta_angles = isnothing(theta_idxs) ? nothing : y[theta_idxs]
        @views delta_angles = isnothing(delta_idxs) ? nothing : y[delta_idxs]
        last_theta = last_theta_ref[]
        last_delta = last_delta_ref[]

        if !isnothing(theta_angles) && isnothing(delta_angles)
            if !all(theta_angles .== last_theta)
                VortexStepMethod.unrefined_deform!(wing, theta_angles, nothing; smooth=false)
                VortexStepMethod.reinit!(body_aero; init_aero=false)
                last_theta .= theta_angles
            end
        elseif !isnothing(theta_angles) && !isnothing(delta_angles)
            if !all(delta_angles .== last_delta) || !all(theta_angles .== last_theta)
                VortexStepMethod.unrefined_deform!(wing, theta_angles, delta_angles; smooth=false)
                VortexStepMethod.reinit!(body_aero; init_aero=false)
                last_theta .= theta_angles
                last_delta .= delta_angles
            end
        elseif isnothing(theta_angles) && !isnothing(delta_angles)
            if !all(delta_angles .== last_delta)
                VortexStepMethod.unrefined_deform!(wing, nothing, delta_angles; smooth=false)
                VortexStepMethod.reinit!(body_aero; init_aero=false)
                last_delta .= delta_angles
            end
        end

        va = isnothing(va_idxs) ? init_va : y[va_idxs]
        om = isnothing(omega_idxs) ? init_omega :
            y[omega_idxs]
        set_va!(body_aero, va, om)

        solve!(solver, body_aero; kwargs...)
        if !aero_coeffs
            results[1:3] .= solver.sol.force
            results[4:6] .= solver.sol.moment
            results[7:end] .= solver.sol.moment_unrefined_dist
        else
            results[1:3] .= solver.sol.force_coeffs
            results[4:6] .= solver.sol.moment_coeffs
            results[7:end] .= solver.sol.cm_unrefined_dist
        end
        return nothing
    end

    results = zeros(3+3+length(solver.sol.moment_unrefined_dist))
    jac = zeros(length(results), length(y))
    backend = AutoFiniteDiff(absstep=1e2solver.atol, relstep=1e2solver.rtol)
    prep = prepare_jacobian(calc_results!, results, backend, y)
    jacobian!(calc_results!, results, jac, prep, backend, y)

    calc_results!(results, y)
    return jac, results
end

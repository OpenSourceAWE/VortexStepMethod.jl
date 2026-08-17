
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
- `cm_unrefined_dist`::MVector{U, Float64}: Averaged airfoil moment coefficients for unrefined sections [-]
- `moment_coeff_unrefined_dist`::MVector{U, Float64}: Summed `moment_frac`-referenced pitching-moment coefficient per unrefined section [-]
- `alpha_unrefined_dist`::MVector{U, Float64}: Averaged angles of attack for unrefined sections [rad]
- `solver_status`::SolverStatus: enum, see [SolverStatus](@ref)
"""
@with_kw mutable struct VSMSolution{P, U, T}
    ### private vectors of solve_base!
    _x_airf_dist::Matrix{T} = zeros(T, P, 3)
    _y_airf_dist::Matrix{T} = zeros(T, P, 3)
    _z_airf_dist::Matrix{T} = zeros(T, P, 3)
    _va_dist::Matrix{T} = zeros(T, P, 3)
    _chord_dist::Vector{T} = zeros(T, P)
    ### end of private vectors
    width_dist::Vector{T} = zeros(T, P)
    panel_area_dist::Vector{T} = zeros(T, P)
    alpha_dist::Vector{T} = zeros(T, P)
    alpha_geometric_dist::Vector{T} = zeros(T, P)
    cl_dist::Vector{T} = zeros(T, P)
    cd_dist::Vector{T} = zeros(T, P)
    cm_dist::Vector{T} = zeros(T, P)
    lift_dist::Vector{T} = zeros(T, P)
    drag_dist::Vector{T} = zeros(T, P)
    panel_moment_dist::Vector{T} = zeros(T, P)
    f_body_3D::Matrix{T} = zeros(T, 3, P)
    m_body_3D::Matrix{T} = zeros(T, 3, P)
    gamma_distribution::Union{Nothing, Vector{T}} = nothing
    force::MVector{3, T} = zeros(MVector{3, T})
    moment::MVector{3, T} = zeros(MVector{3, T})
    force_coeffs::MVector{3, T} = zeros(MVector{3, T})
    moment_coeffs::MVector{3, T} = zeros(MVector{3, T})
    center_of_pressure::Union{Nothing, MVector{3, T}} = nothing
    panel_cp_locations::Vector{MVector{3, T}} = MVector{3, T}[]
    moment_dist::MVector{P, T} = zeros(MVector{P, T})
    moment_coeff_dist::MVector{P, T} = zeros(MVector{P, T})
    moment_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    cl_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    cd_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    cm_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    moment_coeff_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    alpha_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    x_airf_unrefined_dist::Vector{MVector{3, T}} = [zeros(MVector{3, T}) for _ in 1:U]
    y_airf_unrefined_dist::Vector{MVector{3, T}} = [zeros(MVector{3, T}) for _ in 1:U]
    z_airf_unrefined_dist::Vector{MVector{3, T}} = [zeros(MVector{3, T}) for _ in 1:U]
    va_unrefined_dist::Vector{MVector{3, T}} = [zeros(MVector{3, T}) for _ in 1:U]
    chord_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    width_unrefined_dist::MVector{U, T} = zeros(MVector{U, T})
    unrefined_count_dist::Vector{Int} = zeros(Int, U)
    solver_status::SolverStatus = FAILURE
end

# Output of the function gamma_loop!
@with_kw mutable struct LoopResult{P, T}
    converged::Bool                  = false
    gamma_new::MVector{P, T}         = zeros(MVector{P, T})
    alpha_dist::MVector{P, T}        = zeros(MVector{P, T})
    v_a_dist::MVector{P, T}          = zeros(MVector{P, T})
end

@with_kw struct BaseResult{P, T}
    va_norm_dist::MVector{P, T} = zeros(MVector{P, T})
    va_unit_dist::Matrix{T} = zeros(T, P, 3)
end

@inline function check_reference_point(reference_point, ::Type{T}=Float64) where {T}
    msg = "reference_point must be a list/array with 3 numbers."
    reference_point isa AbstractVector || throw(ArgumentError(msg))
    length(reference_point) == 3 || throw(ArgumentError(msg))
    all(x -> x isa Number, reference_point) || throw(ArgumentError(msg))
    try
        return MVector{3, T}(reference_point[1], reference_point[2], reference_point[3])
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

## Artificial viscosity settings
- `is_with_artificial_viscosity`::Bool = false: Enable the Li/Gaunaa spanwise artificial
    viscosity (TORQUE 2026) for post-stall stabilization in the LOOP solver
- `artificial_viscosity_factor`::Float64 = 0.035: Coefficient k in the viscosity scaling
    (the conservative envelope from the paper)

## Additional settings
- `type_initial_gamma_distribution`::InitialGammaDistribution = ZEROS: see: [InitialGammaDistribution](@ref)
- `use_gamma_prev`::Bool = true: reuse provided previous gamma as initial guess when available
- `core_radius_fraction`::Float64 = 0.05: Bound vortex core cut-off, as a fraction of the
    filament length, following Damiani et al. (2019)
- mu::Float64 = 1.81e-5: Dynamic viscosity [N·s/m²]
- `is_only_f_and_gamma_output`::Bool = false: Whether to only output f and gamma
- `reference_point`::MVec3 = [0.0, 0.0, 0.0]: Moment reference point in body frame

## Solution
sol::VSMSolution = VSMSolution(): The result of calling [solve!](@ref) 
"""
@with_kw mutable struct Solver{P, U, T}
    # General settings
    solver_type::SolverType = LOOP
    aerodynamic_model_type::Model = VSM
    density::T = T(1.225)
    max_iterations::Int64 = 1500
    rtol::T = T(1e-5)
    tol_reference_error::T = T(0.001)
    relaxation_factor::T = T(0.03)

    # Nonlin solver fields
    atol::T = T(1e-5)
    nonlin_jac::Matrix{T} = zeros(T, P, P)
    nonlin_residual::MVector{P, T} = zeros(MVector{P, T})
    nonlin_residual_perturbed::MVector{P, T} = zeros(MVector{P, T})
    nonlin_gamma_perturbed::MVector{P, T} = zeros(MVector{P, T})
    nonlin_ipiv::Vector{LinearAlgebra.BlasInt} = zeros(LinearAlgebra.BlasInt, P)

    # Damping settings
    is_with_artificial_damping::Bool = false
    artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}} =(k2=0.1, k4=0.0)

    # Li/Gaunaa spanwise artificial viscosity (TORQUE 2026) settings
    is_with_artificial_viscosity::Bool = false
    artificial_viscosity_factor::T = T(0.035)

    # Additional settings
    type_initial_gamma_distribution::InitialGammaDistribution = ZEROS
    use_gamma_prev::Bool = true
    core_radius_fraction::T = T(0.05)
    mu::T = T(1.81e-5)
    is_only_f_and_gamma_output::Bool = false
    correct_aoa::Bool = false
    reference_point::MVector{3, T} = zeros(MVector{3, T})

    # Intermediate results
    lr::LoopResult{P, T} = LoopResult{P, T}()
    br::BaseResult{P, T} = BaseResult{P, T}()
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:11]
    cache_base::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}  = [LazyBufferCache()]
    cache_lin::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:4]

    # Solution
    sol::VSMSolution{P, U, T} = VSMSolution{P, U, T}()
end

function Solver(body_aero::BodyAerodynamics{P, W, T}; reference_point=[0.0, 0.0, 0.0], kwargs...) where {P, W, T}
    U = sum([wing.n_unrefined_sections for wing in body_aero.wings])
    reference_point_checked = check_reference_point(reference_point, T)
    return Solver{P, U, T}(; reference_point=reference_point_checked, kwargs...)
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
        is_with_artificial_viscosity=ss.is_with_artificial_viscosity,
        artificial_viscosity_factor=ss.artificial_viscosity_factor,
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
function solve!(solver::Solver{P, U, T}, body_aero::BodyAerodynamics, gamma_distribution=solver.sol.gamma_distribution;
        log=false, reference_point=solver.reference_point, moment_frac=0.1) where {P, U, T}

    # calculate intermediate result
    solve_base!(solver, body_aero, gamma_distribution; log)
    gamma_new = solver.lr.gamma_new
    if !isnothing(solver.sol.gamma_distribution)
        solver.sol.gamma_distribution .= gamma_new
    else
        solver.sol.gamma_distribution = gamma_new
    end

    return calc_forces!(solver, body_aero; reference_point, moment_frac)
end

"""
    calc_forces!(solver::Solver, body_aero::BodyAerodynamics;
                 reference_point=solver.reference_point, moment_frac=0.1)

Assemble aerodynamic forces and moments from the circulation already converged
by [`solve_base!`](@ref) and stored in `solver`. Split out of [`solve!`](@ref).
"""
function calc_forces!(solver::Solver{P, U, T}, body_aero::BodyAerodynamics;
        reference_point=solver.reference_point, moment_frac=0.1) where {P, U, T}
    gamma_new = solver.lr.gamma_new

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
                alpha_geometric_dist[i] = atan(
                    -v_normal, -v_tangential)
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
    area_all_panels = zero(T)
    panel_areas = solver.sol.panel_area_dist

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
        moment_coeff_unrefined_dist = solver.sol.moment_coeff_unrefined_dist
        alpha_unrefined_dist = solver.sol.alpha_unrefined_dist
        x_airf_unrefined_dist = solver.sol.x_airf_unrefined_dist
        y_airf_unrefined_dist = solver.sol.y_airf_unrefined_dist
        z_airf_unrefined_dist = solver.sol.z_airf_unrefined_dist
        va_unrefined_dist = solver.sol.va_unrefined_dist
        chord_unrefined_dist = solver.sol.chord_unrefined_dist
        width_unrefined_dist = solver.sol.width_unrefined_dist
        unrefined_count_dist = solver.sol.unrefined_count_dist

        # Zero all unrefined arrays
        moment_unrefined_dist .= 0.0
        cl_unrefined_dist .= 0.0
        cd_unrefined_dist .= 0.0
        cm_unrefined_dist .= 0.0
        moment_coeff_unrefined_dist .= 0.0
        alpha_unrefined_dist .= 0.0
        for i in eachindex(x_airf_unrefined_dist)
            x_airf_unrefined_dist[i] .= 0.0
            y_airf_unrefined_dist[i] .= 0.0
            z_airf_unrefined_dist[i] .= 0.0
            va_unrefined_dist[i] .= 0.0
        end
        chord_unrefined_dist .= 0.0
        width_unrefined_dist .= 0.0
        fill!(unrefined_count_dist, 0)

        panel_idx = 1
        unrefined_idx = 1
        for wing in body_aero.wings
            if wing.n_unrefined_sections > 0
                for local_panel_idx in 1:wing.n_panels
                    panel = body_aero.panels[panel_idx]
                    original_section_idx = wing.refined_panel_mapping[local_panel_idx]
                    target_unrefined_idx = unrefined_idx + original_section_idx - 1

                    # Accumulate coefficients and moments
                    moment_unrefined_dist[target_unrefined_idx] += moment_dist[panel_idx]
                    cl_unrefined_dist[target_unrefined_idx] += solver.sol.cl_dist[panel_idx]
                    cd_unrefined_dist[target_unrefined_idx] += solver.sol.cd_dist[panel_idx]
                    cm_unrefined_dist[target_unrefined_idx] += solver.sol.cm_dist[panel_idx]
                    moment_coeff_unrefined_dist[target_unrefined_idx] += solver.sol.moment_coeff_dist[panel_idx]
                    alpha_unrefined_dist[target_unrefined_idx] += solver.sol.alpha_dist[panel_idx]

                    # Accumulate geometry
                    x_airf_unrefined_dist[target_unrefined_idx] .+= panel.x_airf
                    y_airf_unrefined_dist[target_unrefined_idx] .+= panel.y_airf
                    z_airf_unrefined_dist[target_unrefined_idx] .+= panel.z_airf
                    va_unrefined_dist[target_unrefined_idx] .+= panel.va
                    chord_unrefined_dist[target_unrefined_idx] += panel.chord
                    width_unrefined_dist[target_unrefined_idx] += panel.width

                    unrefined_count_dist[target_unrefined_idx] += 1
                    panel_idx += 1
                end

                # Average coefficients and geometry. width and
                # moment_coeff_unrefined_dist stay summed (extensive).
                for i in 1:wing.n_unrefined_sections
                    target_unrefined_idx = unrefined_idx + i - 1
                    if unrefined_count_dist[target_unrefined_idx] > 0
                        count = unrefined_count_dist[target_unrefined_idx]
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

"""
    solve_base!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing;
                log=false)

Converge the circulation distribution and leave it in `solver.lr.gamma_new`, without
turning it into forces. Fills the solver's panel arrays, builds the AIC matrices,
starts from `gamma_distribution` (or an elliptical/zero guess when it is `nothing` or
`solver.use_gamma_prev` is false) and iterates; a `LOOP` solver that fails to converge
retries once with half the relaxation factor. The circulation half of [`solve!`](@ref),
paired with [`calc_forces!`](@ref). Returns `nothing`.
"""
function solve_base!(solver::Solver{P, U, T}, body_aero::BodyAerodynamics, gamma_distribution=nothing;
               log=false) where {P, U, T}
    
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

@inline smooth_sqrt(x) = sqrt(x + 1e-30)

@inline function update_gamma_candidate!(
    gamma_out,
    gamma_in,
    solver::Solver,
    panels::AbstractVector{<:Panel},
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
        solver.lr.v_a_dist[i] = smooth_sqrt(
            relative_velocity_crossz[i,1]^2 +
            relative_velocity_crossz[i,2]^2 +
            relative_velocity_crossz[i,3]^2)
        va_magw_array[i] = smooth_sqrt(
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

"""
    build_spanwise_laplacian!(laplacian, n_panels)

Fill the `n_panels × n_panels` matrix `laplacian` with the discrete spanwise
Laplacian used by the Li/Gaunaa post-stall artificial-viscosity regularization.
Interior rows use the three-point stencil `gamma[i-1] - 2 gamma[i] + gamma[i+1]`;
the tip rows use the second-order closures of Li, Gaunaa, Pirrung & Lønbæk
(TORQUE 2026, Eq. 15) that enforce `gamma -> 0` at the wing tips. Panels are
assumed ordered consecutively along the span and approximately uniformly spaced.
Returns the matrix unchanged (all zeros) when `n_panels < 3`.
"""
function build_spanwise_laplacian!(laplacian, n_panels)
    laplacian .= 0
    n_panels < 3 && return laplacian
    @inbounds for i in 2:n_panels-1
        laplacian[i, i-1] = 1.0
        laplacian[i, i]   = -2.0
        laplacian[i, i+1] = 1.0
    end
    laplacian[1, 1] = -4.0
    laplacian[1, 2] = 4.0 / 3.0
    laplacian[n_panels, n_panels]   = -4.0
    laplacian[n_panels, n_panels-1] = 4.0 / 3.0
    return laplacian
end

"""
    local_lift_slope!(slopes, panels, alpha_dist, delta=deg2rad(0.5))

Per-panel local lift-curve slope `dCl/dalpha`, evaluated from each panel's own
2-D polar at its current effective angle of attack by central differences. The
slope turns negative in post-stall, which is what activates the
artificial-viscosity regularization in [`gamma_loop!`](@ref).
"""
function local_lift_slope!(slopes, panels, alpha_dist, delta=deg2rad(0.5))
    inv_2delta = 1.0 / (2delta)
    @inbounds for i in eachindex(panels)
        cl_plus  = calculate_cl(panels[i], alpha_dist[i] + delta)
        cl_minus = calculate_cl(panels[i], alpha_dist[i] - delta)
        slopes[i] = (cl_plus - cl_minus) * inv_2delta
    end
    return slopes
end

"""
    apply_artificial_viscosity!(gamma, panels, alpha_dist, laplacian, viscosity_matrix,
                                lift_slope, mu_array, gamma_target, planform_area, factor)

Apply one implicit Li/Gaunaa artificial-viscosity step to `gamma` in place and return
`true` when it fired. The per-panel viscosity is
`mu_i = max(0, -factor * planform_area * Cl'_i / width_i^2)`, with the local lift slope
`Cl'_i` from [`local_lift_slope!`](@ref) evaluated at `alpha_dist`. When any panel is
post-stall (`mu_i > 0`), `gamma` is replaced by the solution of
`(I - diag(mu) L) gamma = gamma`, with `L` the spanwise Laplacian in `laplacian`
(see [`build_spanwise_laplacian!`](@ref)); otherwise `gamma` is left unchanged. The
remaining arguments are preallocated work buffers reused across iterations.
"""
function apply_artificial_viscosity!(gamma, panels, alpha_dist, laplacian, viscosity_matrix,
        lift_slope, mu_array, gamma_target, planform_area, factor)
    n_panels = length(panels)
    local_lift_slope!(lift_slope, panels, alpha_dist)
    any_stalled = false
    @inbounds for i in 1:n_panels
        m = -factor * planform_area * lift_slope[i] / panels[i].width^2
        mu_array[i] = max(zero(eltype(mu_array)), m)
        mu_array[i] > 0 && (any_stalled = true)
    end
    any_stalled || return false
    gamma_target .= gamma
    one_t, zero_t = one(eltype(viscosity_matrix)), zero(eltype(viscosity_matrix))
    @inbounds for col in 1:n_panels, row in 1:n_panels
        viscosity_matrix[row, col] =
            (row == col ? one_t : zero_t) - mu_array[row] * laplacian[row, col]
    end
    ldiv!(gamma, lu!(viscosity_matrix), gamma_target)
    return true
end

"""
    gamma_loop!(solver::Solver, AIC_x::Matrix{Float64},
              AIC_y::Matrix{Float64}, AIC_z::Matrix{Float64},
              panels::AbstractVector{<:Panel}, relaxation_factor::Float64; log=true)

Main iteration loop for calculating circulation distribution.

When `solver.is_with_artificial_viscosity` is set, the LOOP solver replaces the
explicit target `F(gamma)` with the implicit Li/Gaunaa solution
`(I - diag(mu) L) gamma = F(gamma)` before relaxation, stabilizing post-stall
distributions. The Python reference gates this on precomputed per-panel stall
angles; here the polar is an interpolation object, so the expensive linear solve
is gated on `any(mu > 0)` instead, which is behaviourally identical.
"""
function gamma_loop!(
    solver::Solver{P, U, T},
    body_aero::BodyAerodynamics,
    panels::AbstractVector{<:Panel},
    relaxation_factor;
    log::Bool = true
) where {P, U, T}
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
    damp                    .= zero(T)
    v_normal_array           = solver.cache[10][solver.lr.gamma_new]
    v_tangential_array       = solver.cache[11][solver.lr.gamma_new]

    AIC_x = @view body_aero.AIC[:, :, 1]
    AIC_y = @view body_aero.AIC[:, :, 2]
    AIC_z = @view body_aero.AIC[:, :, 3]

    velocity_view_x = @view induced_velocity_all[:, 1]
    velocity_view_y = @view induced_velocity_all[:, 2]
    velocity_view_z = @view induced_velocity_all[:, 3]

    if solver.solver_type == NONLIN
        T <: ForwardDiff.Dual && error(
            "NONLIN solver does not support ForwardDiff in this release. " *
            "Use `solver_type=LOOP` for differentiable solves.")
        gamma_iter = solver.lr.gamma_new
        residual = solver.nonlin_residual
        residual_perturbed = solver.nonlin_residual_perturbed
        gamma_perturbed = solver.nonlin_gamma_perturbed
        jac = solver.nonlin_jac
        ipiv = solver.nonlin_ipiv
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

            _, _, info = LinearAlgebra.LAPACK.getrf!(jac, ipiv; check=false)
            info == 0 || break
            LinearAlgebra.LAPACK.getrs!('N', jac, ipiv, residual)

            max_step = 0.0
            ref = solver.tol_reference_error
            @inbounds for i in 1:n_panels
                s = abs(residual[i])
                s > max_step && (max_step = s)
                gamma_iter[i] -= residual[i]
                g = abs(gamma_iter[i])
                g > ref && (ref = g)
            end
            if max_step < solver.atol + solver.rtol * ref
                solver.lr.converged = true
                break
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
            @inbounds for i in 1:n_panels
                residual[i] -= gamma_iter[i]
            end
        end

        gamma .= gamma_iter
        return nothing
    end

    if solver.solver_type == LOOP
        # Work buffers allocated once: the geometry is frozen during iteration.
        use_viscosity = solver.is_with_artificial_viscosity
        laplacian        = use_viscosity ? zeros(T, n_panels, n_panels) : zeros(T, 0, 0)
        viscosity_matrix = use_viscosity ? zeros(T, n_panels, n_panels) : zeros(T, 0, 0)
        lift_slope       = use_viscosity ? zeros(T, n_panels) : zeros(T, 0)
        mu_array         = use_viscosity ? zeros(T, n_panels) : zeros(T, 0)
        gamma_target     = use_viscosity ? zeros(T, n_panels) : zeros(T, 0)
        planform_area    = zero(T)
        if use_viscosity
            build_spanwise_laplacian!(laplacian, n_panels)
            @inbounds for i in 1:n_panels
                planform_area += panels[i].width * chord_array[i]
            end
        end

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
            # Gate the linear solve on any(mu > 0): fires exactly in post-stall.
            if use_viscosity
                apply_artificial_viscosity!(
                    gamma_new, panels, solver.lr.alpha_dist, laplacian,
                    viscosity_matrix, lift_slope, mu_array, gamma_target,
                    planform_area, solver.artificial_viscosity_factor,
                )
            end

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

function _section_with_eltype(section::Section, ::Type{TD}) where TD
    return Section{TD}(
        MVector{3, TD}(section.LE_point),
        MVector{3, TD}(section.TE_point),
        section.aero_model,
        section.aero_data,
        section.section_aero,
    )
end

function _wing_with_eltype(wing::Wing{P, Float64}, ::Type{TD}) where {P, TD}
    return Wing{P, TD}(
        wing.n_panels,
        wing.n_unrefined_sections,
        wing.spanwise_distribution,
        PanelProperties{P, TD}(),
        MVector{3, TD}(wing.spanwise_direction),
        Section{TD}[_section_with_eltype(s, TD) for s in wing.unrefined_sections],
        Section{TD}[_section_with_eltype(s, TD) for s in wing.refined_sections],
        wing.remove_nan,
        wing.use_prior_polar,
        wing.billowing_percentage,
        TD(wing.crease_frac),
        copy(wing.refined_panel_mapping),
        copy(wing.refined_section_left_idx),
        Vector{TD}(wing.refined_section_weight),
        Section{TD}[_section_with_eltype(s, TD) for s in wing.non_deformed_sections],
        Vector{TD}(wing.theta_dist),
        Vector{TD}(wing.delta_dist),
        TD(wing.mass),
        TD(wing.gamma_tip),
        Matrix{TD}(wing.inertia_tensor),
        MVector{3, TD}(wing.T_cad_body),
        MMatrix{3, 3, TD, 9}(wing.R_cad_body),
        TD(wing.radius),
        wing.le_interp,
        wing.te_interp,
        wing.area_interp,
        [PreallocationTools.LazyBufferCache()],
    )
end

"""
    make_dual_shadow(solver, body_aero, ::Type{TD}) where TD

Build a `BodyAerodynamics`/`Solver` pair whose mutating buffers carry element type
`TD` (typically a `ForwardDiff.Dual`). Wing geometry, polar coefficients and
interpolators are reused from the Float64 originals; circulation/AIC/scratch
buffers are freshly allocated as `TD`-typed.
"""
function make_dual_shadow(solver::Solver{P, U, Float64},
                          body_aero::BodyAerodynamics{P, W, Float64},
                          ::Type{TD}) where {P, U, W, TD}
    length(body_aero.wings) == 1 || throw(ArgumentError(
        "make_dual_shadow currently supports body_aero with one wing"))
    wing_d = _wing_with_eltype(body_aero.wings[1], TD)
    body_aero_d = BodyAerodynamics([wing_d];
        va = MVector{3, TD}(body_aero._va),
        omega = MVector{3, TD}(body_aero.omega),
    )
    solver_d = Solver(body_aero_d;
        solver_type = solver.solver_type,
        aerodynamic_model_type = solver.aerodynamic_model_type,
        density = TD(solver.density),
        max_iterations = solver.max_iterations,
        rtol = TD(solver.rtol),
        tol_reference_error = TD(solver.tol_reference_error),
        relaxation_factor = TD(solver.relaxation_factor),
        atol = TD(solver.atol),
        is_with_artificial_damping = solver.is_with_artificial_damping,
        artificial_damping = solver.artificial_damping,
        is_with_artificial_viscosity = solver.is_with_artificial_viscosity,
        artificial_viscosity_factor = TD(solver.artificial_viscosity_factor),
        type_initial_gamma_distribution = solver.type_initial_gamma_distribution,
        use_gamma_prev = false,
        core_radius_fraction = TD(solver.core_radius_fraction),
        mu = TD(solver.mu),
        is_only_f_and_gamma_output = solver.is_only_f_and_gamma_output,
        correct_aoa = solver.correct_aoa,
        reference_point = MVector{3, TD}(solver.reference_point),
    )
    return body_aero_d, solver_d
end

"""
    linearize(solver, body_aero, y; theta_idxs=1:4, delta_idxs=nothing,
              va_idxs=nothing, omega_idxs=nothing, aero_coeffs=false,
              backend=AutoForwardDiff(), kwargs...)

Jacobian of aerodynamic outputs w.r.t. control and kinematic inputs at `y`. Each `*_idxs`
selects which entries of `y` map to twist angles (one per unrefined section), trailing-edge
deflections (one per unrefined section), apparent wind `(vx, vy, vz)`, and angular rate
`(ωx, ωy, ωz)` respectively.

`backend` accepts any `DifferentiationInterface` backend; `AutoForwardDiff()` (the default)
requires `solver_type=LOOP`. `fd_absstep`/`fd_relstep` are forwarded only when the backend is
`AutoFiniteDiff`.

Returns `(jac, results, converged)` where `results` is `(F, M, moment_unrefined_dist...)` —
or the corresponding coefficients when `aero_coeffs=true` — and `converged` is `false` (with a
warning) if any internal solve missed the solver's tolerances.
"""
function linearize(solver::Solver, body_aero::BodyAerodynamics, y::Vector{T};
        theta_idxs=1:4,
        delta_idxs=nothing,
        va_idxs=nothing,
        omega_idxs=nothing,
        aero_coeffs=false,
        backend = AutoForwardDiff(),
        fd_absstep::Float64=1e-3,
        fd_relstep::Float64=1e-3,
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

    n_failed = Ref(0)
    shadow_ref = Ref{Any}(nothing)

    function calc_results!(results, y_in)
        TI = eltype(y_in)
        if TI === Float64
            body_aero_c = body_aero
            solver_c = solver
            wing_c = wing
        else
            shadow = shadow_ref[]
            if shadow === nothing || eltype(shadow[1]._va) !== TI
                shadow_ref[] = make_dual_shadow(solver, body_aero, TI)
            end
            body_aero_c, solver_c = shadow_ref[]
            wing_c = body_aero_c.wings[1]
        end

        @views theta_angles = isnothing(theta_idxs) ? nothing : y_in[theta_idxs]
        @views delta_angles = isnothing(delta_idxs) ? nothing : y_in[delta_idxs]

        if !isnothing(theta_angles) || !isnothing(delta_angles)
            VortexStepMethod.unrefined_deform!(wing_c, theta_angles, delta_angles; smooth=false)
            VortexStepMethod.reinit!(body_aero_c; init_aero=false)
        end

        va = isnothing(va_idxs) ? MVector{3, TI}(body_aero_c._va) : y_in[va_idxs]
        om = isnothing(omega_idxs) ? MVector{3, TI}(body_aero_c.omega) : y_in[omega_idxs]
        set_va!(body_aero_c, va, om)

        solve!(solver_c, body_aero_c; kwargs...)
        solver_c.lr.converged || (n_failed[] += 1)
        if !aero_coeffs
            results[1:3] .= solver_c.sol.force
            results[4:6] .= solver_c.sol.moment
            results[7:end] .= solver_c.sol.moment_unrefined_dist
        else
            results[1:3] .= solver_c.sol.force_coeffs
            results[4:6] .= solver_c.sol.moment_coeffs
            results[7:end] .= solver_c.sol.moment_coeff_unrefined_dist
        end
        return nothing
    end

    n_results = 3 + 3 + length(solver.sol.moment_unrefined_dist)
    jac = zeros(n_results, length(y))
    results = zeros(n_results)
    be = backend === nothing ?
        AutoFiniteDiff(absstep=fd_absstep, relstep=fd_relstep) : backend
    prep = prepare_jacobian(calc_results!, results, be, y)
    jacobian!(calc_results!, results, jac, prep, be, y)
    calc_results!(results, y)
    converged = n_failed[] == 0
    if !converged
        @warn "linearize: $(n_failed[]) internal solve(s) did not converge; \
               jacobian columns may be unreliable"
    end
    return jac, results, converged
end

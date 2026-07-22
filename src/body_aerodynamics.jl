"""
    @with_kw mutable struct BodyAerodynamics{P,W<:AbstractWing}

Main structure for calculating aerodynamic properties of bodies. Use the constructor to initialize.

# Fields
- panels::Vector{<:Panel}: Vector of refined [Panel](@ref) structs
- wings::Vector{W}: A vector of wings of type `W <: AbstractWing`; a body can have multiple wings
- `va::MVec3` = zeros(MVec3):   A vector of the apparent wind speed, see: [MVec3](@ref)
- `omega`::MVec3 = zeros(MVec3): A vector of the turn rates around the kite body axes
- `gamma_distribution`=zeros(Float64, P): A vector of the circulation
                        of the velocity field; Length: Number of segments. [m²/s]
- `alpha_uncorrected`=zeros(Float64, P): angles of attack per panel
- `alpha_corrected`=zeros(Float64, P):   corrected angles of attack per panel
- `stall_angle_list`=zeros(Float64, P):  stall angle per panel
- `alpha_dist::MVector{P, Float64}` = zeros(Float64, P)
- `v_a_dist::MVector{P, Float64}` = zeros(Float64, P)
- `work_vectors`::NTuple{10, MVec3} = ntuple(_ -> zeros(MVec3), 10)
- `AIC::Array{Float64, 3}` = zeros(3, P, P)
- `projected_area::Float64` = 1.0: The area projected onto the xy-plane of the kite body reference frame [m²]
- `c_ref::Float64` = 1.0: Reference chord length (max panel chord) [m]
- `y::MVector{P, Float64}` = MVector{P,Float64}(zeros(P))
- `cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}` = [LazyBufferCache() for _ in 1:15]
"""
@with_kw mutable struct BodyAerodynamics{P, W<:AbstractWing, T, PN<:Panel{T}}
    panels::Vector{PN}
    wings::Vector{W}
    _va::MVector{3, T} = zeros(MVector{3, T})
    has_distributed_va::Bool = false
    omega::MVector{3, T} = zeros(MVector{3, T})
    gamma_distribution::MVector{P, T} = zeros(MVector{P, T})
    alpha_uncorrected::MVector{P, T} = zeros(MVector{P, T})
    alpha_corrected::MVector{P, T} = zeros(MVector{P, T})
    stall_angle_list::MVector{P, T} = zeros(MVector{P, T})
    alpha_dist::MVector{P, T} = zeros(MVector{P, T})
    v_a_dist::MVector{P, T} = zeros(MVector{P, T})
    work_vectors::NTuple{10, MVector{3, T}} = ntuple(_ -> zeros(MVector{3, T}), 10)
    AIC::Array{T, 3} = zeros(T, 3, P, P)
    projected_area::T = one(T)
    c_ref::T = one(T)
    y::MVector{P, T} = zeros(MVector{P, T})
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:15]
end

"""
    BodyAerodynamics(wings::Vector{T}; 
                     kite_body_origin=zeros(MVec3)) where T <: AbstractWing

Construct a [BodyAerodynamics](@ref) object for aerodynamic calculations.

This constructor handles initialization of panels, coordinate transformations, and
aerodynamic properties, returning a fully initialized structure ready for simulation.

# Arguments
- `wings::Vector{T}`: Vector of wings to analyze, where T is an AbstractWing type

# Keyword Arguments
- `kite_body_origin=zeros(MVec3)`: Origin point of kite body reference frame in CAD reference frame
- `va=[15.0, 0.0, 0.0]`: Apparent wind vector
- `omega=zeros(3)`: Turn rate in kite body frame x y and z

# Returns
- [BodyAerodynamics](@ref) object initialized with panels and wings

# Example
```julia
wing = Wing("wing.yaml"; n_panels=40); refine!(wing)
body_aero = BodyAerodynamics([wing], va=[15.0, 0.0, 0.0], omega=zeros(3))
```
"""
function BodyAerodynamics(
    wings::Vector{W};
    kite_body_origin=zeros(MVec3),
    va=[15.0, 0.0, 0.0],
    omega=zeros(MVec3)
) where {T, W <: AbstractWing{T}}
    # Validate all wings are refined
    for (i, wing) in enumerate(wings)
        if isempty(wing.refined_sections) ||
           length(wing.refined_sections) != wing.n_panels + 1
            throw(ArgumentError(
                "Wing $i has not been refined. " *
                "Call refine!(wing) before creating BodyAerodynamics.\n\n" *
                "Expected workflow:\n" *
                "  wing = Wing(...)\n" *
                "  refine!(wing)\n" *
                "  body_aero = BodyAerodynamics([wing])"
            ))
        end

        if isempty(wing.non_deformed_sections)
            @warn "Wing $i has no non_deformed_sections. " *
                  "Deformation (unrefined_deform!) will not work. " *
                  "This should have been created by refine!." maxlog=1
        end
    end

    # All panels share one concrete type from the wings' (uniform) aero model, so the
    # panel vector stays concretely typed — see panel_interp_types.
    sec0 = first(first(wings).unrefined_sections)
    CL, CD, CM, CP = panel_interp_types(sec0, first(wings).remove_nan)
    panels = Panel{T, CL, CD, CM, CP}[]
    for wing in wings
        for section in wing.unrefined_sections
            section.LE_point .-= kite_body_origin
            section.TE_point .-= kite_body_origin
        end

        # Create panels
        for _ in 1:wing.n_panels
            push!(panels, Panel{T, CL, CD, CM, CP}())
        end
    end

    body_aero = BodyAerodynamics{length(panels), W, T, eltype(panels)}(; panels, wings)
    reinit!(body_aero; va, omega)
    return body_aero
end

function Base.getproperty(obj::BodyAerodynamics, sym::Symbol)
    if sym === :va
        if getfield(obj, :has_distributed_va)
            throw(ArgumentError(
                "body_aero.va is undefined after set_va! with distributed inflow. " *
                "Use panel.va or solver.sol._va_dist for per-panel inflow data."
            ))
        end
        return getfield(obj, :_va)
    end
    return getfield(obj, sym)
end

function Base.setproperty!(obj::BodyAerodynamics, sym::Symbol, val)
    if sym === :va
        set_va!(obj, val)
    elseif sym === :omega
        set_va!(obj, obj._va, val)
    else
        setfield!(obj, sym, val)
    end
end

@inline function _can_skip_panel_aero_reinit(wing::Wing, panels, panel_idx_start::Int)
    wing.use_prior_polar || return false
    wing.n_panels > 0 || return false
    length(wing.refined_sections) == wing.n_panels + 1 || return false

    # Only skip when panel interpolators are already initialized for polar models.
    if isempty(wing.refined_sections)
        return false
    end
    model = wing.refined_sections[1].aero_model
    if !(model in (POLAR_VECTORS, POLAR_MATRICES))
        return false
    end

    for i in 0:(wing.n_panels - 1)
        panel = panels[panel_idx_start + i]
        if panel.cl_interp === nothing || panel.cd_interp === nothing || panel.cm_interp === nothing
            return false
        end
    end
    return true
end

"""
    calculate_stall_angle_list(panels::Vector{<:Panel};
                             begin_aoa=9.0,
                             end_aoa=22.0,
                             step_aoa=1.0,
                             stall_angle_if_none_detected=50.0,
                             cl_initial=-10.0)

Calculate stall angles for each panel.

Returns:
    Vector{Float64}: Stall angles in radians
"""
function calculate_stall_angle_list(panels::Vector{<:Panel};
                                  begin_aoa=9.0,
                                  end_aoa=22.0,
                                  step_aoa=1.0,
                                  stall_angle_if_none_detected=50.0,
                                  cl_initial=-10.0)
    stall_angles = Vector{Float64}(undef, length(panels))
    calculate_stall_angle_list!(stall_angles, panels;
                                begin_aoa, end_aoa, step_aoa,
                                stall_angle_if_none_detected, cl_initial)
    return stall_angles
end

function calculate_stall_angle_list!(stall_angles::AbstractVector,
                                     panels::Vector{<:Panel};
                                     begin_aoa=9.0,
                                     end_aoa=22.0,
                                     step_aoa=1.0,
                                     stall_angle_if_none_detected=50.0,
                                     cl_initial=-10.0)

    # Pre-compute range values to avoid allocation
    n_steps = Int(floor((end_aoa - begin_aoa) / step_aoa)) + 1

    for (idx, panel) in enumerate(panels)
        # Default stall angle if none found
        panel_stall = stall_angle_if_none_detected

        # Start with minimum cl
        cl_old = cl_initial

        # Find stall angle
        for i in 0:(n_steps-1)
            aoa = deg2rad(begin_aoa + i * step_aoa)
            cl = calculate_cl(panel, aoa)
            if cl < cl_old
                panel_stall = aoa
                break
            end
            cl_old = cl
        end

        stall_angles[idx] = panel_stall
    end

    return nothing
end

"""
    reinit!(body_aero::BodyAerodynamics; init_aero, va, omega, refine_mesh, recompute_mapping, sort_sections)

Initialize a BodyAerodynamics struct in-place by setting up panels and coefficients.

# Arguments
- `body_aero::BodyAerodynamics`: The structure to initialize

# Keyword Arguments
- `init_aero::Bool`: Whether to initialize the aero data or not
- `va=[15.0, 0.0, 0.0]`: Apparent wind vector
- `omega=zeros(3)`: Turn rate in kite body frame x y and z

# Returns
nothing
"""
function reinit!(body_aero::BodyAerodynamics{P, W, T};
    init_aero=true,
    va=[15.0, 0.0, 0.0],
    omega=zeros(MVector{3, T})
) where {P, W, T}
    idx = 1
    vec = zeros(MVector{3, T})
    for wing in body_aero.wings
        reinit!(wing)
        validate_section_aero(wing.refined_sections)
        panel_props = wing.panel_props
        wing_init_aero = init_aero && !_can_skip_panel_aero_reinit(wing, body_aero.panels, idx)
        
        # Create panels
        for i in 1:wing.n_panels
            if length(wing.delta_dist) > 0
                # Panel i gets its delta directly from delta_dist[i]
                delta = wing.delta_dist[i]
            else
                delta = zero(T)
            end
            @views reinit!(
                body_aero.panels[idx], 
                wing.refined_sections[i],
                wing.refined_sections[i+1],
                panel_props.aero_centers[i, :],
                panel_props.control_points[i, :],
                panel_props.bound_points_1[i, :],
                panel_props.bound_points_2[i, :],
                panel_props.x_airf[i, :],
                panel_props.y_airf[i, :],
                panel_props.z_airf[i, :],
                delta,
                vec,
                wing.spanwise_direction;
                remove_nan=wing.remove_nan,
                init_aero=wing_init_aero
            )
            idx += 1
        end
    end

    # Initialize rest of the struct
    body_aero.projected_area = sum(calculate_projected_area, body_aero.wings)
    isempty(body_aero.panels) && throw(ArgumentError("Cannot compute c_ref: body_aero has no panels."))
    body_aero.c_ref = maximum(panel.chord for panel in body_aero.panels)
    calculate_stall_angle_list!(body_aero.stall_angle_list, body_aero.panels)
    body_aero.alpha_dist .= 0.0
    body_aero.v_a_dist .= 0.0
    body_aero.AIC .= 0.0
    set_va!(body_aero, va, omega)
    return nothing
end

"""
    _compute_reference_velocity_from_distribution(va_input, n_panels, panel_areas=nothing)

Return a single reference velocity vector from uniform or distributed inflow.
For distributed inflow, the speed is area-weighted RMS and the direction is
the area-weighted mean direction.
"""
@inline function _compute_reference_velocity_from_distribution(
    va_input::AbstractVector,
    n_panels::Int,
    panel_areas::Union{Nothing, AbstractVector}=nothing
)
    length(va_input) == 3 ||
        throw(ArgumentError("'va' must be shape (3,) or ($(n_panels), 3); got length $(length(va_input))"))
    T = eltype(va_input)
    return MVector{3, T}(va_input[1], va_input[2], va_input[3])
end

@inline function _compute_reference_velocity_from_distribution(
    va_input::AbstractMatrix,
    n_panels::Int,
    panel_areas::Union{Nothing, AbstractVector}=nothing
)
    size(va_input) == (n_panels, 3) ||
        throw(ArgumentError("'va' must be shape (3,) or ($(n_panels), 3); got $(size(va_input))"))
    if !isnothing(panel_areas)
        length(panel_areas) == n_panels ||
            throw(ArgumentError("panel_areas must be shape ($(n_panels),), got length $(length(panel_areas))"))
    end

    T = promote_type(eltype(va_input),
                     isnothing(panel_areas) ? Float64 : eltype(panel_areas))
    total_area = zero(T)
    weighted_speed_sq = zero(T)
    direction = zeros(MVector{3, T})
    @inbounds for i in 1:n_panels
        area_i = isnothing(panel_areas) ? one(T) : T(panel_areas[i])
        va1 = va_input[i, 1]
        va2 = va_input[i, 2]
        va3 = va_input[i, 3]
        speed_i = sqrt(va1^2 + va2^2 + va3^2)
        total_area += area_i
        weighted_speed_sq += area_i * speed_i^2
        direction[1] += area_i * va1
        direction[2] += area_i * va2
        direction[3] += area_i * va3
    end
    total_area > 0.0 || throw(ArgumentError("Total panel area must be positive."))

    reference_speed = sqrt(weighted_speed_sq / total_area)
    direction_norm = norm(direction)
    if direction_norm <= 0.0
        direction .= (one(T), zero(T), zero(T))
        direction_norm = one(T)
    end

    return direction ./ direction_norm .* reference_speed
end

"""
    calculate_AIC_matrices!(body_aero::BodyAerodynamics, model::Model, 
                         core_radius_fraction,
                         va_norm_array, 
                         va_unit_array)

Calculate Aerodynamic Influence Coefficient matrices.

See also: [BodyAerodynamics](@ref), [Model](@ref)

Returns: nothing
"""
@inline function calculate_AIC_matrices!(body_aero::BodyAerodynamics{P, W, T}, model::Model,
                              core_radius_fraction,
                              va_norm_array::AbstractVector{T},
                              va_unit_array::AbstractMatrix{T}) where {P, W, T}
    # Determine evaluation point based on model
    evaluation_point = model == VSM ? :control_point : :aero_center
    evaluation_point_on_bound = model == LLT

    # Allocate work vectors for this function (separate from those used by child functions)
    velocity_induced = zeros(MVector{3, T})
    tempvel = zeros(MVector{3, T})
    va_unit = zeros(MVector{3, T})
    U_2D = zeros(MVector{3, T})

    # Python parity: one shared area-weighted wake vector for all panels.
    panel_areas = [panel.chord * panel.width for panel in body_aero.panels]
    va_distribution = zeros(T, length(body_aero.panels), 3)
    @inbounds for i in 1:length(body_aero.panels), k in 1:3
        va_distribution[i, k] = va_unit_array[i, k] * va_norm_array[i]
    end
    wake_velocity = _compute_reference_velocity_from_distribution(
        va_distribution,
        length(body_aero.panels),
        panel_areas
    )
    wake_speed = norm(wake_velocity)
    wake_speed > 0.0 || throw(ArgumentError("Wake reference speed must be positive."))
    va_unit .= wake_velocity ./ wake_speed
    va_norm = wake_speed
    
    # Calculate influence coefficients
    for icp in eachindex(body_aero.panels)
        panel_icp = body_aero.panels[icp]
        ep = evaluation_point == :control_point ? panel_icp.control_point : panel_icp.aero_center
        for jring in eachindex(body_aero.panels)
            panel_jring = body_aero.panels[jring]
            filaments = panel_jring.filaments
            calculate_velocity_induced_single_ring_semiinfinite!(
                velocity_induced,
                tempvel,
                filaments,
                ep,
                evaluation_point_on_bound,
                va_norm,
                va_unit,
                one(T),
                core_radius_fraction,
                body_aero.work_vectors
            )
                      
            # Subtract 2D induced velocity for VSM
            if icp == jring && model == VSM
                calculate_velocity_induced_bound_2D!(U_2D, panel_jring, ep, body_aero.work_vectors)
                velocity_induced .-= U_2D
            end
            body_aero.AIC[:, icp, jring] .= velocity_induced
        end
    end
    return nothing
end

"""
    calculate_circulation_distribution_elliptical_wing(body_aero::BodyAerodynamics, gamma_0=1.0)

Calculate circulation distribution for an elliptical wing.

Returns: nothing
"""
function calculate_circulation_distribution_elliptical_wing(gamma_i, body_aero::BodyAerodynamics, gamma_0=1.0)
    length(body_aero.wings) == 1 || throw(ArgumentError("Multiple wings not yet implemented"))
    
    wing_span = body_aero.wings[1].span
    @debug "Wing span: $wing_span"
    
    # Calculate y-coordinates of control points
    y = body_aero.y
    for (i, panel) in pairs(body_aero.panels) 
        y[i] = panel.control_point[2] 
    end
    
    # Calculate elliptical distribution (clamp to avoid sqrt of negative
    # when control points lie outside the nominal span envelope)
    gamma_i .= gamma_0 * sqrt.(max.(0.0, 1 .- (2 .* y ./ wing_span).^2))
    
    @debug "Calculated circulation distribution: $gamma_i"
    nothing
end

"""
    update_effective_angle_of_attack_if_VSM(body_aero::BodyAerodynamics, gamma,
                                          core_radius_fraction,
                                          z_airf_array,
                                          x_airf_array,
                                          va_array,
                                          va_norm_array,
                                          va_unit_array)

Update angle of attack at aerodynamic center for VSM method.

Returns:
    nothing
"""
function update_effective_angle_of_attack!(alpha_corrected,
    body_aero::BodyAerodynamics, 
    gamma,
    core_radius_fraction,
    z_airf_array,
    x_airf_array,
    va_array,
    va_norm_array,
    va_unit_array)

    # Calculate AIC matrices (keep existing optimized view)
    calculate_AIC_matrices!(body_aero, LLT, core_radius_fraction, va_norm_array, va_unit_array)

    # Get dimensions from existing data
    n_rows = size(body_aero.AIC, 2)
    n_cols = size(body_aero.AIC, 3)

    # Preallocate induced velocity array
    induced_velocity = body_aero.cache[1][va_array]

    # Calculate each component with explicit loops
    for j in 1:3  # For each x/y/z component
        for i in 1:n_rows
            acc = zero(eltype(induced_velocity))  # Type-stable accumulator
            for k in 1:n_cols
                acc += body_aero.AIC[j, i, k] * gamma[k]
            end
            induced_velocity[i, j] = acc
        end
    end

    # In-place relative velocity calculation
    relative_velocity = body_aero.cache[2][va_array]
    relative_velocity .= va_array .+ induced_velocity

    # Preallocate and compute dot products manually
    n = size(relative_velocity, 1)
    v_normal     = body_aero.cache[3][relative_velocity]
    v_tangential = body_aero.cache[4][relative_velocity]
    
    @inbounds for i in 1:n
        vn = 0.0
        vt = 0.0
        for j in 1:3
            vn += z_airf_array[i, j] * relative_velocity[i, j]
            vt += x_airf_array[i, j] * relative_velocity[i, j]
        end
        v_normal[i] = vn
        v_tangential[i] = vt
    end

    # Direct angle calculation without temporary arrays
    @inbounds for i in 1:n
        alpha_corrected[i] = atan(v_normal[i], v_tangential[i])
    end

    nothing
end

@inline function intersect_line_with_plane(
    x_cp, f_unit, plane_point, plane_normal; tol=1e-6
)
    numerator = plane_normal[1]*(plane_point[1]-x_cp[1]) +
                plane_normal[2]*(plane_point[2]-x_cp[2]) +
                plane_normal[3]*(plane_point[3]-x_cp[3])
    denominator = dot3(plane_normal, f_unit)
    abs(denominator) < tol && return nothing
    λ = numerator / denominator
    return MVec3(x_cp[1] + λ*f_unit[1],
                 x_cp[2] + λ*f_unit[2],
                 x_cp[3] + λ*f_unit[3])
end

@inline function point_in_triangle(pt, v0, v1, v2; tol=1e-8)
    uu = 0.0; uv = 0.0; vv = 0.0; wu = 0.0; wv = 0.0
    @inbounds for k in 1:3
        uk = v1[k] - v0[k]
        vk = v2[k] - v0[k]
        wk = pt[k] - v0[k]
        uu += uk * uk
        uv += uk * vk
        vv += vk * vk
        wu += wk * uk
        wv += wk * vk
    end

    denom = uv * uv - uu * vv
    abs(denom) < 1e-12 && return false

    s = (uv * wv - vv * wu) / denom
    t = (uv * wu - uu * wv) / denom
    return (s >= -tol) && (t >= -tol) && (s + t <= 1 + tol)
end

@inline function point_in_quad(pt, corners)
    return _point_in_triangle_col(pt, corners, 1, 2, 3) ||
           _point_in_triangle_col(pt, corners, 1, 3, 4)
end

@inline function _point_in_triangle_col(
    pt, corners, c0, c1, c2; tol=1e-8
)
    uu = 0.0; uv = 0.0; vv = 0.0; wu = 0.0; wv = 0.0
    @inbounds for k in 1:3
        uk = corners[k, c1] - corners[k, c0]
        vk = corners[k, c2] - corners[k, c0]
        wk = pt[k] - corners[k, c0]
        uu += uk * uk
        uv += uk * vk
        vv += vk * vk
        wu += wk * uk
        wv += wk * vk
    end

    denom = uv * uv - uu * vv
    abs(denom) < 1e-12 && return false

    s = (uv * wv - vv * wu) / denom
    t = (uv * wu - uu * wv) / denom
    return (s >= -tol) && (t >= -tol) && (s + t <= 1 + tol)
end

function find_center_of_pressure(
    body_aero::BodyAerodynamics,
    force_array,
    moment_array,
    reference_point;
    force_tol::Float64 = 1e-12
)
    F = force_array
    M0 = moment_array
    r0 = reference_point
    F_norm_sq = dot3(F, F)
    # Treat near-zero forces as "CoP undefined"
    if !(isfinite(F_norm_sq)) || F_norm_sq ≤ force_tol^2
        return nothing
    end

    wv = body_aero.work_vectors
    r0_moment = wv[1]
    f_unit = wv[2]
    normal = wv[3]
    cross_tmp = wv[4]

    cross3!(cross_tmp, F, M0)
    F_norm = sqrt(F_norm_sq)
    @inbounds for k in 1:3
        r0_moment[k] = r0[k] + cross_tmp[k] / F_norm_sq
        f_unit[k] = F[k] / F_norm
    end

    for panel in body_aero.panels
        corners = panel.corner_points
        # cross(v1, v2) where v1 = col2-col1, v2 = col3-col1
        v1x = corners[1,2]-corners[1,1]
        v1y = corners[2,2]-corners[2,1]
        v1z = corners[3,2]-corners[3,1]
        v2x = corners[1,3]-corners[1,1]
        v2y = corners[2,3]-corners[2,1]
        v2z = corners[3,3]-corners[3,1]
        normal[1] = v1y*v2z - v1z*v2y
        normal[2] = v1z*v2x - v1x*v2z
        normal[3] = v1x*v2y - v1y*v2x
        normal_norm = norm3(normal)
        if normal_norm != 0
            normal[1] /= normal_norm
            normal[2] /= normal_norm
            normal[3] /= normal_norm
            # Avoid view allocation for plane point
            cross_tmp[1] = corners[1, 1]
            cross_tmp[2] = corners[2, 1]
            cross_tmp[3] = corners[3, 1]
            intersection = intersect_line_with_plane(
                r0_moment, f_unit, cross_tmp, normal)
            if !isnothing(intersection) &&
               point_in_quad(intersection, corners)
                return MVec3(intersection)
            end
        end
    end

    @warn "No intersection found with any panel " *
          "in center-of-pressure calculation."
    return nothing
end

function compute_panel_center_of_pressures(
    body_aero::BodyAerodynamics,
    f_distribution::AbstractMatrix,
    m_distribution::AbstractMatrix,
    reference_point
)
    n = length(body_aero.panels)
    panel_cp_locations = Vector{MVec3}(undef, n)
    for i in 1:n
        panel = body_aero.panels[i]
        @views F = f_distribution[:, i]
        @views M_ref = m_distribution[:, i]
        ac = panel.aero_center
        chord_dir = panel.x_airf
        span_dir = panel.y_airf
        c = panel.chord

        # Guard against non-finite forces and near-zero forces
        if !all(isfinite, F) || dot3(F, F) ≤ 1e-24
            panel_cp_locations[i] = MVec3(ac)
            continue
        end

        # cross(r, F) where r = ac - reference_point
        rx = ac[1]-reference_point[1]
        ry = ac[2]-reference_point[2]
        rz = ac[3]-reference_point[3]
        crx = ry*F[3] - rz*F[2]
        cry = rz*F[1] - rx*F[3]
        crz = rx*F[2] - ry*F[1]

        # m_pitch = dot(M_ref - cross(r,F), span_dir)
        m_pitch = (M_ref[1]-crx)*span_dir[1] +
                  (M_ref[2]-cry)*span_dir[2] +
                  (M_ref[3]-crz)*span_dir[3]

        # F_perp_mag = dot(cross(chord_dir, F), span_dir)
        cx = chord_dir[2]*F[3] - chord_dir[3]*F[2]
        cy = chord_dir[3]*F[1] - chord_dir[1]*F[3]
        cz = chord_dir[1]*F[2] - chord_dir[2]*F[1]
        F_perp_mag = cx*span_dir[1] + cy*span_dir[2] +
                     cz*span_dir[3]

        if abs(F_perp_mag) < 1e-12
            panel_cp_locations[i] = MVec3(ac)
            continue
        end

        lever = clamp(m_pitch / F_perp_mag,
                      -0.25 * c, 0.75 * c)
        panel_cp_locations[i] = MVec3(
            ac[1] + lever*chord_dir[1],
            ac[2] + lever*chord_dir[2],
            ac[3] + lever*chord_dir[3])
    end

    return panel_cp_locations
end

"""
    calculate_results(body_aero::BodyAerodynamics, gamma_new, 
                     density,
                     core_radius_fraction, mu,
                     alpha_dist, v_a_dist,
                     chord_array, x_airf_array,
                     z_airf_array,
                     va_array, va_norm_array,
                     va_unit_array, panels::Vector{<:Panel},
                     is_only_f_and_gamma_output::Bool)

Calculate final aerodynamic results. Reference point is in the kite body (KB) frame.

Returns:
    Dict: Results including forces, coefficients and distributions
"""
function calculate_results(
    body_aero::BodyAerodynamics,
    gamma_new,
    reference_point,
    density,
    core_radius_fraction,
    mu,
    alpha_dist,
    v_a_dist,
    chord_array,
    x_airf_array,
    z_airf_array,
    va_array,
    va_norm_array,
    va_unit_array,
    panels::Vector{<:Panel},
    is_only_f_and_gamma_output::Bool;
    correct_aoa::Bool=false,
)

    n_panels = length(panels)
    if length(body_aero.cache) < 15
        append!(body_aero.cache, [LazyBufferCache() for _ in 1:(15 - length(body_aero.cache))])
    end

    cl_array = body_aero.cache[5][alpha_dist]
    cd_array = body_aero.cache[6][alpha_dist]
    cm_array = body_aero.cache[7][alpha_dist]
    panel_width_array = body_aero.cache[8][alpha_dist]
    alpha_corrected = body_aero.cache[9][alpha_dist]
    cl_prescribed_va = body_aero.cache[10][alpha_dist]
    cd_prescribed_va = body_aero.cache[11][alpha_dist]
    cs_prescribed_va = body_aero.cache[12][alpha_dist]
    panel_view_3xn = @view body_aero.AIC[:, :, 1]
    f_body_3D = body_aero.cache[13][panel_view_3xn]
    m_body_3D = body_aero.cache[14][panel_view_3xn]
    alpha_geometric = body_aero.cache[15][alpha_dist]

    fill!(f_body_3D, 0.0)
    fill!(m_body_3D, 0.0)

    # Calculate coefficients and geometric AoA for each panel
    for (i, panel) in enumerate(panels)
        cl_array[i] = calculate_cl(panel, alpha_dist[i])
        cd_array[i], cm_array[i] = calculate_cd_cm(
            panel, alpha_dist[i])
        panel_width_array[i] = panel.width
        va_norm = va_norm_array[i]
        x_norm = norm3(panel.x_airf)
        z_norm = norm3(panel.z_airf)
        if va_norm == 0.0 || x_norm == 0.0 || z_norm == 0.0
            alpha_geometric[i] = NaN
        else
            inv_va_norm = 1.0 / va_norm
            v_tangential = -dot3(panel.x_airf, panel.va) *
                           inv_va_norm / x_norm
            v_normal = -dot3(panel.z_airf, panel.va) *
                       inv_va_norm / z_norm
            alpha_geometric[i] = atan(-v_normal, -v_tangential)
        end
    end

    # Calculate alpha corrections based on model type
    if correct_aoa
        update_effective_angle_of_attack!(
            alpha_corrected,
            body_aero,
            gamma_new,
            core_radius_fraction,
            z_airf_array,
            x_airf_array,
            va_array,
            va_norm_array,
            va_unit_array
        )
    else
        alpha_corrected .= alpha_dist
    end

    area_all_panels = 0.0
    lift_wing_3D_sum = 0.0
    drag_wing_3D_sum = 0.0
    side_wing_3D_sum = 0.0

    # Get wing properties and reference velocity
    spanwise_direction = body_aero.wings[1].spanwise_direction
    va_ref_vector = MVec3(0.0, 0.0, 0.0)
    weighted_speed_sq = 0.0
    total_area = 0.0
    @inbounds for i in 1:n_panels
        area_i = chord_array[i] * panel_width_array[i]
        total_area += area_i
        speed_i = va_norm_array[i]
        weighted_speed_sq += area_i * speed_i^2
        va_ref_vector[1] += area_i * va_array[i, 1]
        va_ref_vector[2] += area_i * va_array[i, 2]
        va_ref_vector[3] += area_i * va_array[i, 3]
    end
    total_area > 0.0 || throw(ArgumentError(
        "Total panel area must be positive."))
    reference_speed = sqrt(weighted_speed_sq / total_area)
    direction_norm = norm3(va_ref_vector)
    if direction_norm <= 0.0
        va_ref_vector .= (1.0, 0.0, 0.0)
        direction_norm = 1.0
    end
    @inbounds for k in 1:3
        va_ref_vector[k] = va_ref_vector[k] / direction_norm *
                           reference_speed
    end
    va_ref_mag = norm3(va_ref_vector)
    va_ref_mag > 0.0 || throw(ArgumentError(
        "Reference freestream magnitude must be positive."))
    va_ref_unit = body_aero.work_vectors[1]
    inv_va_ref = 1.0 / va_ref_mag
    @inbounds for k in 1:3
        va_ref_unit[k] = va_ref_vector[k] * inv_va_ref
    end
    dir_lift_ref = body_aero.work_vectors[2]
    cross3!(dir_lift_ref, va_ref_vector, spanwise_direction)
    dir_lift_ref_norm = norm3(dir_lift_ref)
    dir_lift_ref_norm > 0.0 || throw(ArgumentError(
        "Reference lift direction is undefined because " *
        "reference flow is parallel to spanwise direction."))
    @inbounds for k in 1:3
        dir_lift_ref[k] /= dir_lift_ref_norm
    end
    dir_side_ref = body_aero.work_vectors[3]
    cross3!(dir_side_ref, dir_lift_ref, va_ref_unit)
    q_ref = 0.5 * density * va_ref_mag^2

    induced_va_airfoil = body_aero.work_vectors[4]
    dir_lift_induced_va = body_aero.work_vectors[5]
    dir_drag_induced_va = body_aero.work_vectors[6]
    lift_induced_va = body_aero.work_vectors[7]
    drag_induced_va = body_aero.work_vectors[8]
    dir_lift_prescribed_va = body_aero.work_vectors[9]
    temp_vec = body_aero.work_vectors[10]

    # Main calculation loop
    for (i, panel) in enumerate(panels)
        panel_area = panel.chord * panel.width
        area_all_panels += panel_area

        alpha_corrected_i = alpha_corrected[i]
        c_alpha = cos(alpha_corrected_i)
        s_alpha = sin(alpha_corrected_i)
        @inbounds for k in 1:3
            induced_va_airfoil[k] = c_alpha * panel.x_airf[k] +
                                    s_alpha * panel.z_airf[k]
        end
        normalize3!(induced_va_airfoil)

        cross3!(dir_lift_induced_va,
                induced_va_airfoil, panel.y_airf)
        normalize3!(dir_lift_induced_va)
        cross3!(dir_drag_induced_va,
                spanwise_direction, dir_lift_induced_va)
        normalize3!(dir_drag_induced_va)

        q_lift = 0.5 * density * v_a_dist[i]^2
        lift_i = cl_array[i] * q_lift * chord_array[i]
        drag_i = cd_array[i] * q_lift * chord_array[i]
        moment_i = cm_array[i] * q_lift * chord_array[i]^2

        @inbounds for k in 1:3
            lift_induced_va[k] = lift_i * dir_lift_induced_va[k]
            drag_induced_va[k] = drag_i * dir_drag_induced_va[k]
        end

        va_panel_mag = va_norm_array[i]
        va_panel_mag > 0.0 || throw(ArgumentError(
            "Panel $i has non-positive apparent " *
            "velocity magnitude."))
        q_panel = 0.5 * density * va_panel_mag^2
        cross3!(dir_lift_prescribed_va,
                panel.va, spanwise_direction)
        normalize3!(dir_lift_prescribed_va)

        lift_prescribed_va =
            dot3(lift_induced_va, dir_lift_prescribed_va) +
            dot3(drag_induced_va, dir_lift_prescribed_va)
        drag_prescribed_va =
            (dot3(lift_induced_va, panel.va) +
             dot3(drag_induced_va, panel.va)) / va_panel_mag
        cross3!(temp_vec, dir_lift_prescribed_va, panel.va)
        inv_vpm = 1.0 / va_panel_mag
        @inbounds for k in 1:3
            temp_vec[k] *= inv_vpm
        end
        side_prescribed_va =
            dot3(lift_induced_va, temp_vec) +
            dot3(drag_induced_va, temp_vec)

        width = panel.width
        @inbounds for k in 1:3
            f_body_3D[k, i] = (lift_induced_va[k] +
                               drag_induced_va[k]) * width
        end

        lift_wing_3D_sum += lift_prescribed_va * width *
            dot3(dir_lift_prescribed_va, dir_lift_ref)
        drag_wing_3D_sum += drag_prescribed_va * width *
            (dot3(panel.va, va_ref_unit) / va_panel_mag)
        side_wing_3D_sum += side_prescribed_va * width *
            dot3(temp_vec, dir_side_ref)

        inv_qc = 1.0 / (q_panel * panel.chord)
        cl_prescribed_va[i] = lift_prescribed_va * inv_qc
        cd_prescribed_va[i] = drag_prescribed_va * inv_qc
        cs_prescribed_va[i] = side_prescribed_va * inv_qc

        ### Moment ###
        # r_vector = panel.aero_center - reference_point
        # M_shift = cross(r_vector, f_body_3D[:,i])
        # m_body_3D[:,i] = moment_i * panel.y_airf * width + M_shift
        @inbounds for k in 1:3
            dir_lift_prescribed_va[k] = panel.aero_center[k] -
                                        reference_point[k]
            drag_induced_va[k] = f_body_3D[k, i]
        end
        cross3!(temp_vec,
                dir_lift_prescribed_va, drag_induced_va)
        local_moment_scale = moment_i * width
        @inbounds for k in 1:3
            m_body_3D[k, i] = local_moment_scale *
                              panel.y_airf[k] + temp_vec[k]
        end
    end

    if is_only_f_and_gamma_output
        return Dict{String,Any}(
            "F_distribution" => copy(f_body_3D),
            "gamma_distribution" => gamma_new
        )
    end

    # Calculate wing geometry properties
    projected_area = body_aero.projected_area
    wing_span = body_aero.wings[1].span
    aspect_ratio_projected = wing_span^2 / projected_area

    # Calculate Reynolds number
    c_ref = body_aero.c_ref
    reynolds_number = density * va_ref_mag * c_ref / mu

    force_total = body_aero.work_vectors[9]
    moment_total = body_aero.work_vectors[10]
    force_total .= 0.0
    moment_total .= 0.0
    @inbounds for i in 1:n_panels
        force_total[1] += f_body_3D[1, i]
        force_total[2] += f_body_3D[2, i]
        force_total[3] += f_body_3D[3, i]
        moment_total[1] += m_body_3D[1, i]
        moment_total[2] += m_body_3D[2, i]
        moment_total[3] += m_body_3D[3, i]
    end
    center_of_pressure = try
        find_center_of_pressure(body_aero, force_total, moment_total, reference_point)
    catch err
        @warn "Center-of-pressure calculation failed: $(err)"
        nothing
    end
    panel_cp_locations = compute_panel_center_of_pressures(
        body_aero,
        f_body_3D,
        m_body_3D,
        reference_point
    )

    # Create results dictionary
    results = Dict{String,Any}(
        "Fx" => force_total[1],
        "Fy" => force_total[2],
        "Fz" => force_total[3],
        "Mx" => moment_total[1],
        "My" => moment_total[2],
        "Mz" => moment_total[3],
        "lift" => lift_wing_3D_sum,
        "drag" => drag_wing_3D_sum,
        "side" => side_wing_3D_sum,
        "cl" => lift_wing_3D_sum / (q_ref * projected_area),
        "cd" => drag_wing_3D_sum / (q_ref * projected_area),
        "cs" => side_wing_3D_sum / (q_ref * projected_area),
        "cmx" => moment_total[1] / (q_ref * projected_area * c_ref),
        "cmy" => moment_total[2] / (q_ref * projected_area * c_ref),
        "cmz" => moment_total[3] / (q_ref * projected_area * c_ref),
        "cl_distribution" => copy(cl_prescribed_va),
        "cd_distribution" => copy(cd_prescribed_va),
        "cs_distribution" => copy(cs_prescribed_va),
        "F_distribution" => copy(f_body_3D),
        "M_distribution" => copy(m_body_3D),
        "cfx" => (force_total[1] / (q_ref * projected_area)),
        "cfy" => (force_total[2] / (q_ref * projected_area)),
        "cfz" => (force_total[3] / (q_ref * projected_area)),
        "alpha_at_ac" => copy(alpha_corrected),
        "alpha_uncorrected" => alpha_dist,
        "alpha_geometric" => copy(alpha_geometric),
        "gamma_distribution" => gamma_new,
        "area_all_panels" => area_all_panels,
        "projected_area" => projected_area,
        "wing_span" => wing_span,
        "aspect_ratio_projected" => aspect_ratio_projected,
        "Rey" => reynolds_number,
        "q_ref" => q_ref,
        "va_ref" => va_ref_vector,
        "center_of_pressure" => center_of_pressure,
        "panel_cp_locations" => panel_cp_locations
    )

    @debug "Results summary:" cl=results["cl"] cd=results["cd"] cs=results["cs"]
    @debug "Forces:" lift=lift_wing_3D_sum drag=drag_wing_3D_sum side=side_wing_3D_sum
    @debug "Areas:" total=area_all_panels projected=projected_area
    @debug "Aspect ratio:" ar=aspect_ratio_projected

    return results
end


"""
    set_va!(body_aero::BodyAerodynamics, va::VelVector, omega=zeros(MVec3))

Set velocity array and update wake filaments.

# Arguments
- body_aero::BodyAerodynamics: The [BodyAerodynamics](@ref) struct to modify
- `va::VelVector`: Velocity vector of the apparent wind speed           [m/s]
- `omega::VelVector`: Turn rate vector around x y and z axis            [rad/s]
"""
function set_va!(body_aero::BodyAerodynamics{P, W, T}, va::AbstractVector, omega=zeros(MVector{3, T})) where {P, W, T}
    n_panels = length(body_aero.panels)
    va_distribution = zeros(T, n_panels, 3)
    body_aero.omega .= omega

    if all(iszero, omega)
        va_distribution .= reshape(va, 1, 3)
    else
        idx = 1
        for wing in body_aero.wings
            panel_end = idx + wing.n_panels - 1

            # Calculate velocities for each panel in this wing slice
            for j in idx:panel_end
                omega_va = -omega × body_aero.panels[j].control_point
                va_distribution[j, :] .= omega_va .+ va
            end
            idx = panel_end + 1
        end
    end

    # Update panel velocities
    for (i, panel) in enumerate(body_aero.panels)
        panel.va .= va_distribution[i,:]
    end

    # Update wake elements
    frozen_wake!(body_aero, va_distribution)
    body_aero._va .= va
    body_aero.has_distributed_va = false
    return nothing
end

function set_va!(body_aero::BodyAerodynamics, va_distribution::AbstractMatrix)
    size(va_distribution, 1) != length(body_aero.panels) &&
        throw(ArgumentError("Number of rows in va distribution should be equal to number of panels."))

    for (i, panel) in enumerate(body_aero.panels)
        panel.va .= va_distribution[i, :]
    end

    # Update wake elements
    frozen_wake!(body_aero, va_distribution)
    body_aero._va .= [mean(va_distribution[:,i]) for i in 1:3]
    body_aero.has_distributed_va = true
    return nothing
end

"""
    set_va!(body_aero::BodyAerodynamics, settings::VSMSettings)

Set velocity array from VSM settings configuration.

This convenience method extracts flight conditions from VSMSettings and 
constructs the velocity vector in the body reference frame based on:
- Wind speed from settings.condition.wind_speed
- Angle of attack from settings.condition.alpha (converted from degrees)
- Sideslip angle from settings.condition.beta (converted from degrees)

The velocity vector is constructed as:
- X_b (forward): wind_speed * cos(α) * cos(β)  
- Y_b (right): wind_speed * sin(β)
- Z_b (down): wind_speed * sin(α) * cos(β)

# Arguments
- `body_aero::BodyAerodynamics`: The aerodynamic body to modify
- `settings::VSMSettings`: Settings object containing flight conditions

# Example
```julia
settings = VSMSettings("path/to/settings.yaml")
body_aero = BodyAerodynamics([wing])
set_va!(body_aero, settings)
```
"""
function set_va!(body_aero::BodyAerodynamics, settings::VSMSettings)
    α = deg2rad(settings.condition.alpha)
    β = deg2rad(settings.condition.beta)
    wind_speed = settings.condition.wind_speed
    
    va = wind_speed * [
        cos(α)*cos(β),  # X_b (forward)
        sin(β),         # Y_b (right)
        sin(α)*cos(β)   # Z_b (down)
    ]
    
    set_va!(body_aero, va)
end

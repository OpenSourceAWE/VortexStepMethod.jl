"""
    @with_kw mutable struct BodyAerodynamics{P}

Main structure for calculating aerodynamic properties of bodies.

# Fields
- panels::Vector{Panel}: Vector of [Panel](@ref) structs
- wings::Union{Vector{Wing}, Vector{RamAirWing}}: A vector of wings; a body can have multiple wings
- `_va`::MVec3 = zeros(MVec3):   A vector of the apparent wind speed, see: [MVec3](@ref)
- `omega`::MVec3 = zeros(MVec3): A vector of the turn rates around the kite body axes
- `gamma_distribution`=zeros(Float64, P): A vector of the circulation 
                        of the velocity field; Length: Number of segments. [m²/s]
- `alpha_uncorrected`=zeros(Float64, P): angles of attack per panel
- `alpha_corrected`=zeros(Float64, P):   corrected angles of attack per panel
- `stall_angle_list`=zeros(Float64, P):  stall angle per panel
- `alpha_array` = zeros(Float64, P)
- `v_a_array` = zeros(Float64, P)
- `work_vectors`::NTuple{10, MVec3} = ntuple(_ -> zeros(MVec3), 10)
- AIC::Array{Float64, 3} = zeros(3, P, P)
- `projected_area`::Float64 = 1.0: The area projected onto the xy-plane of the kite body reference frame [m²]
"""
@with_kw mutable struct BodyAerodynamics{P}
    panels::Vector{Panel}
    wings::Union{Vector{Wing}, Vector{RamAirWing}}
    _va::MVec3 = zeros(MVec3)
    omega::MVec3 = zeros(MVec3)
    gamma_distribution::MVector{P, Float64} = zeros(MVector{P, Float64})
    alpha_uncorrected::MVector{P, Float64} = zeros(MVector{P, Float64})
    alpha_corrected::MVector{P, Float64} = zeros(MVector{P, Float64})
    stall_angle_list::MVector{P, Float64} = zeros(MVector{P, Float64})
    alpha_array::MVector{P, Float64} = zeros(MVector{P, Float64})
    v_a_array::MVector{P, Float64} = zeros(MVector{P, Float64})
    work_vectors::NTuple{10,MVec3} = ntuple(_ -> zeros(MVec3), 10)
    AIC::Array{Float64, 3} = zeros(3, P, P)
    projected_area::Float64 = one(Float64)
    y::MVector{P, Float64} = zeros(MVector{P, Float64})
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:5]
end

"""
    BodyAerodynamics(wings::Vector{T}; 
                     kite_body_origin=zeros(MVec3)) where T <: AbstractWing

Construct a [BodyAerodynamics](@ref) object for aerodynamic calculations.

# Arguments
- `wings::Vector{T}`: Vector of wings to analyze, where T is an AbstractWing type

# Keyword Arguments
- `kite_body_origin=zeros(MVec3)`: Origin point of kite body reference frame in CAD reference frame

# Returns
- [BodyAerodynamics](@ref) object initialized with panels and wings
"""
function BodyAerodynamics(
    wings::Vector{T};
    kite_body_origin=zeros(MVec3)
) where T <: AbstractWing
    # Initialize panels
    panels = Panel[]
    for wing in wings
        for section in wing.sections
            section.LE_point .-= kite_body_origin
            section.TE_point .-= kite_body_origin
        end
        if wing.spanwise_distribution == UNCHANGED
            wing.n_panels = length(wing.sections) - 1
            wing.refined_sections = wing.sections
        else
            wing.refined_sections = Section[Section() for _ in 1:wing.n_panels+1]
        end

        # Create panels
        for _ in 1:wing.n_panels
            panel = Panel()
            push!(panels, panel)
        end
    end

    body_aero = BodyAerodynamics{length(panels)}(; panels, wings)
    init!(body_aero)
    return body_aero
end

function Base.getproperty(obj::BodyAerodynamics, sym::Symbol)
    if sym === :va
        return getfield(obj, :_va)
    end
    return getfield(obj, sym)
end

function Base.setproperty!(obj::BodyAerodynamics, sym::Symbol, val)
    if sym === :va || sym === :omega
        set_va!(obj, val)
    else
        setfield!(obj, sym, val)
    end
end

"""
    init!(body_aero::BodyAerodynamics)

Initialize a BodyAerodynamics struct in-place by setting up panels and coefficients.

# Arguments
- `body_aero::BodyAerodynamics`: The structure to initialize

# Keyword Arguments
- `init_aero::Bool`: Wether to initialize the aero data or not

# Returns
nothing
"""
function init!(body_aero::BodyAerodynamics; init_aero=true)
    idx = 1
    vec = zeros(MVec3)
    for wing in body_aero.wings
        init!(wing)
        panel_props = wing.panel_props
        
        # Create panels
        for i in 1:wing.n_panels
            if wing isa RamAirWing
                delta = wing.delta_dist[i]
            else
                delta = 0.0
            end
            @views init!(
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
                vec;
                remove_nan=wing.remove_nan,
                init_aero
            )
            idx += 1
        end
    end
    
    # Initialize rest of the struct
    body_aero.projected_area = sum(wing -> calculate_projected_area(wing), body_aero.wings)
    body_aero.stall_angle_list .= calculate_stall_angle_list(body_aero.panels)
    body_aero.alpha_array .= 0.0
    body_aero.v_a_array .= 0.0 
    body_aero.AIC .= 0.0
    return nothing
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
@inline function calculate_AIC_matrices!(body_aero::BodyAerodynamics, model::Model,
                              core_radius_fraction,
                              va_norm_array, 
                              va_unit_array)
    # Determine evaluation point based on model
    evaluation_point = model == VSM ? :control_point : :aero_center
    evaluation_point_on_bound = model == LLT
    
    # Initialize AIC matrices
    velocity_induced, tempvel, va_unit, U_2D = zeros(MVec3), zeros(MVec3), zeros(MVec3), zeros(MVec3)
    
    # Calculate influence coefficients
    for icp in eachindex(body_aero.panels)
        ep = getproperty(body_aero.panels[icp], evaluation_point)
        for jring in eachindex(body_aero.panels)
            va_unit .= @views va_unit_array[jring, :]
            filaments = body_aero.panels[jring].filaments
            va_norm = va_norm_array[jring]
            calculate_velocity_induced_single_ring_semiinfinite!(
                velocity_induced,
                tempvel,
                filaments,
                ep,
                evaluation_point_on_bound,
                va_norm,
                va_unit,
                1.0,
                core_radius_fraction,
                body_aero.work_vectors
            )
                      
            # Subtract 2D induced velocity for VSM
            if icp == jring && model == VSM
                calculate_velocity_induced_bound_2D!(U_2D, body_aero.panels[jring], ep, body_aero.work_vectors)
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
    
    # Calculate elliptical distribution
    gamma_i .= gamma_0 * sqrt.(1 .- (2 .* y ./ wing_span).^2)
    
    @debug "Calculated circulation distribution: $gamma_i"
    nothing
end

"""
    calculate_stall_angle_list(panels::Vector{Panel};
                             begin_aoa=9.0,
                             end_aoa=22.0,
                             step_aoa=1.0,
                             stall_angle_if_none_detected=50.0,
                             cl_initial=-10.0)

Calculate stall angles for each panel.

Returns:
    Vector{Float64}: Stall angles in radians
"""
function calculate_stall_angle_list(panels::Vector{Panel};
                                  begin_aoa=9.0,
                                  end_aoa=22.0,
                                  step_aoa=1.0,
                                  stall_angle_if_none_detected=50.0,
                                  cl_initial=-10.0)
    
    aoa_range = deg2rad.(range(begin_aoa, end_aoa, step=step_aoa))
    stall_angles = Float64[]
    
    for panel in panels
        # Default stall angle if none found
        panel_stall = stall_angle_if_none_detected
        
        # Start with minimum cl
        cl_old = cl_initial
        
        # Find stall angle
        for aoa in aoa_range
            cl = calculate_cl(panel, aoa)
            if cl < cl_old
                panel_stall = aoa
                break
            end
            cl_old = cl
        end
        
        push!(stall_angles, panel_stall)
    end
    
    return stall_angles
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
    Vector{Float64}: Updated angles of attack
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

"""
    calculate_results(body_aero::BodyAerodynamics, gamma_new, 
                     density, aerodynamic_model_type::Model,
                     core_radius_fraction, mu,
                     alpha_array, v_a_array,
                     chord_array, x_airf_array,
                     y_airf_array, z_airf_array,
                     va_array, va_norm_array,
                     va_unit_array, panels::Vector{Panel},
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
    aerodynamic_model_type::Model,
    core_radius_fraction,
    mu,
    alpha_array,
    v_a_array,
    chord_array,
    x_airf_array,
    y_airf_array,
    z_airf_array,
    va_array,
    va_norm_array,
    va_unit_array,
    panels::Vector{Panel},
    is_only_f_and_gamma_output::Bool,
)

    # Initialize arrays
    n_panels = length(panels)
    cl_array = zeros(n_panels)
    cd_array = zeros(n_panels)
    cm_array = zeros(n_panels)
    panel_width_array = zeros(n_panels)
    alpha_corrected = zeros(n_panels)

    # Calculate coefficients for each panel
    for (i, panel) in enumerate(panels)
        cl_array[i] = calculate_cl(panel, alpha_array[i])
        cd_array[i], cm_array[i] = calculate_cd_cm(panel, alpha_array[i])
        panel_width_array[i] = panel.width
    end

    # Calculate forces
    lift = reshape((cl_array .* 0.5 .* density .* v_a_array.^2 .* chord_array), :, 1)
    drag = reshape((cd_array .* 0.5 .* density .* v_a_array.^2 .* chord_array), :, 1)
    moment = reshape((cm_array .* 0.5 .* density .* v_a_array.^2 .* chord_array), :, 1)

    # Calculate alpha corrections based on model type
    if aerodynamic_model_type == VSM
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
    elseif aerodynamic_model_type == LLT
        alpha_corrected .= alpha_array
    end

    # Verify va is not distributed
    if length(body_aero.va) != 3
        throw(ArgumentError("calculate_results not ready for va_distributed input"))
    end

    # Initialize result arrays
    cl_prescribed_va = Float64[]
    cd_prescribed_va = Float64[]
    cs_prescribed_va = Float64[]
    f_body_3D = zeros(3, n_panels)
    m_body_3D = zeros(3, n_panels)
    area_all_panels = 0.0
    
    # Initialize force sums
    lift_wing_3D_sum = 0.0
    drag_wing_3D_sum = 0.0
    side_wing_3D_sum = 0.0

    # Get wing properties
    spanwise_direction = body_aero.wings[1].spanwise_direction
    va_mag = norm(body_aero.va)
    va = body_aero.va
    va_unit = va / va_mag
    q_inf = 0.5 * density * va_mag^2

    # Main calculation loop
    for (i, panel) in enumerate(panels)
        ### Lift and Drag ###
        # Panel geometry
        panel_area = panel.chord * panel.width
        area_all_panels += panel_area

        # Calculate induced velocity direction
        alpha_corrected_i = alpha_corrected[i]
        induced_va_airfoil = cos(alpha_corrected_i) * panel.x_airf + 
                            sin(alpha_corrected_i) * panel.z_airf
        dir_induced_va_airfoil = induced_va_airfoil / norm(induced_va_airfoil)

        # Calculate lift and drag directions
        dir_lift_induced_va = cross(dir_induced_va_airfoil, panel.y_airf)
        dir_lift_induced_va = dir_lift_induced_va / norm(dir_lift_induced_va)
        dir_drag_induced_va = cross(spanwise_direction, dir_lift_induced_va)
        dir_drag_induced_va = dir_drag_induced_va / norm(dir_drag_induced_va)

        # Calculate force vectors
        lift_induced_va = lift[i] * dir_lift_induced_va
        drag_induced_va = drag[i] * dir_drag_induced_va
        ftotal_induced_va = lift_induced_va + drag_induced_va

        # Calculate forces in prescribed wing frame
        dir_lift_prescribed_va = cross(va, spanwise_direction)
        dir_lift_prescribed_va = dir_lift_prescribed_va / norm(dir_lift_prescribed_va)

        # Calculate force components
        lift_prescribed_va = dot(lift_induced_va, dir_lift_prescribed_va) + 
                           dot(drag_induced_va, dir_lift_prescribed_va)
        drag_prescribed_va = dot(lift_induced_va, va_unit) + 
                           dot(drag_induced_va, va_unit)
        side_prescribed_va = dot(lift_induced_va, spanwise_direction) + 
                           dot(drag_induced_va, spanwise_direction)

        # Body frame forces
        f_body_3D[:,i] .= ftotal_induced_va .* panel.width

        # Update sums
        lift_wing_3D_sum += lift_prescribed_va * panel.width
        drag_wing_3D_sum += drag_prescribed_va * panel.width  
        side_wing_3D_sum += side_prescribed_va * panel.width

        # Store coefficients
        push!(cl_prescribed_va, lift_prescribed_va / (q_inf * panel.chord))
        push!(cd_prescribed_va, drag_prescribed_va / (q_inf * panel.chord))
        push!(cs_prescribed_va, side_prescribed_va / (q_inf * panel.chord))

        ### Moment ###
        # (1) Panel aerodynamic center in body frame:
        panel_ac_body = panel.aero_center  # 3D [x, y, z]

        # (2) Convert local (2D) pitching moment to a 3D vector in body coords.
        #     Use the axis around which the moment is defined,
        #     which is the y-axis pointing "spanwise"
        moment_axis_body = panel.y_airf

        # Scale by panel width if your 'moment[i]' is 2D moment-per-unit-span:
        M_local_3D = moment[i] * moment_axis_body * panel.width

        # Vector from panel AC to the chosen reference point:
        r_vector = panel_ac_body - reference_point  # e.g. CG, wing root, etc.

        # Cross product to shift the force from panel AC to ref. point:
        M_shift = cross(r_vector, f_body_3D[:,i])

        # Total panel moment about the reference point:
        m_body_3D[:,i] = M_local_3D + M_shift
    end

    if is_only_f_and_gamma_output
        return Dict{String,Any}(
            "F_distribution" => f_body_3D,
            "gamma_distribution" => gamma_new
        )
    end

    # Calculate wing geometry properties
    projected_area = body_aero.projected_area
    wing_span = body_aero.wings[1].span
    aspect_ratio_projected = wing_span^2 / projected_area

    # Calculate geometric angle of attack
    horizontal_direction = [1.0, 0.0, 0.0]
    alpha_geometric = [rad2deg(acos(dot(panel.x_airf, horizontal_direction) /
                     (norm(panel.x_airf) * norm(horizontal_direction))))
                     for panel in panels]

    # Calculate Reynolds number
    max_chord = maximum(panel.chord for panel in panels)
    reynolds_number = density * va_mag * max_chord / mu

    # Create results dictionary
    results = Dict{String,Any}(
        "Fx" => sum(f_body_3D[1,:]),
        "Fy" => sum(f_body_3D[2,:]),
        "Fz" => sum(f_body_3D[3,:]),
        "Mx" => sum(m_body_3D[1,:]),
        "My" => sum(m_body_3D[2,:]),
        "Mz" => sum(m_body_3D[3,:]),
        "lift" => lift_wing_3D_sum,
        "drag" => drag_wing_3D_sum,
        "side" => side_wing_3D_sum,
        "cl" => lift_wing_3D_sum / (q_inf * projected_area),
        "cd" => drag_wing_3D_sum / (q_inf * projected_area),
        "cs" => side_wing_3D_sum / (q_inf * projected_area),
        "cmx" => sum(m_body_3D[1,:]) / (q_inf * projected_area * max_chord),
        "cmy" => sum(m_body_3D[2,:]) / (q_inf * projected_area * max_chord),
        "cmz" => sum(m_body_3D[3,:]) / (q_inf * projected_area * max_chord),
        "cl_distribution" => cl_prescribed_va,
        "cd_distribution" => cd_prescribed_va,
        "cs_distribution" => cs_prescribed_va,
        "F_distribution" => f_body_3D,
        "cfx" => (sum(f_body_3D[1,:]) / (q_inf * projected_area)),
        "cfy" => (sum(f_body_3D[2,:]) / (q_inf * projected_area)),
        "cfz" => (sum(f_body_3D[3,:]) / (q_inf * projected_area)),
        "alpha_at_ac" => alpha_corrected,
        "alpha_uncorrected" => alpha_array,
        "alpha_geometric" => alpha_geometric,
        "gamma_distribution" => gamma_new,
        "area_all_panels" => area_all_panels,
        "projected_area" => projected_area,
        "wing_span" => wing_span,
        "aspect_ratio_projected" => aspect_ratio_projected,
        "Rey" => reynolds_number
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
- `va`: Velocity vector of the apparent wind speed           [m/s]
- `omega`: Turn rate vector around x y and z axis            [rad/s]
"""
function set_va!(body_aero::BodyAerodynamics, va::AbstractVector, omega=zeros(MVec3))
    
    # Calculate va_distribution based on input type
    va_distribution = if all(omega .== 0.0)
        repeat(reshape(va, 1, 3), length(body_aero.panels))
    elseif !all(omega .== 0.0)
        va_dist = zeros(length(body_aero.panels), 3)
        
        for wing in body_aero.wings
            # Get spanwise positions
            spanwise_positions = [panel.control_point for panel in body_aero.panels]
            
            # Calculate velocities for each panel
            for i in 1:wing.n_panels
                omega_va = -omega × spanwise_positions[i]
                va_dist[i, :] .= omega_va .+ va
            end
        end
        va_dist
    end
    
    # Update panel velocities
    for (i, panel) in enumerate(body_aero.panels)
        panel.va .= va_distribution[i,:]
    end
    
    # Update wake elements
    frozen_wake!(body_aero, va_distribution)
    body_aero._va .= va
    return nothing
end

function set_va!(body_aero::BodyAerodynamics, va_distribution::AbstractVector{AbstractVector}, omega=zeros(MVec3))
    length(va) != length(body_aero.panels) && throw(ArgumentError("Length of va distribution should be equal to number of panels."))
    
    for (i, panel) in enumerate(body_aero.panels)
        panel.va = va_distribution[i]
    end
    
    # Update wake elements
    frozen_wake!(body_aero, va_distribution)
    body_aero._va = va
    return nothing
end


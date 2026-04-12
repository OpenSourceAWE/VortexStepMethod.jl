"""
    @with_kw mutable struct Section

Represents a wing section with leading edge, trailing edge, and aerodynamic properties.

# Fields
- `LE_point::MVec3` = zeros(MVec3): Leading edge point coordinates
- `TE_point::MVec3` = zeros(MVec3): Trailing edge point coordinates
- `aero_model`::AeroModel = INVISCID: [AeroModel](@ref)
- `aero_data`::AeroData = nothing: See: [AeroData](@ref)
"""
@with_kw mutable struct Section
    LE_point::MVec3 = zeros(MVec3)
    TE_point::MVec3 = zeros(MVec3)
    aero_model::AeroModel = INVISCID
    aero_data::AeroData = nothing
end

"""
    Section(LE_point::PosVector, TE_point::PosVector, aero_model)

Create a new wing section with the specified leading edge point, trailing edge point, 
and aerodynamic model.

# Arguments
- `LE_point::PosVector`: Leading edge point coordinates
- `TE_point::PosVector`: Trailing edge point coordinates  
- `aero_model::AeroModel`: Aerodynamic model type (e.g., INVISCID, POLAR_VECTORS)

# Returns
- `Section`: A new section with the specified parameters and no aerodynamic data
"""
function Section(LE_point, TE_point, aero_model)
    return Section(LE_point, TE_point, aero_model, nothing)
end

"""
    reinit!(section::Section, LE_point::PosVector, TE_point::PosVector, aero_model=nothing, aero_data=nothing)

Function to update a [Section](@ref) in place.
"""
function reinit!(section::Section, LE_point, TE_point, aero_model=nothing, aero_data=nothing)
    section.LE_point .= LE_point
    section.TE_point .= TE_point
    (!isnothing(aero_model)) && (section.aero_model = aero_model)
    if !isnothing(aero_data)
        # NTuple is immutable, so we must assign directly
        # For mutable types (Vector, Matrix), we can broadcast for efficiency
        if aero_data isa NTuple || aero_data isa Tuple || isnothing(section.aero_data)
            section.aero_data = aero_data
        else
            section.aero_data .= aero_data
        end
    end
    nothing
end

function reinit!(refined_section::Section, section::Section)
    reinit!(
        refined_section,
        section.LE_point,
        section.TE_point,
        section.aero_model,
        section.aero_data,
    )
end

"""
    PanelProperties

Structure to hold calculated panel properties.

# Fields
- `aero_centers`::Matrix{Float64}
- `control_points`::Matrix{Float64}
- `bound_points_1`::Matrix{Float64}
- `bound_points_2`::Matrix{Float64}
- `x_airf`::Matrix{Float64}: Vector of unit vectors tangential to chord line
- `y_airf`::Matrix{Float64}: Vector of unit vectors in spanwise direction
- `z_airf`::Matrix{Float64}: Vector of unit vectors pointing up (cross of x_airf and y_airf)
"""
@with_kw mutable struct PanelProperties{P}
    aero_centers::Matrix{Float64} = zeros(P, 3)
    control_points::Matrix{Float64} = zeros(P, 3)
    bound_points_1::Matrix{Float64} = zeros(P, 3)
    bound_points_2::Matrix{Float64} = zeros(P, 3)
    x_airf::Matrix{Float64} = zeros(P, 3)
    y_airf::Matrix{Float64} = zeros(P, 3)
    z_airf::Matrix{Float64} = zeros(P, 3)
    coords::Matrix{Float64} = zeros(2(P+1), 3)
end

"""
    update_panel_properties!(panel_props::PanelProperties, section_list::Vector{Section}, n_panels::Int)

Update geometric properties for each panel.

# Arguments
- section_list::Vector{Section}: List of [Section](@ref)s
- `n_panels`::Int: Number of [Panel](@ref)s

# Returns:
`nothing`, updates the [PanelProperties](@ref) in-place
"""
function update_panel_properties!(panel_props::PanelProperties, section_list::Vector{Section}, n_panels)
    coords = panel_props.coords
    aero_centers = panel_props.aero_centers
    control_points = panel_props.control_points
    bound_points_1 = panel_props.bound_points_1
    bound_points_2 = panel_props.bound_points_2
    x_airf = panel_props.x_airf
    y_airf = panel_props.y_airf
    z_airf = panel_props.z_airf
    vec = zeros(MVec3)
    vec2 = zeros(MVec3)
    @debug "Shape of coordinates: $(size(coords))"
    
    for i in 1:n_panels
        coords[2i-1, :] .= section_list[i].LE_point
        coords[2i, :]   .= section_list[i].TE_point
        coords[2i+1, :] .= section_list[i+1].LE_point
        coords[2i+2, :] .= section_list[i+1].TE_point
    end
    
    @debug "Coordinates: $coords"
    
    for i in 1:n_panels
        # Define panel points
        @views begin
            LE_1 = coords[2i-1, :]     # LE_1
            LE_2 = coords[2i+1, :]     # LE_2
            TE_2 = coords[2i+2, :]     # TE_2
            TE_1 = coords[2i, :]       # TE_1
        end
        
        # Calculate control point position
        @views @. vec = coords[2i-1, :] * 0.75 + coords[2i, :] * 0.25 - 
            (coords[2i+1, :] * 0.75 + coords[2i+2, :] * 0.25)
        di = norm(vec)
        
        ncp = if i == 1
            @views @. vec = coords[2i+1, :] * 0.75 + coords[2i+2, :] * 0.25 - 
                            (coords[2i+3, :] * 0.75 + coords[2i+4, :] * 0.25)
            diplus = norm(vec)
            di / (di + diplus)
        elseif i == n_panels
            @views @. vec = coords[2i-3, :] * 0.75 + coords[2i-2, :] * 0.25 - 
                            (coords[2i-1, :] * 0.75 + coords[2i, :] * 0.25)
            dimin = norm(vec)
            dimin / (dimin + di)
        else
            @views @. vec = coords[2i-3, :] * 0.75 + coords[2i-2, :] * 0.25 - 
                            (coords[2i-1, :] * 0.75 + coords[2i, :] * 0.25)
            dimin = norm(vec)
            @views @. vec = coords[2i+1, :] * 0.75 + coords[2i+2, :] * 0.25 - 
                            (coords[2i+3, :] * 0.75 + coords[2i+4, :] * 0.25)
            diplus = norm(vec)
            0.25 * (dimin / (dimin + di) + di / (di + diplus) + 1)
        end
        ncp = 1 - ncp
        
        # Calculate points
        @. begin
            aero_centers[i, :] = (LE_2 * (1 - ncp) + LE_1 * ncp) * 0.75 +
                    (TE_2 * (1 - ncp) + TE_1 * ncp) * 0.25        
            control_points[i, :] = (LE_2 * (1 - ncp) + LE_1 * ncp) * 0.25 +
                        (TE_2 * (1 - ncp) + TE_1 * ncp) * 0.75
            
            bound_points_1[i, :] = LE_1 * 0.75 + TE_1 * 0.25
            bound_points_2[i, :] = LE_2 * 0.75 + TE_2 * 0.25
        end
        
        # Calculate reference frame vectors
        @views begin
            @. vec = (control_points[i, :] - aero_centers[i, :])
            @. vec2 = (LE_1 - LE_2)
            vec .= vec × vec2 
            z_airf[i, :] .= normalize(vec)
            @. vec = control_points[i, :] .- aero_centers[i, :]
            x_airf[i, :] .= normalize(vec)
            @. vec = bound_points_1[i, :] - bound_points_2[i, :]
            y_airf[i, :] .= normalize(vec)
        end
    end
    return nothing
end


"""
    Wing

Represents a wing composed of multiple sections with aerodynamic properties.

# Core Fields (all wings)
- `n_panels::Int16`: Number of panels in aerodynamic mesh
- `n_unrefined_sections::Int16`: Number of unrefined sections (sections before mesh refinement)
- `spanwise_distribution`::PanelDistribution: [PanelDistribution](@ref)
- `spanwise_direction::MVec3`: Wing span direction vector
- `sections::Vector{Section}`: Vector of wing sections, see: [Section](@ref)
- `refined_sections::Vector{Section}`: Vector of refined wing sections, see: [Section](@ref)
- `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not
- `use_prior_polar::Bool`: Keep previously-initialized section/panel polar data when refining geometry updates

# Deformation Fields (optional, for deformable wings)
- `non_deformed_sections::Vector{Section}`: Original undeformed sections
- `theta_dist::Vector{Float64}`: Panel twist angle distribution
- `delta_dist::Vector{Float64}`: Trailing edge deflection distribution

# Physical Properties (optional, for OBJ-based wings)
- `mass::Float64`: Total wing mass in kg (0.0 if not applicable)
- `gamma_tip::Float64`: Angular extent from center to wing tip (0.0 if not applicable)
- `inertia_tensor::Matrix{Float64}`: 3x3 inertia tensor (empty if not applicable)
- `T_cad_body::MVec3`: Translation from CAD to body frame (zeros if not applicable)
- `R_cad_body::MMat3`: Rotation from CAD to body frame (identity if not applicable)
- `radius::Float64`: Wing curvature radius (0.0 if not applicable)
- `le_interp::Union{Nothing, NTuple{3, Extrapolation}}`: Leading edge interpolation
- `te_interp::Union{Nothing, NTuple{3, Extrapolation}}`: Trailing edge interpolation
- `area_interp::Union{Nothing, Extrapolation}`: Area interpolation
- `cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}`: Preallocated buffers

"""
mutable struct Wing{P} <: AbstractWing
    n_panels::Int16
    n_unrefined_sections::Int16
    spanwise_distribution::PanelDistribution
    panel_props::PanelProperties{P}
    spanwise_direction::MVec3
    unrefined_sections::Vector{Section}
    refined_sections::Vector{Section}
    remove_nan::Bool
    use_prior_polar::Bool
    billowing_percentage::Float64  # TE billow as percentage of arc length (0=flat)

    # Grouping
    refined_panel_mapping::Vector{Int16}  # Maps each refined panel index to unrefined section index (1 to n_unrefined_sections)

    # Deformation fields
    non_deformed_sections::Vector{Section}
    theta_dist::Vector{Float64}  # Length: n_panels (panel twist angles)
    delta_dist::Vector{Float64}  # Length: n_panels (panel TE deflection angles)

    # Physical properties (OBJ-based wings)
    mass::Float64
    gamma_tip::Float64
    inertia_tensor::Matrix{Float64}
    T_cad_body::MVec3
    R_cad_body::MMat3
    radius::Float64
    le_interp::Union{Nothing, NTuple{3, Interpolations.Extrapolation}}
    te_interp::Union{Nothing, NTuple{3, Interpolations.Extrapolation}}
    area_interp::Union{Nothing, Interpolations.Extrapolation}
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}
end

# Compatibility constructor for full positional initialization with integer panel counts.
# This keeps call sites that pass Int values working after n_panels/n_unrefined_sections
# were tightened to Int16 fields.
function Wing(
        n_panels::Integer,
        n_unrefined_sections::Integer,
        spanwise_distribution::PanelDistribution,
        panel_props::PanelProperties{P},
        spanwise_direction::MVec3,
        unrefined_sections::Vector{Section},
        refined_sections::Vector{Section},
        remove_nan::Bool,
        use_prior_polar::Bool,
        billowing_percentage::Float64,
        refined_panel_mapping::Vector{Int16},
        non_deformed_sections::Vector{Section},
        theta_dist::Vector{Float64},
        delta_dist::Vector{Float64},
        mass::Float64,
        gamma_tip::Float64,
        inertia_tensor::Matrix{Float64},
        T_cad_body::MVec3,
        R_cad_body::MMat3,
        radius::Float64,
        le_interp::Union{Nothing, NTuple{3, Interpolations.Extrapolation}},
        te_interp::Union{Nothing, NTuple{3, Interpolations.Extrapolation}},
        area_interp::Union{Nothing, Interpolations.Extrapolation},
        cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}
    ) where {P}

    n_panels_i16 = Int16(n_panels)
    n_unrefined_sections_i16 = Int16(n_unrefined_sections)

    Int(n_panels_i16) == P || throw(ArgumentError(
        "n_panels ($n_panels) must match PanelProperties{$P}"
    ))

    return Wing{P}(
        n_panels_i16,
        n_unrefined_sections_i16,
        spanwise_distribution,
        panel_props,
        spanwise_direction,
        unrefined_sections,
        refined_sections,
        remove_nan,
        use_prior_polar,
        billowing_percentage,
        refined_panel_mapping,
        non_deformed_sections,
        theta_dist,
        delta_dist,
        mass,
        gamma_tip,
        inertia_tensor,
        T_cad_body,
        R_cad_body,
        radius,
        le_interp,
        te_interp,
        area_interp,
        cache
    )
end

"""
    Wing(n_panels::Int;
         n_unrefined_sections=nothing,
         spanwise_distribution::PanelDistribution=LINEAR,
         spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]),
         remove_nan::Bool=true,
         use_prior_polar::Bool=false,
         billowing_percentage::Float64=0.0)

Constructor for a [Wing](@ref) struct with default values that initializes the sections
and refined sections as empty arrays. Creates a basic wing suitable for YAML-based construction.

# Parameters
- `n_panels::Int`: Number of panels in aerodynamic mesh
- `n_unrefined_sections::Int`: Number of unrefined sections (inferred from added sections for YAML wings)
- `spanwise_distribution`::PanelDistribution = LINEAR: [PanelDistribution](@ref)
- `spanwise_direction::MVec3` = MVec3([0.0, 1.0, 0.0]): Wing span direction vector
- `remove_nan::Bool`: Whether to remove the NaNs from interpolations or not
- `use_prior_polar::Bool`: Reuse prior refined/panel polar mapping during geometry-only updates
- `billowing_percentage::Float64`: TE billow as percentage of arc length (0=flat)
"""
function Wing(n_panels::Int;
    n_unrefined_sections::Union{Nothing, Int}=nothing,
        spanwise_distribution::PanelDistribution=LINEAR,
        spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]),
        remove_nan::Bool=true,
        use_prior_polar::Bool=false,
        billowing_percentage=0.0)

    # For YAML wings, n_unrefined_sections will be set when sections are added
    # Set to 0 as placeholder for now
    n_unrefined_sections_value::Int16 =
        isnothing(n_unrefined_sections) ? Int16(0) : Int16(n_unrefined_sections)

    panel_props::PanelProperties{n_panels} = PanelProperties{n_panels}()
    spanwise_direction_m::MVec3 = MVec3(spanwise_direction)

    # Initialize with default/empty values for optional fields
    Wing{n_panels}(
        Int16(n_panels), n_unrefined_sections_value, spanwise_distribution, panel_props, spanwise_direction_m,
        Section[], Section[], remove_nan, use_prior_polar, Float64(billowing_percentage),
        # Grouping
        Int16[],
        # Deformation fields
        Section[], zeros(max(0, n_panels)), zeros(max(0, n_panels)),
        # Physical properties (defaults for non-OBJ wings)
        0.0, 0.0, zeros(0, 0), zeros(MVec3), MMat3(I),
        0.0, nothing, nothing, nothing,
        PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}[]
    )
end

"""
    reinit!(wing::AbstractWing)

Reinitialize wing panel properties based on current refined_sections geometry.

This function only updates panel properties (chord, area, etc.) from the existing
refined_sections. It does NOT refine the mesh - call refine_aerodynamic_mesh!(wing)
first if needed.

# Note
After deformation via `unrefined_deform!()` or `deform!()`, call `reinit!` to update
panel properties while preserving the deformed geometry.
"""
function reinit!(wing::AbstractWing)
    # Calculate panel properties from refined sections
    update_panel_properties!(
        wing.panel_props,
        wing.refined_sections,
        wing.n_panels
    )
    return nothing
end

"""
    unrefined_deform!(wing::Wing, theta_angles=nothing, delta_angles=nothing; smooth=false)

Apply deformation angles directly to unrefined wing sections.

For wings that support deformation (OBJ-based wings with non_deformed_sections), this
applies theta_angles and delta_angles directly to unrefined sections and then applies deformation.
For wings without deformation support (YAML-based), this is a no-op that only succeeds
if both angle inputs are nothing.

# Arguments
- `wing::Wing`: The wing to deform
- `theta_angles::AbstractVector`: Twist angles in radians for each unrefined section (or nothing).
  Length must be `n_unrefined_sections`
- `delta_angles::AbstractVector`: Trailing edge deflection angles in radians for each unrefined section (or nothing).
  Length must be `n_unrefined_sections`
- `smooth::Bool`: DEPRECATED - no longer used. Apply smoothing to input angles if needed.

# Algorithm
1. Copies theta_angles and delta_angles directly to wing.theta_dist and wing.delta_dist
2. Calls deform! to update wing geometry and propagate to refined sections

# Errors
- Throws `ArgumentError` if wing doesn't support deformation but angles are provided
- Throws `ArgumentError` if angle vectors don't match n_unrefined_sections

# Returns
- `nothing` (modifies wing in-place)
"""
function unrefined_deform!(wing::Wing, theta_angles=nothing, delta_angles=nothing;
                           smooth=false, smooth_window=nothing)
    # Check if deformation is supported
    can_deform = !isempty(wing.non_deformed_sections)

    # If no deformation requested, just return
    isnothing(theta_angles) && isnothing(delta_angles) && return nothing

    # If deformation requested but not supported, throw error
    if !can_deform
        throw(ArgumentError("This Wing does not support deformation. Only OBJ-based wings created with ObjWing() can be deformed."))
    end

    # Validate inputs
    !isnothing(theta_angles) && length(theta_angles) != wing.n_unrefined_sections &&
        throw(ArgumentError("theta_angles must have length n_unrefined_sections = $(wing.n_unrefined_sections), got $(length(theta_angles))"))
    !isnothing(delta_angles) && length(delta_angles) != wing.n_unrefined_sections &&
        throw(ArgumentError("delta_angles must have length n_unrefined_sections = $(wing.n_unrefined_sections), got $(length(delta_angles))"))

    # Map unrefined sections → panels → sections (no smoothing yet)
    if !isnothing(theta_angles)
        map_unrefined_to_sections!(wing.theta_dist, theta_angles, wing.refined_panel_mapping, wing.n_panels)
    end
    if !isnothing(delta_angles)
        map_unrefined_to_sections!(wing.delta_dist, delta_angles, wing.refined_panel_mapping, wing.n_panels)
    end

    # Apply deformation with optional smoothing
    deform!(wing, smooth=smooth, smooth_window=smooth_window)
    return nothing
end

"""
    map_unrefined_to_sections!(panel_dist, unrefined_angles, panel_mapping, n_panels)

Map angles from unrefined sections to panels.
Steps: unrefined[1:n_unrefined] → panels[1:n_panels]

# Arguments
- `panel_dist::Vector{Float64}`: Output panel angles (length n_panels)
- `unrefined_angles::Vector{Float64}`: Input unrefined section angles
- `panel_mapping::Vector{Int16}`: Maps panel index to unrefined section index
- `n_panels::Int`: Number of panels
"""
function map_unrefined_to_sections!(panel_dist,
                                     unrefined_angles,
                                     panel_mapping,
                                     n_panels)
    # Map unrefined sections to panels
    for i in 1:n_panels
        unrefined_idx = panel_mapping[i]
        panel_dist[i] = unrefined_angles[unrefined_idx]
    end

    return nothing
end

"""
    smooth_distribution!(dist, window_size)

Apply moving average smoothing to a distribution in-place.
Uses a centered window of size `window_size` (must be odd).

# Arguments
- `dist::Vector{Float64}`: Distribution to smooth (modified in-place)
- `window_size::Int`: Size of smoothing window (must be odd)
"""
function smooth_distribution!(dist::Vector{Float64}, window_size::Int)
    n = length(dist)
    window_size <= 1 && return nothing
    n <= window_size && return nothing

    # Create temporary copy
    dist_copy = copy(dist)
    half_window = div(window_size, 2)

    # Apply moving average
    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        dist[i] = sum(dist_copy[start_idx:end_idx]) / (end_idx - start_idx + 1)
    end

    return nothing
end

"""
    deform!(wing::Wing, theta_dist::AbstractVector, delta_dist::AbstractVector;
            smooth=false, smooth_window=nothing)

Deform wing by applying theta and delta distributions at the panel level.

# Arguments
- `wing::Wing`: Wing to deform (must support deformation)
- `theta_dist::AbstractVector`: Twist angles for each panel (length = n_panels)
- `delta_dist::AbstractVector`: TE deflections for each panel (length = n_panels)
- `smooth::Bool`: Whether to apply smoothing (default: false)
- `smooth_window::Union{Nothing, Int}`: Smoothing window size (default: auto-calculated)

# Effects
Updates wing.refined_sections with deformed geometry based on wing.non_deformed_sections
"""
function deform!(wing::Wing, theta_dist::AbstractVector, delta_dist::AbstractVector;
                 smooth=false, smooth_window=nothing)
    !isempty(wing.non_deformed_sections) || throw(ArgumentError("Wing does not support deformation"))

    expected_len = wing.n_panels
    !(length(theta_dist) == expected_len) && throw(ArgumentError("theta_dist must have length $(expected_len), got $(length(theta_dist))"))
    !(length(delta_dist) == expected_len) && throw(ArgumentError("delta_dist must have length $(expected_len), got $(length(delta_dist))"))

    wing.theta_dist .= theta_dist
    wing.delta_dist .= delta_dist

    deform!(wing, smooth=smooth, smooth_window=smooth_window)
end

"""
    deform!(wing::Wing; smooth=false, smooth_window=nothing)

Apply stored theta_dist and delta_dist to deform the wing geometry.
Converts panel angles (n_panels) to section angles (n_panels+1) by averaging adjacent panels.

# Arguments
- `wing::Wing`: Wing to deform (must have non_deformed_sections)
- `smooth::Bool`: Whether to apply smoothing to theta_dist and delta_dist (default: false)
- `smooth_window::Union{Nothing, Int}`: Smoothing window size (default: auto-calculated)

# Effects
Updates wing.refined_sections based on wing.non_deformed_sections and stored distributions
"""
function deform!(wing::Wing; smooth=false, smooth_window=nothing)
    !isempty(wing.non_deformed_sections) || return nothing

    # Apply smoothing if requested
    if smooth
        # Default window size based on refinement ratio
        if isnothing(smooth_window)
            smooth_window = max(3, round(Int, wing.n_panels / wing.n_unrefined_sections))
        end

        # Ensure window is odd for symmetric averaging
        if smooth_window % 2 == 0
            smooth_window += 1
        end

        # Apply moving average to theta_dist and delta_dist (panel-level)
        smooth_distribution!(wing.theta_dist, smooth_window)
        smooth_distribution!(wing.delta_dist, smooth_window)
    end

    local_y = zeros(MVec3)
    chord = zeros(MVec3)
    normal = zeros(MVec3)

    # Process all refined sections (n_panels + 1)
    # Convert panel angles to section angles by averaging
    for i in 1:(wing.n_panels + 1)
        # Determine theta for this section by averaging adjacent panels
        if i == 1
            # First section: use panel 1 angle
            theta = wing.theta_dist[1]
        elseif i == wing.n_panels + 1
            # Last section: use last panel angle
            theta = wing.theta_dist[wing.n_panels]
        else
            # Middle sections: average of panels (i-1) and i
            theta = (wing.theta_dist[i-1] + wing.theta_dist[i]) / 2.0
        end

        section = wing.non_deformed_sections[i]

        # Compute local coordinate system
        if i < wing.n_panels + 1
            section2 = wing.non_deformed_sections[i+1]
            local_y .= normalize(section.LE_point - section2.LE_point)
        else
            # For last section, use same local_y as previous
            section_prev = wing.non_deformed_sections[i-1]
            local_y .= normalize(section_prev.LE_point - section.LE_point)
        end

        chord .= section.TE_point .- section.LE_point
        normal .= chord × local_y
        @. wing.refined_sections[i].TE_point = section.LE_point +
            cos(theta) * chord - sin(theta) * normal
    end
    return nothing
end


"""
    remove_vector_nans(aero_data)

Remove the indices from aero_data where a NaN is found.
"""
function remove_vector_nans(aero_data::NTuple{4,AbstractVector})
    alpha_range, cl_vector, cd_vector, cm_vector = aero_data
    alpha_range = collect(alpha_range)
    nan_indices = Set{Int}()
    for vec in (cl_vector, cd_vector, cm_vector)
        union!(nan_indices, findall(isnan, vec))
    end
    if isempty(nan_indices)
        return aero_data
    end
    # Convert to sorted array for consistent removal
    nan_indices = sort(collect(nan_indices))
    # Create mask for valid indices
    valid_mask = trues(length(alpha_range))
    valid_mask[nan_indices] .= false
    # Remove NaN values from all vectors
    clean_alpha = alpha_range[valid_mask]
    clean_cl = cl_vector[valid_mask]
    clean_cd = cd_vector[valid_mask]
    clean_cm = cm_vector[valid_mask]
    @info "Removed $(length(nan_indices)) indices containing NaNs from aero_data."
    return (clean_alpha, clean_cl, clean_cd, clean_cm)
end

@inline remove_vector_nans(::Nothing) = nothing
@inline remove_vector_nans(::Tuple{Float64, Float64}) =
    throw(ArgumentError("POLAR_VECTORS requires aero_data = (alpha, cl, cd, cm) vectors."))
@inline remove_vector_nans(::Tuple{AbstractVector, AbstractVector, AbstractMatrix, AbstractMatrix, AbstractMatrix}) =
    throw(ArgumentError("POLAR_VECTORS requires aero_data = (alpha, cl, cd, cm) vectors."))

function interpolate_polar_matrix_nans!(
    aero_data::Tuple{AbstractVector, AbstractVector, AbstractMatrix, AbstractMatrix, AbstractMatrix}
)
    interpolate_matrix_nans!.(aero_data[3:5])
    return nothing
end

@inline interpolate_polar_matrix_nans!(::Nothing) = nothing
@inline interpolate_polar_matrix_nans!(::Tuple{Float64, Float64}) =
    throw(ArgumentError("POLAR_MATRICES requires aero_data = (alpha, delta, cl, cd, cm)."))
@inline interpolate_polar_matrix_nans!(::NTuple{4,AbstractVector}) =
    throw(ArgumentError("POLAR_MATRICES requires aero_data = (alpha, delta, cl, cd, cm)."))

"""
    add_section!(wing::Wing, LE_point::PosVector, TE_point::PosVector, 
                 aero_model, aero_data::AeroData=nothing)

Add a new section to the wing.

# Arguments:
- wing::Wing: The [Wing](@ref) to which a section shall be added
- LE_point::PosVector: [PosVector](@ref) of the point on the side of the leading edge
- TE_point::PosVector: [PosVector](@ref) of the point on the side of the trailing edge
- `aero_model`::AeroModel: [AeroModel](@ref)
- `aero_data`::AeroData: See [AeroData](@ref)  
"""
function add_section!(wing::Wing, LE_point,
                     TE_point, aero_model::AeroModel, aero_data::AeroData=nothing)
    if aero_model == POLAR_VECTORS && wing.remove_nan
        aero_data = remove_vector_nans(aero_data)
    elseif aero_model == POLAR_MATRICES && wing.remove_nan
        interpolate_polar_matrix_nans!(aero_data)
    end
    push!(wing.unrefined_sections, Section(LE_point, TE_point, aero_model, aero_data))
    wing.n_unrefined_sections = Int16(length(wing.unrefined_sections))
    return nothing
end

"""
    flip_created_coord_in_pairs_if_needed!(coord::Matrix{Float64})

Ensure coordinates are ordered from positive to negative along y-axis.
"""
function flip_created_coord_in_pairs_if_needed!(coord::Matrix{Float64})
    n_pairs = size(coord, 1) ÷ 2
    reshaped = reshape(coord, (n_pairs, 2, size(coord, 2)))
    
    # Check y-values order
    y_values = reshaped[:, 1, 2]  # y-coordinates of leading edge points
    if !all(y_values[1:end-1] .>= y_values[2:end])
        reshaped = reverse(reshaped, dims=1)
    end
    
    return reshape(reshaped, (size(coord, 1), size(coord, 2)))
end


"""
    update_non_deformed_sections!(wing::AbstractWing)

Create non_deformed_sections to match refined_sections.
This enables deformation support for all wings (YAML and OBJ).
Should be called after refined_sections are populated.
Once populated, non_deformed_sections serves as the undeformed reference geometry.
"""
function update_non_deformed_sections!(wing::AbstractWing)
    n_sections = wing.n_panels + 1

    # Populate or update non_deformed_sections
    if isempty(wing.non_deformed_sections)
        # Initial setup
        wing.non_deformed_sections = [Section() for _ in 1:n_sections]
        for i in 1:n_sections
            reinit!(wing.non_deformed_sections[i], wing.refined_sections[i])
        end
    elseif length(wing.non_deformed_sections) != n_sections
        # Size mismatch - error
        throw(ArgumentError(
            "non_deformed_sections has incorrect size. " *
            "Expected $(n_sections) sections (n_panels+1), got $(length(wing.non_deformed_sections)). " *
            "This indicates an inconsistent wing state."
        ))
    else
        # Correct size - update all sections
        for i in 1:n_sections
            reinit!(wing.non_deformed_sections[i], wing.refined_sections[i])
        end
    end
    return nothing
end

@inline function _has_initialized_section_aero_data(section::Section)
    section.aero_model == INVISCID && return true
    return !isnothing(section.aero_data)
end

@inline function _can_reuse_prior_refined_polar_data(wing::AbstractWing, n_sections::Int)
    wing.use_prior_polar || return false
    length(wing.refined_sections) == n_sections || return false
    return all(_has_initialized_section_aero_data, wing.refined_sections)
end

"""
    copy_sections_to_refined!(wing; reuse_aero_data=false)

Copy unrefined sections to refined sections 1:1 (no interpolation).
If `refined_sections` is empty, allocates via `copy`; otherwise
reinitialises in-place.  Warns if billowing was requested but there
are no intermediate sections to billow.
"""
function copy_sections_to_refined!(
    wing::AbstractWing; reuse_aero_data::Bool=false
)
    if wing.spanwise_distribution == BILLOWING &&
            wing.billowing_percentage > 0
        @warn "Billowing requested but n_panels " *
            "($(wing.n_panels)) == n_provided; no " *
            "intermediate sections to billow. " *
            "Increase n_panels."
    end
    if length(wing.refined_sections) == 0
        wing.refined_sections = copy(wing.unrefined_sections)
    else
        for (refined, unrefined) in zip(
                wing.refined_sections, wing.unrefined_sections)
            if reuse_aero_data
                reinit!(refined, unrefined.LE_point,
                    unrefined.TE_point, unrefined.aero_model)
            else
                reinit!(refined, unrefined)
            end
        end
    end
    return nothing
end

"""
    refine!(wing::AbstractWing; recompute_mapping=true, sort_sections=true)

Refine the wing aerodynamic mesh from unrefined sections to refined sections.

This function interpolates the wing geometry from a coarse set of unrefined sections
to a fine mesh of refined sections (n_panels+1 sections) based on the wing's
spanwise_distribution setting. It also populates non_deformed_sections which
enables deformation support via `unrefined_deform!`.

# Required Workflow
Must be called after wing construction and before creating `BodyAerodynamics`:
```julia
wing = Wing("wing.yaml"; n_panels=40)  # or ObjWing(...) or manual Wing
refine!(wing)                          # Refine mesh
body_aero = BodyAerodynamics([wing])   # Create aerodynamics
```

# Distribution Methods
- `LINEAR`: Linear interpolation between sections
- `COSINE`: Cosine spacing (more panels near tips)
- `SPLIT_PROVIDED`: Split each unrefined section into sub-panels
- `UNCHANGED`: 1:1 copy when n_unrefined_sections == n_panels+1

# Keyword Arguments
- `recompute_mapping::Bool=true`: Recompute the mapping from refined panels to unrefined sections
- `sort_sections::Bool=true`: Sort sections by spanwise position using global `LE_point[2]` (Y-axis). Disable for structural ordering.

# Effects
1. Populates `wing.refined_sections` (n_panels+1 sections)
2. Populates `wing.non_deformed_sections` (copy of refined_sections for deformation reference)
3. Computes `wing.refined_panel_mapping` (panel → unrefined section mapping)
4. Resizes `wing.theta_dist` and `wing.delta_dist` to n_panels

# Example
```julia
# YAML wing
wing = Wing("wing.yaml"; n_panels=40)
refine!(wing)
body_aero = BodyAerodynamics([wing])

# After refinement, deformation is supported
unrefined_deform!(wing, theta_angles, delta_angles)
```
"""
function refine!(wing::AbstractWing; recompute_mapping=true, sort_sections=true)
    # Validate unrefined_sections exist
    if isempty(wing.unrefined_sections)
        throw(ArgumentError(
            "Cannot refine mesh: wing has no unrefined_sections. " *
            "Add sections using add_section! or check wing construction."
        ))
    end

    # Only sort sections if requested (skip for REFINE wings with fixed structural order)
    #TODO: only works if can be sorted by global y position
    if sort_sections
        sorted = true
        for i in 1:(length(wing.unrefined_sections) - 1)
            if wing.unrefined_sections[i].LE_point[2] <
               wing.unrefined_sections[i+1].LE_point[2]
                sorted = false
                break
            end
        end
        sorted || sort!(wing.unrefined_sections;
            by=s -> s.LE_point[2], rev=true)
    end
    n_sections = wing.n_panels + 1
    reuse_aero_data = _can_reuse_prior_refined_polar_data(wing, n_sections)

    if length(wing.refined_sections) == 0
        if wing.spanwise_distribution == UNCHANGED ||
               length(wing.unrefined_sections) == n_sections
            copy_sections_to_refined!(wing; reuse_aero_data)
            recompute_mapping && compute_refined_panel_mapping!(wing)
            update_non_deformed_sections!(wing)
            return nothing
        else
            wing.refined_sections = Section[Section() for _ in 1:wing.n_panels+1]
        end
    end
    
    # Handle special cases
    if wing.spanwise_distribution == UNCHANGED || length(wing.unrefined_sections) == n_sections
        copy_sections_to_refined!(wing; reuse_aero_data)
        recompute_mapping && compute_refined_panel_mapping!(wing)
        update_non_deformed_sections!(wing)
        return nothing
    end

    @debug "Refining aerodynamic mesh from $(length(wing.unrefined_sections)) sections to $n_sections sections."

    # Handle two-section case
    if n_sections == 2
        s1 = wing.unrefined_sections[1]
        s2 = wing.unrefined_sections[end]
        reinit!(wing.refined_sections[1],
            s1.LE_point, s1.TE_point, s1.aero_model,
            reuse_aero_data ? nothing : s1.aero_data)
        reinit!(wing.refined_sections[2],
            s2.LE_point, s2.TE_point, s2.aero_model,
            reuse_aero_data ? nothing : s2.aero_data)
        recompute_mapping && compute_refined_panel_mapping!(wing)
        update_non_deformed_sections!(wing)
        return nothing
    end

    # Handle different distribution types
    if wing.spanwise_distribution == SPLIT_PROVIDED
        refine_mesh_by_splitting_provided_sections!(wing; reuse_aero_data)
    elseif wing.spanwise_distribution in (LINEAR, COSINE)
        refine_mesh_for_linear_cosine_distribution!(
            wing, 1, wing.spanwise_distribution,
            n_sections, wing.unrefined_sections;
            reuse_aero_data)
    elseif wing.spanwise_distribution == BILLOWING
        refine_mesh_with_billowing!(wing; reuse_aero_data)
    else
        throw(ArgumentError("Unsupported spanwise panel distribution: $(wing.spanwise_distribution)"))
    end

    # Compute panel mapping by finding closest unrefined section for each refined panel
    recompute_mapping && compute_refined_panel_mapping!(wing)

    # Update n_unrefined_sections based on actual sections
    wing.n_unrefined_sections = Int16(length(wing.unrefined_sections))

    # Resize theta_dist and delta_dist to match n_panels
    target_size = wing.n_panels
    if length(wing.theta_dist) != target_size
        resize!(wing.theta_dist, target_size)
        fill!(wing.theta_dist, 0.0)
    end
    if length(wing.delta_dist) != target_size
        resize!(wing.delta_dist, target_size)
        fill!(wing.delta_dist, 0.0)
    end

    # Create/update non_deformed_sections to match refined_sections
    update_non_deformed_sections!(wing)

    return nothing
end


"""
    compute_refined_panel_mapping!(wing::AbstractWing)

Compute the mapping from refined panels to unrefined sections by finding
the closest unrefined section for each refined panel (based on section center distance).
Maps each refined panel index to its corresponding unrefined section index
(1 to n_unrefined_sections).
Works after refinement is complete.
"""
function compute_refined_panel_mapping!(wing::AbstractWing)
    n_unref = length(wing.unrefined_sections)
    n_panels = wing.n_panels

    # Ensure mapping array is correctly sized
    if length(wing.refined_panel_mapping) != n_panels
        resize!(wing.refined_panel_mapping, n_panels)
    end

    # Handle case where no refinement occurred
    if n_unref == n_panels + 1
        for i in 1:n_panels
            wing.refined_panel_mapping[i] = Int16(i)
        end
        return nothing
    end

    # For each refined panel, find closest unrefined section
    # using scalar arithmetic to avoid MVec3 allocations
    for pi in 1:n_panels
        r1 = wing.refined_sections[pi]
        r2 = wing.refined_sections[pi + 1]
        # refined panel center (scalar)
        rc1 = (r1.LE_point[1]+r1.TE_point[1]+
               r2.LE_point[1]+r2.TE_point[1]) * 0.25
        rc2 = (r1.LE_point[2]+r1.TE_point[2]+
               r2.LE_point[2]+r2.TE_point[2]) * 0.25
        rc3 = (r1.LE_point[3]+r1.TE_point[3]+
               r2.LE_point[3]+r2.TE_point[3]) * 0.25

        min_dist = Inf
        closest = Int16(1)
        for ui in 1:n_unref
            u = wing.unrefined_sections[ui]
            uc1 = (u.LE_point[1]+u.TE_point[1]) * 0.5
            uc2 = (u.LE_point[2]+u.TE_point[2]) * 0.5
            uc3 = (u.LE_point[3]+u.TE_point[3]) * 0.5
            d1 = rc1-uc1; d2 = rc2-uc2; d3 = rc3-uc3
            dist = d1*d1 + d2*d2 + d3*d3
            if dist < min_dist
                min_dist = dist
                closest = Int16(ui)
            end
        end
        wing.refined_panel_mapping[pi] = closest
    end

    return nothing
end


"""
    calculate_new_aero_data(aero_model,
                            aero_data, 
                            section_index::Int,
                            left_weight::Float64,
                            right_weight::Float64)

Interpolate aerodynamic input between two sections.
"""
function calculate_new_aero_data(
    sections::AbstractVector{Section},
    section_index::Int,
    left_weight::Float64,
    right_weight::Float64
)
    s_l = sections[section_index]
    s_r = sections[section_index + 1]
    aero_model = (s_l.aero_model, s_r.aero_model)
    aero_data = (s_l.aero_data, s_r.aero_data)
    return calculate_new_aero_data(
        aero_model, aero_data,
        1, left_weight, right_weight)
end

function calculate_new_aero_data(aero_model,
                                aero_data,
                                section_index::Int,
                                left_weight::Float64,
                                right_weight::Float64)
    
    model_type = aero_model[section_index]
    model_type_2 = aero_model[section_index+1]
    if !(model_type == model_type_2)
        throw(ArgumentError("Different aero models over the span are not supported"))
    end
    
    if model_type == INVISCID
        return nothing
        
    elseif model_type in (POLAR_VECTORS, POLAR_MATRICES)
        polar_left = aero_data[section_index]
        polar_right = aero_data[section_index + 1]
        
        # Unpack polar data
        if model_type == POLAR_VECTORS
            alpha_left, CL_left, CD_left, CM_left = polar_left
            alpha_right, CL_right, CD_right, CM_right = polar_right

            (
                length(alpha_left) == length(alpha_right) &&
                all(isapprox.(diff(alpha_left), diff(alpha_right)))
            ) || throw(ArgumentError("Alpha steps must be identical."))
            isa(CL_right, AbstractVector) || throw(ArgumentError(
                "Provide polar data in the correct format."
            ))
            
            # Weighted interpolation
            CL_data = CL_left .* left_weight .+ CL_right .* right_weight
            CD_data = CD_left .* left_weight .+ CD_right .* right_weight
            CM_data = CM_left .* left_weight .+ CM_right .* right_weight
            
            return (alpha_left, CL_data, CD_data, CM_data)
            
        elseif model_type == POLAR_MATRICES
            alpha_left, delta_left, CL_left, CD_left, CM_left = polar_left
            alpha_right, delta_right, CL_right, CD_right, CM_right = polar_right
            
            (
                length(alpha_left) == length(alpha_right) &&
                all(isapprox.(diff(alpha_left), diff(alpha_right)))
            ) || throw(ArgumentError("Alpha steps must be identical."))
            (
                length(delta_left) == length(delta_right) &&
                all(isapprox.(diff(delta_left), diff(delta_right)))
            ) || throw(ArgumentError("Delta steps must be identical."))
            isa(CL_right, AbstractMatrix) || throw(ArgumentError(
                "Provide polar data in the correct format."
            ))

            # Weighted interpolation
            CL_data = CL_left .* left_weight .+ CL_right .* right_weight
            CD_data = CD_left .* left_weight .+ CD_right .* right_weight
            CM_data = CM_left .* left_weight .+ CM_right .* right_weight
            
            return (alpha_left, delta_left, CL_data, CD_data, CM_data)
        end

    elseif model_type == LEI_AIRFOIL_BREUKELS
        tube_diameter_left = aero_data[section_index][1]
        tube_diameter_right = aero_data[section_index + 1][1]
        tube_diameter_i = tube_diameter_left * left_weight + tube_diameter_right * right_weight
        
        chamber_height_left = aero_data[section_index][2]
        chamber_height_right = aero_data[section_index + 1][2]
        chamber_height_i = chamber_height_left * left_weight + chamber_height_right * right_weight
        
        @debug "Interpolation weights" left_weight right_weight
        @debug "Interpolated parameters" tube_diameter_i chamber_height_i
        
        return (tube_diameter_i, chamber_height_i)
    else
        throw(ArgumentError("Unsupported aero model: $(model_type)"))
    end
end

"""
    refine_mesh_for_linear_cosine_distribution!(wing, idx, dist,
        n_sections, sections; reuse_aero_data)

Zero-allocation variant that reads directly from a `Vector{Section}`.
"""
function refine_mesh_for_linear_cosine_distribution!(
    wing::AbstractWing,
    idx,
    spanwise_distribution::PanelDistribution,
    n_sections::Int,
    sections::AbstractVector{Section};
    endpoints::Bool=true,
    reuse_aero_data::Bool=false)

    n_input = length(sections)

    @inline _le(s, j) = @inbounds s.LE_point[j]
    @inline _te(s, j) = @inbounds s.TE_point[j]

    # Compute total quarter-chord length (scalar only)
    qc_total = 0.0
    for i in 1:(n_input - 1)
        d = 0.0
        s_i = sections[i]; s_ip = sections[i+1]
        for j in 1:3
            qc_j = (_le(s_ip,j) + 0.25*(_te(s_ip,j) - _le(s_ip,j))) -
                    (_le(s_i,j) + 0.25*(_te(s_i,j) - _le(s_i,j)))
            d += qc_j * qc_j
        end
        qc_total += sqrt(d)
    end

    new_le = MVec3(0.0, 0.0, 0.0)
    new_te = MVec3(0.0, 0.0, 0.0)
    dir = MVec3(0.0, 0.0, 0.0)

    for i in 1:n_sections
        target = if spanwise_distribution == LINEAR
            qc_total * (i - 1) / (n_sections - 1)
        else
            qc_total * (1-cos(π*(i-1)/(n_sections-1))) / 2
        end

        # Linear scan to find segment
        cum = 0.0
        si = 1
        for k in 1:(n_input - 1)
            d = 0.0
            s_k = sections[k]; s_kp = sections[k+1]
            for j in 1:3
                qc_j = (_le(s_kp,j) + 0.25*(_te(s_kp,j)-_le(s_kp,j))) -
                        (_le(s_k,j) + 0.25*(_te(s_k,j)-_le(s_k,j)))
                d += qc_j * qc_j
            end
            next_cum = cum + sqrt(d)
            if next_cum >= target || k == n_input - 1
                si = k
                break
            end
            cum = next_cum
        end
        s_l = sections[si]; s_r = sections[si + 1]

        # Segment length for interpolation weight
        seg_d = 0.0
        for j in 1:3
            qc_j = (_le(s_r,j) + 0.25*(_te(s_r,j)-_le(s_r,j))) -
                    (_le(s_l,j) + 0.25*(_te(s_l,j)-_le(s_l,j)))
            seg_d += qc_j * qc_j
        end
        seg_len = sqrt(seg_d)
        t = seg_len > 1e-30 ? (target - cum) / seg_len : 0.0
        wl = 1 - t; wr = t

        lc_len = 0.0; rc_len = 0.0
        for j in 1:3
            lc_len += (_te(s_l,j) - _le(s_l,j))^2
            rc_len += (_te(s_r,j) - _le(s_r,j))^2
        end
        lc_len = sqrt(lc_len); rc_len = sqrt(rc_len)

        for j in 1:3
            dir[j] = wl * (_te(s_l,j)-_le(s_l,j)) / max(lc_len, 1e-12) +
                     wr * (_te(s_r,j)-_le(s_r,j)) / max(rc_len, 1e-12)
        end
        dir_len = sqrt(dir[1]^2+dir[2]^2+dir[3]^2)
        avg_len = wl * lc_len + wr * rc_len

        for j in 1:3
            qc_j = _le(s_l,j) + 0.25*(_te(s_l,j)-_le(s_l,j)) +
                t * ((_le(s_r,j) + 0.25*(_te(s_r,j)-_le(s_r,j))) -
                     (_le(s_l,j) + 0.25*(_te(s_l,j)-_le(s_l,j))))
            chord_j = dir[j] / max(dir_len, 1e-12) * avg_len
            new_le[j] = qc_j - 0.25 * chord_j
            new_te[j] = qc_j + 0.75 * chord_j
        end

        new_data = reuse_aero_data ? nothing :
            calculate_new_aero_data(sections, si, wl, wr)

        if endpoints || (i != 1 && i != n_sections)
            reinit!(wing.refined_sections[idx],
                new_le, new_te, s_l.aero_model, new_data)
            idx += 1
        end
    end

    return idx
end

"""
    refine_mesh_for_linear_cosine_distribution!(
        wing::AbstractWing,
        idx::Int,
        spanwise_distribution::PanelDistribution,
        n_sections::Int,
        LE::Matrix{Float64},
        TE::Matrix{Float64},
        aero_model,
        aero_data)

Refine wing mesh using linear or cosine spacing.

# Arguments
- `wing`: Wing object
- `idx`: Section start index
- `spanwise_distribution`: [PanelDistribution](@ref)
- `n_sections`: Number of sections to generate
- `LE`: Matrix of leading edge points
- `TE`: Matrix of trailing edge points
- `aero_model`: Vector of aerodynamic models for each section
- `aero_data`: Vector of aerodynamic data for each section

# Keyword arguments
- endpoints=true

Returns:
    idx: Last section index
"""
function refine_mesh_for_linear_cosine_distribution!(
    wing::AbstractWing,
    idx,
    spanwise_distribution::PanelDistribution,
    n_sections::Int,
    LE,
    TE,
    aero_model,
    aero_data;
    endpoints=true,
    reuse_aero_data::Bool=false)

    # 1. Compute quarter chord line
    quarter_chord = LE .+ 0.25 .* (TE .- LE)

    # Calculate segment lengths
    qc_lengths = [norm(quarter_chord[i+1,:] - quarter_chord[i,:]) for i in 1:size(quarter_chord,1)-1]
    qc_total_length = sum(qc_lengths)
    qc_cum_length = vcat(0, cumsum(qc_lengths))

    # 2. Define target lengths
    target_lengths = if spanwise_distribution == LINEAR
        range(0, qc_total_length, n_sections)
    elseif spanwise_distribution == COSINE
        theta = range(0, π, n_sections)
        qc_total_length .* (1 .- cos.(theta)) ./ 2
    else
        throw(ArgumentError("Unsupported distribution: $spanwise_distribution"))
    end

    # Initialize arrays
    new_quarter_chord = zeros(Float64, n_sections, 3)
    new_LE = zeros(Float64, n_sections, 3)
    new_TE = zeros(Float64, n_sections, 3)
    # 3. Calculate new points and interpolate
    for i in 1:n_sections
        target_length = target_lengths[i]

        # Find segment index
        section_index = searchsortedlast(qc_cum_length, target_length)
        section_index = clamp(section_index, 1, length(qc_cum_length)-1)

        # 4. Calculate weights
        segment_start = qc_cum_length[section_index]
        segment_end = qc_cum_length[section_index+1]
        t = (target_length - segment_start) / (segment_end - segment_start)
        left_weight = 1 - t
        right_weight = t

        # 5. Calculate quarter chord point
        new_quarter_chord[i,:] = quarter_chord[section_index,:] +
                                t .* (quarter_chord[section_index+1,:] - quarter_chord[section_index,:])

        # 6. Calculate chord vectors
        left_chord = TE[section_index,:] - LE[section_index,:]
        right_chord = TE[section_index+1,:] - LE[section_index+1,:]

        # Normalize chord vectors
        left_chord_norm = left_chord ./ max(norm(left_chord), 1e-12)
        right_chord_norm = right_chord ./ max(norm(right_chord), 1e-12)

        # Interpolate direction
        avg_direction = left_weight .* left_chord_norm .+ right_weight .* right_chord_norm
        avg_direction = avg_direction ./ max(norm(avg_direction), 1e-12)

        # Interpolate length
        left_length = norm(left_chord)
        right_length = norm(right_chord)
        avg_length = left_weight * left_length + right_weight * right_length

        # Final chord vector
        avg_chord = avg_direction .* avg_length

        # Calculate LE and TE points
        new_LE[i,:] = new_quarter_chord[i,:] .- 0.25 .* avg_chord
        new_TE[i,:] = new_quarter_chord[i,:] .+ 0.75 .* avg_chord

        # Interpolate aero properties unless reusing prior refined section aero.
        new_data = reuse_aero_data ? nothing :
            calculate_new_aero_data(aero_model, aero_data, section_index, left_weight, right_weight)

        # Create new section
        if endpoints || (i != 1 && i != n_sections)
            @views reinit!(wing.refined_sections[idx], new_LE[i,:], new_TE[i,:], aero_model[1], new_data)
            idx += 1
        end
    end

    return idx
end

"""
    refine_mesh_by_splitting_provided_sections!(wing; reuse_aero_data,
        billowing_percentage)

Refine mesh by splitting provided sections into desired number of panels.

When `billowing_percentage > 0`, rotates chord vectors around the leading
edge with a sinusoidal profile to simulate fabric billowing between ribs.
"""
function refine_mesh_by_splitting_provided_sections!(
    wing::AbstractWing;
    reuse_aero_data::Bool=false,
    billowing_percentage::Float64=0.0
)
    n_sections_provided = length(wing.unrefined_sections)
    n_panels_provided = n_sections_provided - 1
    n_panels_desired = wing.n_panels
    
    @debug "Panel counts" n_panels_provided n_panels_desired n_sections_provided
    
    # Check if refinement is needed
    if n_panels_provided == n_panels_desired
        copy_sections_to_refined!(wing; reuse_aero_data)
        return nothing
    end

    # Validate panel count relationship
    if n_panels_desired % n_panels_provided != 0
        throw(ArgumentError(
            "Desired panels ($n_panels_desired) must be multiple of provided panels " *
            "($n_panels_provided). Choose: $(n_panels_provided*2), $(n_panels_provided*3), ..."
        ))
    end
    
    # Calculate distribution
    n_new_sections = wing.n_panels + 1 - n_sections_provided
    n_section_pairs = n_sections_provided - 1
    new_sections_per_pair, remaining = divrem(n_new_sections, n_section_pairs)
    
    # Process each section pair (no extraction — read sections directly)
    sections = wing.unrefined_sections
    idx = 1
    for li in 1:n_section_pairs
        s_l = sections[li]; s_r = sections[li + 1]

        # Add left section of pair
        if reuse_aero_data
            reinit!(wing.refined_sections[idx],
                s_l.LE_point, s_l.TE_point, s_l.aero_model)
        else
            reinit!(wing.refined_sections[idx], s_l)
        end
        idx += 1

        # Calculate new sections for this pair
        num_new = new_sections_per_pair +
            (li <= remaining ? 1 : 0)

        if num_new > 0
            start_idx = idx
            pair_sections = @view sections[li:li+1]
            idx = refine_mesh_for_linear_cosine_distribution!(
                wing, idx, LINEAR,
                num_new + 2,
                pair_sections;
                endpoints=false, reuse_aero_data)

            # Apply billowing by rotating chords around LE
            if billowing_percentage > 0 && idx > start_idx
                le_l = s_l.LE_point; le_r = s_r.LE_point
                diff = le_l - le_r
                span_len = norm(diff)
                y_hat = diff / span_len
                apply_billowing_to_pair!(
                    wing.refined_sections,
                    start_idx, idx - 1,
                    y_hat, span_len, le_r,
                    s_l.TE_point, s_r.TE_point,
                    billowing_percentage)
            end
        end
    end
    
    # Add final section
    if reuse_aero_data
        last_section = wing.unrefined_sections[end]
        reinit!(
            wing.refined_sections[idx],
            last_section.LE_point,
            last_section.TE_point,
            last_section.aero_model
        )
    else
        reinit!(wing.refined_sections[idx], wing.unrefined_sections[end])
    end
    
    # Validate result
    if length(wing.refined_sections) != wing.n_panels + 1
        @warn "Number of panels ($(length(wing.refined_sections)-1)) differs from desired ($(wing.n_panels))"
    end
    
    return nothing
end


"""
    rotated_te(le, te, y_hat, θ)

Compute the trailing-edge position after rotating the chord
vector `te - le` around `y_hat` by `θ` radians (Rodrigues).
Returns an `SVector{3}` (stack-allocated, zero heap allocs).
"""
@inline function rotated_te(le, te, y_hat, θ)
    c1 = te[1] - le[1]; c2 = te[2] - le[2]; c3 = te[3] - le[3]
    ct = cos(θ); st = sin(θ)
    cx1 = y_hat[2] * c3 - y_hat[3] * c2
    cx2 = y_hat[3] * c1 - y_hat[1] * c3
    cx3 = y_hat[1] * c2 - y_hat[2] * c1
    d_y = c1 * y_hat[1] + c2 * y_hat[2] + c3 * y_hat[3]
    k = (1 - ct) * d_y
    SVector(le[1] + ct*c1 + st*cx1 + k*y_hat[1],
            le[2] + ct*c2 + st*cx2 + k*y_hat[2],
            le[3] + ct*c3 + st*cx3 + k*y_hat[3])
end

"""
    billowing_arc_length(sections, start_si, end_si, y_hat,
                         span_len, le_ref, te_left, te_right,
                         angle_max)

Compute the TE arc length that would result from rotating each
section's chord around `y_hat` by `angle_max * sin(π t)` radians,
where `t` is the normalised spanwise position within the rib pair.

Non-allocating (uses scalar and SVector arithmetic only).
"""
function billowing_arc_length(
    sections, start_si, end_si,
    y_hat, span_len, le_ref,
    te_left, te_right, angle_max
)
    p1 = te_left[1]; p2 = te_left[2]; p3 = te_left[3]
    arc = 0.0
    for si in start_si:end_si
        sec = sections[si]
        le = sec.LE_point; te = sec.TE_point
        t = ((le[1] - le_ref[1]) * y_hat[1] +
             (le[2] - le_ref[2]) * y_hat[2] +
             (le[3] - le_ref[3]) * y_hat[3]) / span_len
        θ = -angle_max * sin(π * t)
        cur = rotated_te(le, te, y_hat, θ)
        d1 = cur[1] - p1; d2 = cur[2] - p2; d3 = cur[3] - p3
        arc += sqrt(d1*d1 + d2*d2 + d3*d3)
        p1 = cur[1]; p2 = cur[2]; p3 = cur[3]
    end
    d1 = te_right[1]-p1; d2 = te_right[2]-p2; d3 = te_right[3]-p3
    arc += sqrt(d1*d1 + d2*d2 + d3*d3)
    return arc
end

"""
    apply_billowing_to_pair!(sections, start_si, end_si, y_hat,
                             span_len, le_ref, te_left, te_right,
                             percentage)

Apply billowing to refined sections between two ribs by rotating
each section's chord around `y_hat`.  Uses Newton iteration to
find the rotation amplitude that matches the target TE arc-length
`percentage`, then applies the rotations in-place.

The angle profile is `angle_max * sin(π t)` (zero at ribs,
maximum at centre).  All angles are in radians.

Non-allocating: modifies `sections[start_si:end_si].TE_point`
in-place using scalar arithmetic only.
"""
function apply_billowing_to_pair!(
    sections, start_si, end_si,
    y_hat, span_len, le_ref,
    te_left, te_right, percentage
)
    percentage <= 0 && return nothing
    if percentage >= 100
        throw(ArgumentError(
            "billowing_percentage must be < 100, got $percentage"))
    end
    d1 = te_right[1]-te_left[1]
    d2 = te_right[2]-te_left[2]
    d3 = te_right[3]-te_left[3]
    straight = sqrt(d1*d1 + d2*d2 + d3*d3)
    straight < 1e-12 && return nothing
    target_arc = straight / (1 - percentage / 100)

    # Newton iteration on angle_max (rad)
    angle_max = sqrt(4 * percentage / (100 - percentage))
    for _ in 1:50
        arc = billowing_arc_length(
            sections, start_si, end_si,
            y_hat, span_len, le_ref,
            te_left, te_right, angle_max)
        f = arc - target_arc
        abs(f) < 1e-12 * target_arc && break

        δ = max(1e-8, abs(angle_max) * 1e-8)
        arc_p = billowing_arc_length(
            sections, start_si, end_si,
            y_hat, span_len, le_ref,
            te_left, te_right, angle_max + δ)
        df = (arc_p - arc) / δ
        abs(df) < 1e-30 && break
        angle_max = max(0.0, angle_max - f / df)
    end

    # Apply converged rotations in-place
    for si in start_si:end_si
        sec = sections[si]
        le = sec.LE_point; te = sec.TE_point
        t = ((le[1] - le_ref[1]) * y_hat[1] +
             (le[2] - le_ref[2]) * y_hat[2] +
             (le[3] - le_ref[3]) * y_hat[3]) / span_len
        θ = -angle_max * sin(π * t)
        cur = rotated_te(le, te, y_hat, θ)
        te[1] = cur[1]; te[2] = cur[2]; te[3] = cur[3]
    end
    return nothing
end

"""
    refine_mesh_with_billowing!(wing; reuse_aero_data)

Refine wing mesh using SPLIT_PROVIDED spacing with TE billowing.

Between each pair of unrefined (rib) sections, chord vectors are rotated
around the leading edge to simulate fabric billowing.  The rotation
amplitude is found iteratively so that the TE arc length matches the
wing's `billowing_percentage`.

Delegates to [`refine_mesh_by_splitting_provided_sections!`](@ref).
"""
function refine_mesh_with_billowing!(wing; reuse_aero_data::Bool=false)
    refine_mesh_by_splitting_provided_sections!(
        wing; reuse_aero_data,
        billowing_percentage=wing.billowing_percentage
    )
end

"""
    calculate_span(wing::AbstractWing)

Calculate wing span along spanwise direction.

Returns:
    Float64: Wing span
"""
function calculate_span(wing::AbstractWing)
    # Normalize spanwise direction
    vector_axis = wing.spanwise_direction ./ norm(wing.spanwise_direction)
    
    # Get all points
    all_points = reduce(vcat, [[section.LE_point, section.TE_point] 
                              for section in wing.unrefined_sections])
    
    # Project points and calculate span
    projections = [dot(point, vector_axis) for point in all_points]
    return maximum(projections) - minimum(projections)
end

# Project point onto plane
@inline function project_onto_plane!(point_proj, point, normal)
    point_proj .= point .- (point ⋅ normal) .* normal
    return nothing
end

"""
    calculate_projected_area(wing::AbstractWing, z_plane_vector=[0.0, 0.0, 1.0])

Calculate projected wing area onto plane defined by normal vector.

Returns:
    Float64: Projected area
"""
function calculate_projected_area(wing::AbstractWing, 
                                z_plane_vector=[0.0, 0.0, 1.0])
    # Normalize plane normal vector
    z_plane_vector = z_plane_vector ./ norm(z_plane_vector)

    LE_current_proj = zeros(MVec3)
    TE_current_proj = zeros(MVec3)
    LE_next_proj = zeros(MVec3)
    TE_next_proj = zeros(MVec3)
    
    # Calculate area by decomposing each projected panel quadrilateral
    # into two triangles: (A, B, C) and (A, C, D).
    projected_area = 0.0
    for i in 1:(length(wing.unrefined_sections)-1)
        # Get section points
        LE_current = wing.unrefined_sections[i].LE_point
        TE_current = wing.unrefined_sections[i].TE_point
        LE_next = wing.unrefined_sections[i+1].LE_point
        TE_next = wing.unrefined_sections[i+1].TE_point
        
        # Project points
        project_onto_plane!(LE_current_proj, LE_current, z_plane_vector)
        project_onto_plane!(TE_current_proj, TE_current, z_plane_vector)
        project_onto_plane!(LE_next_proj, LE_next, z_plane_vector)
        project_onto_plane!(TE_next_proj, TE_next, z_plane_vector)
        
        # Two triangles: A=LE_i, B=TE_i, C=TE_{i+1}, D=LE_{i+1}
        projected_area += 0.5 * norm(cross(TE_current_proj - LE_current_proj, TE_next_proj - LE_current_proj))
        projected_area += 0.5 * norm(cross(TE_next_proj - LE_current_proj, LE_next_proj - LE_current_proj))
    end
    
    return projected_area
end

# Add span property to Wing struct
Base.propertynames(w::AbstractWing) = (fieldnames(typeof(w))..., :span)
function Base.getproperty(w::AbstractWing, s::Symbol)
    if s === :span
        return calculate_span(w)
    else
        return getfield(w, s)
    end
end

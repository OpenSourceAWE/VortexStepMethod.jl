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
    refined_section.LE_point .= section.LE_point
    refined_section.TE_point .= section.TE_point
    refined_section.aero_model = section.aero_model
    if isnothing(refined_section.aero_data)
        refined_section.aero_data = section.aero_data
    else
        for i in eachindex(section.aero_data)
            copyto!(refined_section.aero_data[i], section.aero_data[i])
        end
    end
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
mutable struct Wing <: AbstractWing
    n_panels::Int16
    n_unrefined_sections::Int16
    spanwise_distribution::PanelDistribution
    panel_props::PanelProperties
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
    le_interp::Union{Nothing, NTuple{3, Extrapolation}}
    te_interp::Union{Nothing, NTuple{3, Extrapolation}}
    area_interp::Union{Nothing, Extrapolation}
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}
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
        n_unrefined_sections=nothing,
        spanwise_distribution::PanelDistribution=LINEAR,
        spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]),
        remove_nan=true,
        use_prior_polar=false,
        billowing_percentage=0.0)

    # For YAML wings, n_unrefined_sections will be set when sections are added
    # Set to 0 as placeholder for now
    n_unrefined_sections_value = isnothing(n_unrefined_sections) ? Int16(0) : Int16(n_unrefined_sections)

    panel_props = PanelProperties{n_panels}()

    # Initialize with default/empty values for optional fields
    Wing(
        n_panels, n_unrefined_sections_value, spanwise_distribution, panel_props, spanwise_direction,
        Section[], Section[], remove_nan, use_prior_polar, Float64(billowing_percentage),
        # Grouping
        Int16[],
        # Deformation fields
        Section[], zeros(max(0, n_panels)), zeros(max(0, n_panels)),
        # Physical properties (defaults for non-OBJ wings)
        0.0, 0.0, zeros(0, 0), zeros(MVec3), Matrix{Float64}(I, 3, 3),
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
function remove_vector_nans(aero_data)
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
        interpolate_matrix_nans!.(aero_data[3:5])
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
    sort_sections && sort!(wing.unrefined_sections, by=s -> s.LE_point[2], rev=true)
    n_sections = wing.n_panels + 1
    reuse_aero_data = _can_reuse_prior_refined_polar_data(wing, n_sections)

    if length(wing.refined_sections) == 0
        if wing.spanwise_distribution == UNCHANGED ||
               length(wing.unrefined_sections) == n_sections
            wing.refined_sections = copy(wing.unrefined_sections)
            recompute_mapping && compute_refined_panel_mapping!(wing)
            update_non_deformed_sections!(wing)
            return nothing
        else
            wing.refined_sections = Section[Section() for _ in 1:wing.n_panels+1]
        end
    end
    
    # Extract geometry data
    n_current = length(wing.unrefined_sections)
    LE = zeros(Float64, n_current, 3)
    TE = zeros(Float64, n_current, 3)
    aero_model = Vector{typeof(wing.unrefined_sections[1].aero_model)}()
    aero_data = Vector{typeof(wing.unrefined_sections[1].aero_data)}()
    
    for (i, section) in enumerate(wing.unrefined_sections)
        LE[i,:] = section.LE_point
        TE[i,:] = section.TE_point
        push!(aero_model, section.aero_model)
        push!(aero_data, section.aero_data)
    end
    
    # Validate input
    if size(LE,1) != size(TE,1) || size(LE,1) != length(aero_model)
        throw(ArgumentError("LE, TE, and aero_model must have the same length"))
    end
    
    # Handle special cases
    if wing.spanwise_distribution == UNCHANGED || length(wing.unrefined_sections) == n_sections
        for i in eachindex(wing.unrefined_sections)
            if reuse_aero_data
                section = wing.unrefined_sections[i]
                reinit!(
                    wing.refined_sections[i],
                    section.LE_point,
                    section.TE_point,
                    section.aero_model
                )
            else
                reinit!(wing.refined_sections[i], wing.unrefined_sections[i])
            end
        end
        recompute_mapping && compute_refined_panel_mapping!(wing)
        update_non_deformed_sections!(wing)
        return nothing
    end

    @debug "Refining aerodynamic mesh from $(length(wing.unrefined_sections)) sections to $n_sections sections."

    # Handle two-section case
    if n_sections == 2
        reinit!(
            wing.refined_sections[1],
            LE[1,:],
            TE[1,:],
            aero_model[1],
            reuse_aero_data ? nothing : aero_data[1]
        )
        reinit!(
            wing.refined_sections[2],
            LE[end,:],
            TE[end,:],
            aero_model[end],
            reuse_aero_data ? nothing : aero_data[end]
        )
        recompute_mapping && compute_refined_panel_mapping!(wing)
        update_non_deformed_sections!(wing)
        return nothing
    end

    # Handle different distribution types
    if wing.spanwise_distribution == SPLIT_PROVIDED
        refine_mesh_by_splitting_provided_sections!(wing; reuse_aero_data)
    elseif wing.spanwise_distribution in (LINEAR, COSINE)
        refine_mesh_for_linear_cosine_distribution!(
            wing,
            1,
            wing.spanwise_distribution,
            n_sections,
            LE,
            TE,
            aero_model,
            aero_data;
            reuse_aero_data
        )
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
    n_unrefined_sections = length(wing.unrefined_sections)
    n_refined_panels = wing.n_panels

    # Handle case where no refinement occurred
    if n_unrefined_sections == n_refined_panels + 1
        wing.refined_panel_mapping = Int16[i for i in 1:n_refined_panels]
        return nothing
    end

    # Ensure mapping array is allocated
    if length(wing.refined_panel_mapping) != n_refined_panels
        wing.refined_panel_mapping = zeros(Int16, n_refined_panels)
    end

    # Compute centers of unrefined sections
    unrefined_centers = Vector{MVec3}(undef, n_unrefined_sections)
    for i in 1:n_unrefined_sections
        le_point = wing.unrefined_sections[i].LE_point
        te_point = wing.unrefined_sections[i].TE_point
        unrefined_centers[i] = MVec3((le_point + te_point) / 2)
    end

    # For each refined panel, find closest unrefined section
    for refined_panel_idx in 1:n_refined_panels
        le_mid = (wing.refined_sections[refined_panel_idx].LE_point +
                  wing.refined_sections[refined_panel_idx+1].LE_point) / 2
        te_mid = (wing.refined_sections[refined_panel_idx].TE_point +
                  wing.refined_sections[refined_panel_idx+1].TE_point) / 2
        refined_center = MVec3((le_mid + te_mid) / 2)

        # Find closest unrefined section
        min_dist = Inf
        closest_idx = 1
        for unrefined_section_idx in 1:n_unrefined_sections
            dist = norm(refined_center - unrefined_centers[unrefined_section_idx])
            if dist < min_dist
                min_dist = dist
                closest_idx = unrefined_section_idx
            end
        end

        wing.refined_panel_mapping[refined_panel_idx] = Int16(closest_idx)
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
    new_sections = Section[]

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

When `billowing_percentage > 0`, applies sinusoidal TE displacement to
intermediate sections within each rib pair (simulating fabric billowing
between ribs).
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
        for (refined_section, section) in zip(
                wing.refined_sections, wing.unrefined_sections)
            if reuse_aero_data
                reinit!(refined_section, section.LE_point,
                    section.TE_point, section.aero_model)
            else
                reinit!(refined_section, section)
            end
        end
        if billowing_percentage > 0
            LE = [s.LE_point for s in wing.unrefined_sections]
            TE = [s.TE_point for s in wing.unrefined_sections]
            for i in 1:(n_panels_provided - 1)
                y_hat = normalize(LE[i] - LE[i + 1])
                span_len = norm(LE[i] - LE[i + 1])
                start_si = i + 1
                end_si = i + 1
                if start_si <= end_si
                    apply_billowing_to_pair!(
                        wing.refined_sections,
                        start_si, end_si,
                        y_hat, span_len,
                        LE[i + 1], TE[i], TE[i + 1],
                        billowing_percentage)
                end
            end
        end
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
    
    # Extract geometry data
    LE = [section.LE_point for section in wing.unrefined_sections]
    TE = [section.TE_point for section in wing.unrefined_sections]
    aero_model = [section.aero_model for section in wing.unrefined_sections]
    aero_data = [section.aero_data for section in wing.unrefined_sections]
    
    # Process each section pair
    idx = 1
    for left_section_index in 1:n_section_pairs
        # Add left section of pair
        if reuse_aero_data
            left_section = wing.unrefined_sections[left_section_index]
            reinit!(
                wing.refined_sections[idx],
                left_section.LE_point,
                left_section.TE_point,
                left_section.aero_model
            )
        else
            reinit!(wing.refined_sections[idx], wing.unrefined_sections[left_section_index])
        end
        idx += 1

        # Calculate new sections for this pair
        num_new_sections = new_sections_per_pair + (left_section_index <= remaining ? 1 : 0)

        if num_new_sections > 0
            # Prepare pair data
            LE_pair = hcat(LE[left_section_index], LE[left_section_index + 1])'
            TE_pair = hcat(TE[left_section_index], TE[left_section_index + 1])'
            aero_model_pair = [
                aero_model[left_section_index],
                aero_model[left_section_index + 1]
            ]
            aero_data_pair = [
                aero_data[left_section_index],
                aero_data[left_section_index + 1]
            ]

            # Generate sections for this pair
            start_idx = idx
            idx = refine_mesh_for_linear_cosine_distribution!(
                wing,
                idx,
                LINEAR,
                num_new_sections + 2,  # +2 for endpoints
                LE_pair,
                TE_pair,
                aero_model_pair,
                aero_data_pair;
                endpoints=false,
                reuse_aero_data
            )

            # Apply billowing by rotating chords around LE
            if billowing_percentage > 0 && idx > start_idx
                y_hat = normalize(
                    LE[left_section_index] -
                    LE[left_section_index + 1])
                span_len = norm(
                    LE[left_section_index] -
                    LE[left_section_index + 1])
                apply_billowing_to_pair!(
                    wing.refined_sections,
                    start_idx, idx - 1,
                    y_hat, span_len,
                    LE[left_section_index + 1],
                    TE[left_section_index],
                    TE[left_section_index + 1],
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
    idx += 1
    
    # Validate result
    if length(wing.refined_sections) != wing.n_panels + 1
        @warn "Number of panels ($(length(wing.refined_sections)-1)) differs from desired ($(wing.n_panels))"
    end
    
    return nothing
end


"""
    rodrigues_rotate(vec, axis, θ)

Rotate `vec` around unit `axis` by angle `θ` (radians) using the
Rodrigues rotation formula. Non-allocating with MVec3 inputs.
"""
@inline function rodrigues_rotate(vec, axis, θ)
    ct = cos(θ); st = sin(θ)
    d = dot(vec, axis)
    return ct * vec + st * cross(axis, vec) + (1 - ct) * d * axis
end

"""
    billowing_arc_length(sections, start_si, end_si, y_hat,
                         span_len, le_ref, te_left, te_right,
                         angle_max)

Compute the TE arc length that would result from rotating each
section's chord around `y_hat` by `angle_max * sin(π t)` radians,
where `t` is the normalised spanwise position within the rib pair.

Non-allocating: uses only stack-allocated MVec3 arithmetic.
"""
function billowing_arc_length(
    sections, start_si, end_si,
    y_hat, span_len, le_ref,
    te_left, te_right, angle_max
)
    prev_te = MVec3(te_left)
    arc = 0.0
    for si in start_si:end_si
        sec = sections[si]
        t = dot(sec.LE_point - le_ref, y_hat) / span_len
        θ = -angle_max * sin(π * t)
        chord_vec = sec.TE_point - sec.LE_point
        rotated = rodrigues_rotate(chord_vec, y_hat, θ)
        current_te = sec.LE_point + rotated
        arc += norm(current_te - prev_te)
        prev_te = current_te
    end
    arc += norm(te_right - prev_te)
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
in-place using only stack-allocated MVec3 arithmetic.
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
    straight = norm(te_right - te_left)
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
        t = dot(sec.LE_point - le_ref, y_hat) / span_len
        θ = -angle_max * sin(π * t)
        chord_vec = sec.TE_point - sec.LE_point
        rotated = rodrigues_rotate(chord_vec, y_hat, θ)
        sec.TE_point .= sec.LE_point + rotated
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

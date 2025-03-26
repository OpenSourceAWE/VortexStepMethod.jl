
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
    Section(LE_point, TE_point, aero_model)

Create a new wing section with the specified leading edge point, trailing edge point, 
and aerodynamic model.

# Arguments
- `LE_point::MVec3`: Leading edge point coordinates
- `TE_point::MVec3`: Trailing edge point coordinates  
- `aero_model::AeroModel`: Aerodynamic model type (e.g., INVISCID, POLAR_VECTORS)

# Returns
- `Section`: A new section with the specified parameters and no aerodynamic data
"""
function Section(LE_point, TE_point, aero_model)
    return Section(LE_point, TE_point, aero_model, nothing)
end

"""
    init!(section::Section, LE_point, TE_point, aero_model=nothing, aero_data=nothing)

Function to update a [Section](@ref) in place.
"""
function init!(section::Section, LE_point::AbstractVector, TE_point::AbstractVector, aero_model=nothing, aero_data=nothing)
    section.LE_point .= LE_point
    section.TE_point .= TE_point
    (!isnothing(aero_model)) && (section.aero_model = aero_model)
    if !isnothing(aero_data)
        if !isnothing(section.aero_data)
            section.aero_data .= aero_data
        else 
            section.aero_data = aero_data
        end
    end
    nothing
end

function init!(refined_section::Section, section::Section)
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
    update_panel_properties!(section_list::Vector{Section}, n_panels::Int)

Calculate geometric properties for each panel.

# Arguments
- section_list::Vector{Section}: List of [Section](@ref)s
- `n_panels`::Int: Number of [Panel](@ref)s

# Returns:
[PanelProperties](@ref) containing vectors for each property
"""
function update_panel_properties!(panel_props::PanelProperties, section_list::Vector{Section}, n_panels::Int)
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

# Fields
- `n_panels::Int16`: Number of panels in aerodynamic mesh
- `n_groups::Int16`: Number of panel groups
- `spanwise_distribution`::PanelDistribution: [PanelDistribution](@ref)
- `spanwise_direction::MVec3`: Wing span direction vector
- `sections::Vector{Section}`: Vector of wing sections, see: [Section](@ref)
- `refined_sections::Vector{Section}`: Vector of refined wing sections, see: [Section](@ref)
- `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not

"""
mutable struct Wing <: AbstractWing
    n_panels::Int16
    n_groups::Int16
    spanwise_distribution::PanelDistribution
    panel_props::PanelProperties
    spanwise_direction::MVec3
    sections::Vector{Section}
    refined_sections::Vector{Section}
    remove_nan::Bool
end

"""
    Wing(n_panels::Int;
         n_groups=n_panels,
         spanwise_distribution::PanelDistribution=LINEAR,
         spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]),
         remove_nan::Bool=true)

Constructor for a [Wing](@ref) struct with default values that initializes the sections 
and refined sections as empty arrays.

# Parameters
- `n_panels::Int`: Number of panels in aerodynamic mesh
- `n_groups::Int`: Number of panel groups in aerodynamic mesh
- `spanwise_distribution`::PanelDistribution = LINEAR: [PanelDistribution](@ref)
- `spanwise_direction::MVec3` = MVec3([0.0, 1.0, 0.0]): Wing span direction vector
- `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not
"""
function Wing(n_panels::Int;
        n_groups = n_panels,
        spanwise_distribution::PanelDistribution=LINEAR,
        spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]),
        remove_nan=true)
    panel_props = PanelProperties{n_panels}()
    Wing(n_panels, n_groups, spanwise_distribution, panel_props, spanwise_direction, Section[], Section[], remove_nan)
end

function init!(wing::AbstractWing)
    refine_aerodynamic_mesh!(wing)
    
    # Calculate panel properties
    update_panel_properties!(
        wing.panel_props,
        wing.refined_sections,
        wing.n_panels
    )
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
    push!(wing.sections, Section(LE_point, TE_point, aero_model, aero_data))
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
    refine_aerodynamic_mesh!(wing::AbstractWing)

Refine the aerodynamic mesh of the wing based on spanwise panel distribution.

Returns:
    Vector{Section}: List of refined sections
"""
function refine_aerodynamic_mesh!(wing::AbstractWing)
    sort!(wing.sections, by=s -> s.LE_point[2], rev=true)
    n_sections = wing.n_panels + 1
    if length(wing.refined_sections) == 0
        if wing.spanwise_distribution == UNCHANGED || length(wing.sections) == n_sections
            wing.refined_sections = wing.sections
            return nothing
        else
            wing.refined_sections = Section[Section() for _ in 1:wing.n_panels+1]
        end
    end
    
    # Extract geometry data
    n_current = length(wing.sections)
    LE = zeros(Float64, n_current, 3)
    TE = zeros(Float64, n_current, 3)
    aero_model = Vector{typeof(wing.sections[1].aero_model)}()
    aero_data = Vector{typeof(wing.sections[1].aero_data)}()
    
    for (i, section) in enumerate(wing.sections)
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
    if wing.spanwise_distribution == UNCHANGED || length(wing.sections) == n_sections
        for i in eachindex(wing.sections)
            init!(wing.refined_sections[i], wing.sections[i])
        end
        return nothing
    end

    @debug "Refining aerodynamic mesh from $(length(wing.sections)) sections to $n_sections sections."
    
    # Handle two-section case
    if n_sections == 2
        init!(wing.refined_sections[1], LE[1,:], TE[1,:], aero_model[1], aero_data[1])
        init!(wing.refined_sections[2], LE[end,:], TE[end,:], aero_model[end], aero_data[end])
        return nothing
    end
    
    # Handle different distribution types
    if wing.spanwise_distribution == SPLIT_PROVIDED
        return refine_mesh_by_splitting_provided_sections!(wing)
    elseif wing.spanwise_distribution in (LINEAR, COSINE, COSINE_VAN_GARREL)
        return refine_mesh_for_linear_cosine_distribution!(
            wing,
            1,
            wing.spanwise_distribution,
            n_sections,
            LE,
            TE,
            aero_model,
            aero_data
        )
    else
        throw(ArgumentError("Unsupported spanwise panel distribution: $(wing.spanwise_distribution)"))
    end
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

            # Create common alpha array
            !all(isapprox.(alpha_left, alpha_right)) && @error "Make sure you use the same alpha range for all your interpolations."
            !isa(CL_right, AbstractVector) && @error "Provide polar data in the correct format: (alpha, cl, cd, cm)"
            
            # Weighted interpolation
            CL_data = CL_left .* left_weight .+ CL_right .* right_weight
            CD_data = CD_left .* left_weight .+ CD_right .* right_weight
            CM_data = CM_left .* left_weight .+ CM_right .* right_weight
            
            return (alpha_left, CL_data, CD_data, CM_data)
            
        elseif model_type == POLAR_MATRICES
            alpha_left, delta_left, CL_left, CD_left, CM_left = polar_left
            alpha_right, delta_right, CL_right, CD_right, CM_right = polar_right
            
            # Create common alpha array
            !all(isapprox.(alpha_left, alpha_right)) && @error "Make sure you use the same alpha range for all your interpolations."
            !all(isapprox.(delta_left, delta_right)) && @error "Make sure you use the same alpha range for all your interpolations."
            !isa(CL_right, AbstractMatrix) && @error "Provide polar data in the correct format: (alpha, delta, cl, cd, cm)"

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
    endpoints=true)

    # 1. Compute quarter chord line
    quarter_chord = LE .+ 0.25 .* (TE .- LE)

    # Calculate segment lengths
    qc_lengths = [norm(quarter_chord[i+1,:] - quarter_chord[i,:]) for i in 1:size(quarter_chord,1)-1]
    qc_total_length = sum(qc_lengths)
    qc_cum_length = vcat(0, cumsum(qc_lengths))

    # 2. Define target lengths
    target_lengths = if spanwise_distribution == LINEAR
        range(0, qc_total_length, n_sections)
    elseif spanwise_distribution in (COSINE, COSINE_VAN_GARREL)
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

        # Interpolate aero properties
        new_data = calculate_new_aero_data(aero_model, aero_data, section_index, left_weight, right_weight)

        # Create new section
        if endpoints || (i != 1 && i != n_sections)
            @views init!(wing.refined_sections[idx], new_LE[i,:], new_TE[i,:], aero_model[1], new_data)
            idx += 1
        end
    end

    # Apply van Garrel distribution if requested
    if spanwise_distribution == COSINE_VAN_GARREL
        idx = calculate_cosine_van_Garrel!(wing, idx)
    end

    return idx
end


"""
    calculate_cosine_van_Garrel!(wing::AbstractWing, idx)

Calculate van Garrel cosine distribution of sections.
Reference: http://dx.doi.org/10.13140/RG.2.1.2773.8000

Returns:
    idx
"""
function calculate_cosine_van_Garrel!(wing::AbstractWing, idx)
    n = length(sections)
    
    # Calculate chords and quarter chords
    chords = [section.TE_point - section.LE_point for section in sections]
    quarter_chords = [section.LE_point + 0.25 * chord for (section, chord) in zip(sections, chords)]
    
    # Calculate widths
    widths = [norm(quarter_chords[i+1] - quarter_chords[i]) for i in 1:n-1]
    
    # Calculate correction factors
    eta_cp = zeros(n-1)
    
    # First panel
    eta_cp[1] = widths[1] / (widths[1] + widths[2])
    
    # Internal panels
    for j in 2:n-2
        eta_cp[j] = 0.25 * (
            widths[j-1] / (widths[j-1] + widths[j]) +
            widths[j] / (widths[j] + widths[j+1]) + 1
        )
    end
    
    # Last panel
    eta_cp[end] = widths[end-1] / (widths[end-1] + widths[end])
    
    @debug "Correction factors" eta_cp
    
    # Calculate control points
    control_points = [
        quarter_chords[i] + eta * (quarter_chords[i+1] - quarter_chords[i])
        for (i, eta) in enumerate(eta_cp)
    ]
    
    # Generate new sections
    for (i, (control_point, chord)) in enumerate(zip(control_points, chords))
        @views init!(wing.refined_sections, 
            control_point - 0.25 * chord, 
            control_point + 0.75 * chord,
            sections[i].aero_model, 
            sections[i].aero_data
        )
        idx += 1
    end
    return idx
end


"""
    refine_mesh_by_splitting_provided_sections!(wing::AbstractWing)

Refine mesh by splitting provided sections into desired number of panels.

Returns:
    Vector{Section}: Refined sections
"""
function refine_mesh_by_splitting_provided_sections!(wing::AbstractWing)
    n_sections_provided = length(wing.sections)
    n_panels_provided = n_sections_provided - 1
    n_panels_desired = wing.n_panels
    
    @debug "Panel counts" n_panels_provided n_panels_desired n_sections_provided
    
    # Check if refinement is needed
    if n_panels_provided == n_panels_desired
        for (refined_section, section) in zip(wing.refined_sections, wing.sections)
            init!(refined_section, section)
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
    LE = [section.LE_point for section in wing.sections]
    TE = [section.TE_point for section in wing.sections]
    aero_model = [section.aero_model for section in wing.sections]
    aero_data = [section.aero_data for section in wing.sections]
    
    # Process each section pair
    idx = 1
    for left_section_index in 1:n_section_pairs
        # Add left section of pair
        init!(wing.refined_sections[idx], wing.sections[left_section_index])
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
            idx = refine_mesh_for_linear_cosine_distribution!(
                wing,
                idx,
                LINEAR,
                num_new_sections + 2,  # +2 for endpoints
                LE_pair,
                TE_pair,
                aero_model_pair,
                aero_data_pair;
                endpoints=false
            )
        end
    end
    
    # Add final section
    init!(wing.refined_sections[idx], wing.sections[end])
    idx += 1
    
    # Validate result
    if length(wing.refined_sections) != wing.n_panels + 1
        @warn "Number of panels ($(length(new_sections)-1)) differs from desired ($(wing.n_panels))"
    end
    
    return nothing
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
                              for section in wing.sections])
    
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
    
    # Calculate area by summing trapezoid areas
    projected_area = 0.0
    for i in 1:(length(wing.sections)-1)
        # Get section points
        LE_current = wing.sections[i].LE_point
        TE_current = wing.sections[i].TE_point
        LE_next = wing.sections[i+1].LE_point
        TE_next = wing.sections[i+1].TE_point
        
        # Project points
        project_onto_plane!(LE_current_proj, LE_current, z_plane_vector)
        project_onto_plane!(TE_current_proj, TE_current, z_plane_vector)
        project_onto_plane!(LE_next_proj, LE_next, z_plane_vector)
        project_onto_plane!(TE_next_proj, TE_next, z_plane_vector)
        
        # Calculate projected dimensions
        chord_current = norm(TE_current_proj - LE_current_proj)
        chord_next = norm(TE_next_proj - LE_next_proj)
        span = norm(LE_next_proj - LE_current_proj)
        
        # Add trapezoid area
        projected_area += 0.5 * (chord_current + chord_next) * span
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
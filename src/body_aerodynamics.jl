"""
    @with_kw mutable struct BodyAerodynamics{P,W<:AbstractWing}

Main structure for calculating aerodynamic properties of bodies. Use the constructor to initialize.

# Fields
- panels::Vector{<:Panel}: Vector of refined [`Panel`](@ref) structs
- wings::Vector{W}: A vector of wings of type `W <: AbstractWing`; a body can have multiple wings
- `va_vec::MVec3` = zeros(MVec3): apparent wind vector [m/s], see: [`MVec3`](@ref)
- `omega`::MVec3 = zeros(MVec3): A vector of the turn rates around the KA body axes
- `reference_point`::MVec3 = zeros(MVec3): The point `omega` turns the body about [m]
- `gamma_distribution`=zeros(Float64, P): A vector of the circulation
                        of the velocity field; Length: Number of segments. [m²/s]
- `alpha_uncorrected`=zeros(Float64, P): angles of attack per panel
- `alpha_corrected`=zeros(Float64, P):   corrected angles of attack per panel
- `stall_angle_list`=zeros(Float64, P):  stall angle per panel
- `alpha_dist::MVector{P, Float64}` = zeros(Float64, P)
- `v_rel_dist::MVector{P, Float64}` = zeros(Float64, P): norm of the relative velocity
    crossed with the panel spanwise axis, |v_rel × y_airf| [m/s]
- `pitch_rate_dist::MVector{P, Float64}` = zeros(Float64, P): rotation rate of each
    panel about its own spanwise axis, positive nose-up [rad/s]; set by
    [`set_va!`](@ref) and read when the solver has `flow_curvature` enabled
- `work_vectors`::NTuple{10, MVec3} = ntuple(_ -> zeros(MVec3), 10)
- `AIC::Array{Float64, 3}` = zeros(P, P, 3): control-point influence coefficients, the
                        matrix the circulation is solved against; component last so that
                        each `AIC[:, :, k]` slice is a contiguous BLAS matrix
- `AIC_curved_correction::Array{Float64, 3}` = zeros(P, P, 3): stored
    [`curved_bound_correction!`](@ref)
- `curved_correction_core_radius` = NaN: core radius fraction `AIC_curved_correction` was
    built for, `NaN` until built for the current shape
- `AIC_aero_center::Array{Float64, 3}` = zeros(P, P, 3): aerodynamic-centre (LLT)
                        influence coefficients, used only for the corrected angle of attack
- `projected_area::Float64` = 1.0: The area projected onto the xy-plane of the KA body frame [m²]
- `c_ref::Float64` = 1.0: Reference chord length (max panel chord) [m]
- `cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}` = [LazyBufferCache()]
"""
@with_kw mutable struct BodyAerodynamics{P, W<:AbstractWing, T, PN<:Panel{T}}
    panels::Vector{PN}
    wings::Vector{W}
    va_vec::MVector{3, T} = zeros(MVector{3, T})
    has_distributed_va::Bool = false
    omega::MVector{3, T} = zeros(MVector{3, T})
    reference_point::MVector{3, T} = zeros(MVector{3, T})
    gamma_distribution::MVector{P, T} = zeros(MVector{P, T})
    alpha_uncorrected::MVector{P, T} = zeros(MVector{P, T})
    alpha_corrected::MVector{P, T} = zeros(MVector{P, T})
    stall_angle_list::MVector{P, T} = zeros(MVector{P, T})
    alpha_dist::MVector{P, T} = zeros(MVector{P, T})
    v_rel_dist::MVector{P, T} = zeros(MVector{P, T})
    pitch_rate_dist::MVector{P, T} = zeros(MVector{P, T})
    work_vectors::NTuple{10, MVector{3, T}} = ntuple(_ -> zeros(MVector{3, T}), 10)
    AIC::Array{T, 3} = zeros(T, P, P, 3)
    AIC_aero_center::Array{T, 3} = zeros(T, P, P, 3)
    AIC_curved_correction::Array{T, 3} = zeros(T, P, P, 3)
    curved_correction_core_radius::T = T(NaN)
    projected_area::T = one(T)
    c_ref::T = one(T)
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache()]
end

"""
    BodyAerodynamics(wings::Vector{T}; 
                     kite_body_origin=zeros(MVec3)) where T <: AbstractWing

Construct a [`BodyAerodynamics`](@ref) object for aerodynamic calculations.

This constructor handles initialization of panels, coordinate transformations, and
aerodynamic properties, returning a fully initialized structure ready for simulation.

# Arguments
- `wings::Vector{T}`: Vector of wings to analyze, where T is an AbstractWing type

# Keyword Arguments
- `kite_body_origin=zeros(MVec3)`: Origin point of the KA body frame in CAD reference frame
- `va_vec=[15.0, 0.0, 0.0]`: Apparent wind vector [m/s]
- `omega=zeros(3)`: Turn rate about the KA body x, y and z axes

# Returns
- [`BodyAerodynamics`](@ref) object initialized with panels and wings

# Example
```julia
wing = Wing("wing.yaml"; n_panels=40); refine!(wing)
body_aero = BodyAerodynamics([wing], va_vec=[15.0, 0.0, 0.0], omega=zeros(3))
```
"""
function BodyAerodynamics(
    wings::Vector{W};
    kite_body_origin=zeros(MVec3),
    va_vec=[15.0, 0.0, 0.0],
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
    reinit!(body_aero; va_vec, omega)
    return body_aero
end

"""
    wing_span_flip(wing) -> Int8

`-1` when `wing`'s sections run against its `spanwise_direction`, `+1` otherwise: the
`flip` every panel of the wing is reinitialized with ([`reinit!`](@ref)). One answer per
wing, so neighbouring panels cannot disagree and invert a single normal by 180°.
"""
wing_span_flip(wing) =
    dot(first(wing.refined_sections).LE_point - last(wing.refined_sections).LE_point,
        wing.spanwise_direction) < 0 ? Int8(-1) : Int8(1)

function Base.getproperty(obj::BodyAerodynamics, sym::Symbol)
    if sym === :va_vec && getfield(obj, :has_distributed_va)
        throw(ArgumentError(
            "body_aero.va_vec is undefined after set_va! with distributed inflow. " *
            "Use panel.va_vec or solver.sol.va_vec_dist for per-panel inflow data."
        ))
    end
    return getfield(obj, sym)
end

function Base.setproperty!(obj::BodyAerodynamics, sym::Symbol, val)
    if sym === :va_vec
        set_va!(obj, val)
    elseif sym === :omega
        set_va!(obj, getfield(obj, :va_vec), val)
    elseif sym === :reference_point
        set_va!(obj, getfield(obj, :va_vec), obj.omega; reference_point=val)
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

        # A table that stops short of the scan says nothing about a stall in it:
        # past its end it is held flat, which is not a peak.
        if panel.alpha_window > 0 &&
           panel.alpha_ref + panel.alpha_window < deg2rad(begin_aoa)
            stall_angles[idx] = panel_stall
            continue
        end

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
    unrefined_section_range(body_aero::BodyAerodynamics, wing_idx)

Indices of the unrefined sections of wing `wing_idx` in a distribution that runs over the
unrefined sections of all wings in order, such as `moment_unrefined_dist`.
"""
function unrefined_section_range(body_aero::BodyAerodynamics, wing_idx)
    wings = body_aero.wings
    offset = sum(wing.n_unrefined_sections for wing in view(wings, 1:wing_idx-1); init=0)
    return offset .+ (1:wings[wing_idx].n_unrefined_sections)
end

"""
    panel_range(body_aero::BodyAerodynamics, wing_idx)

Indices of the panels of wing `wing_idx` in `body_aero.panels`.
"""
function panel_range(body_aero::BodyAerodynamics, wing_idx)
    wings = body_aero.wings
    offset = sum(wing.n_panels for wing in view(wings, 1:wing_idx-1); init=0)
    return offset .+ (1:wings[wing_idx].n_panels)
end

"""
    unrefined_deform!(body_aero::BodyAerodynamics, theta_angles, delta_angles)

Deform each wing of `body_aero` by its entries of `theta_angles` and `delta_angles` [rad],
which run over the unrefined sections of all wings in order; `nothing` leaves that angle
unchanged. Call [`reinit!`](@ref) afterwards to update the panels.
"""
function unrefined_deform!(body_aero::BodyAerodynamics, theta_angles, delta_angles)
    for (wing_idx, wing) in enumerate(body_aero.wings)
        section_range = unrefined_section_range(body_aero, wing_idx)
        unrefined_deform!(wing,
            isnothing(theta_angles) ? nothing : view(theta_angles, section_range),
            isnothing(delta_angles) ? nothing : view(delta_angles, section_range))
    end
    return nothing
end

"""
    reinit!(body_aero::BodyAerodynamics; init_aero, va_vec, omega)

Initialize a BodyAerodynamics struct in-place by setting up panels and coefficients.

# Arguments
- `body_aero::BodyAerodynamics`: The structure to initialize

# Keyword Arguments
- `init_aero::Bool`: Whether to initialize the aero data or not
- `va_vec=[15.0, 0.0, 0.0]`: Apparent wind vector [m/s]
- `omega=zeros(3)`: Turn rate about the KA body x, y and z axes

# Returns
nothing
"""
function reinit!(body_aero::BodyAerodynamics{P, W, T};
    init_aero=true,
    va_vec=[15.0, 0.0, 0.0],
    omega=zeros(MVector{3, T})
) where {P, W, T}
    idx = 1
    vec = zeros(MVector{3, T})
    for (wing_idx, wing) in enumerate(body_aero.wings)
        reinit!(wing)
        validate_section_aero(wing.refined_sections)
        panel_props = wing.panel_props
        wing_init_aero = init_aero && !_can_skip_panel_aero_reinit(wing, body_aero.panels, idx)
        wing_flip = wing_span_flip(wing) == -1

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
                vec;
                remove_nan=wing.remove_nan,
                init_aero=wing_init_aero,
                flip=wing_flip
            )
            body_aero.panels[idx].crease_frac = wing.crease_frac
            idx += 1
        end
    end

    # Initialize rest of the struct
    body_aero.projected_area = sum(calculate_projected_area, body_aero.wings)
    isempty(body_aero.panels) && throw(ArgumentError("Cannot compute c_ref: body_aero has no panels."))
    body_aero.c_ref = maximum(panel.chord for panel in body_aero.panels)
    calculate_stall_angle_list!(body_aero.stall_angle_list, body_aero.panels)
    body_aero.alpha_dist .= 0.0
    body_aero.v_rel_dist .= 0.0
    body_aero.AIC .= 0.0
    body_aero.curved_correction_core_radius = T(NaN)
    set_va!(body_aero, va_vec, omega)
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
        throw(ArgumentError("va_vec must be shape (3,) or va_vec_dist ($(n_panels), 3); " *
                            "got length $(length(va_input))"))
    T = eltype(va_input)
    return MVector{3, T}(va_input[1], va_input[2], va_input[3])
end

@inline function _compute_reference_velocity_from_distribution(
    va_input::AbstractMatrix,
    n_panels::Int,
    panel_areas::Union{Nothing, AbstractVector}=nothing
)
    size(va_input) == (n_panels, 3) ||
        throw(ArgumentError("va_vec must be shape (3,) or va_vec_dist ($(n_panels), 3); " *
                            "got $(size(va_input))"))
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
    calculate_AIC_matrices!(body_aero::BodyAerodynamics, model::Model, core_radius_fraction,
                            va_dist, va_unit_dist, target=body_aero.AIC)

Calculate Aerodynamic Influence Coefficient matrices.

See also: [`BodyAerodynamics`](@ref), [`Model`](@ref)

Returns: nothing
"""
@inline function calculate_AIC_matrices!(body_aero::BodyAerodynamics{P, W, T}, model::Model,
                              core_radius_fraction,
                              va_dist::AbstractVector{T},
                              va_unit_dist::AbstractMatrix{T},
                              target::AbstractArray{T, 3}=body_aero.AIC) where {P, W, T}
    # Determine evaluation point based on model
    evaluation_point = model == VSM ? :control_point : :aero_center
    evaluation_point_on_bound = model == LLT

    # Allocate work vectors for this function (separate from those used by child functions)
    velocity_induced = zeros(MVector{3, T})
    tempvel = zeros(MVector{3, T})
    va_unit = zeros(MVector{3, T})
    U_2D = zeros(MVector{3, T})

    all(>(0), va_dist) || throw(ArgumentError(
        "Wake reference speed must be positive on every panel."))

    # Calculate influence coefficients
    for jring in eachindex(body_aero.panels)
        panel_jring = body_aero.panels[jring]
        filaments = panel_jring.filaments
        va_unit .= matrix_row(va_unit_dist, jring)
        for icp in eachindex(body_aero.panels)
            panel_icp = body_aero.panels[icp]
            ep = evaluation_point == :control_point ? panel_icp.control_point :
                 panel_icp.aero_center
            calculate_velocity_induced_single_ring_semiinfinite!(
                velocity_induced,
                tempvel,
                filaments,
                ep,
                evaluation_point_on_bound,
                va_dist[jring],
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
            @inbounds for k in 1:3
                target[icp, jring, k] = velocity_induced[k]
            end
        end
    end
    return nothing
end

"""
    calculate_circulation_distribution_elliptical_wing(gamma_i, body_aero::BodyAerodynamics,
                                                       gamma_0=1.0)

Write into `gamma_i` an elliptic circulation of peak `gamma_0` over each wing, its control
points measured along that wing's `spanwise_direction` from the wing's mid-span.
"""
function calculate_circulation_distribution_elliptical_wing(gamma_i,
        body_aero::BodyAerodynamics, gamma_0=1.0)
    for (wing_idx, wing) in enumerate(body_aero.wings)
        lo, hi = spanwise_extent(wing)
        axis = normalize(wing.spanwise_direction)
        for i in panel_range(body_aero, wing_idx)
            span_position = dot(body_aero.panels[i].control_point, axis) - (lo + hi) / 2
            # Clamped: a control point can lie outside the span of the unrefined sections
            gamma_i[i] = gamma_0 * sqrt(max(0.0, 1 - (2span_position / (hi - lo))^2))
        end
    end
    return nothing
end

"""
    thin_airfoil_strips(n_strips=24)

Chordwise stations `x/c` and weights of the flat-plate bound vorticity
`gamma(x) ~ sqrt((c - x) / x)`. The weights sum to one with their centroid at the quarter
chord.
"""
function thin_airfoil_strips(n_strips=24)
    phi = ((1:n_strips) .- 0.5) .* (π / 2) ./ n_strips
    weights = cos.(phi) .^ 2
    return sin.(phi) .^ 2, weights ./ sum(weights)
end

"""Velocity at `point` of unit circulation on a trailing filament from `x1` to `x2`."""
function unit_trailing_velocity(filament::BoundFilament{T}, x1, x2, point, va,
        work_vectors) where {T}
    reinit!(filament, x1, x2)
    velocity = work_vectors[8]
    velocity_3D_trailing_vortex!(velocity, filament, point, one(T), va, work_vectors)
    return SVector{3, T}(velocity)
end

"""Velocity at `point` of unit circulation on a bound filament from `x1` to `x2`."""
function unit_bound_velocity(filament::BoundFilament{T}, x1, x2, point,
        core_radius_fraction, work_vectors) where {T}
    reinit!(filament, x1, x2)
    velocity = work_vectors[8]
    velocity_3D_bound_vortex!(velocity, filament, point, one(T), core_radius_fraction,
        work_vectors)
    return SVector{3, T}(velocity)
end

"""
    unrolled_ring_difference(filament, work_vectors, point, panel, unrolled_1, unrolled_2,
                             chord_leg, stations, weights, core_radius_fraction, va)

Velocity at `point` induced by unit circulation on the horseshoe of `panel` less the same
horseshoe unrolled to the bound points `unrolled_1` and `unrolled_2`, whose legs are
`chord_leg`. Both are spread along the chord over the [`thin_airfoil_strips`](@ref)
`stations` and `weights`. `filament` is scratch.
"""
function unrolled_ring_difference(filament::BoundFilament{T}, work_vectors, point, panel,
        unrolled_1, unrolled_2, chord_leg, stations, weights, core_radius_fraction,
        va) where {T}
    bound_1 = SVector{3, T}(panel.bound_point_1)
    bound_2 = SVector{3, T}(panel.bound_point_2)
    te_1 = SVector{3, T}(panel.TE_point_1)
    te_2 = SVector{3, T}(panel.TE_point_2)
    unrolled_te_1 = unrolled_1 + chord_leg
    unrolled_te_2 = unrolled_2 + chord_leg
    velocity = zero(SVector{3, T})
    # Each pair below reads (from, to, sign), legs of the two systems first.
    for (x1, x2, sign) in ((unrolled_1, unrolled_te_1, 1), (unrolled_te_2, unrolled_2, 1),
                           (bound_1, te_1, -1), (te_2, bound_2, -1))
        velocity += sign * unit_trailing_velocity(filament, x1, x2, point, va, work_vectors)
    end
    for (station, weight) in zip(stations, weights)
        along_leg = (station - 0.25) / 0.75
        node_1 = bound_1 + along_leg * (te_1 - bound_1)
        node_2 = bound_2 + along_leg * (te_2 - bound_2)
        unrolled_node_1 = unrolled_1 + along_leg * chord_leg
        unrolled_node_2 = unrolled_2 + along_leg * chord_leg
        strip = unit_bound_velocity(filament, node_2, node_1, point, core_radius_fraction,
            work_vectors) - unit_bound_velocity(filament, unrolled_node_2, unrolled_node_1,
            point, core_radius_fraction, work_vectors)
        for (x1, x2, sign) in ((node_1, te_1, 1), (te_2, node_2, 1),
                               (unrolled_node_1, unrolled_te_1, -1),
                               (unrolled_te_2, unrolled_node_2, -1))
            strip += sign * unit_trailing_velocity(filament, x1, x2, point, va,
                work_vectors)
        end
        velocity += weight * strip
    end
    return velocity
end

"""
    curved_bound_correction!(correction, body_aero::BodyAerodynamics, core_radius_fraction,
                             va_dist)

Write into `correction` the `(i, j, k)` influence, per unit circulation of panel `j`, of
the curvature, sweep and rib skew of a wing at the aerodynamic centre of its panel `i`:
the actual vortex system less the one unrolled onto the tangent of panel `i`, with the
bound vortex spread over the chord. Zero on a straight wing with perpendicular legs, and
between panels of different wings.
"""
function curved_bound_correction!(correction, body_aero::BodyAerodynamics{P, W, T},
        core_radius_fraction, va_dist) where {P, W, T}
    correction .= 0
    stations, weights = thin_airfoil_strips()
    filament = BoundFilament{T}()
    panels = body_aero.panels
    for wing_idx in eachindex(body_aero.wings)
        wing_panels = panel_range(body_aero, wing_idx)
        node_position = cumsum([zero(T); [panels[j].width for j in wing_panels]])
        next_panel = panels[min(first(wing_panels) + 1, last(wing_panels))]
        is_forward = norm(panels[first(wing_panels)].bound_point_1 - next_panel.bound_point_2) <
            norm(panels[first(wing_panels)].bound_point_2 - next_panel.bound_point_1)
        sequence = is_forward ? one(T) : -one(T)
        for (local_i, i) in enumerate(wing_panels)
            centre = SVector{3, T}(panels[i].aero_center)
            tangent = SVector{3, T}(panels[i].bound_point_1 - panels[i].bound_point_2) /
                panels[i].width
            position_i = (node_position[local_i] + node_position[local_i + 1]) / 2
            for (local_j, j) in enumerate(wing_panels)
                node_low = centre + sequence * (node_position[local_j] - position_i) * tangent
                node_high = centre +
                    sequence * (node_position[local_j + 1] - position_i) * tangent
                unrolled_1, unrolled_2 = is_forward ? (node_high, node_low) :
                    (node_low, node_high)
                chord_leg = 0.75 * panels[j].chord * SVector{3, T}(panels[i].x_airf)
                velocity = unrolled_ring_difference(filament, body_aero.work_vectors, centre,
                    panels[j], unrolled_1, unrolled_2, chord_leg, stations, weights,
                    core_radius_fraction, va_dist[i])
                correction[i, j, :] .= velocity
            end
        end
    end
    return nothing
end

"""
    quarter_chord_AIC!(target, body_aero::BodyAerodynamics, core_radius_fraction, va_dist,
                       va_unit_dist)

Write into `target` the influence matrices of the flow at the aerodynamic centres that
turns the panel forces: the [`LLT`](@ref Model) matrices plus
[`curved_bound_correction!`](@ref), kept in `body_aero` for one shape and core radius.
"""
function quarter_chord_AIC!(target, body_aero::BodyAerodynamics, core_radius_fraction,
        va_dist, va_unit_dist)
    calculate_AIC_matrices!(body_aero, LLT, core_radius_fraction, va_dist, va_unit_dist,
        target)
    if body_aero.curved_correction_core_radius != core_radius_fraction
        curved_bound_correction!(body_aero.AIC_curved_correction, body_aero,
            core_radius_fraction, va_dist)
        body_aero.curved_correction_core_radius = core_radius_fraction
    end
    target .+= body_aero.AIC_curved_correction
    return nothing
end

"""
    update_effective_angle_of_attack!(alpha_corrected, body_aero::BodyAerodynamics, gamma,
                                      core_radius_fraction, va_vec_dist, va_dist,
                                      va_unit_dist)

Write each panel's [`inflow_angle`](@ref) at its aerodynamic centre into
`alpha_corrected`, the induced velocity taken from `gamma` on the
[`quarter_chord_AIC!`](@ref) matrices.
"""
function update_effective_angle_of_attack!(alpha_corrected,
    body_aero::BodyAerodynamics,
    gamma,
    core_radius_fraction,
    va_vec_dist,
    va_dist,
    va_unit_dist)

    # Its own buffer: `AIC` holds the control-point matrix the circulation was solved
    # against, so overwriting it here would leave post-solve readers on the other one.
    quarter_chord_AIC!(body_aero.AIC_aero_center, body_aero, core_radius_fraction, va_dist,
        va_unit_dist)

    induced_velocity = body_aero.cache[1][va_vec_dist]
    for k in 1:3
        mul!(view(induced_velocity, :, k), view(body_aero.AIC_aero_center, :, :, k), gamma)
    end

    @inbounds for (i, panel) in enumerate(body_aero.panels)
        v_rel = matrix_row(va_vec_dist, i) .+ matrix_row(induced_velocity, i)
        alpha_corrected[i] = inflow_angle(v_rel, SVector{3}(panel.y_airf),
                                          SVector{3}(panel.z_airf))
    end
    nothing
end

"""
    induced_velocity_at(body_aero::BodyAerodynamics, point, gamma, core_radius_fraction)

Velocity [m/s] induced at `point` by the horseshoe vortices of all panels, panel `j`
carrying circulation `gamma[j]` [m²/s] and its frozen wake.
"""
function induced_velocity_at(body_aero::BodyAerodynamics{P, W, T}, point, gamma,
        core_radius_fraction) where {P, W, T}
    velocity_ring = body_aero.work_vectors[8]
    velocity_filament = body_aero.work_vectors[9]
    velocity = zero(SVector{3, T})
    for (j, panel) in enumerate(body_aero.panels)
        wake = panel.filaments[4]
        calculate_velocity_induced_single_ring_semiinfinite!(velocity_ring,
            velocity_filament, panel.filaments, point, false, wake.va, wake.direction,
            gamma[j], core_radius_fraction, body_aero.work_vectors)
        velocity += SVector{3, T}(velocity_ring)
    end
    return velocity
end

"""
    attached_trailed_loads!(panel_forces, m_body_3D, body_aero::BodyAerodynamics, gamma,
                            density, core_radius_fraction, reference_point)

Write into the `3×n` matrix `panel_forces` [N] the Kutta–Joukowski force on the chordwise
trailed vortex legs between bound vortex and trailing edge, and add its moment [N·m] about
`reference_point` to `m_body_3D`. Each section boundary carries the net circulation of its
two panels and sees, at its three-quarter-chord point, their mean inflow plus the
[`induced_velocity_at`](@ref) that point; its load goes half to each panel, or whole to a
tip panel.
"""
function attached_trailed_loads!(panel_forces, m_body_3D,
        body_aero::BodyAerodynamics{P, W, T}, gamma, density, core_radius_fraction,
        reference_point) where {P, W, T}
    panel_forces .= 0
    three_quarter_chord = (0.75 - 0.25) / (1 - 0.25)
    panels = body_aero.panels
    for wing_idx in eachindex(body_aero.wings)
        wing_panels = panel_range(body_aero, wing_idx)
        for boundary in first(wing_panels):(last(wing_panels) + 1)
            is_first = boundary == first(wing_panels)
            is_last = boundary == last(wing_panels) + 1
            before = is_first ? boundary : boundary - 1
            after = is_last ? before : boundary
            bound_point = is_last ? panels[after].bound_point_2 : panels[after].bound_point_1
            te_point = is_last ? panels[after].TE_point_2 : panels[after].TE_point_1
            leg = SVector{3, T}(te_point) - SVector{3, T}(bound_point)
            point = SVector{3, T}(bound_point) + three_quarter_chord * leg
            velocity = (SVector{3, T}(panels[before].va_vec) +
                        SVector{3, T}(panels[after].va_vec)) / 2 +
                induced_velocity_at(body_aero, point, gamma, core_radius_fraction)
            net_gamma = (is_last ? zero(T) : gamma[after]) -
                (is_first ? zero(T) : gamma[before])
            force = (density * net_gamma) * cross(velocity, leg)
            moment = cross(point - SVector{3, T}(reference_point), force)
            share = before == after ? one(T) : T(0.5)
            @inbounds for k in 1:3
                panel_forces[k, before] += share * force[k]
                m_body_3D[k, before] += share * moment[k]
                before == after && continue
                panel_forces[k, after] += share * force[k]
                m_body_3D[k, after] += share * moment[k]
            end
        end
    end
    return nothing
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
    return SVector{3}(x_cp[1] + λ*f_unit[1],
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

"""
    find_center_of_pressure!(center_of_pressure, body_aero::BodyAerodynamics, force,
                             moment, reference_point; force_tol=1e-12)

Set `center_of_pressure` to the first point where the line of action of `force` and
`moment` about `reference_point` crosses a panel, or to `NaN` where it crosses none.
"""
function find_center_of_pressure!(
    center_of_pressure,
    body_aero::BodyAerodynamics,
    force,
    moment,
    reference_point;
    force_tol::Float64 = 1e-12
)
    F = force
    M0 = moment
    r0 = reference_point
    F_norm_sq = dot3(F, F)
    center_of_pressure .= NaN
    # Treat near-zero forces as "CoP undefined"
    if !(isfinite(F_norm_sq)) || F_norm_sq ≤ force_tol^2
        return center_of_pressure
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
                center_of_pressure .= intersection
                return center_of_pressure
            end
        end
    end
    return center_of_pressure
end

"""
    compute_panel_center_of_pressures!(panel_cp_locations, body_aero::BodyAerodynamics,
                                       f_distribution, m_distribution, reference_point)

Set each entry of `panel_cp_locations` to the point on its panel's chord, clamped between
leading and trailing edge, where the panel's column of `f_distribution` gives its column of
`m_distribution` about `reference_point`; to the aerodynamic center where that force has
no finite component normal to the chord.
"""
function compute_panel_center_of_pressures!(
    panel_cp_locations,
    body_aero::BodyAerodynamics,
    f_distribution::AbstractMatrix,
    m_distribution::AbstractMatrix,
    reference_point
)
    for i in eachindex(body_aero.panels)
        panel = body_aero.panels[i]
        @views F = f_distribution[:, i]
        @views M_ref = m_distribution[:, i]
        ac = panel.aero_center
        chord_dir = panel.x_airf
        span_dir = panel.y_airf
        c = panel.chord

        # Guard against non-finite forces and near-zero forces
        if !all(isfinite, F) || dot3(F, F) ≤ 1e-24
            panel_cp_locations[i] .= ac
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
            panel_cp_locations[i] .= ac
            continue
        end

        lever = clamp(m_pitch / F_perp_mag,
                      -0.25 * c, 0.75 * c)
        panel_cp_locations[i] .= (
            ac[1] + lever*chord_dir[1],
            ac[2] + lever*chord_dir[2],
            ac[3] + lever*chord_dir[3])
    end

    return panel_cp_locations
end

"""
    set_pitch_rate_dist!(body_aero, omega)

Fill `body_aero.pitch_rate_dist` from a rigid-body turn rate by projecting it
onto each panel's own spanwise axis. Panels with different dihedral see
different rates from the same `omega`.
"""
function set_pitch_rate_dist!(body_aero::BodyAerodynamics, omega)
    for (i, panel) in enumerate(body_aero.panels)
        body_aero.pitch_rate_dist[i] = dot3(omega, panel.y_airf)
    end
    return nothing
end

"""
    prescribed_va_directions(va, spanwise)

Lift and side unit vectors, `(; dir_lift, dir_side)`, of inflow `va` on a wing along
`spanwise`: lift normal to both, side normal to lift and `va`.
"""
@inline function prescribed_va_directions(va, spanwise)
    dir_lift = normalize(cross(va, spanwise))
    return (; dir_lift, dir_side=cross(dir_lift, va) / norm(va))
end

"""
    set_va!(body_aero::BodyAerodynamics, va_vec::VelVector, omega=zeros(MVec3);
            reference_point=body_aero.reference_point)

Set a uniform apparent wind and a body turn rate, and update the wake filaments. Each
panel sees `va_vec - omega × (control_point - reference_point)`.

# Arguments
- body_aero::BodyAerodynamics: The [`BodyAerodynamics`](@ref) struct to modify
- `va_vec::VelVector`: Velocity vector of the apparent wind speed       [m/s]
- `omega::VelVector`: Turn rate vector around x y and z axis            [rad/s]
- `reference_point`: Point the body turns about, stored on `body_aero`  [m]

`omega` is also projected onto each panel's spanwise axis into
`pitch_rate_dist`, which the solver reads when `flow_curvature` is enabled.
"""
function set_va!(body_aero::BodyAerodynamics{P, W, T}, va_vec::AbstractVector,
                 omega=zeros(MVector{3, T});
                 reference_point=body_aero.reference_point) where {P, W, T}
    body_aero.omega .= omega
    body_aero.reference_point .= reference_point
    set_pitch_rate_dist!(body_aero, omega)

    va_vec_dist = zeros(T, P, 3)
    for (i, panel) in enumerate(body_aero.panels)
        panel.va_vec .= va_vec .-
            omega × (panel.control_point .- body_aero.reference_point)
        va_vec_dist[i, :] .= panel.va_vec
    end

    # Update wake elements
    frozen_wake!(body_aero, va_vec_dist)
    getfield(body_aero, :va_vec) .= va_vec
    body_aero.has_distributed_va = false
    return nothing
end

"""
    set_va!(body_aero::BodyAerodynamics, va_vec_dist::AbstractMatrix;
            pitch_rate_dist=nothing)

Set a per-panel inflow distribution. `pitch_rate_dist` gives each panel's rotation
rate about its own spanwise axis [rad/s], positive nose-up; build it with
[`section_pitch_rate`](@ref) when the structure deforms, since twist and flapping
rates differ per section and no single body rate describes them. It is reset to
zero when omitted, because this method takes no `omega` and a stale one would
silently feed the `flow_curvature` moment.
"""
function set_va!(body_aero::BodyAerodynamics, va_vec_dist::AbstractMatrix;
                 pitch_rate_dist=nothing)
    size(va_vec_dist, 1) != length(body_aero.panels) &&
        throw(ArgumentError(
            "Number of rows in va_vec_dist should be equal to number of panels."))
    if isnothing(pitch_rate_dist)
        body_aero.pitch_rate_dist .= 0
    else
        length(pitch_rate_dist) != length(body_aero.panels) &&
            throw(ArgumentError("Length of pitch rate distribution should be equal to number of panels."))
        body_aero.pitch_rate_dist .= pitch_rate_dist
    end

    for (i, panel) in enumerate(body_aero.panels)
        panel.va_vec .= va_vec_dist[i, :]
    end

    # Update wake elements
    frozen_wake!(body_aero, va_vec_dist)
    getfield(body_aero, :va_vec) .= [mean(va_vec_dist[:,i]) for i in 1:3]
    body_aero.has_distributed_va = true
    return nothing
end

"""
    apparent_wind(alpha, beta, va)

Apparent wind vector in the body frame [m/s] at angle of attack `alpha` [rad], sideslip
`beta` [rad] and apparent wind speed `va` [m/s].
"""
apparent_wind(alpha, beta, va) =
    va .* [cos(alpha) * cos(beta), sin(beta), sin(alpha) * cos(beta)]

"""
    set_va!(body_aero::BodyAerodynamics, settings::VSMSettings)

Set the uniform inflow of `body_aero` to the [`apparent_wind`](@ref) at the `alpha` and
`beta` [°] and apparent wind speed `va` [m/s] of `settings.condition`, turning the body
about `body_aero.reference_point` at its `yaw_rate` [°/s] about Z_b.

# Example
```julia
settings = VSMSettings("path/to/settings.yaml")
body_aero = BodyAerodynamics([wing])
set_va!(body_aero, settings)
```
"""
function set_va!(body_aero::BodyAerodynamics, settings::VSMSettings)
    condition = settings.condition
    va_vec = apparent_wind(deg2rad(condition.alpha), deg2rad(condition.beta),
        condition.va)
    set_va!(body_aero, va_vec, [0.0, 0.0, deg2rad(condition.yaw_rate)])
end

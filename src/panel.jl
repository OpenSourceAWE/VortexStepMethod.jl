"""
    @with_kw mutable struct Panel

Represents a panel in a vortex step method simulation. All points and vectors are in the kite body (KB) frame.

# Fields
- `TE_point_1`::MVec3=zeros(MVec3): First trailing edge point
- `LE_point_1`::MVec3=zeros(MVec3): First leading edge point
- `TE_point_2`::MVec3=zeros(MVec3): Second trailing edge point
- `LE_point_2`::MVec3=zeros(MVec3): Second leading edge point
- `chord`::Float64=0: Panel chord length
- `va`::MVec3=zeros(MVec3): Panel velocity
- `corner_points`::MMatrix{3, 4, Float64}=zeros(MMatrix{3, 4, Float64}: Panel corner points
- `aero_model`::AeroModel=INVISCID: Aerodynamic model type [AeroModel](@ref)
- `aero_center::Vector{Float64}`: Panel aerodynamic center
- cl_coeffs::Vector{Float64}=zeros(Float64, 3)
- cd_coeffs::Vector{Float64}=zeros(Float64, 3)
- cm_coeffs::Vector{Float64}=zeros(Float64, 3)
- cl_interp::CL = nothing: lift interpolation (its type is a struct parameter)
- cd_interp::CD = nothing: drag interpolation
- cm_interp::CM = nothing: moment interpolation
- section_aero::SA = nothing: optional surface aero table, see [SectionAero](@ref)
- `control_point`::Vector{MVec3}: Panel control point
- `bound_point_1`::Vector{MVec3}: First bound point
- `bound_point_2`::Vector{MVec3}: Second bound point
- `x_airf`::MVec3=zeros(MVec3): Unit vector tangential to chord line
- `y_airf`::MVec3=zeros(MVec3): Unit vector in spanwise direction
- `z_airf`::MVec3=zeros(MVec3): Unit vector, cross product of x_airf and y_airf
- `width`::Float64=0: Panel width
- filaments::Tuple{BoundFilament,BoundFilament,BoundFilament,SemiInfiniteFilament,SemiInfiniteFilament} = (
        BoundFilament(),
        BoundFilament(),
        BoundFilament(),
        SemiInfiniteFilament(),
        SemiInfiniteFilament()
    ): Panel filaments, see: [BoundFilament](@ref)
- `delta`::T=0: flap trailing-edge deflection [rad]
- `crease_frac`::T=0: chordwise flap-hinge fraction (0–1); 0 disables the plate kink
"""
@with_kw mutable struct Panel{T, CL, CD, CM, SA}
    TE_point_1::MVector{3, T} = zeros(MVector{3, T})
    LE_point_1::MVector{3, T} = zeros(MVector{3, T})
    TE_point_2::MVector{3, T} = zeros(MVector{3, T})
    LE_point_2::MVector{3, T} = zeros(MVector{3, T})
    chord::T = zero(T)
    va::MVector{3, T} = zeros(MVector{3, T})
    corner_points::MMatrix{3, 4, T, 12} = zeros(MMatrix{3, 4, T, 12})
    aero_model::AeroModel = INVISCID
    cl_coeffs::Vector{Float64} = zeros(Float64, 3)
    cd_coeffs::Vector{Float64} = zeros(Float64, 3)
    cm_coeffs::Vector{Float64} = zeros(Float64, 3)
    cl_interp::CL = nothing
    cd_interp::CD = nothing
    cm_interp::CM = nothing
    section_aero::SA = nothing
    aero_center::MVector{3, T} = zeros(MVector{3, T})
    control_point::MVector{3, T} = zeros(MVector{3, T})
    bound_point_1::MVector{3, T} = zeros(MVector{3, T})
    bound_point_2::MVector{3, T} = zeros(MVector{3, T})
    x_airf::MVector{3, T} = zeros(MVector{3, T})
    y_airf::MVector{3, T} = zeros(MVector{3, T})
    z_airf::MVector{3, T} = zeros(MVector{3, T})
    width::T = zero(T)
    filaments::Tuple{BoundFilament{T},BoundFilament{T},BoundFilament{T},SemiInfiniteFilament{T},SemiInfiniteFilament{T}} = (
        BoundFilament{T}(),
        BoundFilament{T}(),
        BoundFilament{T}(),
        SemiInfiniteFilament{T}(),
        SemiInfiniteFilament{T}()
    )
    delta::T = zero(T)
    crease_frac::T = zero(T)
end

"""
    Panel{T}(; kwargs...)

Construct an empty (INVISCID) panel: the interpolation fields default to `nothing`, so
their type parameters are `Nothing`. Panels backed by polars are built with concrete
interpolation parameters via [`panel_interp_types`](@ref).
"""
Panel{T}(; kwargs...) where {T} = Panel{T, Nothing, Nothing, Nothing, Nothing}(; kwargs...)

function init_pos!(
    panel::Panel{T},
    section_1::Section,
    section_2::Section,
    aero_center,
    control_point,
    bound_point_1,
    bound_point_2,
    x_airf,
    y_airf,
    z_airf,
    delta,
    vec::AbstractVector{T}
) where T
    # Initialize basic geometry
    panel.TE_point_1 .= section_1.TE_point
    panel.LE_point_1 .= section_1.LE_point
    panel.TE_point_2 .= section_2.TE_point
    panel.LE_point_2 .= section_2.LE_point
    panel.chord = panel_chord(
        SVector{3}(section_1.LE_point), SVector{3}(section_1.TE_point),
        SVector{3}(section_2.LE_point), SVector{3}(section_2.TE_point))
    panel.corner_points[:, 1] = panel.LE_point_1
    panel.corner_points[:, 2] = panel.TE_point_1
    panel.corner_points[:, 3] = panel.TE_point_2
    panel.corner_points[:, 4] = panel.LE_point_2
    vec .= bound_point_2 .- bound_point_1
    panel.width = smooth_norm(vec)
    reinit!(panel.filaments[1], bound_point_2, bound_point_1, vec)
    reinit!(panel.filaments[2], bound_point_1, panel.TE_point_1, vec)
    reinit!(panel.filaments[3], panel.TE_point_2, bound_point_2, vec)

    panel.bound_point_1 .= bound_point_1
    panel.bound_point_2 .= bound_point_2
    panel.aero_center .= aero_center
    panel.control_point .= control_point
    panel.x_airf .= x_airf
    panel.y_airf .= y_airf
    panel.z_airf .= z_airf
    panel.delta = delta
    return nothing
end

"""
    build_interps(section_1, section_2, remove_nan) -> (cl, cd, cm, section_aero)

Build the averaged aerodynamic interpolations for the panel between two sections.
Returns `(cl_interp, cd_interp, cm_interp, section_aero)`, each `nothing` for models
that do not use it (INVISCID, POLY). `cl`/`cm` clamp (`Flat`) past the alpha range but
extrapolate linearly (`Line`) over delta; `cd` extrapolates linearly in both. The
concrete return types parameterise [`Panel`](@ref) — see [`panel_interp_types`](@ref).
"""
function build_interps(section_1::Section, section_2::Section, remove_nan)
    am = section_1.aero_model
    cl_i = cd_i = cm_i = nothing
    if am in (POLAR_VECTORS, POLAR_MATRICES)
        extrap_flat = remove_nan ? Flat() : NaN
        extrap_line = remove_nan ? Line() : NaN
        aero_1 = section_1.aero_data
        aero_2 = section_2.aero_data
        if am == POLAR_VECTORS
            (aero_1 isa Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}} &&
             aero_2 isa Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}) ||
                throw(ArgumentError("POLAR_VECTORS requires aero_data = (alpha, cl, cd, cm) vectors."))
            all(size.(aero_1) .== size.(aero_2)) ||
                throw(ArgumentError("Polar data must have same shape"))
            (length(aero_1[1]) == length(aero_2[1]) &&
             all(isapprox.(diff(aero_1[1]), diff(aero_2[1])))) ||
                throw(ArgumentError("Alpha steps must be identical."))
            alphas = Vector{Float64}(aero_1[1])
            cl = Vector{Float64}((aero_1[2] .+ aero_2[2]) ./ 2)
            cd = Vector{Float64}((aero_1[3] .+ aero_2[3]) ./ 2)
            cm = Vector{Float64}((aero_1[4] .+ aero_2[4]) ./ 2)
            cl_i = linear_interpolation(alphas, cl; extrapolation_bc=extrap_flat)
            cd_i = linear_interpolation(alphas, cd; extrapolation_bc=extrap_line)
            cm_i = linear_interpolation(alphas, cm; extrapolation_bc=extrap_flat)
        else
            (aero_1 isa Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}} &&
             aero_2 isa Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}) ||
                throw(ArgumentError("POLAR_MATRICES requires aero_data = (alpha, delta, cl, cd, cm)."))
            all(size.(aero_1) .== size.(aero_2)) ||
                throw(ArgumentError("Polar data must have same shape"))
            (length(aero_1[1]) == length(aero_2[1]) &&
             all(isapprox.(diff(aero_1[1]), diff(aero_2[1])))) ||
                throw(ArgumentError("Alpha steps must be identical."))
            (length(aero_1[2]) == length(aero_2[2]) &&
             all(isapprox.(diff(aero_1[2]), diff(aero_2[2])))) ||
                throw(ArgumentError("Delta steps must be identical."))
            alphas = Vector{Float64}(aero_1[1])
            deltas = Vector{Float64}(aero_1[2])
            cl = Matrix{Float64}((aero_1[3] .+ aero_2[3]) ./ 2)
            cd = Matrix{Float64}((aero_1[4] .+ aero_2[4]) ./ 2)
            cm = Matrix{Float64}((aero_1[5] .+ aero_2[5]) ./ 2)
            cl_bc = remove_nan ? (extrap_flat, extrap_line) : extrap_flat
            cm_bc = remove_nan ? (extrap_flat, extrap_line) : extrap_flat
            cl_i = linear_interpolation((alphas, deltas), cl; extrapolation_bc=cl_bc)
            cd_i = linear_interpolation((alphas, deltas), cd; extrapolation_bc=extrap_line)
            cm_i = linear_interpolation((alphas, deltas), cm; extrapolation_bc=cm_bc)
        end
    end
    aero = nothing
    if section_1.section_aero !== nothing
        a1, a2 = section_1.section_aero, section_2.section_aero
        aero = SectionAero(a1.alpha_range, a1.delta_range,
            (a1.x .+ a2.x) ./ 2, (a1.y .+ a2.y) ./ 2,
            (a1.cp .+ a2.cp) ./ 2, (a1.cf .+ a2.cf) ./ 2)
    end
    return cl_i, cd_i, cm_i, aero
end

"""
    panel_interp_types(section, remove_nan) -> (CL, CD, CM, CP)

Field types for a [`Panel`](@ref) whose aero model matches `section`. Since a wing has
one aero model, every panel shares these, so `Vector{Panel{...}}` is concretely typed.
"""
function panel_interp_types(section::Section, remove_nan)
    cl, cd, cm, cp = build_interps(section, section, remove_nan)
    return (Union{Nothing, typeof(cl)}, Union{Nothing, typeof(cd)},
            Union{Nothing, typeof(cm)}, Union{Nothing, typeof(cp)})
end

function init_aero!(panel::Panel, section_1::Section, section_2::Section;
                    remove_nan = true)
    panel.aero_model = section_1.aero_model
    section_1.aero_model == section_2.aero_model ||
        throw(ArgumentError("Both sections must have the same aero model, not " *
            "$(section_1.aero_model) and $(section_2.aero_model)"))
    if panel.aero_model == POLY
        c1, c2 = section_1.aero_data, section_2.aero_data
        (c1 isa NTuple{3, Vector{Float64}} && c2 isa NTuple{3, Vector{Float64}}) ||
            throw(ArgumentError("POLY requires aero_data = (cl_coeffs, cd_coeffs, cm_coeffs)."))
        panel.cl_coeffs = (c1[1] .+ c2[1]) ./ 2
        panel.cd_coeffs = (c1[2] .+ c2[2]) ./ 2
        panel.cm_coeffs = (c1[3] .+ c2[3]) ./ 2
    elseif !(panel.aero_model in (POLAR_VECTORS, POLAR_MATRICES, INVISCID))
        throw(ArgumentError("Unsupported aero model: $(panel.aero_model)"))
    end
    panel.cl_interp, panel.cd_interp, panel.cm_interp, panel.section_aero =
        build_interps(section_1, section_2, remove_nan)
    return nothing
end

"""
    reinit!(panel, section_1, section_2, aero_center, control_point, bound_point_1,
            bound_point_2, x_airf, y_airf, z_airf, delta, vec; kwargs...)

Reinitialize a panel's geometry, horseshoe filaments and aerodynamic interpolations.

`flip` reverses the section order so `y_airf` points along `+spanwise_direction` and
`z_airf` to the airfoil upper surface, making the aero independent of section ordering.
The caller owns it and must not derive it from the live geometry, which would invert
normals mid-run; [`reinit!(::BodyAerodynamics)`](@ref) decides it once per wing.
"""
function reinit!(
    panel::Panel,
    section_1::Section,
    section_2::Section,
    aero_center,
    control_point,
    bound_point_1,
    bound_point_2,
    x_airf,
    y_airf,
    z_airf,
    delta,
    vec;
    init_aero = true,
    remove_nan = true,
    flip::Bool = false
)
    if flip
        section_1, section_2 = section_2, section_1
        bound_point_1, bound_point_2 = bound_point_2, bound_point_1
    end
    init_pos!(panel, section_1, section_2, aero_center, control_point, bound_point_1, bound_point_2,
        x_airf, y_airf, z_airf, delta, vec)
    if flip
        panel.y_airf .*= -1
        panel.z_airf .*= -1
    end
    init_aero && init_aero!(panel, section_1, section_2; remove_nan)
    return nothing
end


"""
    panel_axes(panel::Panel)

The panel's stored airfoil frame and size in the shape [`panel_axes`](@ref)
returns, so a panel built by the geometry pass feeds the same functions as one built
straight from section points.
"""
@inline panel_axes(panel::Panel) =
    (x_airf = SVector{3}(panel.x_airf), y_airf = SVector{3}(panel.y_airf),
     z_airf = SVector{3}(panel.z_airf), chord = panel.chord,
     width = panel.width)

"""
    calculate_relative_alpha_and_relative_velocity(panel::Panel, induced_velocity::Vector{Float64})

Calculate the relative angle of attack and relative velocity of the panel.

# Arguments
- `panel::Panel`: The panel object
- `induced_velocity::Vector{Float64}`: Induced velocity at the control point

# Returns
- `Tuple{Float64,Vector{Float64}}`: Tuple containing:
  - alpha: Relative angle of attack of the panel (in radians)
  - relative_velocity: Relative velocity vector of the panel
"""
function calculate_relative_alpha_and_relative_velocity(
    panel::Panel{T},
    induced_velocity::AbstractVector{T}
) where T
    flow = panel_inflow(panel_axes(panel), panel.va, panel.va, induced_velocity)
    return flow.alpha, flow.v_eff
end

"""
    calculate_relative_alpha_and_velocity(panel::Panel, induced_velocity)

Calculate relative angle of attack and relative velocity of the panel.
"""
function calculate_relative_alpha_and_velocity(panel::Panel, induced_velocity)
    flow = panel_inflow(panel_axes(panel), panel.va, panel.va, induced_velocity)
    return flow.alpha, flow.v_eff
end

"""
    calculate_cl(panel::Panel, alpha)
    calculate_cl(panel::Panel, alpha, delta)

Calculate lift coefficient for given angle of attack `alpha` [rad]. The 3-arg
form evaluates the `(α, δ)` polar (`POLAR_MATRICES`) at the passed flap
deflection `delta` [rad] instead of the panel's stored `delta`; the 2-arg form
forwards with `panel.delta`. Other aero models ignore `delta`.

# Returns
- `Float64`: Lift coefficient (Cl)
"""
calculate_cl(panel::Panel, alpha) = calculate_cl(panel, alpha, panel.delta)
function calculate_cl(panel::Panel{Tp}, alpha::Ta, delta::Td) where {Tp, Ta, Td}
    R = promote_type(Tp, Ta, Td)
    isnan(alpha) && return R(NaN)
    if panel.aero_model == POLY
        cl = evalpoly(rad2deg(alpha), panel.cl_coeffs)
        if abs(alpha) > (π/9)
            cl = 2 * cos(alpha) * sin(alpha)^2
        end
        return R(cl)
    elseif panel.aero_model == INVISCID
        return R(2π * alpha)
    end
    interp = panel.cl_interp
    interp === nothing &&
        throw(ArgumentError("cl_interp is not initialized for $(panel.aero_model)."))
    return panel.aero_model == POLAR_VECTORS ? R(interp(alpha)) :
                                               R(interp(alpha, delta))
end


"""
    calculate_cd(panel::Panel, alpha)
    calculate_cd(panel::Panel, alpha, delta)

Calculate the drag coefficient for the given angle of attack. The 3-arg form
evaluates the `(α, δ)` polar at the passed flap deflection `delta`; see
[`calculate_cl`](@ref).
"""
calculate_cd(panel::Panel, alpha) = calculate_cd(panel, alpha, panel.delta)
function calculate_cd(panel::Panel{Tp}, alpha::Ta, delta::Td) where {Tp, Ta, Td}
    R = promote_type(Tp, Ta, Td)
    isnan(alpha) && return R(NaN)
    if panel.aero_model == POLY
        if abs(alpha) > (π/9)  # Outside ±20 degrees
            return R(2 * sin(alpha)^3)
        end
        return R(evalpoly(rad2deg(alpha), panel.cd_coeffs))
    elseif panel.aero_model in (POLAR_VECTORS, POLAR_MATRICES)
        cd_interp = panel.cd_interp
        cd_interp === nothing &&
            throw(ArgumentError("cd_interp is not initialized for $(panel.aero_model)."))
        return panel.aero_model == POLAR_VECTORS ? R(cd_interp(alpha)) :
                                                   R(cd_interp(alpha, delta))
    elseif !(panel.aero_model == INVISCID)
        throw(ArgumentError("Unsupported aero model: $(panel.aero_model)"))
    end
    return zero(R)
end

"""
    calculate_cm(panel::Panel, alpha)
    calculate_cm(panel::Panel, alpha, delta)

Calculate the pitching-moment coefficient for the given angle of attack. The
3-arg form evaluates the `(α, δ)` polar at the passed flap deflection `delta`;
see [`calculate_cl`](@ref).
"""
calculate_cm(panel::Panel, alpha) = calculate_cm(panel, alpha, panel.delta)
function calculate_cm(panel::Panel{Tp}, alpha::Ta, delta::Td) where {Tp, Ta, Td}
    R = promote_type(Tp, Ta, Td)
    isnan(alpha) && return R(NaN)
    if panel.aero_model == POLY
        return R(evalpoly(rad2deg(alpha), panel.cm_coeffs))
    elseif panel.aero_model in (POLAR_VECTORS, POLAR_MATRICES)
        cm_interp = panel.cm_interp
        cm_interp === nothing &&
            throw(ArgumentError("cm_interp is not initialized for $(panel.aero_model)."))
        return panel.aero_model == POLAR_VECTORS ? R(cm_interp(alpha)) :
                                                   R(cm_interp(alpha, delta))
    elseif !(panel.aero_model == INVISCID)
        throw(ArgumentError("Unsupported aero model: $(panel.aero_model)"))
    end
    return zero(R)
end

"""
    calculate_cd_cm(panel::Panel, alpha)

Calculate drag and moment coefficients for the given angle of attack
([`calculate_cd`](@ref) and [`calculate_cm`](@ref)).
"""
calculate_cd_cm(panel::Panel{Tp}, alpha::Ta) where {Tp, Ta} =
    (calculate_cd(panel, alpha), calculate_cm(panel, alpha))

"""
    calculate_filaments_for_plotting(panel::Panel)

Calculate filaments for plotting with their positions and colors.

# Returns
- `Vector{Tuple{Vector{Float64}, Vector{Float64}, String}}`: List of tuples containing:
  - First point (x1)
  - Second point (x2)
  - Color string
"""
function calculate_filaments_for_plotting(panel::Panel)
    filaments_plot = []
    
    for (i, filament) in enumerate(panel.filaments)
        x1 = filament.x1
        
        if isdefined(filament, :x2) && !isnothing(filament.x2)
            x2 = filament.x2
            # Color based on filament type
            color = i == 1 ? "magenta" : "green"  # bound vs trailing
        else
            # For semi-infinite filaments
            x2 = x1 + 2 * panel.chord * (panel.va / norm(panel.va))
            color = "orange"
            
            if filament.filament_direction == -1
                x1, x2 = x2, x1  # swap points
                color = "red"
            end
        end
        
        push!(filaments_plot, (x1, x2, color))
    end
    
    return filaments_plot
end

"""
    calculate_velocity_induced_single_ring_semiinfinite!(
        velind::MVec3,
        tempvel::MVec3,
        filaments,
        evaluation_point::MVec3,
        evaluation_point_on_bound::Bool,
        va_norm::Float64,
        va_unit::MVec3,
        gamma::Float64,
        core_radius_fraction::Float64,
        work_vectors::NTuple{10, MVec3}
    )

Calculate the velocity induced by a vortex ring at a control point.

# Arguments
- velind
- tempvel
- filaments
- `evaluation_point`::MVec3:         Point where induced velocity is evaluated
- `evaluation_point_on_bound`::Bool: Whether evaluation point is on bound vortex
- `va_norm`::Float64:                Norm of apparent velocity
- `va_unit`::MVec3:                  Unit vector of apparent velocity
- `gamma`::Float64:                  Circulation strength
- `core_radius_fraction`::Float64:   Vortex core radius as fraction of panel width
- `work_vectors`::NTuple{10, MVec3}    Pre-allocated temporary variables

# Returns
- nothing
"""
@inline function calculate_velocity_induced_single_ring_semiinfinite!(
    velind::AbstractVector{T},
    tempvel::AbstractVector{T},
    filaments,
    evaluation_point::AbstractVector{T},
    evaluation_point_on_bound::Bool,
    va_norm::T,
    va_unit::AbstractVector{T},
    gamma::T,
    core_radius_fraction::Real,
    work_vectors
) where T
    velind .= 0.0

    # Filament 1: bound filament (BoundFilament — compiler knows type)
    f1 = filaments[1]
    f1.initialized || throw(ArgumentError(
        "Filament not initialized: 1. " *
        "Maybe you forgot to call set_va! before running solve."))
    if evaluation_point_on_bound
        tempvel .= 0.0
    else
        velocity_3D_bound_vortex!(
            tempvel, f1, evaluation_point, gamma,
            core_radius_fraction, work_vectors)
    end
    velind .+= tempvel

    # Filament 2: trailing filament (BoundFilament)
    f2 = filaments[2]
    f2.initialized || throw(ArgumentError(
        "Filament not initialized: 2."))
    velocity_3D_trailing_vortex!(
        tempvel, f2, evaluation_point, gamma,
        va_norm, work_vectors)
    velind .+= tempvel

    # Filament 3: trailing filament (BoundFilament)
    f3 = filaments[3]
    f3.initialized || throw(ArgumentError(
        "Filament not initialized: 3."))
    velocity_3D_trailing_vortex!(
        tempvel, f3, evaluation_point, gamma,
        va_norm, work_vectors)
    velind .+= tempvel

    # Filament 4: semi-infinite trailing (SemiInfiniteFilament)
    f4 = filaments[4]
    f4.initialized || throw(ArgumentError(
        "Filament not initialized: 4."))
    velocity_3D_trailing_vortex_semiinfinite!(
        tempvel, f4, va_unit, evaluation_point, gamma,
        va_norm, work_vectors)
    velind .+= tempvel

    # Filament 5: semi-infinite trailing (SemiInfiniteFilament)
    f5 = filaments[5]
    f5.initialized || throw(ArgumentError(
        "Filament not initialized: 5."))
    velocity_3D_trailing_vortex_semiinfinite!(
        tempvel, f5, va_unit, evaluation_point, gamma,
        va_norm, work_vectors)
    velind .+= tempvel

    return nothing
end

"""
    calculate_velocity_induced_bound_2D!(U2D, panel::Panel, evaluation_point, work_vectors)

Calculate velocity induced by bound vortex filaments at the control point.
Only needed for VSM, as LLT bound and filament align, thus no induced velocity.

# Arguments
- `U_2D`:             Resulting 2D velocity vector
- `panel::Panel`:     Panel object
- `evaluation_point`: Point where induced velocity is evaluated
- `work_vectors`:     Pre-allocated temporary variables
"""
function calculate_velocity_induced_bound_2D!(
    U_2D,
    panel::Panel,
    evaluation_point,
    work_vectors
)
    r3, r0, cross_ = work_vectors
    # r3 perpendicular to the bound vortex
    r3 .= evaluation_point .-
          (panel.bound_point_1 .+ panel.bound_point_2) ./ 2

    # r0 is the direction of the bound vortex
    r0 .= panel.bound_point_1 .- panel.bound_point_2

    # Calculate cross product
    cross3!(cross_, r0, r3)

    # Calculate induced velocity
    coeff = norm3(r0) / (2π * norm3(cross_)^2)
    @inbounds for k in 1:3
        U_2D[k] = cross_[k] * coeff
    end
    return nothing
end

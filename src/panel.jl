# static types for interpolations
const I1 = Interpolations.FilledExtrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Vector{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Float64}
const I2 = Interpolations.Extrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Vector{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Interpolations.Flat{Nothing}}
const I3 = Interpolations.FilledExtrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Float64}
const I4 = Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Interpolations.Flat{Nothing}}
# Line extrapolation types for cd
const I5 = Interpolations.Extrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Vector{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Interpolations.Line{Nothing}}
const I6 = Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Interpolations.Line{Nothing}}

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
- cl_interp::Union{Nothing, I1, I2} = nothing
- cd_interp::Union{Nothing, I1, I2} = nothing
- cm_interp::Union{Nothing, I1, I2} = nothing
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
"""
@with_kw mutable struct Panel{T}
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
    cl_interp::Union{Nothing, I1, I2, I3, I4} = nothing
    cd_interp::Union{Nothing, I1, I2, I3, I4, I5, I6} = nothing
    cm_interp::Union{Nothing, I1, I2, I3, I4} = nothing
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
end

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
    panel.chord = (
        norm(panel.TE_point_1 - panel.LE_point_1) +
        norm(panel.TE_point_2 - panel.LE_point_2)
    ) / 2
    panel.corner_points[:, 1] = panel.LE_point_1
    panel.corner_points[:, 2] = panel.TE_point_1
    panel.corner_points[:, 3] = panel.TE_point_2
    panel.corner_points[:, 4] = panel.LE_point_2
    vec .= bound_point_2 .- bound_point_1
    panel.width = norm(vec)
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

function init_aero!(
    panel::Panel,
    section_1::Section,
    section_2::Section;
    remove_nan = true
)
    # Validate aero model consistency
    panel.aero_model = section_1.aero_model
    aero_model_2 = section_2.aero_model
    if !(panel.aero_model == aero_model_2)
        throw(ArgumentError("Both sections must have the same aero model, not $(panel.aero_model) and $aero_model_2"))
    end
    
    if panel.aero_model == LEI_AIRFOIL_BREUKELS
        panel.cl_coeffs, panel.cd_coeffs, panel.cm_coeffs = compute_lei_coeffs(section_1, section_2)

    elseif panel.aero_model in (POLAR_VECTORS, POLAR_MATRICES)
        if remove_nan
            extrap_flat = Flat()
            extrap_line = Line()
        else
            extrap_flat = NaN
            extrap_line = NaN
        end

        if panel.aero_model == POLAR_VECTORS
            aero_1 = section_1.aero_data
            aero_2 = section_2.aero_data
            aero_1 isa Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}} ||
                throw(ArgumentError("POLAR_VECTORS requires aero_data = (alpha, cl, cd, cm) vectors."))
            aero_2 isa Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}} ||
                throw(ArgumentError("POLAR_VECTORS requires aero_data = (alpha, cl, cd, cm) vectors."))
            aero_1_t = aero_1::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
            aero_2_t = aero_2::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

            if !all(size.(aero_1_t) .== size.(aero_2_t))
                throw(ArgumentError("Polar data must have same shape"))
            end

            alphas_1 = aero_1_t[1]
            alphas_2 = aero_2_t[1]
            (
                length(alphas_1) == length(alphas_2) &&
                all(isapprox.(diff(alphas_1), diff(alphas_2)))
            ) || throw(ArgumentError("Alpha steps must be identical."))

            polar_data = (
                Vector{Float64}((aero_1_t[2] + aero_2_t[2]) / 2),
                Vector{Float64}((aero_1_t[3] + aero_2_t[3]) / 2),
                Vector{Float64}((aero_1_t[4] + aero_2_t[4]) / 2)
            )
            alphas = Vector{Float64}(alphas_1)

            panel.cl_interp = linear_interpolation(alphas, polar_data[1]; extrapolation_bc=extrap_flat)
            panel.cd_interp = linear_interpolation(alphas, polar_data[2]; extrapolation_bc=extrap_line)
            panel.cm_interp = linear_interpolation(alphas, polar_data[3]; extrapolation_bc=extrap_flat)

        elseif panel.aero_model == POLAR_MATRICES
            aero_1 = section_1.aero_data
            aero_2 = section_2.aero_data
            aero_1 isa Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}} ||
                throw(ArgumentError("POLAR_MATRICES requires aero_data = (alpha, delta, cl, cd, cm)."))
            aero_2 isa Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}} ||
                throw(ArgumentError("POLAR_MATRICES requires aero_data = (alpha, delta, cl, cd, cm)."))
            aero_1_t = aero_1::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
            aero_2_t = aero_2::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}

            if !all(size.(aero_1_t) .== size.(aero_2_t))
                throw(ArgumentError("Polar data must have same shape"))
            end

            alphas_1 = aero_1_t[1]
            alphas_2 = aero_2_t[1]
            deltas_1 = aero_1_t[2]
            deltas_2 = aero_2_t[2]
            (
                length(alphas_1) == length(alphas_2) &&
                all(isapprox.(diff(alphas_1), diff(alphas_2)))
            ) || throw(ArgumentError("Alpha steps must be identical."))
            (
                length(deltas_1) == length(deltas_2) &&
                all(isapprox.(diff(deltas_1), diff(deltas_2)))
            ) || throw(ArgumentError("Delta steps must be identical."))

            polar_data = (
                Matrix{Float64}((aero_1_t[3] + aero_2_t[3]) / 2),
                Matrix{Float64}((aero_1_t[4] + aero_2_t[4]) / 2),
                Matrix{Float64}((aero_1_t[5] + aero_2_t[5]) / 2)
            )
            alphas = Vector{Float64}(alphas_1)
            deltas = Vector{Float64}(deltas_1)

            panel.cl_interp = linear_interpolation((alphas, deltas), polar_data[1]; extrapolation_bc=extrap_flat)
            panel.cd_interp = linear_interpolation((alphas, deltas), polar_data[2]; extrapolation_bc=extrap_line)
            panel.cm_interp = linear_interpolation((alphas, deltas), polar_data[3]; extrapolation_bc=extrap_flat)
        else
            throw(ArgumentError("Polar data in wrong format for model $(panel.aero_model)."))
        end

    elseif !(panel.aero_model == INVISCID)
        throw(ArgumentError("Unsupported aero model: $(panel.aero_model)"))
    end
end

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
    remove_nan = true
)
    init_pos!(panel, section_1, section_2, aero_center, control_point, bound_point_1, bound_point_2,
        x_airf, y_airf, z_airf, delta, vec)
    init_aero && init_aero!(panel, section_1, section_2; remove_nan)
    return nothing
end


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
    # Calculate relative velocity and angle of attack
    # Constants throughout iterations: panel.va, panel.x_airf, panel.y_airf
    relative_velocity = panel.va .+ induced_velocity
    v_normal = dot(panel.z_airf, relative_velocity)
    v_tangential = dot(panel.x_airf, relative_velocity)
    alpha = atan(v_normal, v_tangential)
    
    return alpha, relative_velocity
end

"""
    compute_lei_coeffs(section_1::Section, section_2::Section)

Compute lift, drag and moment coefficients for Lei airfoil using Breukels model.
"""
function compute_lei_coeffs(section_1::Section, section_2::Section)
    section_1.aero_data isa NTuple{2, Float64} ||
        throw(ArgumentError("LEI_AIRFOIL_BREUKELS requires aero_data = (tube_diameter, camber)."))
    section_2.aero_data isa NTuple{2, Float64} ||
        throw(ArgumentError("LEI_AIRFOIL_BREUKELS requires aero_data = (tube_diameter, camber)."))

    # Average tube diameter and camber from both sections
    t1, k1 = section_1.aero_data::NTuple{2, Float64}
    t2, k2 = section_2.aero_data::NTuple{2, Float64}
    t = (t1 + t2) / 2
    k = (k1 + k2) / 2

    # Lift coefficient constants
    C = Dict(
        20 => -0.008011, 21 => -0.000336, 22 => 0.000992,
        23 => 0.013936, 24 => -0.003838, 25 => -0.000161,
        26 => 0.001243, 27 => -0.009288, 28 => -0.002124,
        29 => 0.012267, 30 => -0.002398, 31 => -0.000274,
        32 => 0.0, 33 => 0.0, 34 => 0.0,
        35 => -3.371000, 36 => 0.858039, 37 => 0.141600,
        38 => 7.201140, 39 => -0.676007, 40 => 0.806629,
        41 => 0.170454, 42 => -0.390563, 43 => 0.101966
    )

    # Compute S values
    S = Dict{Int64,Float64}()
    S[9] = C[20]*t^2 + C[21]*t + C[22]
    S[10] = C[23]*t^2 + C[24]*t + C[25]
    S[11] = C[26]*t^2 + C[27]*t + C[28]
    S[12] = C[29]*t^2 + C[30]*t + C[31]
    S[13] = C[32]*t^2 + C[33]*t + C[34]
    S[14] = C[35]*t^2 + C[36]*t + C[37]
    S[15] = C[38]*t^2 + C[39]*t + C[40]
    S[16] = C[41]*t^2 + C[42]*t + C[43]

    # Compute lambda values for cl
    λ = [
        S[9]*k + S[10],
        S[11]*k + S[12],
        S[13]*k + S[14],
        S[15]*k + S[16]
    ]

    # Drag coefficient constants and computation
    cd_coeffs = [
        ((0.546094*t + 0.022247)*k^2 + 
         (-0.071462*t - 0.006527)*k + 
         (0.002733*t + 0.000686)),
        0.0,
        ((0.123685*t + 0.143755)*k + 
         (0.495159*t^2 - 0.105362*t + 0.033468))
    ]

    # Moment coefficient constants and computation
    cm_coeffs = [
        ((-0.284793*t - 0.026199)*k + 
         (-0.024060*t + 0.000559)),
        0.0,
        ((-1.787703*t + 0.352443)*k + 
         (-0.839323*t + 0.137932))
    ]

    return λ, cd_coeffs, cm_coeffs
end

"""
    calculate_relative_alpha_and_velocity(panel::Panel, induced_velocity)

Calculate relative angle of attack and relative velocity of the panel.
"""
function calculate_relative_alpha_and_velocity(panel::Panel, induced_velocity)
    relative_velocity = panel.va + induced_velocity
    v_normal = dot(panel.z_airf, relative_velocity)
    v_tangential = dot(panel.x_airf, relative_velocity)
    alpha = atan(v_normal, v_tangential)
    return alpha, relative_velocity
end

"""
    calculate_cl(panel::Panel, alpha::Float64)

Calculate lift coefficient for given angle of attack.

# Arguments
- `panel::Panel`: Panel object
- `alpha::Float64`: Angle of attack in radians

# Returns
- `Float64`: Lift coefficient (Cl)
"""
function calculate_cl(panel::Panel{Tp}, alpha::Ta) where {Tp, Ta}
    R = promote_type(Tp, Ta)
    isnan(alpha) && return R(NaN)
    if panel.aero_model == LEI_AIRFOIL_BREUKELS
        cl = evalpoly(rad2deg(alpha), reverse(panel.cl_coeffs))
        if abs(alpha) > (π/9)
            cl = 2 * cos(alpha) * sin(alpha)^2
        end
        return R(cl)
    elseif panel.aero_model == INVISCID
        return R(2π * alpha)
    elseif panel.aero_model == POLAR_VECTORS
        interp = panel.cl_interp
        interp isa Union{I1, I2} || throw(ArgumentError("cl_interp is not initialized for POLAR_VECTORS."))
        return R((interp::Union{I1, I2})(alpha))
    elseif panel.aero_model == POLAR_MATRICES
        interp = panel.cl_interp
        interp isa Union{I3, I4} || throw(ArgumentError("cl_interp is not initialized for POLAR_MATRICES."))
        return R((interp::Union{I3, I4})(alpha, panel.delta))
    else
        throw(ArgumentError("Unsupported aero model: $(panel.aero_model)"))
    end
end


"""
    calculate_cd_cm(panel::Panel, alpha::Float64)

Calculate drag and moment coefficients for given angle of attack.
"""
function calculate_cd_cm(panel::Panel{Tp}, alpha::Ta) where {Tp, Ta}
    R = promote_type(Tp, Ta)
    isnan(alpha) && return R(NaN), R(NaN)
    if panel.aero_model == LEI_AIRFOIL_BREUKELS
        cd = evalpoly(rad2deg(alpha), reverse(panel.cd_coeffs))
        cm = evalpoly(rad2deg(alpha), reverse(panel.cm_coeffs))
        if abs(alpha) > (π/9)  # Outside ±20 degrees
            cd = 2 * sin(alpha)^3
        end
        return R(cd), R(cm)
    elseif panel.aero_model == POLAR_VECTORS
        cd_interp = panel.cd_interp
        cm_interp = panel.cm_interp
        cd_interp isa Union{I1, I2, I5} || throw(ArgumentError("cd_interp is not initialized for POLAR_VECTORS."))
        cm_interp isa Union{I1, I2} || throw(ArgumentError("cm_interp is not initialized for POLAR_VECTORS."))
        return R((cd_interp::Union{I1, I2, I5})(alpha)),
               R((cm_interp::Union{I1, I2})(alpha))
    elseif panel.aero_model == POLAR_MATRICES
        cd_interp = panel.cd_interp
        cm_interp = panel.cm_interp
        cd_interp isa Union{I3, I4, I6} || throw(ArgumentError("cd_interp is not initialized for POLAR_MATRICES."))
        cm_interp isa Union{I3, I4} || throw(ArgumentError("cm_interp is not initialized for POLAR_MATRICES."))
        return R((cd_interp::Union{I3, I4, I6})(alpha, panel.delta)),
               R((cm_interp::Union{I3, I4})(alpha, panel.delta))
    elseif !(panel.aero_model == INVISCID)
        throw(ArgumentError("Unsupported aero model: $(panel.aero_model)"))
    end
    return zero(R), zero(R)
end

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
    coeff = smooth_norm3(r0) / (2π * smooth_norm3(cross_)^2)
    @inbounds for k in 1:3
        U_2D[k] = cross_[k] * coeff
    end
    return nothing
end

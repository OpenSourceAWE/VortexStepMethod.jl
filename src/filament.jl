
"""
Abstract type for vortex filaments
"""
abstract type Filament end

# Constants for all filament types
const ALPHA0 = 1.25643  # Oseen parameter
const NU = 1.48e-5     # Kinematic viscosity of air

"""
    BoundFilament

Represents a bound vortex filament defined by two points.

# Fields
- x1::MVec3=zeros(MVec3): First point
- x2::MVec3=zeros(MVec3): Second point
- length=zero(Float64):   Filament length
- r0::MVec3=zeros(MVec3): Vector from x1 to x2
- initialized::Bool = false
"""
@with_kw mutable struct BoundFilament <: Filament
    x1::MVec3         = zeros(MVec3)
    x2::MVec3         = zeros(MVec3)
    length::Float64   = zero(Float64)
    r0::MVec3         = zeros(MVec3)
    initialized::Bool = false
end

function init!(filament::BoundFilament, x1::PosVector, x2::PosVector, vec=zeros(MVec3))
    filament.x1 .= x1
    filament.x2 .= x2
    vec .= x2 .- x1
    filament.length = norm(vec)
    filament.r0 .= x2 .- x1
    filament.initialized = true
    return nothing
end

"""
    velocity_3D_bound_vortex(vel, filament::BoundFilament, XVP, 
                           gamma, core_radius_fraction, work_vectors)

Calculate induced velocity by a bound vortex filament at a point in space.
"""
function velocity_3D_bound_vortex!(
    vel,
    filament::BoundFilament,
    XVP,
    gamma,
    core_radius_fraction,
    work_vectors
)
    r1, r2, r1Xr2, r1Xr0, r2Xr0, r1r2norm, r1_proj, r2_proj, r1_projXr2_proj, vel_ind_proj = work_vectors
    r0 = filament.r0
    r1 .= XVP .- filament.x1
    r2 .= XVP .- filament.x2

    # Cut-off radius
    epsilon = core_radius_fraction * norm(r0)

    cross3!(r1Xr2, r1, r2)
    cross3!(r1Xr0, r1, r0)
    r1r2norm .= r1./norm(r1) .- r2./norm(r2)
    
    # Check point location relative to filament
    if norm(r1Xr0) / norm(r0) > epsilon
        vel .= (gamma / (4π)) .* r1Xr2 ./ (norm(r1Xr2)^2) .* 
            dot(r0, r1r2norm)
    elseif norm(r1Xr0) / norm(r0) == 0
        vel .= 0.0
    else
        @debug "inside core radius"
        @debug "distance from control point to filament: $(norm(r1Xr0) / norm(r0))"
        
        # Project onto core radius
        cross3!(r2Xr0, r2, r0)
        r1_proj .= dot(r1, r0) .* r0 ./ (norm(r0)^2) .+ 
                  epsilon .* r1Xr0 ./ norm(r1Xr0)
        r2_proj .= dot(r2, r0) .* r0 ./ (norm(r0)^2) .+ 
                  epsilon .* r2Xr0 ./ norm(r2Xr0)
        cross3!(r1_projXr2_proj, r1_proj, r2_proj)

        vel_ind_proj .= (gamma / (4π)) .* r1_projXr2_proj ./ (norm(r1_projXr2_proj)^2) .* 
                      dot(r0, r1_proj/norm(r1_proj) .- r2_proj/norm(r2_proj))

        vel .= norm(r1Xr0) ./ (norm(r0) * epsilon) .* vel_ind_proj
    end
    nothing
end

"""
    velocity_3D_trailing_vortex(vel, filament::BoundFilament, 
                              XVP, gamma, v_a, work_vectors)

Calculate induced velocity by a trailing vortex filament.

# Arguments
- `XVP`: Control point coordinates
- `gamma`: Vortex strength
- `v_a`: Inflow velocity magnitude
- work_vectors: preallocated array of intermediate variables

Reference: Rick Damiani et al. "A vortex step method for nonlinear airfoil polar data 
as implemented in KiteAeroDyn".
"""
@inline function velocity_3D_trailing_vortex!(
    vel,
    filament::BoundFilament,
    XVP,
    gamma,
    v_a,
    work_vectors
)
    r0, r1, r2, r_perp, r1Xr2, r1Xr0, r2Xr0, normr1r2 = work_vectors[1:8]
    r0 .= filament.x2 .- filament.x1  # Vortex filament
    r1 .= XVP .- filament.x1          # Control point to first end
    r2 .= XVP .- filament.x2          # Control point to second end

    # Vector perpendicular to core radius
    r_perp .= dot(r1, r0) .* r0 ./ (norm(r0)^2)
    
    # Cut-off radius
    epsilon = sqrt(4 * ALPHA0 * NU * norm(r_perp) / v_a)

    cross3!(r1Xr2, r1, r2)
    cross3!(r1Xr0, r1, r0)
    cross3!(r2Xr0, r2, r0)

    normr1r2 .= (r1./norm(r1)) .- (r2./norm(r2))

    # Check point location relative to filament
    if norm(r1Xr0) / norm(r0) > epsilon
        vel .= (gamma / (4π)) .* r1Xr2 ./ (norm(r1Xr2)^2) .* 
            dot(r0, normr1r2)
    elseif norm(r1Xr0) / norm(r0) == 0
        vel .= 0.0
    else
        # Project onto core radius
        r1_proj = dot(r1, r0) * r0 / (norm(r0)^2) + 
                  epsilon * r1Xr0 / norm(r1Xr0)
        r2_proj = dot(r2, r0) * r0 / (norm(r0)^2) + 
                  epsilon * r2Xr0 / norm(r2Xr0)
        
        cross3!(r1Xr2, r1_proj, r2_proj)

        vel_ind_proj = (gamma / (4π)) * r1Xr2 / (norm(r1Xr2)^2) * 
                      dot(r0, r1_proj/norm(r1_proj) - r2_proj/norm(r2_proj))
        
        vel .= norm(r1Xr0) / (norm(r0) * epsilon) * vel_ind_proj
    end
    nothing
end

"""
    SemiInfiniteFilament

Represents a semi-infinite vortex filament.

# Fields
- x1::MVec3=zeros(MVec3):           Starting point
- direction::MVec3=zeros(MVec3):    Direction vector
- `vel_mag`::Float64=zero(Float64): Velocity magnitude
- `filament_direction`::Int64=0   : Direction indicator (-1 or 1)
- initialized::Bool=false
"""
@with_kw mutable struct SemiInfiniteFilament <: Filament
    x1::MVec3 = zeros(MVec3)
    direction::MVec3 = zeros(MVec3)
    vel_mag::Float64 = zero(Float64)
    filament_direction::Int64 = zero(Int64)
    initialized::Bool = false
end

function init!(filament::SemiInfiniteFilament, x1::PosVector, direction::PosVector, vel_mag::Real, filament_direction::Real)
    filament.x1 .= x1
    filament.direction .= direction
    filament.vel_mag = vel_mag
    filament.filament_direction = filament_direction
    filament.initialized = true
    return nothing
end

"""
    velocity_3D_trailing_vortex_semiinfinite(filament::SemiInfiniteFilament, 
                                             Vf, XVP, GAMMA, v_a, work_vectors)

Calculate induced velocity by a semi-infinite trailing vortex filament.
"""
function velocity_3D_trailing_vortex_semiinfinite!(
    vel,
    filament::SemiInfiniteFilament,
    Vf,
    XVP,
    GAMMA,
    v_a,
    work_vectors
)
    r1, r_perp, r1XVf = work_vectors[1:3]
    GAMMA = -GAMMA * filament.filament_direction
    r1 .= XVP .- filament.x1

    # Calculate core radius
    r_perp .= dot(r1, Vf) .* Vf
    epsilon = sqrt(4 * ALPHA0 * NU * norm(r_perp) / v_a)

    cross3!(r1XVf, r1, Vf)

    if norm(r1XVf) / norm(Vf) > epsilon
        K = GAMMA / (4π) / norm(r1XVf)^2 * (1 + dot(r1, Vf) / norm(r1))
        vel .= K .* r1XVf
    elseif norm(r1XVf) / norm(Vf) == 0
        vel .= 0.0
    else
        r1_proj = dot(r1, Vf) * Vf + 
                  epsilon * (r1/norm(r1) - Vf) / norm(r1/norm(r1) - Vf)
        K = GAMMA / (4π) / norm(cross(r1_proj, Vf))^2 * 
            (1 + dot(r1_proj, Vf) / norm(r1_proj))
        vel .= K * cross(r1_proj, Vf)
    end
    nothing
end


"""
    cross3!(result::AbstractVector{T}, a::AbstractVector{T}, b::AbstractVector{T}) where T

Compute cross product of 3D vectors in-place.
"""
@inline function cross3!(result::AbstractVector{T}, a::AbstractVector{T}, b::AbstractVector{T}) where T
    x = a[2]*b[3] - a[3]*b[2]
    y = a[3]*b[1] - a[1]*b[3]
    z = a[1]*b[2] - a[2]*b[1]
    result[1] = x
    result[2] = y
    result[3] = z
    nothing
end
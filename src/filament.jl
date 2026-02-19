
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

function reinit!(filament::BoundFilament, x1, x2, vec=zeros(MVec3))
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
    r1, r2, r1Xr2, r1Xr0, r2Xr0, r1r2norm, r1_proj, r2_proj,
        r1_projXr2_proj, vel_ind_proj = work_vectors
    r0 = filament.r0
    r1 .= XVP .- filament.x1
    r2 .= XVP .- filament.x2

    # Cut-off radius
    nr0 = norm3(r0)
    epsilon = core_radius_fraction * nr0

    cross3!(r1Xr2, r1, r2)
    cross3!(r1Xr0, r1, r0)
    nr1 = norm3(r1)
    nr2 = norm3(r2)
    @inbounds for k in 1:3
        r1r2norm[k] = r1[k]/nr1 - r2[k]/nr2
    end

    # Check point location relative to filament
    nr1Xr0 = norm3(r1Xr0)
    if nr1Xr0 / nr0 > epsilon
        nr1Xr2 = norm3(r1Xr2)
        coeff = (gamma / (4π)) / (nr1Xr2^2) * dot3(r0, r1r2norm)
        @inbounds for k in 1:3
            vel[k] = coeff * r1Xr2[k]
        end
    elseif nr1Xr0 / nr0 < 1e-12 * epsilon
        vel .= 0.0
    else
        @debug "inside core radius"
        @debug "distance from control point to filament: $(nr1Xr0 / nr0)"

        # Project onto core radius
        cross3!(r2Xr0, r2, r0)
        nr0sq = nr0 * nr0
        nr2Xr0 = norm3(r2Xr0)
        d_r1_r0 = dot3(r1, r0)
        d_r2_r0 = dot3(r2, r0)
        @inbounds for k in 1:3
            r1_proj[k] = d_r1_r0 * r0[k] / nr0sq +
                         epsilon * r1Xr0[k] / nr1Xr0
            r2_proj[k] = d_r2_r0 * r0[k] / nr0sq +
                         epsilon * r2Xr0[k] / nr2Xr0
        end
        cross3!(r1_projXr2_proj, r1_proj, r2_proj)

        nr1pXr2p = norm3(r1_projXr2_proj)
        nr1_proj = norm3(r1_proj)
        nr2_proj = norm3(r2_proj)
        d_sum = 0.0
        @inbounds for k in 1:3
            d_sum += r0[k] * (r1_proj[k]/nr1_proj -
                              r2_proj[k]/nr2_proj)
        end
        coeff = (gamma / (4π)) / (nr1pXr2p^2) * d_sum
        @inbounds for k in 1:3
            vel_ind_proj[k] = coeff * r1_projXr2_proj[k]
        end

        scale = nr1Xr0 / (nr0 * epsilon)
        @inbounds for k in 1:3
            vel[k] = scale * vel_ind_proj[k]
        end
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
    r0 = work_vectors[1]
    r1 = work_vectors[2]
    r2 = work_vectors[3]
    r_perp = work_vectors[4]
    r1Xr2 = work_vectors[5]
    r1Xr0 = work_vectors[6]
    r2Xr0 = work_vectors[7]
    normr1r2 = work_vectors[8]

    r0 .= filament.x2 .- filament.x1
    r1 .= XVP .- filament.x1
    r2 .= XVP .- filament.x2

    # Vector perpendicular to core radius
    nr0 = norm3(r0)
    nr0sq = nr0 * nr0
    d_r1_r0 = dot3(r1, r0)
    @inbounds for k in 1:3
        r_perp[k] = d_r1_r0 * r0[k] / nr0sq
    end

    # Cut-off radius
    epsilon = sqrt(4 * ALPHA0 * NU * norm3(r_perp) / v_a)

    cross3!(r1Xr2, r1, r2)
    cross3!(r1Xr0, r1, r0)
    cross3!(r2Xr0, r2, r0)

    nr1 = norm3(r1)
    nr2 = norm3(r2)
    @inbounds for k in 1:3
        normr1r2[k] = r1[k]/nr1 - r2[k]/nr2
    end

    # Check point location relative to filament
    nr1Xr0 = norm3(r1Xr0)
    if nr1Xr0 / nr0 > epsilon
        nr1Xr2 = norm3(r1Xr2)
        coeff = (gamma / (4π)) / (nr1Xr2^2) * dot3(r0, normr1r2)
        @inbounds for k in 1:3
            vel[k] = coeff * r1Xr2[k]
        end
    elseif nr1Xr0 / nr0 < 1e-12 * epsilon
        vel .= 0.0
    else
        # Project onto core radius — reuse r_perp, normr1r2
        r1_proj = r_perp
        r2_proj = normr1r2
        nr2Xr0 = norm3(r2Xr0)
        d_r2_r0 = dot3(r2, r0)
        @inbounds for k in 1:3
            r1_proj[k] = d_r1_r0 * r0[k] / nr0sq +
                         epsilon * r1Xr0[k] / nr1Xr0
            r2_proj[k] = d_r2_r0 * r0[k] / nr0sq +
                         epsilon * r2Xr0[k] / nr2Xr0
        end

        cross3!(r1Xr2, r1_proj, r2_proj)
        nr1Xr2_val = norm3(r1Xr2)
        nr1_proj = norm3(r1_proj)
        nr2_proj = norm3(r2_proj)
        d_sum = 0.0
        @inbounds for k in 1:3
            d_sum += r0[k] * (r1_proj[k]/nr1_proj -
                              r2_proj[k]/nr2_proj)
        end
        coeff = (gamma / (4π)) / (nr1Xr2_val^2) * d_sum
        scale = nr1Xr0 / (nr0 * epsilon)
        @inbounds for k in 1:3
            vel[k] = scale * coeff * r1Xr2[k]
        end
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

function reinit!(filament::SemiInfiniteFilament, x1::PosVector, direction::PosVector, vel_mag::Real, filament_direction::Real)
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
    r1 = work_vectors[1]
    r_perp = work_vectors[2]
    r1XVf = work_vectors[3]
    GAMMA = -GAMMA * filament.filament_direction
    r1 .= XVP .- filament.x1

    # Calculate core radius
    d_r1_Vf = dot3(r1, Vf)
    @inbounds for k in 1:3
        r_perp[k] = d_r1_Vf * Vf[k]
    end
    epsilon = sqrt(4 * ALPHA0 * NU * norm3(r_perp) / v_a)

    cross3!(r1XVf, r1, Vf)

    nr1XVf = norm3(r1XVf)
    nVf = norm3(Vf)
    nr1 = norm3(r1)
    if nr1XVf / nVf > epsilon
        K = GAMMA / (4π) / (nr1XVf^2) * (1 + d_r1_Vf / nr1)
        @inbounds for k in 1:3
            vel[k] = K * r1XVf[k]
        end
    elseif nr1XVf / nVf < 1e-12 * epsilon
        vel .= 0.0
    else
        r1_proj = work_vectors[4]
        cross_tmp = work_vectors[5]
        # temp = r1/nr1 - Vf
        @inbounds for k in 1:3
            cross_tmp[k] = r1[k]/nr1 - Vf[k]
        end
        n_tmp = norm3(cross_tmp)
        @inbounds for k in 1:3
            r1_proj[k] = d_r1_Vf * Vf[k] +
                         epsilon * cross_tmp[k] / n_tmp
        end
        cross3!(cross_tmp, r1_proj, Vf)
        K = GAMMA / (4π) / (norm3(cross_tmp)^2) *
            (1 + dot3(r1_proj, Vf) / norm3(r1_proj))
        @inbounds for k in 1:3
            vel[k] = K * cross_tmp[k]
        end
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

@inline norm3(a) = sqrt(a[1]*a[1] + a[2]*a[2] + a[3]*a[3])
@inline dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline function normalize3!(v)
    n = norm3(v)
    n > 0 && (v[1] /= n; v[2] /= n; v[3] /= n)
    nothing
end

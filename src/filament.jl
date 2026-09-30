
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
@with_kw mutable struct BoundFilament{T} <: Filament
    x1::MVector{3, T}   = zeros(MVector{3, T})
    x2::MVector{3, T}   = zeros(MVector{3, T})
    length::T           = zero(T)
    r0::MVector{3, T}   = zeros(MVector{3, T})
    initialized::Bool   = false
end

function reinit!(filament::BoundFilament{T}, x1, x2, vec=zeros(MVector{3, T})) where {T}
    filament.x1 .= x1
    filament.x2 .= x2
    vec .= x2 .- x1
    filament.length = norm(vec)
    filament.r0 .= x2 .- x1
    filament.initialized = true
    return nothing
end

"""
    velocity_3D_bound_vortex!(vel, filament::BoundFilament, XVP,
                              gamma, core_radius_fraction, work_vectors)

Calculate induced velocity by a bound vortex filament at a point in space, with a
core radius of `core_radius_fraction` times the filament length.
"""
function velocity_3D_bound_vortex!(
    vel,
    filament::BoundFilament,
    XVP,
    gamma,
    core_radius_fraction,
    work_vectors
)
    epsilon = core_radius_fraction * filament.length
    velocity_3D_vortex_segment!(vel, filament, XVP, gamma, epsilon, work_vectors)
end

"""
    velocity_3D_trailing_vortex!(vel, filament::BoundFilament,
                                 XVP, gamma, va, work_vectors)

Calculate induced velocity by a trailing vortex filament, with a Lamb–Oseen core
radius grown over the axial distance of `XVP` from the filament start.

# Arguments
- `XVP`: Control point coordinates
- `gamma`: Vortex strength
- `va`: Inflow velocity magnitude
- work_vectors: preallocated array of intermediate variables

Reference: Rick Damiani et al. "A vortex step method for nonlinear airfoil polar data
as implemented in KiteAeroDyn".
"""
@inline function velocity_3D_trailing_vortex!(
    vel,
    filament::BoundFilament,
    XVP,
    gamma,
    va,
    work_vectors
)
    r1 = work_vectors[1]
    r1 .= XVP .- filament.x1
    axial_distance = abs(dot3(r1, filament.r0)) / filament.length
    epsilon = lamb_oseen_core_radius(axial_distance, va)
    velocity_3D_vortex_segment!(vel, filament, XVP, gamma, epsilon, work_vectors)
end

"""
    velocity_3D_vortex_segment!(vel, filament::BoundFilament, XVP,
                                gamma, epsilon, work_vectors)

Calculate the Biot–Savart velocity induced by a straight vortex segment at `XVP`.
Inside the core radius `epsilon` the velocity is evaluated on the core boundary and
scaled linearly with the distance to the axis. Without a core, it is zero within
1e-12 of the axis, relative to the distance from `x1`.
"""
@inline function velocity_3D_vortex_segment!(
    vel,
    filament::BoundFilament,
    XVP,
    gamma,
    epsilon,
    work_vectors
)
    r1, r2, r1Xr2, r1Xr0, r1r2norm = work_vectors
    r0 = filament.r0
    nr0 = filament.length
    r1 .= XVP .- filament.x1
    r2 .= XVP .- filament.x2

    cross3!(r1Xr0, r1, r0)
    axis_distance = norm3(r1Xr0) / nr0
    if on_axis_without_core(axis_distance, epsilon, norm3(r1))
        vel .= 0.0
    elseif axis_distance > epsilon
        cross3!(r1Xr2, r1, r2)
        nr1 = norm3(r1)
        nr2 = norm3(r2)
        @inbounds for k in 1:3
            r1r2norm[k] = r1[k]/nr1 - r2[k]/nr2
        end
        nr1Xr2 = norm3(r1Xr2)
        coeff = (gamma / (4π)) / (nr1Xr2^2) * dot3(r0, r1r2norm)
        @inbounds for k in 1:3
            vel[k] = coeff * r1Xr2[k]
        end
    else
        nr0sq = nr0 * nr0
        end_terms = core_end_term(dot3(r1, r0), nr0sq, epsilon) -
                    core_end_term(dot3(r2, r0), nr0sq, epsilon)
        coeff = -core_coefficient(gamma, end_terms, nr0sq, epsilon)
        @inbounds for k in 1:3
            vel[k] = coeff * r1Xr0[k]
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
- `va`::Float64=zero(Float64): apparent wind speed [m/s]
- `filament_direction`::Int64=0   : Direction indicator (-1 or 1)
- initialized::Bool=false
"""
@with_kw mutable struct SemiInfiniteFilament{T} <: Filament
    x1::MVector{3, T} = zeros(MVector{3, T})
    direction::MVector{3, T} = zeros(MVector{3, T})
    va::T = zero(T)
    filament_direction::Int64 = zero(Int64)
    initialized::Bool = false
end

function reinit!(filament::SemiInfiniteFilament{T}, x1::AbstractVector,
                 direction::AbstractVector, va::Real, filament_direction::Real) where T
    filament.x1 .= x1
    filament.direction .= direction
    filament.va = va
    filament.filament_direction = filament_direction
    filament.initialized = true
    return nothing
end

"""
    velocity_3D_trailing_vortex_semiinfinite!(vel, filament::SemiInfiniteFilament,
                                              Vf, XVP, GAMMA, va, work_vectors)

Calculate the velocity induced at `XVP` by a semi-infinite trailing vortex filament along
`Vf`, with a Lamb–Oseen core radius grown over the axial distance of `XVP` from `x1`.
Inside the core the velocity is scaled linearly with the distance to the axis. Without a
core, it is zero within 1e-12 of the axis, relative to the distance from `x1`.

# Arguments
- `Vf`: unit direction of the filament [-]
- `XVP`: evaluation point [m]
- `GAMMA`: vortex strength [m²/s]
- `va`: apparent wind speed [m/s]
- `work_vectors`: preallocated 3-vectors for intermediate results
"""
function velocity_3D_trailing_vortex_semiinfinite!(
    vel,
    filament::SemiInfiniteFilament,
    Vf,
    XVP,
    GAMMA,
    va,
    work_vectors
)
    r1 = work_vectors[1]
    r1XVf = work_vectors[3]
    GAMMA = -GAMMA * filament.filament_direction
    r1 .= XVP .- filament.x1

    d_r1_Vf = dot3(r1, Vf)
    nVf = norm3(Vf)
    epsilon = lamb_oseen_core_radius(abs(d_r1_Vf) / nVf, va)

    cross3!(r1XVf, r1, Vf)
    nr1XVf = norm3(r1XVf)
    axis_distance = nr1XVf / nVf
    nr1 = norm3(r1)
    if on_axis_without_core(axis_distance, epsilon, nr1)
        vel .= 0.0
        return nothing
    elseif axis_distance > epsilon
        K = GAMMA / (4π) / (nr1XVf^2) * (1 + d_r1_Vf / nr1)
    else
        nVfsq = nVf * nVf
        K = core_coefficient(GAMMA, 1 + core_end_term(d_r1_Vf, nVfsq, epsilon), nVfsq,
                             epsilon)
    end
    @inbounds for k in 1:3
        vel[k] = K * r1XVf[k]
    end
    nothing
end

"""
    lamb_oseen_core_radius(axial_distance, va)

Lamb–Oseen core radius [m] of a trailing vortex `axial_distance` [m] downstream of its
start, in an apparent wind of `va` [m/s].
"""
lamb_oseen_core_radius(axial_distance, va) = sqrt(4 * ALPHA0 * NU * axial_distance / va)

"""
    on_axis_without_core(axis_distance, epsilon, nr1)

True when both the distance to the filament axis and the core radius `epsilon` are
within 1e-12 of `nr1`, the distance from the filament start.
"""
@inline on_axis_without_core(axis_distance, epsilon, nr1) =
    max(axis_distance, epsilon) <= 1e-12 * nr1

"""
    core_end_term(d, nr0sq, epsilon)

Contribution of one filament end to the core velocity, for `d` the dot product of the
end-to-point vector with the axis vector of squared length `nr0sq`.
"""
@inline core_end_term(d, nr0sq, epsilon) = d / sqrt(d^2 / nr0sq + epsilon^2)

"""
    core_coefficient(gamma, end_terms, nr0sq, epsilon)

Factor on the cross product of the point vector and the axis vector that gives the
velocity inside a core of radius `epsilon`, linear in the distance to the axis.
"""
@inline core_coefficient(gamma, end_terms, nr0sq, epsilon) =
    (gamma / (4π)) * end_terms / (epsilon^2 * nr0sq)

"""
    cross3!(result::AbstractVector{T}, a::AbstractVector{T}, b::AbstractVector{T}) where T

Compute cross product of 3D vectors in-place.
"""
@inline function cross3!(result::AbstractVector{T}, a::AbstractVector{T},
                         b::AbstractVector{T}) where T
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

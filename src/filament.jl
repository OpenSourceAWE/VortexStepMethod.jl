
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
"""
struct BoundFilament <: Filament
    x1::Vector{Float64}    # First point
    x2::Vector{Float64}    # Second point
    length::Float64        # Filament length
    r0::Vector{Float64}    # Vector from x1 to x2

    function BoundFilament(x1::Vector{Float64}, x2::Vector{Float64})
        new(x1, x2, norm(x2 - x1), x2 - x1)
    end
end

"""
    velocity_3D_bound_vortex(filament::BoundFilament, XVP::Vector{Float64}, 
                           gamma::Float64, core_radius_fraction::Float64)

Calculate induced velocity by a bound vortex filament at a point in space.
"""
function velocity_3D_bound_vortex(
    filament::BoundFilament,
    XVP::Vector{Float64},
    gamma::Float64,
    core_radius_fraction::Float64
)
    r0 = filament.r0
    r1 = XVP - filament.x1
    r2 = XVP - filament.x2

    # Cross products
    r1Xr0 = cross(r1, r0)
    r2Xr0 = cross(r2, r0)

    # Cut-off radius
    epsilon = core_radius_fraction * norm(r0)
    
    # Check point location relative to filament
    if norm(r1Xr0) / norm(r0) > epsilon
        r1Xr2 = cross(r1, r2)
        return (gamma / (4π)) * r1Xr2 / (norm(r1Xr2)^2) * 
               dot(r0, r1/norm(r1) - r2/norm(r2))
    elseif norm(r1Xr0) / norm(r0) == 0
        return zeros(3)
    else
        @info "inside core radius"
        @info "distance from control point to filament: $(norm(r1Xr0) / norm(r0))"
        
        # Project onto core radius
        r1_proj = dot(r1, r0) * r0 / (norm(r0)^2) + 
                  epsilon * r1Xr0 / norm(r1Xr0)
        r2_proj = dot(r2, r0) * r0 / (norm(r0)^2) + 
                  epsilon * r2Xr0 / norm(r2Xr0)
        r1Xr2_proj = cross(r1_proj, r2_proj)
        
        vel_ind_proj = (gamma / (4π)) * r1Xr2_proj / (norm(r1Xr2_proj)^2) * 
                      dot(r0, r1_proj/norm(r1_proj) - r2_proj/norm(r2_proj))
        
        return norm(r1Xr0) / (norm(r0) * epsilon) * vel_ind_proj
    end
end

"""
    velocity_3D_trailing_vortex(filament::BoundFilament, 
                              XVP::Vector{Float64}, 
                              gamma::Float64, 
                              Uinf::Float64)

Calculate induced velocity by a trailing vortex filament.

# Arguments
- `XVP`: Control point coordinates
- `gamma`: Vortex strength
- `Uinf`: Inflow velocity magnitude

Reference: Rick Damiani et al. "A vortex step method for nonlinear airfoil polar data 
as implemented in KiteAeroDyn".
"""
function velocity_3D_trailing_vortex(
    filament::BoundFilament,
    XVP::Vector{Float64},
    gamma::Float64,
    Uinf::Float64
)
    r0 = filament.x2 - filament.x1  # Vortex filament
    r1 = XVP - filament.x1          # Control point to first end
    r2 = XVP - filament.x2          # Control point to second end

    # Vector perpendicular to core radius
    r_perp = dot(r1, r0) * r0 / (norm(r0)^2)
    
    # Cut-off radius
    epsilon = sqrt(4 * ALPHA0 * NU * norm(r_perp) / Uinf)

    # Cross products
    r1Xr0 = cross(r1, r0)
    r2Xr0 = cross(r2, r0)

    # Check point location relative to filament
    if norm(r1Xr0) / norm(r0) > epsilon
        r1Xr2 = cross(r1, r2)
        return (gamma / (4π)) * r1Xr2 / (norm(r1Xr2)^2) * 
               dot(r0, r1/norm(r1) - r2/norm(r2))
    elseif norm(r1Xr0) / norm(r0) == 0
        return zeros(3)
    else
        # Project onto core radius
        r1_proj = dot(r1, r0) * r0 / (norm(r0)^2) + 
                  epsilon * r1Xr0 / norm(r1Xr0)
        r2_proj = dot(r2, r0) * r0 / (norm(r0)^2) + 
                  epsilon * r2Xr0 / norm(r2Xr0)
        r1Xr2_proj = cross(r1_proj, r2_proj)
        
        vel_ind_proj = (gamma / (4π)) * r1Xr2_proj / (norm(r1Xr2_proj)^2) * 
                      dot(r0, r1_proj/norm(r1_proj) - r2_proj/norm(r2_proj))
        
        return norm(r1Xr0) / (norm(r0) * epsilon) * vel_ind_proj
    end
end

"""
    SemiInfiniteFilament

Represents a semi-infinite vortex filament.
"""
struct SemiInfiniteFilament <: Filament
    x1::Vector{Float64}         # Starting point
    direction::Vector{Float64}  # Direction vector
    vel_mag::Float64           # Velocity magnitude
    filament_direction::Int    # Direction indicator (-1 or 1)
end

"""
    velocity_3D_trailing_vortex_semiinfinite(filament::SemiInfiniteFilament, 
                                           Vf::Vector{Float64}, XVP::Vector{Float64},
                                           GAMMA::Float64, Uinf::Float64)

Calculate induced velocity by a semi-infinite trailing vortex filament.
"""
function velocity_3D_trailing_vortex_semiinfinite(
    filament::SemiInfiniteFilament,
    Vf::Vector{Float64},
    XVP::Vector{Float64},
    GAMMA::Float64,
    Uinf::Float64
)
    GAMMA = -GAMMA * filament.filament_direction
    r1 = XVP - filament.x1
    r1XVf = cross(r1, Vf)
    
    # Calculate core radius
    r_perp = dot(r1, Vf) * Vf
    epsilon = sqrt(4 * ALPHA0 * NU * norm(r_perp) / Uinf)
    
    if norm(r1XVf) / norm(Vf) > epsilon
        K = GAMMA / (4π) / norm(r1XVf)^2 * (1 + dot(r1, Vf) / norm(r1))
        return K * r1XVf
    elseif norm(r1XVf) / norm(Vf) == 0
        return zeros(3)
    else
        r1_proj = dot(r1, Vf) * Vf + 
                  epsilon * (r1/norm(r1) - Vf) / norm(r1/norm(r1) - Vf)
        r1XVf_proj = cross(r1_proj, Vf)
        K = GAMMA / (4π) / norm(r1XVf_proj)^2 * 
            (1 + dot(r1_proj, Vf) / norm(r1_proj))
        return K * r1XVf_proj
    end
end
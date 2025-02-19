using LinearAlgebra

"""
    Panel

Represents a panel in a vortex step method simulation.

# Fields
- `TE_point_1::Vector{Float64}`: First trailing edge point
- `LE_point_1::Vector{Float64}`: First leading edge point
- `TE_point_2::Vector{Float64}`: Second trailing edge point
- `LE_point_2::Vector{Float64}`: Second leading edge point
- `chord::Float64`: Panel chord length
- `va::Union{Nothing,Vector{Float64}}`: Panel velocity
- `corner_points::Matrix{Float64}`: Panel corner points
- `panel_aero_model::String`: Aerodynamic model type
- `aerodynamic_center::Vector{Float64}`: Panel aerodynamic center
- `control_point::Vector{Float64}`: Panel control point
- `bound_point_1::Vector{Float64}`: First bound point
- `bound_point_2::Vector{Float64}`: Second bound point
- `x_airf::Vector{Float64}`: Unit vector perpendicular to chord line
- `y_airf::Vector{Float64}`: Unit vector parallel to chord line
- `z_airf::Vector{Float64}`: Unit vector in spanwise direction
- `width::Float64`: Panel width
- `filaments::Vector{BoundFilament}`: Panel filaments
"""
mutable struct Panel
    TE_point_1::MVec3
    LE_point_1::MVec3
    TE_point_2::MVec3
    LE_point_2::MVec3
    chord::Float64
    va::Union{Nothing,Vector{Float64}}
    corner_points::Matrix{Float64}
    panel_aero_model::String
    cl_coefficients::Union{Nothing,Vector{Float64}}
    cd_coefficients::Union{Nothing,Vector{Float64}}
    cm_coefficients::Union{Nothing,Vector{Float64}}
    panel_polar_data::Union{Nothing,Matrix{Float64}}
    aerodynamic_center::MVec3
    control_point::MVec3
    bound_point_1::MVec3
    bound_point_2::MVec3
    x_airf::MVec3
    y_airf::MVec3
    z_airf::MVec3
    width::Float64
    filaments::Vector{Union{BoundFilament, SemiInfiniteFilament}}

    function Panel(
        section_1::Section,
        section_2::Section,
        aerodynamic_center::PosVector,
        control_point::PosVector,
        bound_point_1::PosVector,
        bound_point_2::PosVector,
        x_airf::PosVector,
        y_airf::PosVector,
        z_airf::PosVector
    )
        # Initialize basic geometry
        TE_point_1 = copy(section_1.TE_point)
        LE_point_1 = copy(section_1.LE_point)
        TE_point_2 = copy(section_2.TE_point)
        LE_point_2 = copy(section_2.LE_point)
        
        chord = mean([
            norm(TE_point_1 - LE_point_1),
            norm(TE_point_2 - LE_point_2)
        ])
        
        corner_points = hcat(LE_point_1, TE_point_1, TE_point_2, LE_point_2)
        
        # Validate aero model consistency
        if section_1.aero_input[1] != section_2.aero_input[1]
            throw(ArgumentError("Both sections must have the same aero_input"))
        end
        panel_aero_model = isa(section_1.aero_input, String) ? section_1.aero_input : section_1.aero_input[1]
        
        # Initialize aerodynamic properties
        cl_coeffs = nothing
        cd_coeffs = nothing
        cm_coeffs = nothing
        polar_data = nothing
        
        if panel_aero_model == "lei_airfoil_breukels"
            cl_coeffs, cd_coeffs, cm_coeffs = compute_lei_coefficients(section_1, section_2)
        elseif panel_aero_model == "polar_data"
            aero_1 = section_1.aero_input[2]
            aero_2 = section_2.aero_input[2]
            if size(aero_1) != size(aero_2)
                throw(ArgumentError("Polar data must have same shape"))
            end
            polar_data = (aero_1 + aero_2) / 2
        elseif panel_aero_model != "inviscid"
            throw(ArgumentError("Unsupported aero model: $panel_aero_model"))
        end
        
        # Calculate width
        width = norm(bound_point_2 - bound_point_1)
        
        # Initialize filaments
        filaments = [
            BoundFilament(bound_point_2, bound_point_1),
            BoundFilament(bound_point_1, TE_point_1),
            BoundFilament(TE_point_2, bound_point_2)
        ]

        new(
            TE_point_1, LE_point_1, TE_point_2, LE_point_2,
            chord, nothing, corner_points, panel_aero_model,
            cl_coeffs, cd_coeffs, cm_coeffs, polar_data,
            aerodynamic_center, control_point,
            bound_point_1, bound_point_2,
            x_airf, y_airf, z_airf,
            width, filaments
        )
    end
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
    panel::Panel, 
    induced_velocity::Vector{Float64}
)
    # Calculate relative velocity and angle of attack
    # Constants throughout iterations: panel.va, panel.x_airf, panel.y_airf
    relative_velocity = panel.va .+ induced_velocity
    v_normal = dot(panel.x_airf, relative_velocity)
    v_tangential = dot(panel.y_airf, relative_velocity)
    alpha = atan(v_normal / v_tangential)
    
    return alpha, relative_velocity
end

"""
    compute_lei_coefficients(section_1::Section, section_2::Section)

Compute lift, drag and moment coefficients for Lei airfoil using Breukels model.
"""
function compute_lei_coefficients(section_1::Section, section_2::Section)
    # Average tube diameter and camber from both sections
    t1, k1 = section_1.aero_input[2]
    t2, k2 = section_2.aero_input[2]
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
    S = Dict{Int,Float64}()
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
    calculate_relative_alpha_and_velocity(panel::Panel, induced_velocity::Vector{Float64})

Calculate relative angle of attack and relative velocity of the panel.
"""
function calculate_relative_alpha_and_velocity(panel::Panel, induced_velocity::Vector{Float64})
    relative_velocity = panel.va + induced_velocity
    v_normal = dot(panel.x_airf, relative_velocity)
    v_tangential = dot(panel.y_airf, relative_velocity)
    alpha = atan(v_normal / v_tangential)
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
function calculate_cl(panel::Panel, alpha::Float64)
    if panel.panel_aero_model == "lei_airfoil_breukels"
        cl = evalpoly(rad2deg(alpha), reverse(panel.cl_coefficients))
        if abs(alpha) > (π/9)
            cl = 2 * cos(alpha) * sin(alpha)^2
        end
        return cl
    elseif panel.panel_aero_model == "inviscid"
        return 2π * alpha
    elseif panel.panel_aero_model == "polar_data"
        return linear_interpolation(
            panel.panel_polar_data[:,1],
            panel.panel_polar_data[:,2];
            extrapolation_bc=Line()
        )(alpha)
    else
        throw(ArgumentError("Unsupported aero model: $(panel.panel_aero_model)"))
    end
end

"""
    calculate_cd_cm(panel::Panel, alpha::Float64)

Calculate drag and moment coefficients for given angle of attack.
"""
function calculate_cd_cm(panel::Panel, alpha::Float64)
    if panel.panel_aero_model == "lei_airfoil_breukels"
        cd = evalpoly(rad2deg(alpha), reverse(panel.cd_coefficients))
        cm = evalpoly(rad2deg(alpha), reverse(panel.cm_coefficients))
        if abs(alpha) > (π/9)  # Outside ±20 degrees
            cd = 2 * sin(alpha)^3
        end
        return cd, cm
    elseif panel.panel_aero_model == "inviscid"
        return 0.0, 0.0
    elseif panel.panel_aero_model == "polar_data"
        cd = linear_interpolation(
            panel.panel_polar_data[:,1],
            panel.panel_polar_data[:,3];
            extrapolation_bc=Line()
        )(alpha)
        cm = linear_interpolation(
            panel.panel_polar_data[:,1],
            panel.panel_polar_data[:,4];
            extrapolation_bc=Line()
        )(alpha)
        return cd, cm
    else
        throw(ArgumentError("Unsupported aero model: $(panel.panel_aero_model)"))
    end
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
    calculate_velocity_induced_single_ring_semiinfinite(
        panel::Panel,
        evaluation_point::Vector{Float64},
        evaluation_point_on_bound::Bool,
        va_norm::Float64,
        va_unit::Vector{Float64},
        gamma::Float64,
        core_radius_fraction::Float64
    )

Calculate the velocity induced by a vortex ring at a control point.

# Arguments
- `evaluation_point`: Point where induced velocity is evaluated
- `evaluation_point_on_bound`: Whether evaluation point is on bound vortex
- `va_norm`: Norm of apparent velocity
- `va_unit`: Unit vector of apparent velocity
- `gamma`: Circulation strength
- `core_radius_fraction`: Vortex core radius as fraction of panel width

# Returns
- `Vector{Float64}`: Induced velocity vector
"""
function calculate_velocity_induced_single_ring_semiinfinite(
    panel::Panel,
    evaluation_point::PosVector,
    evaluation_point_on_bound::Bool,
    va_norm::Float64,
    va_unit::VelVector,
    gamma::Float64,
    core_radius_fraction::Float64,
    work_vectors::NTuple{10,Vector{Float64}}
)
    velind = zeros(Float64, 3)
    tempvel = zeros(3)

    # Process each filament
    for (i, filament) in enumerate(panel.filaments)
        if i == 1  # bound filament
            if evaluation_point_on_bound
                tempvel .= zeros(Float64, 3)
            else
                velocity_3D_bound_vortex!(
                    tempvel,
                    filament,
                    evaluation_point,
                    gamma,
                    core_radius_fraction,
                    work_vectors
                )
            end
        elseif i ∈ (2, 3)  # trailing filaments
            velocity_3D_trailing_vortex!(
                tempvel,
                filament,
                evaluation_point,
                gamma,
                va_norm,
                work_vectors
            )
        elseif i ∈ (4, 5)  # semi-infinite trailing filaments
            velocity_3D_trailing_vortex_semiinfinite!(
                tempvel,
                filament,
                va_unit,
                evaluation_point,
                gamma,
                va_norm,
                work_vectors
            )
        else
            tempvel .= zeros(Float64, 3)
        end
        velind .+= tempvel
    end

    return velind
end

"""
    calculate_velocity_induced_bound_2D(panel::Panel, evaluation_point::Vector{Float64})

Calculate velocity induced by bound vortex filaments at the control point.
Only needed for VSM, as LLT bound and filament align, thus no induced velocity.

# Arguments
- `panel::Panel`: Panel object
- `evaluation_point::Vector{Float64}`: Point where induced velocity is evaluated

# Returns
- `Vector{Float64}`: Induced velocity at the control point
"""
function calculate_velocity_induced_bound_2D(
    panel::Panel, 
    evaluation_point::PosVector
)
    # r3 perpendicular to the bound vortex
    r3 = evaluation_point - (panel.bound_point_1 + panel.bound_point_2) / 2
    
    # r0 is the direction of the bound vortex
    r0 = panel.bound_point_1 - panel.bound_point_2
    
    # Calculate cross product
    cross_ = cross(r0, r3)
    
    # Calculate induced velocity
    return (cross_ / sum(cross_.^2) / 2π) * norm(r0)
end
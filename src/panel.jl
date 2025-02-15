
"""
    Panel

Structure representing a wing panel with geometric and aerodynamic properties.
"""
mutable struct Panel
    # Geometric properties
    TE_point_1::Vector{Float64}
    LE_point_1::Vector{Float64}
    TE_point_2::Vector{Float64}
    LE_point_2::Vector{Float64}
    chord::Float64
    corner_points::Matrix{Float64}
    width::Float64
    
    # Aerodynamic properties
    panel_aero_model::String
    aerodynamic_center::Vector{Float64}
    control_point::Vector{Float64}
    bound_point_1::Vector{Float64}
    bound_point_2::Vector{Float64}
    x_airf::Vector{Float64}
    y_airf::Vector{Float64}
    z_airf::Vector{Float64}
    va::Union{Vector{Float64}, Nothing}
    
    # Coefficients
    cl_coefficients::Union{Vector{Float64}, Nothing}
    cd_coefficients::Union{Vector{Float64}, Nothing}
    cm_coefficients::Union{Vector{Float64}, Nothing}
    panel_polar_data::Union{Matrix{Float64}, Nothing}
    
    # Filaments
    filaments::Vector{Filament}

    function Panel(
        section_1::Section,
        section_2::Section,
        aerodynamic_center::Vector{Float64},
        control_point::Vector{Float64},
        bound_point_1::Vector{Float64},
        bound_point_2::Vector{Float64},
        x_airf::Vector{Float64},
        y_airf::Vector{Float64},
        z_airf::Vector{Float64}
    )
        # Validate aero model consistency
        section_1.aero_input[1] == section_2.aero_input[1] || 
            throw(ArgumentError("Both sections must have same aero_input"))
        
        # Initialize geometric properties
        TE_point_1 = copy(section_1.TE_point)
        LE_point_1 = copy(section_1.LE_point)
        TE_point_2 = copy(section_2.TE_point)
        LE_point_2 = copy(section_2.LE_point)
        
        chord = mean([
            norm(TE_point_1 - LE_point_1),
            norm(TE_point_2 - LE_point_2)
        ])
        
        corner_points = hcat(LE_point_1, TE_point_1, TE_point_2, LE_point_2)
        width = norm(bound_point_2 - bound_point_1)
        
        # Initialize aerodynamic properties
        panel_aero_model = section_1.aero_input[1]
        
        # Initialize coefficients based on aero model
        cl_coeff, cd_coeff, cm_coeff, polar_data = 
            initialize_aero_coefficients(panel_aero_model, section_1, section_2)
        
        # Initialize filaments
        filaments = [
            BoundFilament(bound_point_2, bound_point_1),
            BoundFilament(bound_point_1, TE_point_1),
            BoundFilament(TE_point_2, bound_point_2)
        ]
        
        new(
            TE_point_1, LE_point_1, TE_point_2, LE_point_2,
            chord, corner_points, width,
            panel_aero_model,
            aerodynamic_center, control_point,
            bound_point_1, bound_point_2,
            x_airf, y_airf, z_airf, nothing,
            cl_coeff, cd_coeff, cm_coeff, polar_data,
            filaments
        )
    end
end

"""
    calculate_relative_alpha_and_relative_velocity(panel::Panel, 
                                                induced_velocity::Vector{Float64})

Calculate relative angle of attack and velocity for panel.
"""
function calculate_relative_alpha_and_relative_velocity(
    panel::Panel,
    induced_velocity::Vector{Float64}
)
    isnothing(panel.va) && throw(ArgumentError("Panel velocity not set"))
    
    relative_velocity = panel.va + induced_velocity
    v_normal = dot(panel.x_airf, relative_velocity)
    v_tangential = dot(panel.y_airf, relative_velocity)
    alpha = atan(v_normal, v_tangential)
    
    return alpha, relative_velocity
end

"""
    calculate_cl(panel::Panel, alpha::Float64)

Calculate lift coefficient for given angle of attack.
"""
function calculate_cl(panel::Panel, alpha::Float64)
    if panel.panel_aero_model == "lei_airfoil_breukels"
        cl = Polynomial(panel.cl_coefficients)(rad2deg(alpha))
        if abs(alpha) > π/9  # Outside ±20 degrees
            cl = 2 * cos(alpha) * sin(alpha)^2
        end
        return cl
    elseif panel.panel_aero_model == "inviscid"
        return 2π * alpha
    elseif panel.panel_aero_model == "polar_data"
        return linear_interpolation(
            panel.panel_polar_data[:, 1],
            panel.panel_polar_data[:, 2],
            alpha
        )
    else
        throw(ArgumentError("Unsupported aerodynamic model"))
    end
end

"""
    calculate_cd_cm(panel::Panel, alpha::Float64)

Calculate drag and moment coefficients for given angle of attack.
"""
function calculate_cd_cm(panel::Panel, alpha::Float64)
    if panel.panel_aero_model == "lei_airfoil_breukels"
        cd = Polynomial(panel.cd_coefficients)(rad2deg(alpha))
        cm = Polynomial(panel.cm_coefficients)(rad2deg(alpha))
        
        # Outside ±20 degrees
        if abs(alpha) > π/9
            cd = 2 * sin(alpha)^3
        end
        
        return cd, cm
        
    elseif panel.panel_aero_model == "inviscid"
        return 0.0, 0.0
        
    elseif panel.panel_aero_model == "polar_data"
        cd = linear_interpolation(
            panel.panel_polar_data[:, 1],
            panel.panel_polar_data[:, 2],
            alpha
        )
        cm = linear_interpolation(
            panel.panel_polar_data[:, 1],
            panel.panel_polar_data[:, 3],
            alpha
        )
        return cd, cm
    else
        throw(ArgumentError("Unsupported aerodynamic model"))
    end
end

"""
    calculate_velocity_induced_bound_2D(panel::Panel, 
                                     evaluation_point::Vector{Float64})

Calculate induced velocity by bound vortex filaments at control point.
Only needed for VSM, as LLT bound and filament align.
"""
function calculate_velocity_induced_bound_2D(
    panel::Panel,
    evaluation_point::Vector{Float64}
)
    # Calculate r3 perpendicular to bound vortex
    r3 = evaluation_point - (panel.bound_point_1 + panel.bound_point_2) / 2
    
    # Calculate r0 as direction of bound vortex
    r0 = panel.bound_point_1 - panel.bound_point_2
    
    cross_product = cross(r0, r3)
    return cross_product / sum(cross_product.^2) / 2π * norm(r0)
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

Calculate velocity induced by a ring at specified control point.
"""
function calculate_velocity_induced_single_ring_semiinfinite(
    panel::Panel,
    evaluation_point::Vector{Float64},
    evaluation_point_on_bound::Bool,
    va_norm::Float64,
    va_unit::Vector{Float64},
    gamma::Float64,
    core_radius_fraction::Float64
)
    velind = zeros(3)
    
    for (i, filament) in enumerate(panel.filaments)
        tempvel = if i == 1  # bound
            evaluation_point_on_bound ? 
                zeros(3) : 
                velocity_3D_bound_vortex(filament, evaluation_point, gamma, core_radius_fraction)
        elseif i ∈ (2, 3)  # trailing1 or trailing2
            velocity_3D_trailing_vortex(filament, evaluation_point, gamma, va_norm)
        elseif i ∈ (4, 5)  # trailing_semi_inf1 or trailing_semi_inf2
            velocity_3D_trailing_vortex_semiinfinite(filament, va_unit, evaluation_point, gamma, va_norm)
        else
            zeros(3)
        end
        
        velind .+= tempvel
    end
    
    return velind
end

"""
    calculate_filaments_for_plotting(panel::Panel)

Calculate filaments for visualization with proper directions and colors.
"""
function calculate_filaments_for_plotting(panel::Panel)
    isnothing(panel.va) && throw(ArgumentError("Panel velocity not set"))
    
    filaments = Vector{Vector{Any}}()
    
    for (i, filament) in enumerate(panel.filaments)
        x1 = filament.x1
        
        if isdefined(filament, :x2) && !isnothing(filament.x2)
            x2 = filament.x2
            color = i == 0 ? "magenta" : "green"
        else
            # Handle semi-infinite filaments
            x2 = x1 + 2 * panel.chord * (panel.va / norm(panel.va))
            color = "orange"
            
            if filament.filament_direction == -1
                x1, x2 = x2, x1
                color = "red"
            end
        end
        
        push!(filaments, [x1, x2, color])
    end
    
    return filaments
end

"""
    initialize_aero_coefficients(
        model::String,
        section_1::Section,
        section_2::Section
    )

Initialize aerodynamic coefficients based on model type.
"""
function initialize_aero_coefficients(
    model::String,
    section_1::Section,
    section_2::Section
)
    if model == "lei_airfoil_breukels"
        return initialize_lei_airfoil_coefficients(section_1, section_2)
    elseif model == "inviscid"
        return nothing, nothing, nothing, nothing
    elseif model == "polar_data"
        return initialize_polar_data(section_1, section_2)
    else
        throw(ArgumentError("Unsupported aerodynamic model: $model"))
    end
end
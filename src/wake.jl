
"""
    frozen_wake(body_aero::BodyAerodynamics, va_distribution)

Update the filaments of the panels with frozen wake model.
Uses one shared wake vector computed from area-weighted distributed inflow.

Replaces older filaments if present by checking length of filaments.

# Arguments
- `body_aero`::BodyAerodynamics: see: [`BodyAerodynamics`](@ref) 
- `va_distribution::Matrix{Float64}`: Array of velocity vectors at each panel

# Returns
- nothing
"""
function frozen_wake!(body_aero::BodyAerodynamics, va_distribution)
    n_panels = length(body_aero.panels)
    size(va_distribution) == (n_panels, 3) ||
        throw(ArgumentError("va_distribution must be shape ($(n_panels), 3), got $(size(va_distribution))"))

    panel_areas = [panel.chord * panel.width for panel in body_aero.panels]
    wake_velocity = _compute_reference_velocity_from_distribution(
        va_distribution,
        n_panels,
        panel_areas
    )
    wake_speed = norm(wake_velocity)
    if wake_speed <= 0.0
        # Zero apparent flow: keep existing wake state and avoid NaN/Inf updates.
        return nothing
    end
    wake_direction = wake_velocity / wake_speed

    for panel in body_aero.panels
        reinit!(
            panel.filaments[4],
            panel.TE_point_1, 
            wake_direction, 
            wake_speed, 
            1
        )
        reinit!(
            panel.filaments[5], 
            panel.TE_point_2, 
            wake_direction, 
            wake_speed, 
            -1
        )
    end
    return nothing
end

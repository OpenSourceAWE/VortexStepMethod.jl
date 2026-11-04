
"""
    frozen_wake(body_aero::BodyAerodynamics, va_vec_dist)

Update the filaments of the panels with frozen wake model.
Each panel's wake follows its own row of `va_vec_dist`.

Replaces older filaments if present by checking length of filaments.

# Arguments
- `body_aero`::BodyAerodynamics: see: [`BodyAerodynamics`](@ref) 
- `va_vec_dist::Matrix{Float64}`: Array of velocity vectors at each panel

# Returns
- nothing
"""
function frozen_wake!(body_aero::BodyAerodynamics, va_vec_dist)
    n_panels = length(body_aero.panels)
    size(va_vec_dist) == (n_panels, 3) ||
        throw(ArgumentError(
            "va_vec_dist must be shape ($(n_panels), 3), got $(size(va_vec_dist))"))

    for (i, panel) in enumerate(body_aero.panels)
        wake_velocity = matrix_row(va_vec_dist, i)
        wake_speed = norm(wake_velocity)
        # Zero apparent flow: keep existing wake state and avoid NaN/Inf updates.
        wake_speed <= 0.0 && continue
        wake_direction = wake_velocity / wake_speed
        reinit!(panel.filaments[4], panel.TE_point_1, wake_direction, wake_speed, 1)
        reinit!(panel.filaments[5], panel.TE_point_2, wake_direction, wake_speed, -1)
    end
    return nothing
end


"""
    frozen_wake(body_aero::BodyAerodynamics, va_distribution)

Update the filaments of the panels with frozen wake model.

Replaces older filaments if present by checking length of filaments.

# Arguments
- `body_aero`::BodyAerodynamics: see: [BodyAerodynamics](@ref) 
- `va_distribution::Matrix{Float64}`: Array of velocity vectors at each panel

# Returns
- nothing
"""
function frozen_wake!(body_aero::BodyAerodynamics, va_distribution)
    for (i, panel) in enumerate(body_aero.panels)
        va_i = va_distribution[i, :]
        vel_mag = norm(va_i)
        direction = va_i / vel_mag
        init!(
            panel.filaments[4],
            panel.TE_point_1, 
            direction, 
            vel_mag, 
            1
        )
        init!(
            panel.filaments[5], 
            panel.TE_point_2, 
            direction, 
            vel_mag, 
            -1
        )
    end
    return nothing
end
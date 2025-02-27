
"""
    frozen_wake(va_distribution::Matrix{Float64}, panels::Vector{Panel})

Update the filaments of the panels with frozen wake model.

Replaces older filaments if present by checking length of filaments.

# Arguments
- `va_distribution::Matrix{Float64}`: Array of velocity vectors at each panel
- `panels::Vector{Panel}`: List of panels

# Returns
- `Vector{Panel}`: List of panels with updated filaments
"""
function frozen_wake(va_distribution::Matrix{Float64}, panels::Vector{Panel})
    for (i, panel) in enumerate(panels)
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
    
    return panels
end
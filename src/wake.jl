
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
        
        # Check and update filaments
        if length(panel.filaments) == 3
            # Add new wake filaments
            push!(panel.filaments, 
                  SemiInfiniteFilament(
                      panel.TE_point_1, 
                      direction, 
                      vel_mag, 
                      1))
            
            push!(panel.filaments, 
                  SemiInfiniteFilament(
                      panel.TE_point_2, 
                      direction, 
                      vel_mag, 
                      -1))
                      
        elseif length(panel.filaments) == 5
            # Replace existing wake filaments
            panel.filaments[4] = SemiInfiniteFilament(
                panel.TE_point_1, 
                direction, 
                vel_mag, 
                1)
                
            panel.filaments[5] = SemiInfiniteFilament(
                panel.TE_point_2, 
                direction, 
                vel_mag, 
                -1)
                
        else
            throw(ArgumentError("Panel has unexpected number of filaments: $(length(panel.filaments))"))
        end
    end
    
    return panels
end
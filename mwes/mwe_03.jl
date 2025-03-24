# Fix allocations of this array comprehension:
# y = [panel.control_point[2] for panel in body_aero.panels]
# This mwe shows that the line above is allocating a lot,
# and provides the function test_new() that shows a better way of implementing.

using VortexStepMethod, PreallocationTools

# Step 1: Define wing parameters
n_panels = 20          # Number of panels
span = 20.0            # Wing span [m]
chord = 1.0            # Chord length [m]
v_a = 20.0             # Magnitude of inflow velocity [m/s]
density = 1.225        # Air density [kg/m³]
alpha_deg = 30.0       # Angle of attack [degrees]
alpha = deg2rad(alpha_deg)

# Step 2: Create wing geometry with linear panel distribution
wing = Wing(n_panels, spanwise_distribution=LINEAR)

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],    # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    INVISCID)
add_section!(wing, 
    [0.0, -span/2, 0.0],   # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    INVISCID)

# Step 3: Initialize aerodynamics
body_aero::BodyAerodynamics = BodyAerodynamics([wing])

y = [panel.control_point[2] for panel in body_aero.panels]
n = @allocated y = [panel.control_point[2] for panel in body_aero.panels]

function test_old(body_aero, gamma_i)
    y = [panel.control_point[2] for panel in body_aero.panels]
    y
end

function test_new(body_aero, gamma_i)
    y = body_aero.y
    for (i, panel) in pairs(body_aero.panels) 
        y[i] = panel.control_point[2] 
    end
    y
end


gamma_i=zeros(length(body_aero.panels))

# make sure the two test functions get compiled
test_old(body_aero, gamma_i)
test_new(body_aero, gamma_i)

println("Old function, using an array comprehension:")
n = @allocated test_old(body_aero, gamma_i)
println("Allocations: ", n)

println("New function, using a pre-allocated array and a for loop:")
m = @allocated test_new(body_aero, gamma_i)
println("Allocations: ", m)

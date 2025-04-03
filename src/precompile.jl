@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    # list = [OtherType("hello"), OtherType("world!")]
    path = dirname(pathof(@__MODULE__))

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        vss = vs("vsm_settings_dual.yaml")
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

        gamma_initial = zeros(length(body_aero.panels))
        calculate_circulation_distribution_elliptical_wing(gamma_initial, body_aero)

        nothing
    end
end

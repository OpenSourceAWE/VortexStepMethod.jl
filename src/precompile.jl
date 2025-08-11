@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    
    @compile_workload begin
        # Minimal precompilation to avoid segfault
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        
        # Step 1: Define wing parameters
        n_panels = 2           # Reduced number of panels for faster precompilation
        span = 10.0            # Wing span [m]
        chord = 1.0            # Chord length [m]

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

        # Step 3: Initialize aerodynamics (simplified)
        body_aero = BodyAerodynamics([wing])

        nothing
    end
end

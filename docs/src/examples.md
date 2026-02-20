```@meta
CurrentModule = VortexStepMethod
```

## Rectangular Wing
##### Wing Aerodynamics Analysis using Vortex-Step Method (VSM) and Lift Line Theory (LLT)
This example demonstrates the process of setting up and analyzing a wing's
aerodynamics using VSM and LLT. We'll cover the following steps:
1. Importing necessary libraries.
2. Define the wing parameters
3. Create wing geometry with linear panel distribution and add the wing sections
4. Initializing the wing aerodynamics and set the inflow conditions.
5. Initialize solvers for both LLT and VSM methods
6. Running a simulation with both methods
7. Plotting combined analysis

First, install Julia and set up the examples environment as explained in the
section [Running the examples as developer](@ref). Then launch Julia:
```bash
julia --project=examples
```

#### Step 1: Importing the necessary libraries:
```julia
julia> using LinearAlgebra
julia> using GLMakie
julia> using VortexStepMethod
```

#### Step 2: Define wing parameters
```julia
julia> n_panels = 20          # Number of panels
julia> span = 20.0            # Wing span [m]
julia> chord = 1.0            # Chord length [m]
julia> v_a = 20.0             # Magnitude of inflow velocity [m/s]
julia> alpha_deg = 30.0       # Angle of attack [degrees]
julia> alpha = deg2rad(alpha_deg)
```

#### Step 3: Create wing geometry with linear panel distribution
```julia
julia> wing = Wing(n_panels, spanwise_distribution=LINEAR)
```

##### Add wing sections - defining only tip sections with inviscid airfoil model
```julia
julia> add_section!(wing,
           [0.0, span/2, 0.0],   # Left tip LE
           [chord, span/2, 0.0],  # Left tip TE
           INVISCID)
julia> add_section!(wing,
           [0.0, -span/2, 0.0],  # Right tip LE
           [chord, -span/2, 0.0], # Right tip TE
           INVISCID)
```

##### Refine the mesh
```julia
julia> refine!(wing)
```

#### Step 4: Initialize aerodynamics
```julia
julia> body_aero = BodyAerodynamics([wing])
```
We need to pass here an array of wing objects, because a body can have
multiple wings.

###### Set inflow conditions
```julia
julia> vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
julia> set_va!(body_aero, vel_app, [0, 0, 0.1])
```

#### Step 5: Initialize solvers for both LLT and VSM methods
```julia
julia> llt_solver = Solver(body_aero; aerodynamic_model_type=LLT)
julia> vsm_solver = Solver(body_aero; aerodynamic_model_type=VSM)
```

#### Step 6: Solve using both methods
```julia
julia> results_llt = solve(llt_solver, body_aero)
julia> results_vsm = solve(vsm_solver, body_aero)
```

##### Print results comparison
```julia
julia> println("\nLifting Line Theory Results:")
julia> println("CL = $(round(results_llt["cl"], digits=4))")
julia> println("CD = $(round(results_llt["cd"], digits=4))")
julia> println("\nVortex Step Method Results:")
julia> println("CL = $(round(results_vsm["cl"], digits=4))")
julia> println("CD = $(round(results_vsm["cd"], digits=4))")
julia> println("Projected area = $(round(results_vsm["projected_area"], digits=4)) m²")
```

#### Step 7: Plot combined analysis
```julia
julia> angle_range = range(0, 20, 20)
julia> plot_combined_analysis(
           [llt_solver, vsm_solver],
           [body_aero, body_aero],
           [results_llt, results_vsm];
           solver_label=["LLT", "VSM"],
           angle_range=angle_range,
           angle_type="angle_of_attack",
           v_a=v_a,
           title="Rectangular Wing",
           is_show=true,
       )
```

## More examples
From the examples environment (`julia --project=examples`), you can execute
more examples by typing:
```julia
julia> include("examples/menu.jl")
```
You should see the following menu:
```
Choose function to execute or `q` to quit:
 > V3_kite = include("V3_kite.jl")
   pyramid_model = include("pyramid_model.jl")
   rectangular_wing = include("rectangular_wing.jl")
   ram_air_kite = include("ram_air_kite.jl")
   stall_model = include("stall_model.jl")
   bench = include("bench.jl")
   cleanup = include("cleanup.jl")
   quit
```
You can select one of the examples using the `<UP>` and `<DOWN>` keys.
Press `<ENTER>` to run the selected example.
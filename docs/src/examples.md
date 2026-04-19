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
After launching Julia with (`jl`), you can execute
more examples by typing:

```julia
julia> include("examples/menu.jl")
```

or, you prefer to use the ControlPlots library:

```julia
julia> include("examples_cp/menu_cp.jl")
```

You should see the following menu:

```text
Choose function to execute or `q` to quit:
 > V3_kite = include("V3_kite.jl")
   billowing = include("billowing.jl")
   pyramid_model = include("pyramid_model.jl")
   rectangular_wing = include("rectangular_wing.jl")
   ram_air_kite = include("ram_air_kite.jl")
   stall_model = include("stall_model.jl")
   bench = include("bench.jl")
   cleanup = include("cleanup.jl")
   GLMakie.activate!()
   CairoMakie.activate!()
   help_me = VortexStepMethod.help("https://opensourceawe.github.io/VortexStepMethod.jl/dev")
   quit
```

You can select one of the examples using the `<UP>` and `<DOWN>` keys.
Press `<ENTER>` to run the selected example.

## Plotting Backends

The examples in this package support three plotting backends. Here is a comparison to help you choose:

### GLMakie
**Advantages:**
- Interactive plots: zoom, pan, rotate 3D scenes in a native window.
- Hardware-accelerated rendering via OpenGL — fast for large datasets.
- Supports animations and live-updating plots.

**Disadvantages:**
- Requires a display server (does not work in headless/server environments without a virtual framebuffer).
- Heavier dependency: needs OpenGL drivers and a GPU.
- Longer initial load time compared to the other backends.

### CairoMakie
**Advantages:**
- Fully software-rendered — works in headless environments (CI, servers, SSH sessions).
- Produces high-quality vector output (SVG, PDF) suitable for publication.
- Lighter dependency than GLMakie (no GPU required).

**Disadvantages:**
- Plots are static — no interactive zoom or pan.
- Slower for very large or complex scenes because rendering is done in software.
- 3D support is limited compared to GLMakie.

### ControlPlots (based on PyPlot / Matplotlib)
**Advantages:**
- Simple API, easy to learn for students
- In addition, the Matplotlib API for users coming from Python/Matplotlib is supported.
- Works in headless environments; can save to PNG, SVG, PDF, etc.
- Very lightweight Julia-side dependency (delegates work to Python).

**Disadvantages:**
- Requires a working Python installation with Matplotlib (via `PyCall`).
- Can cause issues when multithreading is enabled.
- No native Makie ecosystem integration (e.g. cannot use `Makie.Observable` for live updates).
- Interactivity is limited and depends on the Matplotlib backend in use.
- Extra setup complexity when Python or Matplotlib are not already installed.

| Feature | GLMakie | CairoMakie | ControlPlots |
|---|---|---|---|
| Interactive (zoom/pan) | yes | no | yes |
| Headless / server | no* | yes | yes |
| Vector output (PDF/SVG) | no | yes | yes |
| GPU required | yes | no | no |
| 3D support | full | limited | limited |
| Load time | slow | medium | fast |

\* GLMakie can run headless with a virtual framebuffer (e.g. `Xvfb`), but this requires additional setup.
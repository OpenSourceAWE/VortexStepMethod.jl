# Tips and tricks

## What can this model simulate
The following bodies can be simulated:

- conventional bodies, consisting of one or more wings
- leading edge inflatable (LEI) kites
- RAM-air kites

To build the geometry of a RAM-air kite, a 3D .obj file can be used as input. In addition a `.dat` file is needed.
It should have two columns, one for the `x` and one for the `y` coordinate of the 2D polar that is used.

## Unrefined Section Distribution
When creating a wing, panel forces and moments are automatically computed for each unrefined section. The unrefined sections correspond to the original geometry sections you define, while panels represent the refined mesh used for aerodynamic calculations.

### How It Works
The solver automatically tracks which panels belong to which original unrefined section. After refinement (e.g., splitting each section into multiple panels), the aerodynamic forces and moments are aggregated back to the unrefined section level.

```julia
# Create wing with 4 sections, refined to 40 panels
wing = Wing(40)
add_section!(wing, [0, 5, 0], [1, 5, 0], INVISCID)    # Section 1
add_section!(wing, [0, 2.5, 0], [1, 2.5, 0], INVISCID) # Section 2
add_section!(wing, [0, 0, 0], [1, 0, 0], INVISCID)     # Section 3
add_section!(wing, [0, -5, 0], [1, -5, 0], INVISCID)   # Section 4
```

The 40 panels are distributed across the 3 unrefined panels (sections 1-2, 2-3, 3-4). Forces and moments are computed per panel during the solve, then automatically aggregated to the 3 unrefined panels for output.

This approach is useful for:
- LEI kites where you want loads per rib
- Wings with discrete control surfaces
- Cases where physical structure doesn't align with uniform panel distribution
- Dynamic simulations where you have fewer structural segments than panels needed for accurate VSM aerodynamics. For example, a 6-segment structural model can be combined with 40-panel aerodynamics, with loads automatically mapped back to the 6 structural segments.

## RAM-air kite model
If running the example `ram_air_kite.jl` fails, try to run the `cleanup.jl` script and then try again. Background: this example caches the calculated polars. Reading cached polars can fail after an update.

## Output formats
Currently, the `solve!()` function returns the results as [VSMSolution](@ref) struct. The function solve() returns a dictionary with the results. The `solve!()` function is faster, and the `solve()` contains many more entries, therefore the first function is good for integration in dynamic models and the second one better suited for aerodynamic analysis.

## Performance
Calling `reinit!(body_aero; init_aero=false)` is very fast. After calling `unrefined_deform!(wing, theta_angles, delta_angles)`, you have to run `reinit!(body_aero; init_aero=false)` to apply the deformed wing to the body aerodynamics. This is in turn necessary for the linearization from deformation to aerodynamic coefficients for RAM-air kites.

## Contributing
Please, read [CONTRIBUTING.md](https://github.com/OpenSourceAWE/VortexStepMethod.jl/blob/main/CONTRIBUTING.md)

## Building the documentation locally
You can build and serve the documentation locally after checking out the
source code with git:
```bash
julia --project=docs scripts/build_docu.jl
```
A browser window should pop up automatically with live-reloading enabled.
# Tips and tricks

## What can this model simulate
The following bodies can be simulated:

- conventional bodies, consisting of one or more wings
- leading edge inflatable (LEI) kites
- RAM-air kites

To build the geometry of a RAM-air kite, a 3D .obj file can be used as input. In addition a `.dat` file is needed.
It should have two columns, one for the `x` and one for the `y` coordinate of the 2D polar that is used.

## Panel Grouping Methods
When creating a wing, you can specify how panels should be grouped for moment and force calculations using the `grouping_method` parameter. Two methods are available:

### EQUAL_SIZE (Default)
Divides refined panels into equally-sized sequential groups. This is the original behavior.

```julia
wing = Wing(40; n_groups=4, grouping_method=EQUAL_SIZE)
```

In this example, with 40 panels and 4 groups, each group will contain 10 consecutive panels (panels 1-10, 11-20, 21-30, 31-40).

### REFINE
Groups refined panels back to their original unrefined section. This is useful when you want group moments and forces to represent the original wing structure, regardless of panel refinement.

```julia
# Create wing with 4 unrefined sections (3 panels)
wing = Wing(40; n_groups=3, grouping_method=REFINE)
add_section!(wing, [0, 5, 0], [1, 5, 0], INVISCID)    # Section 1
add_section!(wing, [0, 2.5, 0], [1, 2.5, 0], INVISCID) # Section 2
add_section!(wing, [0, 0, 0], [1, 0, 0], INVISCID)     # Section 3
add_section!(wing, [0, -5, 0], [1, -5, 0], INVISCID)   # Section 4
```

**Important:** When using `REFINE`, `n_groups` must equal the number of unrefined panels (number of sections - 1). The solver will automatically map each refined panel to its closest original unrefined panel and sum their moments and forces accordingly.

This is particularly useful for:
- LEI kites where you want loads per rib
- Wings with discrete control surfaces
- Cases where physical structure doesn't align with uniform panel distribution
- Dynamic simulations where you have fewer structural segments than panels needed for accurate VSM aerodynamics. For example, a 6-segment structural model can be combined with 40-panel aerodynamics by using `n_groups=6` and `grouping_method=REFINE` to map aerodynamic loads back to the structural segments.

## RAM-air kite model
If running the example `ram_air_kite.jl` fails, try to run the `cleanup.jl` script and then try again. Background: this example caches the calculated polars. Reading cached polars can fail after an update.

## Output formats
Currently, the `solve!()` function returns the results as [VSMSolution](@ref) struct. The function solve() returns a dictionary with the results. The `solve!()` function is faster, and the `solve()` contains many more entries, therefore the first function is good for integration in dynamic models and the second one better suited for aerodynamic analysis.

## Performance
Calling `reinit!(body_aero; init_aero=false)` is very fast. After calling `deform!(wing)`, you have to run `reinit!(body_aero; init_aero=false)` to apply the deformed wing to the body aerodynamics. This is in turn necessary for the linearization from deformation to aerodynamic coefficients for RAM-air kites.

## Contributing
Please, read [CONTRIBUTING.md](https://github.com/OpenSourceAWE/VortexStepMethod.jl/blob/main/CONTRIBUTING.md)

## Building the documentation locally
You can build the documentation locally after checking out the source code with git, launching Julia and executing:
```julia
include("scripts/build_docu.jl")
```
A browser window should pop up automatically.
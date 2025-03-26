## Tips and tricks

### What can this model simulate
The following bodies can be simulated:

- conventional bodies, consisting of one or more wings
- leading edge inflatable (LEI) kites
- RAM-air kites

To build the geometry of a RAM-air kite, a 3D .obj file can be used as input. In addition a `.dat` file is needed.
It should have two columns, one for the `x` and one for the `y` coordinate of the 2D polar that is used.

### Output formats
Currently, the `solve!()` function returns the results as [VSMSolution](@ref) struct. The function solve() returns a
dictionary with the results. The `solve!()` function is faster, and the `solve()` contains many more entries, therefore
the first function is good for integration in dynamic models and the second one better suited for aerodynamic analysis.

### Performance
Calling `init!(body_aero; init_aero=false)` is very fast. After calling `deform!(wing)`, you have to run `init!(body_aero; init_aero=false)` to apply the deformed wing to the body aerodynamics. This is in turn necessary for the linearization from deformation to aerodynamic coefficients for RAM-air kites.

### Building the documentation locally
You can build the documentation locally after checking out the source code with git, launching Julia and executing:
```
include("scripts/build_docu.jl")
```
A browser window should pop up automatically.
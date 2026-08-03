```@meta
CurrentModule = VortexStepMethod
```

# Aerodynamic models of 3D wings using the Vortex Step Method

The Vortex Step Method (VSM) is an enhanced lifting line method that improves upon the classic approach by solving the circulation system at the three-quarter chord position, among the most important details. This adjustment allows for more accurate calculations of lift and drag forces, particularly addressing the shortcomings in induced drag prediction. 
VSM is further refined by coupling it with 2D viscous airfoil polars, making it well-suited for complex geometries, 
including low aspect ratio wings, as well as configurations with sweep, dihedral, and anhedral angles.

The software presented here includes a couple of examples: a rectangular wing, a leading-edge inflatable kite and a ram-air kite.

This package was translated from the Python code version 1.0.0 available at [https://github.com/ocayon/Vortex-Step-Method](https://github.com/ocayon/Vortex-Step-Method) with some extensions as documented in [News.md](https://github.com/OpenSourceAWE/VortexStepMethod.jl/blob/main/NEWS.md).

## Installation
Install [Julia 1.11](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html)
or later, if you haven't already.

Before installing this software it is suggested to create a new project, for
example like this:

```bash
mkdir vsm
cd vsm
julia --project=.
```

Then add VortexStepMethod from Julia's package manager, by typing `]` to enter
the package manager (`pkg>` prompt), then:

```julia
(vsm) pkg> add VortexStepMethod
```

Press backspace to return to the `julia>` prompt. You can run the unit tests
from the `pkg>` prompt with:

```julia
(vsm) pkg> test VortexStepMethod
```

You can also type `?` to enter help mode (`help?>` prompt) to look up
documentation for any function.

To run the examples, type at the `julia>` prompt:

```julia
julia> using VortexStepMethod
julia> VortexStepMethod.install_examples()
julia> include("examples/menu.jl")
```

This copies the example scripts to an `examples/` folder and installs the
required packages ([GLMakie](https://docs.makie.org/stable/), CSV, etc.).

## Running the examples as developer
If you have git installed, check out this repo because it makes it easier to understand the code:

```bash
mkdir repos
cd repos
git clone https://github.com/OpenSourceAWE/VortexStepMethod.jl
cd VortexStepMethod.jl/bin
./install
```

You can launch Julia with:

```bash
jl
```

or with:

```bash
./bin/run_julia
```

Then you can display a menu with the available examples using the GLMakie library:

```julia
menu()
```

To browse the code, it is suggested to use
[VSCode](https://code.visualstudio.com/) with the Julia plugin.

## Input
Three kinds of input data is needed:

- The wing geometry, defined by section:
  - for the rectangular wing two sections, two points in CAD reference frame + polars  
    (three different options to provide them) per section
  - kite wing: model of polars included, n sections to define

- The airflow and turn rate:
  - `v_app` vector and `omega` (turn rate) vector in Kite Body (KB) reference frame

- The configuration:
  - how many panels  
    --> two sections make a panel.

Wing geometry can also be loaded from YAML files or `.obj` files. See the examples for details.

### Example for defining the required input:

```julia

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
    [0.0, span/2, 0.0],   # Left tip LE
    [chord, span/2, 0.0],  # Left tip TE
    INVISCID)
add_section!(wing,
    [0.0, -span/2, 0.0],  # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    INVISCID)

# Refine the mesh
refine!(wing)

# Step 3: Initialize aerodynamics
body_aero = BodyAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
set_va!(body_aero, vel_app)
```

It is possible to import the wing geometry using an `.obj` file as shown in the example `ram_air_kite.jl`. During the import the polars are calculated automatically, using NeuralFoil by default or XFoil as a viscous cross-check. This approach is valid for rigid wings and ram-air kites, but not for leading edge inflatable kites. See [From CAD mesh to aerodynamic model](@ref) for the full pipeline.

Surfplan geometries are handled by the built-in `SurfplanAdapter` submodule: the upstream (Python) [SurfplanAdapter](https://github.com/jellepoland/SurfplanAdapter) exports a Surfplan file to an adapter directory, and `VortexStepMethod.SurfplanAdapter.surfplan_to_aero_yaml` converts that into a solver-ready `geometry.yaml`.

## Output
- the aerodynamic forces Fx, Fy, Fz
- the aerodynamic moments Mx, My, Mz
- the force coefficients CL, CD, CS (side force coefficient)
- the status of the solver (is the result valid)

In addition, the spanwise distribution of these and additional values are available.

See also the [documentation](https://OpenSourceAWE.github.io/VortexStepMethod.jl/dev/).

## Citation
If you use this project in your research, please consider citing it. 
Citation details can be found in the [CITATION.cff](https://github.com/OpenSourceAWE/VortexStepMethod.jl/blob/main/CITATION.cff) file included in this repository.

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/OpenSourceAWE/VortexStepMethod.jl/blob/main/LICENSE) file for details.

## WAIVER
Technische Universiteit Delft hereby disclaims all copyright interest in the package “VortexStepMethod.jl” written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering

### Copyright
Copyright (c) 2022 Oriol Cayon

Copyright (c) 2024 Oriol Cayon, Jelle Poland, TU Delft

Copyright (c) 2025, 2026 Oriol Cayon, Jelle Poland, Bart van de Lint, Uwe Fechner

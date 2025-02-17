[![Build Status](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Simulation of a Kite using the Vortex Step Method

The Vortex Step Method (VSM) is an enhanced lifting line method that improves upon the classic approach by solving the 
circulation system at the three-quarter chord position, among the most important details. This adjustment allows for 
more accurate calculations of lift and drag forces, particularly addressing the shortcomings in induced drag prediction. 
VSM is further refined by coupling it with 2D viscous airfoil polars, making it well-suited for complex geometries, 
including low aspect ratio wings, as well as configurations with sweep, dihedral, and anhedral angles.

The software presented here includes a couple of examples: a rectangular wing and a leading-edge inflatable kite.

This package was translated from the Python code available at https://github.com/awegroup/Vortex-Step-Method .

## Installation
Install [Julia 1.10](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) or later, 
if you haven't already. On Linux, make sure that Python3 and Matplotlib are installed:
```
sudo apt install python3-matplotlib
```

Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project=.
```
Then add VortexStepMethod from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl"
``` 
at the Julia prompt. You can run the unit tests with the command:
```julia
pkg"test VortexStepMethod"
```

## Input
Three kinds of input data is needed:

- The wing geometry, defined by section:
  - rec wing two section, two point + polars
  - kite: model of polars included, n sections to define

- The airflow:
  - v_app vector

- The configuration:
  - how many panels  
    --> two sections make a panel.

There is no input file yet, the input has to be defined in the code.

### Example for defining the required input:
```julia

# Step 1: Define wing parameters
n_panels = 20          # Number of panels
span = 20.0            # Wing span [m]
chord = 1.0            # Chord length [m]
Umag = 20.0            # Magnitude of inflow velocity [m/s]
density = 1.225        # Air density [kg/m³]
alpha_deg = 30.0       # Angle of attack [degrees]
alpha = deg2rad(alpha_deg)

# Step 2: Create wing geometry with linear panel distribution
wing = Wing(n_panels, spanwise_panel_distribution="linear")

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],    # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    "inviscid")
add_section!(wing, 
    [0.0, -span/2, 0.0],   # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    "inviscid")

# Step 3: Initialize aerodynamics
wa = WingAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* Umag
set_va!(wa, (vel_app, 0.0))  # Second parameter is yaw rate
```

Surfplan output file can be converted to an input for the vortex step method with a tool that is in this repo.

## Output
- CL, CD, CS (side force coefficient)
- the spanwise distribution of forces  
  --> moment coefficients (not yet implemented)

## Citation
If you use this project in your research, please consider citing it. 
Citation details can be found in the [CITATION.cff](CITATION.cff) file included in this repository.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Copyright
Copyright (c) 2022 Oriol Cayon

Copyright (c) 2024 Oriol Cayon, Jelle Poland, TU Delft

Copyright (c) 2025 Oriol Cayon, Jelle Poland, Bart van de Lint
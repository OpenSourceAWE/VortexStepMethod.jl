[![Build Status](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Simulation of a Kite using the Vortex Step Method

The Vortex Step Method (VSM) is an enhanced lifting line method that improves upon the classic approach by solving the circulation system at the three-quarter chord position, among the most important details. This adjustment allows for more accurate calculations of lift and drag forces, particularly addressing the shortcomings in induced drag prediction. VSM is further refined by coupling it with 2D viscous airfoil polars, making it well-suited for complex geometries, including low aspect ratio wings, as well as configurations with sweep, dihedral, and anhedral angles.

The software presented here includes a couple of examples: a rectangular wing and a leading-edge inflatable kite.

This package was translated from the Python code available at https://github.com/awegroup/Vortex-Step-Method .

## Installation
Install [Julia 1.10](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) or later, if you haven't already. On Linux, make sure that Python3 and Matplotlib are installed:
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
pkg"add VortexStepMethod"
``` 
at the Julia prompt. You can run the unit tests with the command:
```julia
pkg"test VortexStepMethod"
```

## Input
- geometry, defined by section
  - rec wing two section, two point + polars
  - kite: model of polars included, n sections to define

- flow
  - v_app vector

- config
  - how many panels
    -> two sections make a panel

Surfplan output file can be converted to an input for the vortex step method with a tool that is in this repo.

## Output
- cl, cd, cs (side force coefficient)
- spanwise distribution of forces
  -> moment coefficients (not yet implemented)

## Citation
If you use this project in your research, please consider citing it. 
Citation details can be found in the [CITATION.cff](CITATION.cff) file included in this repository.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Copyright
Copyright (c) 2022 Oriol Cayon

Copyright (c) 2024 Oriol Cayon, Jelle Poland, TU Delft

Copyright (c) 2025 Oriol Cayon, Jelle Poland, Bart van de Lint
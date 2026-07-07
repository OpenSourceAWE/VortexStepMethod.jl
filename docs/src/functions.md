```@meta
CurrentModule = VortexStepMethod
```
## Functions for creating the geometry
```@docs
add_section!
refine!
calculate_span
calculate_projected_area
load_polar_data
```

## Surface-pressure (Cp) tables
The Cp table types and IO live in the core package; the generation of Cp tables
(from XFoil/NeuralFoil) lives in `AirfoilAero`.
```@docs
CpData
CpPolar
read_cp_data
write_cp_data
cp_distribution
delta_cp
```

## Airfoil aerodynamics (AirfoilAero)
These live in the `AirfoilAero` sub-package (`lib/AirfoilAero`), which converts
airfoil coordinates to polars and Cp tables. Load it with `using AirfoilAero`.
```@meta
CurrentModule = AirfoilAero
```
```@docs
fit_kulfan_parameters
KulfanFitMethod
LeastSquaresFit
EnvelopeFit
AbstractAirfoilSolver
XFoilSolver
NeuralFoilSolver
SectionSolution
DeformedSection
deform_section
analyze_section
analyze_sweep
generate_cp_polar
```

## OBJ mesh conversion (ObjAdapter)
These live in the `ObjAdapter` sub-package (`lib/ObjAdapter`), which converts a
3D wing `.obj` mesh to the native YAML/CSV geometry format. Load it with
`using ObjAdapter`.
```@meta
CurrentModule = ObjAdapter
```
```@docs
obj_to_yaml
perpendicular_sections
auto_rotation
plot_airfoils
plot_slices_3d
```
```@meta
CurrentModule = VortexStepMethod
```

## Setting the inflow conditions and solving
```@docs
set_va!
solve
solve!
reinit!(body_aero::BodyAerodynamics)
linearize
calculate_results
```

## Main Plotting Functions
The plotting functions are implemented as [package extensions](https://pkgdocs.julialang.org/v1.11/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)). They are available when `GLMakie` (or `ControlPlots`) is loaded before `VortexStepMethod`. The examples use `GLMakie`.
```@docs
plot_geometry
plot_distribution
plot_polars
plot_polar_data
plot_combined_analysis
plot_section_polars
```

## Helper Functions
```@docs
save_plot
show_plot
```

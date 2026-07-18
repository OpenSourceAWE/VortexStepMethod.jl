```@meta
CurrentModule = VortexStepMethod
```
## Functions for creating the geometry
```@docs
add_section!
refine!
calculate_span
calculate_projected_area
```

## Surface aero (contour + Cp + cf) tables
The per-section surface aero table type and its npz IO live in the core package; the
generation (from XFoil/NeuralFoil) lives in `AirfoilAero`.
```@docs
SectionAero
section_surface
read_section_aero
write_section_aero
```

## Airfoil aerodynamics (AirfoilAero)
These live in the `AirfoilAero` submodule of `VortexStepMethod`, which converts
airfoil coordinates to polars and Cp tables. Load it with
`using VortexStepMethod.AirfoilAero`.
```@meta
CurrentModule = VortexStepMethod.AirfoilAero
```
```@docs
fit_kulfan_parameters
kulfan_to_coordinates
KulfanFitMethod
LeastSquaresFit
ShrinkWrap
shrink_wrap
AbstractAirfoilSolver
XFoilSolver
NeuralFoilSolver
SectionSolution
DeformedSection
deform_section
analyze_section
analyze_sweep
neuralfoil_aero
generate_aero_matrices
generate_polar_from_coordinates
generate_polar_from_dat
generate_section_aero
```

## OBJ mesh conversion (ObjAdapter)
These live in the `ObjAdapter` submodule of `VortexStepMethod`, which converts a
3D wing `.obj` mesh to the native YAML/CSV geometry format. Load it with
`using VortexStepMethod.ObjAdapter`.
```@meta
CurrentModule = VortexStepMethod.ObjAdapter
```
```@docs
obj_to_yaml
generate_section_polars
write_yaml
perpendicular_sections
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
reinit!(body_aero::BodyAerodynamics{P, W, T}) where {P, W, T}
linearize
calculate_results
```

## Main Plotting Functions
The plotting functions are implemented as [package extensions](https://pkgdocs.julialang.org/v1.11/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)). They are available once a Makie backend (`GLMakie` or `CairoMakie`) and `MakieControlPlots` are loaded before `VortexStepMethod`. The examples use `GLMakie`.
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

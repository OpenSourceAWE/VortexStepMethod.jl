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
obj_to_yaml
fit_kulfan_parameters
KulfanFitMethod
LeastSquaresFit
EnvelopeFit
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
plot_airfoils
plot_section_polars
```

## Helper Functions
```@docs
save_plot
show_plot
```

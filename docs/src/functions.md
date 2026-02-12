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
The plotting functions are implemented as [package extension](https://pkgdocs.julialang.org/v1.11/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)). This means that they are only available if the package `ControlPlots.jl` was loaded BEFORE loading `VortexStepMethod.jl`.
```@docs
plot_geometry
plot_distribution
plot_circulation_distribution
plot_polars
plot_polar_data
plot_combined_analysis
```

## Helper Functions
```@docs
save_plot
show_plot
```

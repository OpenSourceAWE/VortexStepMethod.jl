```@meta
CurrentModule = VortexStepMethod
```
## Functions for creating the geometry
```@docs
add_section!
```

## Setting the inflow conditions and solving
```@docs
set_va!
solve
```

## Main Plotting Functions
The plotting functions are implemented as [package extension](https://pkgdocs.julialang.org/v1.11/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)). This means that they are only available if the package `ControlPlots.jl` was loaded BEFORE loading `VortexStepMethod.jl`.
```@docs
plot_geometry
plot_distribution
plot_polars
```

## Helper Functions
```@docs
save_plot
show_plot
```

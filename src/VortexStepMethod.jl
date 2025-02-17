module VortexStepMethod

using LinearAlgebra
using StaticArrays
using Logging
using Statistics
using Colors
using DelimitedFiles
using Plots
using Measures
using LaTeXStrings

# Export public interface
export Wing, Section
export WingAerodynamics
export Solver, solve
export calculate_results, solve_circulation_distribution
export add_section!, set_va!
export calculate_span, calculate_projected_area
export plot_wing, plot_circulation_distribution, plot_geometry, plot_distribution, plot_polars

"""
   const MVec3    = MVector{3, Float64}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const MVec3    = MVector{3, Float64}

# Include core functionality
include("wing_geometry.jl")
include("filament.jl")
include("panel.jl")
include("wake.jl")
include("wing_aerodynamics.jl")
include("solver.jl")

# include plotting
include("color_palette.jl")
include("plotting.jl")

end # module
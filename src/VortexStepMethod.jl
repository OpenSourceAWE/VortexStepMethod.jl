module VortexStepMethod

using LinearAlgebra
using StaticArrays
using Logging
using Statistics
using Colors
using DelimitedFiles
using ControlPlots
using Measures
using LaTeXStrings
using NonlinearSolve
using Interpolations: linear_interpolation, Line

# Export public interface
export Wing, Section, KiteWing
export WingAerodynamics
export Solver, solve
export calculate_results, solve_circulation_distribution
export add_section!, set_va!
export calculate_span, calculate_projected_area
export plot_wing, plot_circulation_distribution, plot_geometry, plot_distribution, plot_polars
export show_plot, save_plot

"""
   const MVec3    = MVector{3, Float64}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const MVec3    = MVector{3, Float64}

"""
   const PosVector=Union{MVec3, Vector}

Position vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const PosVector=Union{MVec3, Vector, SizedVector{3, Float64, Vector{Float64}}}
const VelVector=Union{MVec3, Vector, SizedVector{3, Float64, Vector{Float64}}}

abstract type AbstractWing end

# Include core functionality
include("wing_geometry.jl")
include("kite_geometry.jl")
include("filament.jl")
include("panel.jl")
include("wake.jl")
include("wing_aerodynamics.jl")
include("solver.jl")

# include plotting
include("plotting.jl")

end # module
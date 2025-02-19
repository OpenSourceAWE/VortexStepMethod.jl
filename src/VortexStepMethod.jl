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

"""
   const MVec3    = MVector{3, Float64}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const MVec3    = MVector{3, Float64}

"""
   const PosVector=Union{MVec3, Vector}

Position vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const PosVector=Union{MVec3, Vector}
const VelVector=Union{MVec3, Vector}

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
include("color_palette.jl")
include("plotting.jl")

function create_wing_from_obj(obj_path, wing_mass, circle_center, radius, α_tip)
   if !isfile(obj_path)
      error("OBJ file not found: $obj_path")
   end
   data_dir = dirname(obj_path)
   wing_path = joinpath(data_dir, "measurements.bin")
   if isfile(wing_path)
      wing = deserialize(wing_path)
      new_wing = wing.mass != wing_mass ||
               wing.radius != radius ||
               wing.α_tip != α_tip
      @info "Using cached wing struct. Delete $wing_path to create a new wing."
      return deserialize()
   end
   if new_wing
      vertices = read_vertices(obj_path)
      faces = read_faces(obj_path)
      com = calculate_com(vertices, faces)
      inertia_tensor = calculate_inertia_tensor(vertices, faces, wing_mass, com)
      gamma_range = -abs(deg2rad(α_tip)):0.01:abs(deg2rad(α_tip))
      (interp_max, interp_min, gammas, max_xs, min_xs, inertia_tensor) = create_interpolations(vertices, circle_center, gamma_range)
   end
   (interp_max, interp_min, gammas, max_xs, min_xs, com) = deserialize(wing_path)
end

end # module
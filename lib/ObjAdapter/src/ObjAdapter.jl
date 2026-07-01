module ObjAdapter

using LinearAlgebra
using Interpolations
using Statistics
using NonlinearSolve
using Printf
using DelimitedFiles
import YAML
using StaticArrays
using AirfoilAero

include("obj_geometry.jl")
include("obj_slice.jl")
include("polar_generation.jl")
include("obj_to_yaml.jl")

export obj_to_yaml, write_geometry_yaml, airfoils_from_yaml
export slice_obj_wing, slice_obj_at_positions, perpendicular_sections
export generate_neuralfoil_polars
export read_faces, center_to_com!, calculate_inertia_tensor, calc_inertia_y_rotation

end

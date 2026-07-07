module ObjAdapter

using LinearAlgebra
using Interpolations
using Statistics
using NonlinearSolve
using Printf
using DelimitedFiles
import YAML
using StaticArrays
using ..AirfoilAero
import ..VortexStepMethod

include("obj_geometry.jl")
include("obj_slice.jl")
include("polar_generation.jl")
include("obj_to_yaml.jl")

"""
    plot_airfoil_fit(x, y; kwargs...)

Plot a single airfoil with its Kulfan CST fit. Implemented in `ObjAdapterMakieExt`;
load `Makie` (or `GLMakie`) to use it.
"""
function plot_airfoil_fit end

"""
    plot_airfoils(geometry_file::String; kwargs...)

Plot every airfoil of a YAML wing geometry. The airfoil shapes are read from the
`.dat` files referenced by the geometry's `wing_airfoils`. Routes to the active
plotting backend; implemented in `ObjAdapterMakieExt` (load `Makie`/`GLMakie`).
"""
function plot_airfoils end
function plot_airfoils(geometry_file::String; kwargs...)
    plot_airfoils(geometry_file, VortexStepMethod._active_backend(); kwargs...)
end

"""
    plot_slices_3d(obj_path; n_slices=10, rotation=I, is_show=true)

3D diagnostic: draw the OBJ surface mesh with the interpolated leading-/trailing-edge
curves and each station's sliced airfoil contour (+ chord line) overlaid, to inspect
how the slicer cuts a mesh. Implemented in `ObjAdapterMakieExt` (load `Makie`/`GLMakie`).
"""
function plot_slices_3d end

export obj_to_yaml, obj_to_matrix_yaml, resolve_aero_geometry
export write_geometry_yaml, write_yaml, airfoils_from_yaml
export perpendicular_sections, auto_rotation
export generate_neuralfoil_polars
export read_faces, center_to_com!, calculate_inertia_tensor, calc_inertia_y_rotation
export plot_airfoil_fit, plot_airfoils, plot_slices_3d

end

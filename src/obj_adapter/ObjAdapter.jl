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
`.dat` files referenced by the geometry's `wing_airfoils`. Implemented in the
`VortexStepMethodMakieExt` extension (load a Makie backend, e.g. `GLMakie`).
"""
function plot_airfoils end

"""
    plot_slices_3d(path; n_slices=10, rotation=I, delta=0.0, is_show=true)

3D slice diagnostic with a hover 2D airfoil panel. `path` is either a mesh `.obj`
(live preview: slice and wrap here) or a generated [`obj_to_yaml`](@ref) output
directory (audit: plot the written `.dat` airfoils the polar pipeline analysed).
Implemented in the Makie extension (load `Makie`/`GLMakie`); see it for all options.
"""
function plot_slices_3d end

export obj_to_yaml, resolve_aero_geometry
export write_geometry_yaml, write_yaml, airfoils_from_yaml
export perpendicular_sections
export generate_section_polars
export read_faces, center_to_com!, calculate_inertia_tensor, calc_inertia_y_rotation
export plot_airfoil_fit, plot_airfoils, plot_slices_3d

end

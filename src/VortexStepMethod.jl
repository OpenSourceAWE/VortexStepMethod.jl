module VortexStepMethod

using Arrow
using LinearAlgebra
using StaticArrays
using Logging
using Statistics
using Colors
using DelimitedFiles
using DefaultApplication
using Measures
using LaTeXStrings
using NonlinearSolve
using SciMLBase
import NonlinearSolve: solve, solve!
using Interpolations
using Parameters
using Serialization
using Timers
using PreallocationTools
using PrecompileTools
using Pkg
using DifferentiationInterface
using ForwardDiff
import YAML
using StructMapping

# Export public interface
export SolverSettings, VSMSettings, WingSettings
export ObjWing, Section, Wing, refine!, reinit!
export BodyAerodynamics
export Solver, VSMSolution, linearize, solve, solve!, solve_base!, calc_forces!
export calculate_results
export add_section!, set_va!, section_pitch_rate
export calculate_projected_area, calculate_span
export MVec3

export LLT, Model, VSM
export AeroModel, INVISCID, POLY, LEI_AIRFOIL_BREUKELS, POLAR_MATRICES, POLAR_VECTORS
export KulfanParameters
export BILLOWING, COSINE, LINEAR, PanelDistribution, SPLIT_PROVIDED, UNCHANGED
export ELLIPTIC, InitialGammaDistribution, ZEROS
export FAILURE, FEASIBLE, INFEASIBLE, SolverStatus
export LOOP, NONLIN, SolverType
export load_polar_data

# Per-section surface aero (contour + Cp + cf). Human-readable CSV/dat IO: reading here,
# writing in AirfoilAero.
export SectionAero, section_surface, read_section_aero

export plot_combined_analysis, plot_distribution, plot_geometry, plot_polar_data,
    plot_polars, plot_section_polars, save_plot, show_plot

# Plotting functions live in the `VortexStepMethodMakieExt` extension, loaded once a
# Makie backend (`GLMakie` or `CairoMakie`) and `MakieControlPlots` are available. The
# declarations below carry the public docstrings; the extension provides the methods.

"""
    plot_geometry(body_aero::BodyAerodynamics, title; kwargs...)

Plot wing geometry from different viewpoints and optionally save/show plots.

# Arguments
- `body_aero`: the [`BodyAerodynamics`](@ref) to plot
- `title`: plot title

# Keyword arguments
- `data_type`: file extension for saving (default: `".png"`)
- `save_path`: path for saving the graphic (default: `nothing`)
- `is_save`: whether to save the graphic (default: `false`)
- `is_show`: whether to display the graphic (default: `false`)
- `view_elevation`: initial view elevation angle in degrees (default: `15`)
- `view_azimuth`: initial view azimuth angle in degrees (default: `-120`)
- `use_tex`: use external `pdflatex` for rendering (default: `false`; ignored by Makie)
"""
function plot_geometry end

"""
    plot_distribution(y_coordinates_list, results_list, label_list; kwargs...)

Plot spanwise distributions of aerodynamic properties.

# Arguments
- `y_coordinates_list`: list of spanwise coordinate arrays
- `results_list`: list of result dictionaries from [`solve!`](@ref)
- `label_list`: list of labels for each result

# Keyword arguments
- `title`: plot title (default: `"spanwise_distribution"`)
- `data_type`: file extension for saving (default: `".png"`)
- `save_path`: path to save plots (default: `nothing`)
- `is_save`: whether to save (default: `false`)
- `is_show`: whether to display (default: `true`)
- `use_tex`: use external `pdflatex` for rendering (default: `false`; ignored by Makie)
"""
function plot_distribution end

"""
    plot_polars(solver_list, body_aero_list, label_list; kwargs...)

Plot polar data comparing different solvers and configurations.

# Arguments
- `solver_list`: list of aerodynamic solvers
- `body_aero_list`: list of [`BodyAerodynamics`](@ref) objects
- `label_list`: list of labels for each configuration

# Keyword arguments
- `literature_path_list`: optional paths to literature data CSV files (default: `String[]`)
- `angle_range`: range of angles to analyze in degrees (default: `range(0, 20, 2)`)
- `angle_type`: `"angle_of_attack"` or `"side_slip"` (default: `"angle_of_attack"`)
- `angle_of_attack`: AoA for the polar sweep (default: `0.0`) [°]
- `side_slip`: side slip angle (default: `0.0`) [°]
- `v_a`: apparent wind speed magnitude (default: `10.0`) [m/s]
- `title`: plot title (default: `"polar"`)
- `data_type`: file extension for saving (default: `".png"`)
- `save_path`: path to save plots (default: `nothing`)
- `is_save`: whether to save (default: `true`)
- `is_show`: whether to display (default: `true`)
- `use_tex`: use external `pdflatex` for rendering (default: `false`; ignored by Makie)
- `cl_over_cd`: plot CL/CD vs angle instead of CL vs CD (default: `true`)
"""
function plot_polars end

"""
    plot_polar_data(body_aero::BodyAerodynamics; kwargs...)

Plot polar data (Cl, Cd, Cm) as 3-D surfaces against angle of attack and trailing edge deflection.

# Arguments
- `body_aero`: the [`BodyAerodynamics`](@ref) to plot (must use `POLAR_MATRICES` aero model)

# Keyword arguments
- `alphas`: AoA values in radians (default: `deg2rad.(-5:0.3:25)`)
- `delta_tes`: trailing edge deflection angles in radians (default: `deg2rad.(-5:0.3:25)`)
- `is_show`: whether to display (default: `true`)
- `use_tex`: use external `pdflatex` for rendering (default: `false`; ignored by Makie)
"""
function plot_polar_data end

"""
    plot_combined_analysis(solver, body_aero, results; kwargs...)

Create a combined analysis by calling `plot_geometry`, `plot_distribution`, and `plot_polars`
in sequence.

# Arguments
- `solver`: solver or vector of solvers
- `body_aero`: [`BodyAerodynamics`](@ref) object or vector thereof
- `results`: results dictionary (or vector) from [`solve!`](@ref)

# Keyword arguments
- `solver_label`: label string for the solver (backward-compatible alias for `labels`)
- `labels`: optional label string or vector
- `angle_range`: range of angles for polar plots (default: `range(0, 20, length=20)`)
- `angle_type`: `"angle_of_attack"` or `"side_slip"` (default: `"angle_of_attack"`)
- `angle_of_attack`: AoA in degrees (default: `0.0`)
- `side_slip`: side slip angle in degrees (default: `0.0`)
- `v_a`: wind speed in m/s (default: `10.0`)
- `title`: overall figure title (default: `"Combined Analysis"`)
- `view_elevation`: geometry view elevation in degrees (default: `15`)
- `view_azimuth`: geometry view azimuth in degrees (default: `-120`)
- `is_show`: whether to display (default: `true`)
- `use_tex`: use external `pdflatex` for rendering (default: `false`; ignored by Makie)
- `literature_path_list`: paths to literature CSV files (default: `String[]`)
- `data_type`: file extension for saving (default: `".png"`)
- `save_path`: directory to save files (default: `nothing`)
- `is_save`: whether to save (default: `false`)
- `cl_over_cd`: plot CL/CD vs angle (default: `true`)
"""
function plot_combined_analysis end

"""
    plot_section_polars(body_aero::BodyAerodynamics, coefficient=:cl; kwargs...)

Plot one polar coefficient (`:cl`, `:cd`, or `:cm`) against angle of attack for
every section of a wing using stored `POLAR_VECTORS` data. Rendered through
`MakieControlPlots`.

# Arguments
- `body_aero`: the [`BodyAerodynamics`](@ref) to plot
- `coefficient`: `:cl`, `:cd`, or `:cm` (default: `:cl`)

# Keyword arguments
- `is_show`: whether to display (default: `true`)
- `is_save`: whether to save (default: `false`)
- `save_path`: directory to save the figure (default: `nothing`)
- `data_type`: file extension for saving (default: `".png"`)
"""
function plot_section_polars end

function save_plot end
function show_plot end

"""
   const MVec3    = MVector{3, Float64}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const MVec3    = MVector{3, Float64}
const MMat3    = MMatrix{3, 3, Float64, 9}

"""
   const PosVector=Union{MVec3, Vector}

Position vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const PosVector=Union{MVec3, Vector}

"""
   const VelVector=Union{MVec3, Vector}

Velocity vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const VelVector=Union{MVec3, Vector}

"""
   Model `VSM` `LLT`

Enumeration of the implemented model types.

# Elements
- VSM: Vortex Step Method
- LLT: Lifting Line Theory
"""
@enum Model VSM LLT

"""
    WingType `RECTANGULAR` `CURVED` `ELLIPTICAL`

Enumeration of the implemented wing types.

# Elements:
- RECTANGULAR
- CURVED
- ELLIPTICAL
"""
@enum WingType  RECTANGULAR CURVED ELLIPTICAL

"""
   AeroModel `POLY` `POLAR_VECTORS` `POLAR_MATRICES` `INVISCID`

Enumeration of the implemented aerodynamic models. See also: [`AeroData`](@ref)

# Elements
- `POLY`: α-polynomial coefficients for cl/cd/cm (e.g. Breukels LEI coeffs, generated
  by the `AirfoilAero` package). Core only evaluates the polynomial.
- `POLAR_VECTORS`: Polar vectors as function of alpha (lookup tables with interpolation).
  A panel's table may be rewritten at run time by [`set_polar!`](@ref
  VortexStepMethod.set_polar!), which is how a live polar source regenerates it from the
  panel's deformed shape each solve.
- `POLAR_MATRICES`: Polar matrices as function of alpha and delta (lookup tables with interpolation)
- INVISCID

`LEI_AIRFOIL_BREUKELS` is a deprecated alias of `POLY`.

where `alpha` is the angle of attack, `delta` is trailing edge angle.
"""
@enum AeroModel begin
   POLY
   POLAR_VECTORS
   POLAR_MATRICES
   INVISCID
end

"""
    LEI_AIRFOIL_BREUKELS

Deprecated alias of `POLY` (see [`AeroModel`](@ref)). The Breukels `(tube_diameter, camber)` → coeff
derivation now lives in `AirfoilAero.lei_poly_coeffs`; sections carry the resulting
`(cl_coeffs, cd_coeffs, cm_coeffs)`.
"""
const LEI_AIRFOIL_BREUKELS = POLY

"""
   PanelDistribution `LINEAR` `COSINE` `SPLIT_PROVIDED` `UNCHANGED` `BILLOWING`

Enumeration of the implemented panel distributions.

# Elements
- LINEAR               # Linear distribution
- COSINE               # Cosine distribution
- `SPLIT_PROVIDED`     # Split provided sections
- `UNCHANGED`          # 1:1 copy of unrefined to refined sections (no interpolation)
- `BILLOWING`          # Split provided + sinusoidal TE billowing between ribs
"""
@enum PanelDistribution begin
   LINEAR             # Linear distribution
   COSINE             # Cosine distribution
   SPLIT_PROVIDED     # Split provided sections
   UNCHANGED          # 1:1 copy of unrefined to refined sections
   BILLOWING          # Split provided + sinusoidal TE billowing
end

"""
   InitialGammaDistribution ELLIPTIC ZEROS

Enumeration of the implemented initial gamma distributions.

# Elements
- ELLIPTIC
- ZEROS
"""
@enum InitialGammaDistribution ELLIPTIC ZEROS

"""
   SolverStatus FEASIBLE INFEASIBLE FAILURE

Enumeration to report back the validity of the result of the solve! function.
Used in the [`VSMSolution`](@ref) struct.

# Elements
- FEASIBLE: The gamma distribution is physically feasible
- INFEASIBLE: The gamma distribution is physically infeasible
- FAILURE: The result did not converge within the maximal number of iterations
"""
@enum SolverStatus FEASIBLE INFEASIBLE FAILURE

"""
    SolverType

Enumeration specifying the method used to solve for circulation distribution.

# Values
- `LOOP`: Converging gamma loop - iterative approach that repeatedly updates circulation values until convergence
- `NONLIN`: Nonlinear solver - uses a numerical solver to solve the circulation equations
"""
@enum SolverType LOOP NONLIN

abstract type AbstractWing{T} end

"""
    AeroData= Union{
        Nothing,
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    }

Union of different definitions of the aerodynamic properties of a wing section. See also: [`AeroModel`](@ref)
  - nothing for INVISCID
  - (`cl_coeffs`, `cd_coeffs`, `cm_coeffs`) α-polynomial coefficients for `POLY`
  - (`alpha_range`, `cl_vector`, `cd_vector`, `cm_vector`) for `POLAR_VECTORS`
  - (`alpha_range`, `delta_range`, `cl_matrix`, `cd_matrix`, `cm_matrix`) for `POLAR_MATRICES`

where `alpha` is the angle of attack [rad], `delta` is trailing edge angle [rad], `cl` the lift coefficient,
`cd` the drag coefficient and `cm` the pitching moment coefficient. The camber of a kite refers to 
the curvature of its airfoil shape. The camber is typically measured as the maximum distance 
between the mean camber line (the line equidistant from the upper and lower surfaces) 
and the chord line of the airfoil.
"""
const AeroData = Union{
        Nothing,
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    }

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))

function menu()
   # Load the examples menu using a portable path
    ex = joinpath(PACKAGE_ROOT, "examples", "menu.jl")
   Base.include(Main, normpath(ex))
end

"""
    copy_examples()

Copy all example scripts to the folder "examples"
(it will be created if it doesn't exist).
"""
function copy_examples()
    PATH = "examples"
    if ! isdir(PATH)
        mkdir(PATH)
    end
    src_path = joinpath(PACKAGE_ROOT, PATH)
    copy_files(PATH, readdir(src_path))
end

function install_examples(add_packages=true)
    copy_examples()
    pkg_root = PACKAGE_ROOT
    src = joinpath(pkg_root, "data")
    isdir(src) && cp(src, "data"; force=true)
    if add_packages
        for pkg in ("GLMakie", "LaTeXStrings", "Xfoil", "CSV", "DataFrames")
            if !(pkg ∈ keys(Pkg.project().dependencies))
                Pkg.add(pkg)
            end
        end
    end
end

function copy_files(relpath, files)
    if ! isdir(relpath) 
        mkdir(relpath)
    end
    src_path = joinpath(PACKAGE_ROOT, relpath)
    for file in files
        cp(joinpath(src_path, file), joinpath(relpath, file), force=true)
        chmod(joinpath(relpath, file), 0o774)
    end
    files
end

function help(url) 
    if Sys.islinux()
        io = IOBuffer()
        run(pipeline(`xdg-open $url`, stderr = io))
        # ignore any error messages
        _ = take!(io)
    else
        DefaultApplication.open(url)
    end
    nothing
end

# Include core functionality
include("panel_aerodynamics.jl")
include("settings.jl")
include("section_aero.jl")
include("wing_geometry.jl")
include("polars.jl")
include("yaml_geometry.jl")
include("filament.jl")
include("panel.jl")
include("body_aerodynamics.jl")
include("wake.jl")
include("solver.jl")

include("plotting_helpers.jl")

# Airfoil-polar generation and OBJ-mesh conversion, folded in as internal submodules
# (formerly the AirfoilAero and ObjAdapter packages). Access as
# `VortexStepMethod.AirfoilAero` / `.ObjAdapter`, or `using VortexStepMethod.AirfoilAero`.
include("airfoil_aero/AirfoilAero.jl")
include("obj_adapter/ObjAdapter.jl")
include("surfplan_adapter/SurfplanAdapter.jl")

include("precompile.jl")


end # module

module VortexStepMethod

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
import NonlinearSolve: solve, solve!
using Interpolations
using Parameters
using Serialization
using Timers
using PreallocationTools
using PrecompileTools
using Pkg
using DifferentiationInterface
import YAML
using StructMapping
using Xfoil

# Export public interface
export SolverSettings, VSMSettings, WingSettings
export ObjWing, Section, Wing, refine!, reinit!
export BodyAerodynamics
export Solver, VSMSolution, linearize, solve, solve!, solve_base!
export calculate_results
export add_section!, set_va!
export calculate_projected_area, calculate_span
export MVec3
export LLT, Model, VSM
export AeroModel, INVISCID, LEI_AIRFOIL_BREUKELS, POLAR_MATRICES, POLAR_VECTORS
export COSINE, LINEAR, PanelDistribution, SPLIT_PROVIDED, UNCHANGED
export ELLIPTIC, InitialGammaDistribution, ZEROS
export FAILURE, FEASIBLE, INFEASIBLE, SolverStatus
export LOOP, NONLIN, SolverType
export load_polar_data

export plot_circulation_distribution, plot_combined_analysis, plot_distribution, plot_geometry,
    plot_polar_data, plot_polars, save_plot, show_plot

# the following functions are defined in ext/VortexStepMethodExt.jl
function plot_geometry end
function plot_distribution end
function plot_circulation_distribution end
function plot_polars end
function save_plot end
function show_plot end
function plot_polar_data end
function plot_combined_analysis end

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
   AeroModel `LEI_AIRFOIL_BREUKELS` `POLAR_VECTORS` `POLAR_MATRICES` `INVISCID`

Enumeration of the implemented aerodynamic models. See also: [AeroData](@ref)

# Elements
- `LEI_AIRFOIL_BREUKELS`: Polynom approximation for leading edge inflatable kites
- `POLAR_VECTORS`: Polar vectors as function of alpha (lookup tables with interpolation)
- `POLAR_MATRICES`: Polar matrices as function of alpha and delta (lookup tables with interpolation)
- INVISCID

where `alpha` is the angle of attack, `delta` is trailing edge angle.
"""
@enum AeroModel begin
   LEI_AIRFOIL_BREUKELS
   POLAR_VECTORS
   POLAR_MATRICES
   INVISCID
end

"""
   PanelDistribution `LINEAR` `COSINE` `SPLIT_PROVIDED` `UNCHANGED`

Enumeration of the implemented panel distributions.

# Elements
- LINEAR               # Linear distribution
- COSINE               # Cosine distribution
- `SPLIT_PROVIDED`     # Split provided sections
- `UNCHANGED`          # 1:1 copy of unrefined to refined sections (no interpolation)
"""
@enum PanelDistribution begin
   LINEAR             # Linear distribution
   COSINE             # Cosine distribution
   SPLIT_PROVIDED     # Split provided sections
   UNCHANGED          # 1:1 copy of unrefined to refined sections
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
Used in the [VSMSolution](@ref) struct.

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

abstract type AbstractWing end

"""
    AeroData= Union{
        Nothing,
        NTuple{2, Float64},
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    }

Union of different definitions of the aerodynamic properties of a wing section. See also: [AeroModel](@ref)
  - nothing for INVISCID
  - (`tube_diameter`, camber) for `LEI_AIRFOIL_BREUKELS`
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
        NTuple{2, Float64},
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
        String(take!(io)) 
    else
        DefaultApplication.open(url)
    end
    nothing
end

# Include core functionality
include("settings.jl")
include("wing_geometry.jl")
include("polars.jl")
include("obj_geometry.jl")
include("yaml_geometry.jl")
include("filament.jl")
include("panel.jl")
include("body_aerodynamics.jl")
include("wake.jl")
include("solver.jl")
include("precompile.jl")


end # module

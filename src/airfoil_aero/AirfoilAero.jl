module AirfoilAero

using LinearAlgebra
using Parameters
using Statistics
using Interpolations
using NPZ
using Xfoil
using Printf: @sprintf
using ..VortexStepMethod: SectionAero, interpolate_matrix_nans!

include("kulfan.jl")
include("shrink_wrap.jl")
include("neuralfoil.jl")
include("poly.jl")
include("airfoil_solvers/common.jl")
include("airfoil_solvers/xfoil_solver.jl")
include("airfoil_solvers/neuralfoil_solver.jl")
include("airfoil_io.jl")
include("polar_gen.jl")
include("polar_export.jl")
include("section_aero_gen.jl")

export KulfanParameters, KulfanFitMethod, LeastSquaresFit
export ShrinkWrap, shrink_wrap
export fit_kulfan_parameters, kulfan_to_coordinates
export NeuralFoilModel, NeuralFoilResult, load_neuralfoil_model
export neuralfoil_aero, neuralfoil_section
export AbstractAirfoilSolver, XFoilSolver, NeuralFoilSolver
export SectionSolution, DeformedSection, deform_section, analyze_section, analyze_sweep
export create_2d_polars, generate_aero_matrices, generate_section_aero, lei_poly_coeffs
export write_section_aero
export turn_trailing_edge!, get_lower_upper
export read_dat_coordinates, write_dat, write_polar_csv
export generate_polar_from_coordinates, generate_polar_from_dat, resolve_airfoil

end

module AirfoilAero

using LinearAlgebra
using Parameters
using Statistics
using Interpolations
using NPZ
using Xfoil
using VortexStepMethod: CpData, CpPolar, interpolate_matrix_nans!, write_aero_matrix

include("kulfan.jl")
include("neuralfoil.jl")
include("polar_gen.jl")
include("airfoil_solvers/common.jl")
include("airfoil_solvers/xfoil_solver.jl")
include("airfoil_solvers/neuralfoil_solver.jl")
include("cp_gen.jl")

export KulfanParameters, KulfanFitMethod, LeastSquaresFit, EnvelopeFit
export fit_kulfan_parameters, kulfan_to_coordinates
export NeuralFoilModel, NeuralFoilResult, load_neuralfoil_model
export neuralfoil_aero, neuralfoil_section
export AbstractAirfoilSolver, XFoilSolver, NeuralFoilSolver
export SectionSolution, DeformedSection, deform_section, analyze_section, analyze_sweep
export create_polars, generate_cp_polar
export turn_trailing_edge!, get_lower_upper

end

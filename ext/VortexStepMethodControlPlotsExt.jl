module VortexStepMethodControlPlotsExt
using ControlPlots, LaTeXStrings, VortexStepMethod, LinearAlgebra, Statistics, DelimitedFiles
import ControlPlots: plt
import VortexStepMethod: calculate_filaments_for_plotting

export plot_wing, plot_circulation_distribution, plot_geometry, plot_distribution, plot_polars, save_plot, show_plot, plot_matrix_data

include("../src/plotting.jl")

end
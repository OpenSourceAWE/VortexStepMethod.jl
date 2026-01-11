using GLMakie
using VortexStepMethod

# Configuration
OBJ_PATH = joinpath("data", "ram_air_kite", "ram_air_kite_body.obj")
OUTPUT_DIR = joinpath("data", "ram_air_kite", "polars_neuralfoil")
N_SLICES = 10
RE = 1e6  # Reynolds number

# Plot airfoil slices with Kulfan CST fits
println("Plotting airfoil slices from OBJ mesh...")
fig = plot_airfoil_slices(OBJ_PATH; n_slices=N_SLICES)
display(fig)

# Generate polars using NeuralFoil
println("\nGenerating NeuralFoil polars...")
generate_neuralfoil_polars(
    OBJ_PATH,
    OUTPUT_DIR;
    n_slices=N_SLICES,
    Re=RE,
    alpha_range=-10:0.5:20,
    model_size="xlarge",
    verbose=true
)

println("\nPolars written to: $OUTPUT_DIR")
nothing

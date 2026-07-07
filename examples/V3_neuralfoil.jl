using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using VortexStepMethod
using ObjAdapter
using AirfoilAero: EnvelopeFit
using LinearAlgebra

# Configuration
OBJ_PATH = joinpath("data", "TUDELFT_V3_KITE", "V3_25.obj")
POLARS_DIR = joinpath("data", "TUDELFT_V3_KITE", "polars_neuralfoil")
N_SLICES = 19  # Match number of CFD polars
RE = 1e6

# Reorient the mesh to the slicer's chord/span/up convention from its extents.
ROTATION = :auto

# Marched leading-edge stations across the span (finer => smoother edge trace).
N_BINS = 100

# The V3 canopy has an inflatable LE tube whose lower surface recirculates. The
# wrapping-distance constraint always encloses the tube; loosen tightness over it
# (x ≈ 0.05–0.35c) so the fit fairs smoothly over the tube instead of hugging the
# concave recirculation. The arc-length term keeps the rest of the airfoil smooth.
n_grid = 100
tightness_lower = ones(n_grid)
tightness_lower[5:35] .= 0.05
FIT_METHOD = EnvelopeFit(min_distance=0.005, tightness=10.0,
                         tightness_lower=tightness_lower)

# 3D slice diagnostic: mesh + LE/TE curves + section contours and their fits.
# Hover a slice to inspect its 2D airfoil (raw points + the envelope fit).
fig_slices_3d = plot_slices_3d(OBJ_PATH; n_slices=N_SLICES, rotation=ROTATION,
                               n_bins=N_BINS, fit_method=FIT_METHOD)
save("V3_slices_3d.png", fig_slices_3d)

# # Generate NeuralFoil polars if they don't exist
# if !isdir(POLARS_DIR) || length(readdir(POLARS_DIR)) < N_SLICES
    println("Generating NeuralFoil polars from OBJ mesh...")
    generate_neuralfoil_polars(
        OBJ_PATH,
        POLARS_DIR;
        n_slices=N_SLICES,
        Re=RE,
        rotation=ROTATION,
        n_bins=N_BINS,
        fit_method=FIT_METHOD,
        verbose=true
    )
# end

# Flight conditions
v_a = 10.0
angle_range = range(-5, 25, length=31)

# Load settings and create wing with CFD polars
println("\nCreating wing with CFD polars...")
settings_cfd = VSMSettings("TUDELFT_V3_KITE/vsm_settings.yaml")
wing_cfd = Wing(settings_cfd)
refine!(wing_cfd)
body_cfd = BodyAerodynamics([wing_cfd])
VortexStepMethod.reinit!(body_cfd)
solver_cfd = Solver(body_cfd, settings_cfd)

# Create modified settings for NeuralFoil polars
println("Creating wing with NeuralFoil polars...")

# Read and modify the YAML to use neuralfoil polars
using YAML
yaml_path = joinpath("data", "TUDELFT_V3_KITE", "aero_geometry.yaml")
geom = YAML.load_file(yaml_path)

# Update airfoil entries to use neuralfoil polars
for (i, airfoil) in enumerate(geom["wing_airfoils"]["data"])
    airfoil[2] = "polars"
    airfoil[3] = Dict("csv_file_path" => "polars_neuralfoil/$i.csv")
end

# Write the NeuralFoil geometry (every airfoil now points at a NeuralFoil polar)
# and load it as a normal YAML wing — a uniform POLAR_VECTORS geometry, so no
# manual section surgery is needed.
temp_yaml = joinpath("data", "TUDELFT_V3_KITE", "aero_geometry_neuralfoil.yaml")
YAML.write_file(temp_yaml, geom)

wing_nf = Wing(temp_yaml; n_panels=50, spanwise_distribution=LINEAR)
refine!(wing_nf)
body_nf = BodyAerodynamics([wing_nf])
VortexStepMethod.reinit!(body_nf)
solver_nf = Solver(body_nf, settings_cfd)

# Compute polars
function compute_polar(solver, body_aero, angle_range; v_a=10.0)
    cl = zeros(length(angle_range))
    cd = zeros(length(angle_range))

    for (i, aoa_deg) in enumerate(angle_range)
        α = deg2rad(aoa_deg)
        vel = [cos(α), 0.0, sin(α)] * v_a
        set_va!(body_aero, vel)
        results = VortexStepMethod.solve(solver, body_aero; log=false)
        cl[i] = results["cl"]
        cd[i] = results["cd"]
    end
    return cl, cd
end

println("\nComputing CFD polar...")
cl_cfd, cd_cfd = compute_polar(solver_cfd, body_cfd, angle_range; v_a)

println("Computing NeuralFoil polar...")
cl_nf, cd_nf = compute_polar(solver_nf, body_nf, angle_range; v_a)

# Plot comparison
fig = Figure(size=(1200, 500))

ax1 = Axis(fig[1, 1]; xlabel="α [°]", ylabel="CL [-]", title="Lift Coefficient")
lines!(ax1, collect(angle_range), cl_cfd; label="CFD polars", linewidth=2)
lines!(ax1, collect(angle_range), cl_nf; label="NeuralFoil polars", linewidth=2)
axislegend(ax1; position=:rb)

ax2 = Axis(fig[1, 2]; xlabel="α [°]", ylabel="CD [-]", title="Drag Coefficient")
lines!(ax2, collect(angle_range), cd_cfd; label="CFD polars", linewidth=2)
lines!(ax2, collect(angle_range), cd_nf; label="NeuralFoil polars", linewidth=2)
axislegend(ax2; position=:lt)

ax3 = Axis(fig[1, 3]; xlabel="CD [-]", ylabel="CL [-]", title="Drag Polar")
lines!(ax3, cd_cfd, cl_cfd; label="CFD polars", linewidth=2)
lines!(ax3, cd_nf, cl_nf; label="NeuralFoil polars", linewidth=2)
axislegend(ax3; position=:rb)

Label(fig[0, :], "TU Delft V3 Kite: CFD vs NeuralFoil (v=$v_a m/s, Re=$RE)"; fontsize=20)

display(fig)
nothing

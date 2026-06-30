using Pkg
if Base.active_project() \!= joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using MakieControlPlots
using VortexStepMethod
using LinearAlgebra

# Configuration
OBJ_PATH = joinpath("data", "TUDELFT_V3_KITE", "V3_25.obj")
POLARS_DIR = joinpath("data", "TUDELFT_V3_KITE", "polars_neuralfoil")
N_SLICES = 19  # Match number of CFD polars
RE = 1e6

# Generate NeuralFoil polars if they don't exist
if !isdir(POLARS_DIR) || length(readdir(POLARS_DIR)) < N_SLICES
    println("Generating NeuralFoil polars from OBJ mesh...")

    # Plot airfoil slices to verify
    fig_slices = plot_airfoil_slices(OBJ_PATH; n_slices=N_SLICES, is_show=false)
    save("V3_slices.png", fig_slices)
    println("Saved V3_slices.png")

    generate_neuralfoil_polars(
        OBJ_PATH,
        POLARS_DIR;
        n_slices=N_SLICES,
        Re=RE,
        verbose=true
    )
end

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

# Write temporary modified geometry
temp_yaml = joinpath("data", "TUDELFT_V3_KITE", "aero_geometry_neuralfoil.yaml")
YAML.write_file(temp_yaml, geom)

# Create settings pointing to modified geometry
settings_nf = VSMSettings("TUDELFT_V3_KITE/vsm_settings.yaml")
# Manually update the geometry file reference in the wing
wing_nf = Wing(settings_nf)
# Replace sections with NeuralFoil polars
for (i, section) in enumerate(wing_nf.unrefined_sections)
    polar_path = joinpath(POLARS_DIR, "$i.csv")
    if isfile(polar_path)
        aero_data, _ = load_polar_data(polar_path)
        wing_nf.unrefined_sections[i] = Section(
            section.LE_point, section.TE_point, POLAR_VECTORS, aero_data
        )
    end
end
refine!(wing_nf)
body_nf = BodyAerodynamics([wing_nf])
VortexStepMethod.reinit!(body_nf)
solver_nf = Solver(body_nf, settings_nf)

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
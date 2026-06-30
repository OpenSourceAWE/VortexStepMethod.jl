"""
Convert a 3D wing `.obj` mesh to the native YAML geometry route with the
NeuralFoil adapter, then load and solve it like any other YAML wing.

The adapter slices the mesh, fits each section to Kulfan parameters, evaluates
NeuralFoil, and writes `airfoils/*.dat`, `polars/*.csv`, and a `geometry.yaml`.

This example fits each section with `EnvelopeFit`, which wraps tightly around the
slice points with clearance `min_distance` instead of least-squaring through
them. The envelope is robust to the noisy interior-structure points (ribs, spars)
of a ram-air kite slice that otherwise pull the default `LeastSquaresFit` inward.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using VortexStepMethod
using LinearAlgebra

PLOT = true
USE_TEX = false
v_a = 15.0
project_dir = dirname(@__DIR__)
obj_path = joinpath(project_dir, "data", "ram_air_kite", "ram_air_kite_body.obj")
output_dir = joinpath(project_dir, "output", "ram_air_kite_converted")

# --- Convert the .obj mesh to the YAML route ---
println("Converting $(basename(obj_path)) to YAML...")
geometry_file = obj_to_yaml(obj_path, output_dir;
    n_sections=10, Re=5e5, alpha_range=-20:1:20, model_size="xlarge",
    fit_method=EnvelopeFit(min_distance=0.0001))

# --- Everything below is the standard YAML wing route ---
wing = Wing(geometry_file; n_panels=20, spanwise_distribution=LINEAR)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)

solver = Solver(body_aero; aerodynamic_model_type=VSM, rtol=1e-5, solver_type=NONLIN)
set_va!(body_aero, [cos(deg2rad(8)) * v_a, 0.0, sin(deg2rad(8)) * v_a])
results = VortexStepMethod.solve(solver, body_aero; log=true)

if PLOT
    plot_geometry(body_aero, "Ram air kite (converted from .obj)"; is_show=true,
        view_elevation=15, view_azimuth=-120, use_tex=USE_TEX)

    # Airfoils and per-section polars recovered from the converted geometry
    plot_airfoils(geometry_file; is_show=true)
    plot_section_polars(body_aero, :cl; is_show=true)
    plot_section_polars(body_aero, :cd; is_show=true)

    plot_polars([solver], [body_aero], ["VSM (NeuralFoil polars from .obj)"];
        angle_range=range(-5, 20, length=26), v_a=v_a,
        title="Ram air kite: obj_to_yaml route", is_save=false, use_tex=USE_TEX)
end

nothing

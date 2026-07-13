"""
Convert a 3D wing `.obj` mesh to the native YAML geometry route with the
NeuralFoil adapter, then load and solve it like any other YAML wing.

The adapter slices the mesh, shrink-wraps each section into a clean airfoil,
evaluates NeuralFoil, and writes `airfoils/*.dat`, `polars/*.csv`, and a
`geometry.yaml`.

`ShrinkWrap` wraps each slice's point cloud into a clean closed airfoil with a
`clearance` gap, robust both to the noisy interior-structure points (ribs, spars)
of a ram-air kite slice that otherwise pull a plain least-squares fit inward, and
to the complex, harsh corners of an LEI kite slice.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using GLMakie
using MakieControlPlots
using VortexStepMethod
using VortexStepMethod.ObjAdapter
using VortexStepMethod.AirfoilAero: ShrinkWrap, NeuralFoilSolver
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
    n_sections=36, Re=5e5, alpha_range=-20:1:20,
    aero_solver=NeuralFoilSolver(model_size="xlarge"), wrap_method=ShrinkWrap())

# --- Report per-section fit error: how far any raw point sticks OUTSIDE the
# fitted envelope (interior structure points are ignored by the envelope fit). ---
function max_outside(x_fit, y_fit, x_raw, y_raw)
    isempty(x_raw) && return 0.0
    le = argmin(x_fit)
    xu, yu = reverse(x_fit[1:le]), reverse(y_fit[1:le])
    xl, yl = x_fit[le:end], y_fit[le:end]
    interp(xs, ys, px) = begin
        px <= xs[1] && return ys[1]
        px >= xs[end] && return ys[end]
        k = searchsortedfirst(xs, px)
        t = (px - xs[k-1]) / (xs[k] - xs[k-1])
        ys[k-1] + t * (ys[k] - ys[k-1])
    end
    v = 0.0
    for (px, py) in zip(x_raw, y_raw)
        uy = interp(xu, yu, px)
        ly = interp(xl, yl, px)
        v = max(v, py > uy ? py - uy : (py < ly ? ly - py : 0.0))
    end
    return v
end

println("Per-section fit error (max raw point outside the fitted envelope):")
for af in airfoils_from_yaml(geometry_file)
    err = max_outside(af.x, af.y, af.x_raw, af.y_raw)
    println("  Airfoil $(af.id): max_outside=$(round(err, digits=4))" *
            (err > 0.01 ? "  <-- poor fit" : ""))
end

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
    plot_airfoils(geometry_file; symmetric=true, is_show=true)
    plot_section_polars(body_aero, :cl; is_show=true)
    plot_section_polars(body_aero, :cd; is_show=true)

    plot_polars([solver], [body_aero], ["VSM (NeuralFoil polars from .obj)"];
        angle_range=range(-5, 20, length=26), v_a=v_a,
        title="Ram air kite: obj_to_yaml route", is_save=false, use_tex=USE_TEX)
end

nothing

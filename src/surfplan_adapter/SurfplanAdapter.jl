module SurfplanAdapter

using LinearAlgebra
import YAML
using ..AirfoilAero: shrink_wrap, ShrinkWrap, read_dat_coordinates,
    generate_airfoils, write_geometry_yaml, NeuralFoilSolver, AbstractAirfoilSolver

"""
    surfplan_to_aero_yaml(adapter_dir, output_dir; aero_solver=NeuralFoilSolver(),
        wrap_method=ShrinkWrap(), alpha_range=-180:1:180, delta_range=nothing,
        Re=nothing, crease_frac=0.75, force=false, verbose=true) -> geometry_yaml_path

Generate a pressure-ready `geometry.yaml` (per-node surface `cp`/`cf` tables plus
polars) from a SurfplanAdapter aero export, so the wing can be flown with the
`AeroPressure` continuous coupling instead of only integrated polars. The Surfplan
counterpart of [`obj_to_yaml`](@ref): where `.obj` slices a mesh, this reads the
already-clean per-rib airfoil `.dat` profiles the Python adapter exported.

Reads `adapter_dir/aero_geometry.yaml` for each section's leading/trailing-edge
placement (`wing_sections`) and each airfoil's `.dat` path (`wing_airfoils`
`info_dict.dat_file_path`), shrink-wraps every unique profile, and runs the shared
[`generate_airfoils`](@ref) core with `aero_solver` ([`NeuralFoilSolver`](@ref) by
default, [`XFoilSolver`](@ref) opt-in). Sections that share an airfoil generate its
tables once. `Re` defaults to the export's `wing_airfoils.reynolds`; `alpha_range`
defaults to the full `-180:1:180` sweep rather than the export's narrow polar range.
`table_format` writes the per-node surface tables as `:csv` (default, readable) or
`:arrow` (binary, an order of magnitude faster to load).

The `.txt` → adapter-YAML step (the upstream Python `SurfplanAdapter`) is the
documented prerequisite. Load the result with `Wing(geometry_yaml_path)`.

An existing `geometry.yaml` in `output_dir` is reused as-is, skipping the expensive
polar generation; set `force=true` to regenerate the polars.
"""
function surfplan_to_aero_yaml(adapter_dir::AbstractString, output_dir::AbstractString;
        aero_solver::AbstractAirfoilSolver=NeuralFoilSolver(),
        wrap_method::ShrinkWrap=ShrinkWrap(),
        alpha_range=-180:1:180, delta_range=nothing, Re=nothing,
        crease_frac=0.75, force::Bool=false, verbose::Bool=true,
        table_format::Symbol=:csv)
    yaml_path = joinpath(output_dir, "geometry.yaml")
    if !force && isfile(yaml_path)
        verbose && @info "Reusing existing $yaml_path (pass force=true to regenerate)"
        return yaml_path
    end
    data = YAML.load_file(joinpath(adapter_dir, "aero_geometry.yaml"))
    wing_sections = data["wing_sections"]["data"]
    wing_airfoils = data["wing_airfoils"]
    Re = something(Re, Float64(get(wing_airfoils, "reynolds", 1.0e6)))

    dat_of = Dict{Int, String}()
    for row in wing_airfoils["data"]
        id, info = Int(row[1]), row[3]
        path = get(info, "dat_file_path", "")
        isempty(path) &&
            error("wing_airfoil $id has no dat_file_path (need airfoil profiles)")
        dat_of[id] = isabspath(path) ? path : joinpath(adapter_dir, path)
    end

    sections = [(; airfoil_id = Int(row[1]),
                 LE_point = Float64.(row[2:4]), TE_point = Float64.(row[5:7]))
                for row in wing_sections]

    airfoils = NamedTuple[]
    for id in sort(unique(s.airfoil_id for s in sections))
        haskey(dat_of, id) ||
            error("section airfoil_id $id is absent from wing_airfoils")
        x_raw, y_raw = read_dat_coordinates(dat_of[id])
        x_fit, y_fit = shrink_wrap(x_raw, y_raw, wrap_method)
        push!(airfoils, (; id, x_fit, y_fit, x_raw, y_raw))
    end

    airfoil_rows, ok = generate_airfoils(airfoils, output_dir; Re, alpha_range,
        delta_range, aero_solver, crease_frac, verbose, table_format)
    isempty(ok) && error("No airfoil produced a valid polar from $adapter_dir")

    remap(id) = id in ok ? id : ok[argmin(abs.(ok .- id))]
    section_rows = Vector{Any}[]
    for s in sections
        fid = remap(s.airfoil_id)
        push!(section_rows, Any[fid, s.LE_point[1], s.LE_point[2], s.LE_point[3],
                                s.TE_point[1], s.TE_point[2], s.TE_point[3]])
    end
    sort!(section_rows; by = row -> row[3])
    sort!(airfoil_rows; by = row -> row[1])
    write_geometry_yaml(yaml_path, section_rows, airfoil_rows)
    verbose && @info "Wrote pressure geometry to $yaml_path " *
        "($(length(section_rows)) sections, $(length(ok)) airfoils)"
    return yaml_path
end

export surfplan_to_aero_yaml

end

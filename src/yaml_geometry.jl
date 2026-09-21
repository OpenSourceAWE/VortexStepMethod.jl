# Data structures for YAML wing geometry
@with_kw struct WingAirfoilInfo
    csv_file_path::String
    dat_file::String = ""
    cp_file::String = ""
    cf_file::String = ""
    cl_file_path::String = ""
    cd_file_path::String = ""
    cm_file_path::String = ""
end

@with_kw struct WingSectionData
    airfoil_id::Int64
    LE_x::Float64
    LE_y::Float64
    LE_z::Float64
    TE_x::Float64
    TE_y::Float64
    TE_z::Float64
end

span_order_key(section::WingSectionData) = section.LE_y

@with_kw struct WingAirfoilData
    airfoil_id::Int64
    type::String
    info_dict::WingAirfoilInfo
end

"""
    assemble_polar_matrix(alpha, delta, coefficients) -> (alphas, deltas, grids...)

Scatter long-format polar rows (radian `alpha`/`delta` per row, one coefficient per
column of `coefficients`) onto the sorted unique `alphas × deltas` grid, one matrix per
coefficient. Grid points without a row stay `NaN`.
"""
function assemble_polar_matrix(alpha, delta, coefficients::AbstractMatrix)
    alphas = sort(unique(alpha))
    deltas = sort(unique(delta))
    ai = Dict(a => i for (i, a) in enumerate(alphas))
    di = Dict(d => j for (j, d) in enumerate(deltas))
    grids = [fill(NaN, length(alphas), length(deltas)) for _ in axes(coefficients, 2)]
    for k in eachindex(alpha), (c, grid) in enumerate(grids)
        grid[ai[alpha[k]], di[delta[k]]] = coefficients[k, c]
    end
    return (alphas, deltas, grids...)
end

"""
    load_polar_data(path::String) -> (aero_data, model)

Load an airfoil polar table, CSV or Arrow as the suffix of `path` says (see
[`read_node_table`](@ref)), with columns `alpha`, `cl`, `cd` and `cm` (case-insensitive,
any order, alpha in degrees). A `delta` column makes it a long-format `POLAR_MATRICES`
table, returned as `(alpha, delta, cl, cd, cm)` on its sorted grid; without one it is
`POLAR_VECTORS` `(alpha, cl, cd, cm)`. Angles come back in radians. A missing, empty or
malformed file warns and returns `(nothing, INVISCID)`.
"""
function load_polar_data(path::String)
    if !isfile(path)
        @warn "Polar file not found: \"$path\""
        return (nothing, INVISCID)
    end
    alpha, delta, values, columns = try
        read_node_table(path)
    catch e
        @warn "Error reading polar file $path: $e"
        return (nothing, INVISCID)
    end
    wanted = [findfirst(==(name), lowercase.(columns)) for name in ("cl", "cd", "cm")]
    if any(isnothing, wanted) || isempty(alpha)
        @warn "Polar file $path needs alpha, cl, cd and cm columns and a data row"
        return (nothing, INVISCID)
    end
    coefficients = values[:, wanted]
    delta === nothing ||
        return (assemble_polar_matrix(alpha, delta, coefficients), POLAR_MATRICES)
    return (alpha, coefficients[:, 1], coefficients[:, 2], coefficients[:, 3]),
        POLAR_VECTORS
end

"""
    load_matrix_polar_data(cl_path, cd_path, cm_path) -> (aero_data, POLAR_MATRICES)

Read the three `(alpha × delta)` coefficient matrices (see [`read_aero_matrix`](@ref))
and assemble the `POLAR_MATRICES` `aero_data = (alpha, delta, cl, cd, cm)`.
"""
function load_matrix_polar_data(cl_path::String, cd_path::String, cm_path::String)
    cl_matrix, alpha, delta = read_aero_matrix(cl_path)
    cd_matrix, _, _ = read_aero_matrix(cd_path)
    cm_matrix, _, _ = read_aero_matrix(cm_path)
    return (collect(alpha), collect(delta), cl_matrix, cd_matrix, cm_matrix),
        POLAR_MATRICES
end

"""
    Wing(geometry_file::String; n_panels=20, spanwise_distribution=LINEAR,
         spanwise_direction=[0.0, 1.0, 0.0], remove_nan=true,
         use_prior_polar=false, billowing_percentage=0.0, prn=false)

Constructs a `Wing` object from a YAML geometry file.

# Arguments
- `geometry_file::String`: Path to a YAML file describing the wing geometry and airfoils.

# Keyword Arguments
- `n_panels::Int`: Number of spanwise panels (default: 20).
- `spanwise_distribution`: Spanwise panel distribution type (default: `LINEAR`).
- `spanwise_direction::Vector{Float64}`: Direction of the spanwise axis (default: `[0.0, 1.0, 0.0]`). Must be the global Y axis.
- `remove_nan::Bool`: Remove NaN values from the geometry (default: `true`).
- `use_prior_polar::Bool`: Reuse prior refined/panel polar mapping on geometry updates (default: `false`).
- `billowing_percentage::Float64`: TE billow as percentage of arc length (default: `0.0`).
- `prn::Bool`: Print informational messages during construction (default: `false`).

# Returns
- `Wing`: A fully constructed and initialized `Wing` object.

# Description
This function reads a YAML configuration file to define the geometry and airfoil data for a multi-section wing.
Each section and corresponding airfoil is parsed from the YAML file, polar data is loaded, and each section is added
to the wing. The geometry logic currently assumes the spanwise direction is `[0.0, 1.0, 0.0]` (aligned with the global Y axis).
The number of unrefined sections is automatically inferred from the sections in the geometry file.

# Errors
- Throws an `ArgumentError` if `spanwise_direction` is not `[0.0, 1.0, 0.0]`.

# Example
```julia
wing = Wing("wing_geometry.yaml"; n_panels=30, prn=true)
```
"""
function Wing(
    geometry_file::String;
    n_panels=20,
    spanwise_distribution=LINEAR,
    spanwise_direction=[0.0, 1.0, 0.0],
    remove_nan=true,
    use_prior_polar=false,
    billowing_percentage=0.0,
    crease_frac=0.75,
    prn=false,
    sort_sections=true
)

    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))

    prn && @info "Reading YAML wing configuration from $geometry_file"
    
    # Load YAML file following Uwe's suggestion
    data = YAML.load_file(geometry_file)

    # Convert wing sections
    wing_sections_data = data["wing_sections"]
    sections = WingSectionData[]
    for row in wing_sections_data["data"]
        section_dict = Dict(zip(wing_sections_data["headers"], row))
        push!(sections, WingSectionData(
            airfoil_id = section_dict["airfoil_id"],
            LE_x = section_dict["LE_x"],
            LE_y = section_dict["LE_y"], 
            LE_z = section_dict["LE_z"],
            TE_x = section_dict["TE_x"],
            TE_y = section_dict["TE_y"],
            TE_z = section_dict["TE_z"]
        ))
    end
    normalize_span_order!(sections)

    # Convert wing airfoils
    wing_airfoils_data = data["wing_airfoils"]
    airfoils = WingAirfoilData[]
    airfoil_poly_map = Dict{Int64, NTuple{3, Vector{Float64}}}()
    airfoil_type_map = Dict{Int64, String}()
    for row in wing_airfoils_data["data"]
        airfoil_dict = Dict(zip(wing_airfoils_data["headers"], row))
        airfoil_type_map[airfoil_dict["airfoil_id"]] = String(airfoil_dict["type"])
        info = airfoil_dict["info_dict"]
        if haskey(info, "cl_coeffs")
            airfoil_poly_map[airfoil_dict["airfoil_id"]] = (
                Float64.(info["cl_coeffs"]), Float64.(info["cd_coeffs"]),
                Float64.(info["cm_coeffs"]))
        end
        push!(airfoils, WingAirfoilData(
            airfoil_id = airfoil_dict["airfoil_id"],
            type = airfoil_dict["type"],
            info_dict = WingAirfoilInfo(
                csv_file_path = get(airfoil_dict["info_dict"], "csv_file_path", ""),
                dat_file = get(airfoil_dict["info_dict"], "dat_file", ""),
                cp_file = get(airfoil_dict["info_dict"], "cp_file", ""),
                cf_file = get(airfoil_dict["info_dict"], "cf_file", ""),
                cl_file_path = get(airfoil_dict["info_dict"], "cl_file_path", ""),
                cd_file_path = get(airfoil_dict["info_dict"], "cd_file_path", ""),
                cm_file_path = get(airfoil_dict["info_dict"], "cm_file_path", ""))
        ))
    end

    # Create CSV file mapping from airfoils
    airfoil_csv_map = Dict{Int64, String}()
    airfoil_surface_map = Dict{Int64, NTuple{3, String}}()
    airfoil_matrix_map = Dict{Int64, NTuple{3, String}}()
    for airfoil in airfoils
        if !isempty(airfoil.info_dict.csv_file_path)
            airfoil_csv_map[airfoil.airfoil_id] = airfoil.info_dict.csv_file_path
        end
        if !isempty(airfoil.info_dict.cp_file) && !isempty(airfoil.info_dict.cf_file)
            airfoil_surface_map[airfoil.airfoil_id] = (airfoil.info_dict.dat_file,
                airfoil.info_dict.cp_file, airfoil.info_dict.cf_file)
        end
        if !isempty(airfoil.info_dict.cl_file_path)
            airfoil_matrix_map[airfoil.airfoil_id] = (airfoil.info_dict.cl_file_path,
                airfoil.info_dict.cd_file_path, airfoil.info_dict.cm_file_path)
        end
    end
    
    # n_unrefined_sections is set automatically as sections are added
    wing = Wing(n_panels;
        spanwise_distribution=spanwise_distribution,
        spanwise_direction=MVec3(spanwise_direction),
        remove_nan=remove_nan,
        use_prior_polar=use_prior_polar,
        billowing_percentage=Float64(billowing_percentage),
        crease_frac=Float64(crease_frac)
    )

    # Parse sections and populate wing
    for section in sections
        # Get coordinates directly from struct fields
        le_coord = [section.LE_x, section.LE_y, section.LE_z]
        te_coord = [section.TE_x, section.TE_y, section.TE_z]

        base_dir = dirname(geometry_file)
        resolve(p) = (!isempty(p) && !isabspath(p)) ? joinpath(base_dir, p) : p

        # Core accepts only resolved forms (poly/polars/inviscid); rich types resolve in AirfoilAero.
        airfoil_type = get(airfoil_type_map, section.airfoil_id, "")
        if haskey(airfoil_poly_map, section.airfoil_id)
            aero_data = airfoil_poly_map[section.airfoil_id]
            aero_model = POLY
        elseif haskey(airfoil_matrix_map, section.airfoil_id)
            cl_p, cd_p, cm_p = airfoil_matrix_map[section.airfoil_id]
            aero_data, aero_model = load_matrix_polar_data(
                resolve(cl_p), resolve(cd_p), resolve(cm_p))
        elseif airfoil_type in ("neuralfoil", "breukels_regression", "masure_regression")
            throw(ArgumentError("airfoil_id $(section.airfoil_id) has unresolved " *
                "type \"$airfoil_type\"; resolve it with AirfoilAero " *
                "(resolve_aero_geometry) into polars/poly before loading."))
        else
            csv_file_path = resolve(get(airfoil_csv_map, section.airfoil_id, ""))
            aero_data, aero_model = load_polar_data(csv_file_path)
        end

        surface = get(airfoil_surface_map, section.airfoil_id, nothing)
        section_aero = isnothing(surface) ? nothing :
            read_section_aero(resolve(surface[1]), resolve(surface[2]), resolve(surface[3]))

        prn && println("Section airfoil_id $(section.airfoil_id): Using $aero_model model")

        add_section!(wing, le_coord, te_coord, aero_model, aero_data, section_aero)
    end

    refine!(wing; sort_sections)
    return wing
end

"""
    Wing(settings::VSMSettings)

Create a wing model from VSM settings configuration.

This constructor is a convenience wrapper that extracts wing configuration
from VSMSettings and creates a Wing using either:
- YAML geometry file (geometry_file field), or
- OBJ + DAT files (obj_file and dat_file fields)

The constructor automatically determines which path to use based on which
fields are populated in the settings.

# Arguments
- `settings`: VSMSettings object containing wing configuration

# Returns
A fully initialized `Wing` instance ready for aerodynamic simulation.

# Example
```julia
# Using YAML geometry
settings = VSMSettings("path/to/vsm_settings.yaml")
wing = Wing(settings)

# Settings can specify either:
# - geometry_file: "path/to/wing.yaml"  # YAML-based
# - obj_file + dat_file                  # OBJ-based
```
"""
function Wing(settings::VSMSettings; sort_sections::Bool=true)
    wing_settings = settings.wings[1]

    # Check which geometry format to use
    has_yaml = !isempty(wing_settings.geometry_file)
    has_obj = !isempty(wing_settings.obj_file)
    has_dat = !isempty(wing_settings.dat_file)

    if has_yaml && (has_obj || has_dat)
        throw(ArgumentError(
            "Cannot specify both geometry_file and obj_file/dat_file"
        ))
    end

    if has_obj && !has_dat
        throw(ArgumentError(
            "obj_file requires dat_file to be specified"
        ))
    end

    if has_dat && !has_obj
        throw(ArgumentError(
            "dat_file requires obj_file to be specified"
        ))
    end

    if has_yaml
        # Use YAML geometry constructor
        Wing(wing_settings.geometry_file;
            n_panels=wing_settings.n_panels,
            spanwise_distribution=wing_settings.spanwise_panel_distribution,
            remove_nan=wing_settings.remove_nan,
            use_prior_polar=wing_settings.use_prior_polar,
            billowing_percentage=wing_settings.billowing_percentage,
            crease_frac=wing_settings.crease_frac,
            sort_sections
        )
    elseif has_obj && has_dat
        throw(ArgumentError(
            "OBJ/DAT geometry is handled by the ObjAdapter package: convert to a " *
            "standard YAML with ObjAdapter first, then load it via geometry_file."
        ))
    else
        throw(ArgumentError(
            "WingSettings must specify either geometry_file or " *
            "both obj_file and dat_file"
        ))
    end
end

"""
    n_unrefined_sections(wing_settings::WingSettings) -> Int

Number of unrefined sections the [`Wing`](@ref) built from `wing_settings` carries: the
rows of its `geometry_file`'s `wing_sections`.
"""
function n_unrefined_sections(wing_settings::WingSettings)
    isempty(wing_settings.geometry_file) && throw(ArgumentError(
        "wing \"$(wing_settings.name)\" has no geometry_file to count its sections from"))
    return length(YAML.load_file(wing_settings.geometry_file)["wing_sections"]["data"])
end

"""
    ObjWing(obj_path[, dat_path]; n_panels, Re, alpha_range, delta_range,
            n_sections, spanwise_direction, aero_solver, remake, output_dir,
            crease_frac, verbose) → Wing

Convenience constructor retained for backward compatibility. Converts an OBJ mesh
to a YAML wing geometry via [`ObjAdapter.obj_to_yaml`](@ref) and returns a [`Wing`](@ref).

`dat_path` is accepted but ignored — airfoil shapes are extracted directly from the
OBJ geometry. Use `aero_solver=AirfoilAero.XFoilSolver()` to reproduce old XFoil-based
polars; the default is `AirfoilAero.NeuralFoilSolver()`.

`alpha_range` and `delta_range` are in degrees (matching [`ObjAdapter.obj_to_yaml`](@ref)).

By default (`remake=false`) an existing `geometry.yaml` in `output_dir` is reused,
skipping the expensive polar generation. Set `remake=true` to force regeneration.
"""
function ObjWing(obj_path, dat_path=nothing;
                 n_panels::Int=56,
                 Re::Real=1e6,
                 alpha_range=-5:1:20,
                 delta_range=-5:1:20,
                 n_sections::Union{Nothing, Int}=nothing,
                 spanwise_direction=[0.0, 1.0, 0.0],
                 spanwise_distribution=UNCHANGED,
                 remove_nan::Bool=true,
                 aero_solver=AirfoilAero.NeuralFoilSolver(),
                 remake::Bool=false,
                 output_dir::String=mktempdir(),
                 crease_frac=0.75,
                 verbose::Bool=false)
    yaml_path = joinpath(output_dir, "geometry.yaml")
    if remake || !isfile(yaml_path)
        n_sec = isnothing(n_sections) ? n_panels + 1 : n_sections
        yaml_path = ObjAdapter.obj_to_yaml(obj_path, output_dir;
                                           n_sections=n_sec,
                                           Re,
                                           alpha_range,
                                           delta_range,
                                           aero_solver,
                                           spanwise_direction,
                                           crease_frac,
                                           verbose)
    end
    return Wing(yaml_path; n_panels, spanwise_distribution, spanwise_direction, remove_nan)
end

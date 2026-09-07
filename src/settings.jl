function parse_enum(::Type{T}, s::String) where T <: Enum
    for instance in instances(T)
        if string(instance) == s
            return instance
        end
    end
    throw(ArgumentError("Invalid $(T) value: \"$s\". " *
          "Valid values: $(join(string.(instances(T)), ", "))"))
end

@with_kw mutable struct ConditionSettings
    wind_speed::Float64 = 10.0      # wind speed [m/s]
    alpha::Float64 = 5.0            # angle of attack [°]
    beta::Float64 = 0.0             # sideslip angle [°]
    yaw_rate::Float64 = 0.0         # yaw rate [°/s]
end

"""
    MeshSettings

How a wing's `.obj` mesh becomes airfoil slices, used within [`WingSettings`](@ref):
which mesh, and how to cut it.

`obj_file` here is the mesh the sections are *generated from*, an input to the
`geometry_file` a wing then flies. That is not the wing's own `obj_file`, which is
an alternative to `geometry_file` — a wing built straight from an obj and a dat,
with no polar generation in between — and which cannot be given alongside one.

Apart from `n_sections`, which `obj_to_yaml` requires, every default here is the one
that call and [`ShrinkWrap`](@ref) already apply, so a wing naming no `mesh:` block
slices as an unconfigured call.

# Fields
- `obj_file`: Mesh the sections are sliced from, relative to the data directory
    (default `""`, no mesh).
- `n_sections`: Sections sliced from the mesh (default `45`).
- `n_bins`: Leading-edge stations marched across the span; more gives a smoother
    trace (default `60`).
- `rotation`: Rows of the mesh-to-slicer rotation, which brings the mesh into the
    slicer's convention of x = chord, y = span, z = up (default the identity).
- `wingtip_distance`: Arc length [m] the outermost sections stop short of each tip
    (default `0.0`).
- `clearance`: Shrink-wrap offset [chord fraction] the contour holds outside every
    cloud point, and the radius its convex corners are rounded at (default
    `0.006`).
- `min_concave_radius`: Shrink-wrap rolling-ball radius [chord fraction], bridging
    gaps and crevices in the point cloud (default `0.02`).
"""
@with_kw mutable struct MeshSettings
    obj_file::String = ""
    n_sections::Int64 = 45
    n_bins::Int64 = 60
    rotation::Vector{Vector{Float64}} = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                         [0.0, 0.0, 1.0]]
    wingtip_distance::Float64 = 0.0
    clearance::Float64 = 0.006
    min_concave_radius::Float64 = 0.02
end

"""
    AirfoilSettings

The 2D section backend and the polars it tabulates, used within
[`WingSettings`](@ref). One block answers for both the tables a mesh is sliced into
and the live polars a deformed section is re-solved on, so the two cannot be
generated at different transition settings or off different networks.

# Fields
- `solver`: Section backend, `"neuralfoil"` or `"xfoil"` (default `"neuralfoil"`).
- `model_size`: NeuralFoil network size (default `"large"`).
- `n_crit`: e^N transition criticality; lower is dirtier, so transition is earlier
    (default `9.0`, the standard clean-tunnel value).
- `xtr_upper` / `xtr_lower`: Forced transition as a chord fraction, applied to
    whichever backend `solver` names (default `0.05`).
- `alpha_range`: Angle-of-attack sweep [deg] as `[first, step, last]`
    (default `[-180, 1, 180]`).
- `delta_range`: Flap-deflection sweep [deg] as `[first, step, last]`; `nothing`
    for a dataset generated without a flap sweep (default `nothing`).
- `live_offsets`: Angles [deg] off the reference angle a live polar is sampled at
    (default `-12:3:12`).
- `v_app`: Apparent wind [m/s] the polars' Reynolds number is taken at
    (default `25.0`).
- `chord_ref`: Reference (maximum panel) chord [m], which Reynolds is defined
    against (default `1.0`).
- `table_format`: Per-node table format, `:csv` (readable) or `:arrow` (about ten
    times faster to load; default `:arrow`). Written in YAML as a plain string,
    and carried as the `Symbol` the generator takes.
"""
@with_kw mutable struct AirfoilSettings
    solver::String = "neuralfoil"
    model_size::String = "large"
    n_crit::Float64 = 9.0
    xtr_upper::Float64 = 0.05
    xtr_lower::Float64 = 0.05
    alpha_range::Vector{Float64} = [-180.0, 1.0, 180.0]
    delta_range::Union{Nothing, Vector{Float64}} = nothing
    live_offsets::Vector{Float64} = collect(-12.0:3.0:12.0)
    v_app::Float64 = 25.0
    chord_ref::Float64 = 1.0
    table_format::Symbol = :arrow
end

"""
    WingSettings

Settings for a single wing, used within [`VSMSettings`](@ref).

# Fields
- `name`: Wing identifier (default `"main_wing"`)
- `geometry_file`: Path to wing geometry YAML file
- `obj_file`: Path to `.obj` geometry file
- `dat_file`: Path to `.dat` airfoil file
- `n_panels`: Number of panels (default `40`)
- `spanwise_panel_distribution`: Panel distribution type
    (default [`LINEAR`](@ref PanelDistribution))
- `spanwise_direction`: Spanwise direction vector
    (default `[0, 1, 0]`)
- `remove_nan`: Whether to remove NaN values from polar data
    (default `true`)
- `use_prior_polar`: Reuse prior refined/panel polar mapping on
    reinit/refine updates (default `false`)
- `billowing_percentage`: TE billow as percentage of arc length
    (default `0.0`; only used with `BILLOWING` distribution).
- `crease_frac`: Chordwise flap-hinge fraction (0–1) the polars are deflected
    about and the δ-deflected plate/skin is drawn with (default `0.75`).
- `mesh`: [`MeshSettings`](@ref) — how `obj_file` is sliced into sections.
- `airfoil`: [`AirfoilSettings`](@ref) — the 2D backend and the polars it writes.
"""
@with_kw mutable struct WingSettings
    name::String = "main_wing"
    geometry_file::String = ""          # path to wing geometry YAML file
    obj_file::String = ""               # path to .obj geometry file
    dat_file::String = ""               # path to .dat airfoil file
    n_panels::Int64 = 40
    spanwise_panel_distribution::PanelDistribution = LINEAR
    spanwise_direction::MVec3 = [0.0, 1.0, 0.0]
    remove_nan::Bool = true
    use_prior_polar::Bool = false
    billowing_percentage::Float64 = 0.0 # TE billow as % of arc length
    crease_frac::Float64 = 0.75
    mesh::MeshSettings = MeshSettings()
    airfoil::AirfoilSettings = AirfoilSettings()
end

"""
    SolverSettings

Solver configuration, used within [`VSMSettings`](@ref).

# Fields
- `n_panels`: Total number of panels (default `40`)
- `aerodynamic_model_type`: [`VSM`](@ref Model) or
    [`LLT`](@ref Model) (default `VSM`)
- `solver_type`: `"LOOP"` or `"NONLIN"` (default `"LOOP"`)
- `density`: Air density (kg/m^3) (default `1.225`)
- `max_iterations`: Maximum solver iterations (default `1500`)
- `rtol`: Relative tolerance (default `1e-5`)
- `tol_reference_error`: Reference error tolerance
    (default `0.001`)
- `relaxation_factor`: Convergence relaxation factor
    (default `0.03`)
- `artificial_damping`: Enable artificial damping
    (default `false`)
- `k2`, `k4`: Artificial damping parameters
- `is_with_artificial_viscosity`: Enable Li/Gaunaa post-stall
    artificial viscosity (default `false`)
- `artificial_viscosity_factor`: Viscosity scaling coefficient k
    (default `0.035`)
- `type_initial_gamma_distribution`:
    [`ELLIPTIC`](@ref InitialGammaDistribution) or `ZEROS`
    (default `ZEROS`)
- `use_gamma_prev`: Reuse provided previous gamma as initial guess when
    available (default `true`)
- `core_radius_fraction`: Bound vortex core cut-off, as a fraction of the filament
    length, following Damiani et al. (2019) (default `0.05`)
- `mu`: Dynamic viscosity (N*s/m^2) (default `1.81e-5`)
- `calc_only_f_and_gamma`: Only output forces and circulation
    (default `false`)
- `correct_aoa`: Perform angle of attack correction
    (default `false`)
- `flow_curvature`: Add the thin-airfoil pitch-rate moment increment to each
    section (default `false`)
"""
@with_kw mutable struct SolverSettings
    n_panels::Int64 = 40
    aerodynamic_model_type::Model = VSM
    solver_type::String = "LOOP"    # type of solver
    density::Float64 = 1.225                # air density  [kg/m³] 
    max_iterations::Int64 = 1500
    rtol::Float64 = 1e-5                    # relative error   [-]
    tol_reference_error::Float64 = 0.001
    relaxation_factor::Float64 = 0.03       # relaxation factor for convergence
    artificial_damping::Bool = false        # whether to apply artificial damping
    k2::Float64 = 0.1                       # artificial damping parameter
    k4::Float64 = 0.0                       # artificial damping parameter
    is_with_artificial_viscosity::Bool = false  # Li/Gaunaa post-stall artificial viscosity
    artificial_viscosity_factor::Float64 = 0.035 # viscosity scaling coefficient k
    type_initial_gamma_distribution::InitialGammaDistribution = ZEROS # see: [`InitialGammaDistribution`](@ref)
    use_gamma_prev::Bool = true             # if false, always reinitialize gamma from type_initial_gamma_distribution
    core_radius_fraction::Float64 = 0.05
    mu::Float64 = 1.81e-5                   # dynamic viscosity [N·s/m²]
    calc_only_f_and_gamma::Bool=false       # whether to only output f and gamma
    correct_aoa::Bool=false                 # perform aoa correction
    flow_curvature::Bool=false              # thin-airfoil pitch-rate moment increment
end

"""
    VSMSettings

Top-level settings container for a VortexStepMethod simulation.
Can be constructed from keyword arguments or loaded from a YAML
file with `VSMSettings(filename)`.

# Fields
- `condition`: `ConditionSettings` (wind speed, alpha, beta,
    yaw rate)
- `wings`: Vector of [`WingSettings`](@ref)
- `solver_settings`: [`SolverSettings`](@ref)

# Example
```julia
settings = VSMSettings("vsm_settings.yaml")
wing = Wing(settings)
```
"""
@Base.kwdef mutable struct VSMSettings
    condition::ConditionSettings = ConditionSettings()
    wings::Vector{WingSettings} = []
    solver_settings::SolverSettings = SolverSettings()
end

function VSMSettings(filename; data_prefix=true)
    # Uwe's suggested 3-line approach using StructMapping.jl (adapted)
    if data_prefix
        data = YAML.load_file(joinpath("data", filename))
    else
        data = YAML.load_file(filename)
    end
    
    # Use StructMapping for basic structure conversion
    # But handle special fields manually due to enum conversion needs
    vsm_settings = VSMSettings()
    
    # Convert condition settings using StructMapping (if present)
    if haskey(data, "condition")
        vsm_settings.condition = convertdict(ConditionSettings, data["condition"])
    end
    
    # Convert wing settings manually due to enum conversions
    n_panels = 0

    if haskey(data, "wings")
        for wing_data in data["wings"]
            wing = WingSettings()
            wing.name = wing_data["name"]
            if haskey(wing_data, "geometry_file")
                wing.geometry_file = wing_data["geometry_file"]
            end
            if haskey(wing_data, "obj_file")
                wing.obj_file = wing_data["obj_file"]
            end
            if haskey(wing_data, "dat_file")
                wing.dat_file = wing_data["dat_file"]
            end
            wing.n_panels = wing_data["n_panels"]
            # Handle deprecated n_groups parameter
            if haskey(wing_data, "n_groups")
                @warn "n_groups in settings file is deprecated and ignored. Use n_unrefined_sections or let it be inferred automatically." maxlog=1
            end
            wing.spanwise_panel_distribution = parse_enum(PanelDistribution, wing_data["spanwise_panel_distribution"])
            wing.spanwise_direction = MVec3(wing_data["spanwise_direction"])
            if haskey(wing_data, "grouping_method")
                @warn "grouping_method in settings file is deprecated and ignored." maxlog=1
            end
            wing.remove_nan = wing_data["remove_nan"]
            wing.use_prior_polar = get(wing_data, "use_prior_polar", false)

            if haskey(wing_data, "billowing_percentage")
                wing.billowing_percentage =
                    Float64(wing_data["billowing_percentage"])
            end
            if haskey(wing_data, "crease_frac")
                wing.crease_frac = Float64(wing_data["crease_frac"])
            end
            haskey(wing_data, "mesh") &&
                (wing.mesh = convertdict(MeshSettings, wing_data["mesh"]))
            haskey(wing_data, "airfoil") &&
                (wing.airfoil = airfoil_settings(wing_data["airfoil"]))

            push!(vsm_settings.wings, wing)
            n_panels += wing.n_panels
        end
    end
    
    # Convert solver settings using StructMapping base, then override special fields
    if haskey(data, "solver_settings")
        solver_data = data["solver_settings"]
        
        # Create a copy of solver_data with string fields for enums removed
        solver_data_clean = copy(solver_data)
        # Backward-compatible alias from user-facing YAML key.
        if haskey(solver_data_clean, "use_gamme_prev") && !haskey(solver_data_clean, "use_gamma_prev")
            solver_data_clean["use_gamma_prev"] = solver_data_clean["use_gamme_prev"]
        end
        haskey(solver_data_clean, "use_gamme_prev") && delete!(solver_data_clean, "use_gamme_prev")
        delete!(solver_data_clean, "aerodynamic_model_type")
        delete!(solver_data_clean, "type_initial_gamma_distribution")
        
        # Use StructMapping for the basic fields
        vsm_settings.solver_settings = convertdict(SolverSettings, solver_data_clean)
        
        # Handle enum conversions manually
        vsm_settings.solver_settings.aerodynamic_model_type = parse_enum(Model, solver_data["aerodynamic_model_type"])
        vsm_settings.solver_settings.type_initial_gamma_distribution = parse_enum(InitialGammaDistribution, solver_data["type_initial_gamma_distribution"])

        # Override with calculated totals
        vsm_settings.solver_settings.n_panels = n_panels
    end
    
    return vsm_settings
end

"""
    airfoil_settings(data) -> AirfoilSettings

Read a wing's `airfoil:` block, checking that `solver` and `table_format` name
backends that exist. Every key is optional and falls back to the default;
`delta_range: null` (or an absent key) means no flap sweep.
"""
function airfoil_settings(data)
    airfoil = convertdict(AirfoilSettings, data)
    airfoil.solver in ("neuralfoil", "xfoil") || throw(ArgumentError(
        "Invalid airfoil solver \"$(airfoil.solver)\". " *
        "Valid values: neuralfoil, xfoil"))
    airfoil.table_format in (:csv, :arrow) || throw(ArgumentError(
        "Invalid table_format \"$(airfoil.table_format)\". " *
        "Valid values: csv, arrow"))
    return airfoil
end

function Base.show(io::IO, vsm_settings::VSMSettings)
    println(io, "VSMSettings:")
    print(io, replace(repr(vsm_settings.condition), "\n" => "\n    "))
    for (i, wing) in pairs(vsm_settings.wings)
        if i==1
            print(io, "    ")
        end
        print(io, replace(repr(wing), "\n" => "\n    "))
    end
    print(io, replace(repr(vsm_settings.solver_settings), "\n" => "\n    "))
end

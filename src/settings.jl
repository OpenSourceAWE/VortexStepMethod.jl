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
- `crease_frac`: Chordwise flap-hinge fraction (0–1) for drawing the
    δ-deflected plate/skin (default `0.75`).
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

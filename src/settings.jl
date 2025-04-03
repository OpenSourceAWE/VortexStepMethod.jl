@with_kw mutable struct WingSettings
    name::String = "main_wing"
    n_panels::Int64 = 40
    n_groups::Int64 = 40 
    spanwise_panel_distribution::PanelDistribution = LINEAR
    spanwise_direction::MVec3 = [0.0, 1.0, 0.0]
    remove_nan = true
end

@with_kw mutable struct SolverSettings
    n_panels::Int64 = 40
    n_groups::Int64 = 40 
    aerodynamic_model_type::Model = VSM
    density::Float64 = 1.225                # air density  [kg/m³] 
    max_iterations::Int64 = 1500
    rtol::Float64 = 1e-5                    # relative error   [-]
    tol_reference_error::Float64 = 0.001
    relaxation_factor::Float64 = 0.03       # relaxation factor for convergence
    artificial_damping::Bool = false        # whether to apply artificial damping
    k2::Float64 = 0.1                       # artificial damping parameter
    k4::Float64 = 0.0                       # artificial damping parameter
    type_initial_gamma_distribution::InitialGammaDistribution = ELLIPTIC # see: [InitialGammaDistribution](@ref)
    core_radius_fraction::Float64 = 1e-20
    mu::Float64 = 1.81e-5                   # dynamic viscosity [N·s/m²]
    calc_only_f_and_gamma::Bool=false       # whether to only output f and gamma
end

@Base.kwdef mutable struct VSMSettings
    wings::Vector{WingSettings} = []
    solver_settings::SolverSettings = SolverSettings()
end

function vs(filename)
    res = VSMSettings()
    data = YAML.load_file(joinpath("data", filename))
    # update solver settings
    solver = data["solver_settings"]
    res.solver_settings.n_panels = solver["n_panels"]
    res.solver_settings.n_groups = solver["n_groups"]
    res.solver_settings.aerodynamic_model_type = eval(Symbol(solver["aerodynamic_model_type"]))
    res.solver_settings.density = solver["density"]
    res.solver_settings.max_iterations = solver["max_iterations"]
    res.solver_settings.rtol = solver["rtol"]
    res.solver_settings.tol_reference_error = solver["tol_reference_error"]
    res.solver_settings.relaxation_factor = solver["relaxation_factor"]
    res.solver_settings.artificial_damping = solver["artificial_damping"]
    res.solver_settings.k2 = solver["k2"]
    res.solver_settings.k4 = solver["k4"]
    res.solver_settings.type_initial_gamma_distribution = eval(Symbol(solver["type_initial_gamma_distribution"]))
    res.solver_settings.core_radius_fraction = solver["core_radius_fraction"]
    res.solver_settings.mu = solver["mu"]
    res.solver_settings.calc_only_f_and_gamma = solver["calc_only_f_and_gamma"]

    # add and update wing settings
    for (i, wing) in pairs(data["wings"])
        push!(res.wings, WingSettings())
        res.wings[i].name = wing["name"]
        res.wings[i].n_panels = wing["n_panels"]
        res.wings[i].n_groups = wing["n_groups"]
        res.wings[i].spanwise_panel_distribution = eval(Symbol(wing["spanwise_panel_distribution"]))
        res.wings[i].spanwise_direction = MVec3(wing["spanwise_direction"])
        res.wings[i].remove_nan = wing["remove_nan"]
    end
    res
end

function Base.show(io::IO, vs::VSMSettings)
    println(io, "VSMSettings:")
    for (i, wing) in pairs(vs.wings)
        if i==1
            print(io, "    ")
        end
        print(io, replace(repr(wing), "\n" => "\n    "))
    end
    print(io, replace(repr(vs.solver_settings), "\n" => "\n    "))
end

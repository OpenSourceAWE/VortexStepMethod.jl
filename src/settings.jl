filename = "vsm_settings.yaml"
data = YAML.load_file(joinpath("data", filename))

@with_kw mutable struct WingSettings
    name::String = "main_wing"
    n_panels::Int64 = 40
    n_groups::Int64 = 40 
    spanwise_panel_distribution::PanelDistribution = LINEAR
    spanwise_direction::MVec3 = [0.0, 1.0, 0.0]
    remove_nan = true
end

@with_kw mutable struct SolverSettings
    aerodynamic_model_type::Model = VSM
    max_iterations::Int64 = 1000
end

@Base.kwdef mutable struct VSMSettings
    wings::Vector{WingSettings} = [WingSettings(), WingSettings()]
    solver_settings::SolverSettings = SolverSettings()
end

const VSM_SETTINGS = VSMSettings()

function vs(filename=filename)
    res = VSM_SETTINGS
    data = YAML.load_file(joinpath("data", filename))
    res.solver_settings.max_iterations = data["solver_settings"]["max_iterations"]
    res.solver_settings.aerodynamic_model_type = eval(Symbol(data["solver_settings"]["aerodynamic_model_type"]))
    for (i, wing) in pairs(data["wings"])
        println(i)
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
    for wing in vs.wings
        print(io, "    ", replace(repr(wing), "\n" => "\n    "))
    end
    print(io, replace(repr(vs.solver_settings), "\n" => "\n    "))
end

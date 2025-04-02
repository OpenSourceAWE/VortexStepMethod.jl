filename = "vsm_settings.yaml"
data = YAML.load_file(joinpath(dirname(dirname(pathof(VortexStepMethod))), "data", filename))

@Base.kwdef mutable struct WingSettings
    name::String = "main_wing"
    n_panels::Int64 = 40
    n_groups::Int64 = 40 
    spanwise_panel_distribution::String = "LINEAR"
    spanwise_direction::MVec3 = [0.0, 1.0, 0.0]
    remove_nan = true
end

@Base.kwdef mutable struct SolverSettings
    aerodynamic_model_type::String = "VSM"
    max_iterations::Int64 = 1000
end

@Base.kwdef mutable struct VSMSettings
    wings::Vector{WingSettings} = [WingSettings()]
    solver_settings::SolverSettings = SolverSettings()
end

function Base.show(io::IO, vs::VSMSettings)
    println(io, "VSMSettings:")
    for wing in vs.wings
        print(io, "    ", replace(repr(wing), "\n" => "\n    "))
    end
    print(io, replace(repr(vs.solver_settings), "\n" => "\n    "))
end

# read the vsm_settings.yaml file

import YAML
using Parameters

filename = "vsm_settings.yaml"
data = YAML.load_file(joinpath("data", filename))

@with_kw mutable struct WingSettings
    name::String = "main_wing"
    n_groups::Int64 = 40 
    n_panels::Int64 = 40
    spanwise_panel_distribution::String = "LINEAR"
end

@with_kw mutable struct SolverSettings
    aerodynamic_model_type::String = "VSM"
    max_iterations::Int64 = 1000
end

@with_kw mutable struct VSMSettings
    wings::Vector{WingSettings} = [WingSettings()]
    solver_settings::SolverSettings = SolverSettings()
end


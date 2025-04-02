# read the vsm_settings.yaml file

import YAML
using Parameters, StaticArrays

const MVec3    = MVector{3, Float64}

filename = "vsm_settings.yaml"
data = YAML.load_file(joinpath("data", filename))

@with_kw mutable struct WingSettings
    name::String = "main_wing"
    n_panels::Int64 = 40
    n_groups::Int64 = 40 
    spanwise_panel_distribution::String = "LINEAR"
    spanwise_direction::MVec3 = [0.0, 1.0, 0.0]
    remove_nan = true
end

@with_kw mutable struct SolverSettings
    aerodynamic_model_type::String = "VSM"
    max_iterations::Int64 = 1000
end

@with_kw mutable struct VSMSettings
    wings::Vector{WingSettings} = [WingSettings()]
    solver_settings::SolverSettings = SolverSettings()
end


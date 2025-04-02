using Parameters

@with_kw mutable struct WingSettings
    n_panels::Int
    spanwise_panel_distribution::String = "uniform"
end

@with_kw mutable struct VSMSettings
    aerodynamic_model_type::String = "VSM"
    max_iterations::Int = 1000
end

@with_kw mutable struct Configuration
    wings::Vector{Dict{String, WingSettings}}
    solver_settings::VSMSettings
end
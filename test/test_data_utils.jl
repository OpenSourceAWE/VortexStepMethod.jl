# Test data utilities for managing shared test configurations
# Organized by source module being tested

using YAML
using Random: randstring
using Logging
using VortexStepMethod.ObjAdapter: obj_to_matrix_yaml

"""
    ram_air_matrix_wing(; n_panels, n_sections=4, wind_vel=15.0,
                        alpha_range=deg2rad.(-5:5:15), delta_range=deg2rad.(-3:3:3))

Build a ram-air-kite `Wing` via convert-then-load: the obj mesh is converted to a
standard geometry YAML backed by shared NeuralFoil `(alpha, delta)` `POLAR_MATRICES`
(`ObjAdapter.obj_to_matrix_yaml` auto-slices a representative section from the mesh),
then loaded with `Wing(yaml)`. Replaces the removed live `ObjWing` route in tests.

The expensive conversion (mesh slicing + NeuralFoil matrix generation) is written once
into `test/generated/` (gitignored), keyed by `(n_sections, wind_vel, alpha_range,
delta_range)`; if that geometry already exists it is reused, so it is generated only
once per configuration across all runs. Each call still reloads a fresh `Wing`.
"""
function ram_air_matrix_wing(; n_panels, n_sections=4, wind_vel=15.0,
        alpha_range=deg2rad.(-5:5:15), delta_range=deg2rad.(-3:3:3))
    data_dir = joinpath(dirname(@__DIR__), "data", "ram_air_kite")
    obj = joinpath(data_dir, "ram_air_kite_body.obj")
    key = (n_sections, wind_vel, collect(alpha_range), collect(delta_range))
    gen_dir = joinpath(@__DIR__, "generated", "ram_matrix_$(string(hash(key); base=16))")
    yaml = joinpath(gen_dir, "geometry.yaml")
    isfile(yaml) || obj_to_matrix_yaml(obj, gen_dir; n_sections, wind_vel,
        alpha_range, delta_range, verbose=false)
    return Wing(yaml; n_panels)
end

"""
    suppress_warnings(f)

Temporarily suppress all warnings while executing function `f`.
Useful for tests that expect and want to ignore certain warnings.

# Example
```julia
result = suppress_warnings(() -> some_function_that_warns())
```
"""
function suppress_warnings(f)
    old_logger = global_logger()
    try
        global_logger(NullLogger())
        return f()
    finally
        global_logger(old_logger)
    end
end

"""
    test_data_path(module_name, relative_path...)

Get the absolute path to test data files relative to the test/module_name directory.

# Arguments
- `module_name`: Name of the source module being tested (e.g., "yaml_geometry", "body_aerodynamics", "solver")
- `relative_path...`: Path components relative to the module directory

# Examples
```julia
wing_file = test_data_path("yaml_geometry", "wings", "simple_wing.yaml")
polar_file = test_data_path("yaml_geometry", "polars", "standard_airfoil.csv")
body_aero_wing = test_data_path("body_aerodynamics", "wings", "test_wing.yaml")
```
"""
test_data_path(module_name, relative_path...) =
    joinpath(@__DIR__, string(module_name), map(string, relative_path)...)

"""
    create_temp_wing_settings(module_name, wing_file; kwargs...)

Create a temporary VSMSettings file that references the given wing geometry file.
Useful for tests that need to modify settings while using standard wing geometries.

# Arguments
- `module_name`: Name of the source module being tested
- `wing_file`: Wing geometry file name (relative to test/data/module_name/wings/)
- `kwargs...`: Additional settings to override defaults

# Returns
- Path to temporary settings file (caller should clean up)

# Example
```julia
settings_file = create_temp_wing_settings("body_aerodynamics", "test_wing.yaml"; alpha=15.0, wind_speed=25.0)
# Use settings_file...
rm(settings_file)
```
"""
function create_temp_wing_settings(module_name, wing_file;
    name="test_wing",
    n_panels=4,
    spanwise_panel_distribution="COSINE",
    spanwise_direction=[0.0, 1.0, 0.0],
    remove_nan=true,
    use_prior_polar=false,
    aerodynamic_model_type="VSM",
    density=1.225,
    type_initial_gamma_distribution="ZEROS",
    alpha=10.0,
    beta=5.0,
    wind_speed=15.0,
)
    wing_file_path = isabspath(wing_file) ? wing_file : test_data_path(module_name, wing_file)
    wing_file_path = replace(normpath(wing_file_path), '\\' => '/')

    doc = Dict(
        "wings" => [Dict(
            "name" => name,
            "geometry_file" => wing_file_path,
            "n_panels" => n_panels,
            "spanwise_panel_distribution" => spanwise_panel_distribution,
            "spanwise_direction" => spanwise_direction,
            "remove_nan" => remove_nan,
            "use_prior_polar" => use_prior_polar,
        )],
        "solver_settings" => Dict(
            "aerodynamic_model_type" => aerodynamic_model_type,
            "density" => density,
            "type_initial_gamma_distribution" => type_initial_gamma_distribution,
        ),
        "condition" => Dict(
            "alpha" => alpha,
            "beta" => beta,
            "wind_speed" => wind_speed,
        ),
    )

    temp_file = joinpath(tempdir(), "temp_wing_settings_$(randstring(8)).yaml")
    open(temp_file, "w") do io
        YAML.write(io, doc)
    end
    return temp_file
end

"""
    get_standard_wing_file(module_name)

Get the path to the standard simple wing file for a given module.

# Example
```julia
wing_file = get_standard_wing_file("body_aerodynamics")  # Returns path to test_wing.yaml
```
"""
function get_standard_wing_file(module_name)
    if module_name == "yaml_geometry"
        return test_data_path("yaml_geometry", "wings", "simple_wing.yaml")
    elseif module_name == "body_aerodynamics"
        return test_data_path("body_aerodynamics", "wings", "test_wing.yaml")
    elseif module_name == "solver"
        return test_data_path("solver", "wings", "solver_test_wing.yaml")
    else
        error("Unknown module: $module_name")
    end
end

"""
    get_complete_settings_file(module_name)

Get the path to a complete settings file for a given module.

# Example
```julia
settings_file = get_complete_settings_file("body_aerodynamics")
```
"""
function get_complete_settings_file(module_name)
    if module_name == "body_aerodynamics"
        return test_data_path("body_aerodynamics", "complete_settings.yaml")
    elseif module_name == "solver"
        return test_data_path("solver", "solver_settings.yaml")
    else
        # Fall back to creating a temporary file
        return create_temp_wing_settings(module_name, basename(get_standard_wing_file(module_name)))
    end
end

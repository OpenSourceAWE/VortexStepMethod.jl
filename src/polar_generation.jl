"""
Polar generation pipeline using NeuralFoil.

This module provides functions to generate aerodynamic polars from OBJ meshes
using the NeuralFoil neural network.
"""

using DelimitedFiles
using Printf: @sprintf

"""
    generate_neuralfoil_polars(obj_path::String, output_dir::String;
                               n_slices::Int, Re::Real,
                               alpha_range=-180:1:180,
                               model_size::String="xlarge",
                               weights_dir=nothing,
                               n_crit=9.0,
                               verbose=true)

Generate aerodynamic polars for wing sections using NeuralFoil.

# Arguments
- `obj_path`: Path to OBJ file containing 3D wing mesh
- `output_dir`: Directory to write polar CSV files

# Keyword Arguments
- `n_slices`: Number of spanwise slices (typically n_panels+1)
- `Re`: Reynolds number for polar evaluation
- `alpha_range`: Angle of attack range in degrees (default -180:1:180)
- `model_size`: NeuralFoil model size (default "xlarge")
- `weights_dir`: Directory containing NeuralFoil weights
- `n_crit`: Critical amplification factor (default 9.0)
- `verbose`: Print progress information (default true)

# Output
Writes CSV files to `output_dir/{n}.csv` in format:
```
alpha,Cd,Cs,Cl,Cm
-180.0,0.0,0.0,0.0,0.0
...
```
"""
function generate_neuralfoil_polars(obj_path::String, output_dir::String;
                                    n_slices::Int, Re::Real,
                                    alpha_range=-180:1:180,
                                    model_size::String="xlarge",
                                    weights_dir=nothing,
                                    n_crit=9.0,
                                    verbose=true)
    # Validate inputs
    if !isfile(obj_path)
        error("OBJ file not found: $obj_path")
    end

    # Create output directory
    mkpath(output_dir)

    verbose && @info "Slicing OBJ mesh at $n_slices positions..."

    # Slice the mesh
    airfoils, y_positions = slice_obj_wing(obj_path, n_slices)

    verbose && @info "Processing $(length(airfoils)) airfoil sections..."

    # Convert alpha range to vector
    alphas = collect(Float64, alpha_range)

    # Load NeuralFoil model once
    model = load_neuralfoil_model(model_size; weights_dir)

    # Process each slice
    for (i, (x, y)) in enumerate(airfoils)
        if isempty(x) || length(x) < 5
            @warn "Slice $i: Insufficient points, skipping"
            continue
        end

        verbose && print("  Slice $i (y=$(round(y_positions[i], digits=3)))... ")

        try
            # Fit Kulfan parameters
            params = fit_kulfan_parameters(x, y)

            # Evaluate NeuralFoil
            result = neuralfoil_aero(params, alphas, Float64(Re);
                                     model_size, weights_dir, n_crit)

            # Write polar to CSV
            write_polar_csv(joinpath(output_dir, "$i.csv"), result)

            verbose && println("done (CL_max=$(round(maximum(result.CL), digits=2)))")

        catch e
            @warn "Slice $i failed: $e"
            continue
        end
    end

    verbose && @info "Polars written to $output_dir"
    return nothing
end

"""
    write_polar_csv(filepath::String, result::NeuralFoilResult)

Write NeuralFoil results to CSV file in standard polar format.
"""
function write_polar_csv(filepath::String, result::NeuralFoilResult)
    n = length(result.alpha)

    # Build data matrix
    # Format: alpha, Cd, Cs, Cl, Cm
    data = zeros(n, 5)
    data[:, 1] = result.alpha
    data[:, 2] = result.CD
    data[:, 3] = zeros(n)  # Cs = 0 for 2D airfoils
    data[:, 4] = result.CL
    data[:, 5] = result.CM

    # Write with header
    open(filepath, "w") do io
        println(io, "alpha,Cd,Cs,Cl,Cm")
        for i in 1:n
            println(io, join([@sprintf("%.16g", v) for v in data[i, :]], ","))
        end
    end
end

"""
    generate_polar_from_coordinates(x::Vector, y::Vector, output_path::String;
                                    Re::Real, alpha_range=-180:1:180, kwargs...)

Generate a single polar from airfoil coordinates.

# Arguments
- `x, y`: Airfoil coordinates (will be normalized)
- `output_path`: Path to write CSV file
- `Re`: Reynolds number

# Keyword Arguments
- `alpha_range`: Angle of attack range in degrees (default -180:1:180)
- Additional kwargs passed to neuralfoil_aero
"""
function generate_polar_from_coordinates(x::Vector, y::Vector, output_path::String;
                                         Re::Real, alpha_range=-180:1:180, kwargs...)
    alphas = collect(Float64, alpha_range)

    params = fit_kulfan_parameters(x, y)
    result = neuralfoil_aero(params, alphas, Float64(Re); kwargs...)

    write_polar_csv(output_path, result)
    return result
end

"""
    generate_polar_from_dat(dat_path::String, output_path::String;
                            Re::Real, kwargs...)

Generate polar from a .dat airfoil coordinate file.
"""
function generate_polar_from_dat(dat_path::String, output_path::String;
                                 Re::Real, kwargs...)
    if !isfile(dat_path)
        error("DAT file not found: $dat_path")
    end

    x, y = read_dat_coordinates(dat_path)

    if isempty(x)
        error("No valid coordinates found in $dat_path")
    end

    return generate_polar_from_coordinates(x, y, output_path; Re, kwargs...)
end

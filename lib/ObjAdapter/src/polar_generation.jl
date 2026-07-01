"""
Obj-based polar generation: slice a mesh and drive AirfoilAero per section. The
airfoil-general `.dat`/coordinate → CSV entry points live in AirfoilAero.
"""

"""
    generate_neuralfoil_polars(obj_path, output_dir; n_slices, Re,
                               alpha_range=-180:1:180, model_size="xlarge",
                               weights_dir=nothing, n_crit=9.0, verbose=true)

Slice `obj_path` into `n_slices` airfoils and write a `POLAR_VECTORS` CSV per slice
to `output_dir/{i}.csv` via [`generate_polar_from_coordinates`](@ref).
"""
function generate_neuralfoil_polars(obj_path::String, output_dir::String;
                                    n_slices::Int, Re::Real,
                                    alpha_range=-180:1:180,
                                    model_size::String="xlarge",
                                    weights_dir=nothing,
                                    n_crit=9.0,
                                    verbose=true)
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    mkpath(output_dir)
    verbose && @info "Slicing OBJ mesh at $n_slices positions..."
    airfoils, y_positions = slice_obj_wing(obj_path, n_slices)
    verbose && @info "Processing $(length(airfoils)) airfoil sections..."
    load_neuralfoil_model(model_size; weights_dir)
    for (i, (x, y)) in enumerate(airfoils)
        if isempty(x) || length(x) < 5
            @warn "Slice $i: Insufficient points, skipping"
            continue
        end
        verbose && print("  Slice $i (y=$(round(y_positions[i], digits=3)))... ")
        try
            result = generate_polar_from_coordinates(x, y,
                joinpath(output_dir, "$i.csv"); Re=Float64(Re), alpha_range,
                model_size, weights_dir, n_crit)
            verbose && println("done (CL_max=$(round(maximum(result.CL), digits=2)))")
        catch e
            @warn "Slice $i failed: $e"
            continue
        end
    end
    verbose && @info "Polars written to $output_dir"
    return nothing
end

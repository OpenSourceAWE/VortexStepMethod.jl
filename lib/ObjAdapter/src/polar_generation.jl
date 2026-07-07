"""
Obj-based polar generation: slice a mesh and drive AirfoilAero per section. The
airfoil-general `.dat`/coordinate → CSV entry points live in AirfoilAero.
"""

"""
    generate_neuralfoil_polars(obj_path, output_dir; n_slices, Re,
                               alpha_range=-180:1:180, model_size="xlarge",
                               weights_dir=nothing, n_crit=9.0, rotation=I,
                               fit_method=EnvelopeFit(min_distance=0.01), verbose=true)

Slice `obj_path` into `n_slices` airfoils ([`perpendicular_sections`](@ref)) and write
a `POLAR_VECTORS` CSV per slice to `output_dir/{i}.csv` via
[`generate_polar_from_coordinates`](@ref). The default `fit_method` wraps each slice
with a thin [`EnvelopeFit`](@ref) envelope, so an open single-membrane canopy slice
(no closed airfoil) still yields a valid thin cambered airfoil. Pass `rotation` (or
`:auto`) to reorient the mesh first.
"""
function generate_neuralfoil_polars(obj_path::String, output_dir::String;
                                    n_slices::Int, Re::Real,
                                    alpha_range=-180:1:180,
                                    model_size::String="xlarge",
                                    weights_dir=nothing,
                                    n_crit=9.0,
                                    rotation=I,
                                    n_bins=60,
                                    fit_method=EnvelopeFit(min_distance=0.01),
                                    verbose=true)
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    mkpath(output_dir)
    verbose && @info "Slicing OBJ mesh at $n_slices positions..."
    vertices, faces = read_faces(obj_path)
    sections = perpendicular_sections(vertices, faces, n_slices; n_bins, rotation)
    verbose && @info "Processing $(length(sections)) airfoil sections..."
    load_neuralfoil_model(model_size; weights_dir)
    for (i, sec) in enumerate(sections)
        x, y = sec.x_airfoil, sec.y_airfoil
        if isempty(x) || length(x) < 5
            @warn "Slice $i: Insufficient points, skipping"
            continue
        end
        verbose && print("  Slice $i (y=$(round(sec.LE_point[2], digits=3)))... ")
        try
            result = generate_polar_from_coordinates(x, y,
                joinpath(output_dir, "$i.csv"); Re=Float64(Re), alpha_range,
                fit_method, model_size, weights_dir, n_crit)
            verbose && println("done (CL_max=$(round(maximum(result.CL), digits=2)))")
        catch e
            @warn "Slice $i failed: $e"
            continue
        end
    end
    verbose && @info "Polars written to $output_dir"
    return nothing
end

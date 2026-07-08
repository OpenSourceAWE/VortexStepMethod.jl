"""
Obj-based polar generation: slice a mesh and drive AirfoilAero per section. The
airfoil-general `.dat`/coordinate → CSV entry points live in AirfoilAero.
"""

"""
    generate_section_polars(obj_path, output_dir; n_slices, Re,
                            alpha_range=-180:1:180, delta_range=nothing, rotation=I,
                            n_bins=60, fit_method=EnvelopeFit(min_distance=0.01),
                            solver=NeuralFoilSolver(), verbose=true)

Slice `obj_path` into `n_slices` airfoils ([`perpendicular_sections`](@ref)) and write
a polar CSV per slice to `output_dir/{i}.csv` via
[`generate_polar_from_coordinates`](@ref). The default `fit_method` wraps each slice
with a thin [`EnvelopeFit`](@ref) envelope, so an open single-membrane canopy slice
(no closed airfoil) still yields a valid thin cambered airfoil. `solver` picks the
backend — [`NeuralFoilSolver`](@ref) (full ±180° range) or [`XFoilSolver`](@ref) (use a
narrow `alpha_range` where XFoil converges). Pass a `delta_range` of trailing-edge
deflections to write long-format `(alpha, delta)` `POLAR_MATRICES` CSVs instead of
`POLAR_VECTORS`. Pass `rotation` (or `:auto`) to reorient the mesh first.
"""
function generate_section_polars(obj_path::String, output_dir::String;
                                 n_slices::Int, Re::Real,
                                 alpha_range=-180:1:180,
                                 delta_range=nothing,
                                 rotation=I,
                                 n_bins=60,
                                 wingtip_distance=0.0,
                                 fit_method=EnvelopeFit(min_distance=0.01),
                                 solver::AbstractAirfoilSolver=NeuralFoilSolver(),
                                 verbose=true)
    isfile(obj_path) || error("OBJ file not found: $obj_path")
    mkpath(output_dir)
    verbose && @info "Slicing OBJ mesh at $n_slices positions..."
    vertices, faces = read_faces(obj_path)
    sections = perpendicular_sections(vertices, faces, n_slices;
                                      n_bins, rotation, wingtip_distance)
    verbose && @info "Processing $(length(sections)) airfoil sections..."
    for (i, sec) in enumerate(sections)
        x, y = sec.x_airfoil, sec.y_airfoil
        if isempty(x) || length(x) < 5
            @warn "Slice $i: Insufficient points, skipping"
            continue
        end
        verbose && print("  Slice $i (y=$(round(sec.LE_point[2], digits=3)))... ")
        try
            res = generate_polar_from_coordinates(x, y,
                joinpath(output_dir, "$i.csv"); Re=Float64(Re), alpha_range,
                fit_method, solver, delta_range)
            clvals = res isa AbstractVector ? (s.cl for s in res) : vec(res[1])
            finite_cl = [v for v in clvals if !isnan(v)]
            cl_max = isempty(finite_cl) ? NaN : maximum(finite_cl)
            verbose && println("done (CL_max=$(round(cl_max, digits=2)))")
        catch e
            @warn "Slice $i failed: $e"
            continue
        end
    end
    verbose && @info "Polars written to $output_dir"
    return nothing
end

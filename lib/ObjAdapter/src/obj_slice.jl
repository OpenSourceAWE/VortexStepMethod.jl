"""
OBJ mesh slicing for extracting 2D airfoil cross-sections.

This module provides functions to slice 3D wing meshes at specific spanwise
positions to extract 2D airfoil profiles.
"""

"""
    slice_mesh_at_y(vertices::Vector, faces::Vector, y_pos::Real; tol=1e-6)

Slice a triangulated mesh at a constant y-position.

# Arguments
- `vertices`: Vector of 3D vertex coordinates [x, y, z]
- `faces`: Vector of face indices (triangles)
- `y_pos`: Y-coordinate of the slice plane

# Returns
- Vector of line segments [(p1, p2), ...] where p1, p2 are [x, z] coordinates
"""
function slice_mesh_at_y(vertices::Vector, faces::Vector, y_pos::Real; tol=1e-6)
    segments = Tuple{Vector{Float64}, Vector{Float64}}[]

    for face in faces
        v1 = vertices[face[1]]
        v2 = vertices[face[2]]
        v3 = vertices[face[3]]

        # Find intersection points of triangle with y = y_pos plane
        points = find_triangle_plane_intersection(v1, v2, v3, y_pos; tol)

        if length(points) == 2
            # Valid intersection - two points form a line segment
            push!(segments, (points[1], points[2]))
        end
    end

    return segments
end

"""
    find_triangle_plane_intersection(v1, v2, v3, y_pos; tol=1e-6)

Find intersection points of a triangle with the plane y = y_pos.

Returns 0, 1, or 2 points in [x, z] format.
"""
function find_triangle_plane_intersection(v1, v2, v3, y_pos; tol=1e-6)
    points = Vector{Float64}[]

    # Check each edge
    for (va, vb) in [(v1, v2), (v2, v3), (v3, v1)]
        y1, y2 = va[2], vb[2]

        # Check if edge crosses the plane
        if (y1 - y_pos) * (y2 - y_pos) < 0
            # Linear interpolation to find crossing point
            t = (y_pos - y1) / (y2 - y1)
            x = va[1] + t * (vb[1] - va[1])
            z = va[3] + t * (vb[3] - va[3])
            push!(points, [x, z])
        elseif abs(y1 - y_pos) < tol
            # Vertex on plane
            push!(points, [va[1], va[3]])
        end
    end

    # Remove duplicates
    unique_points = Vector{Float64}[]
    for p in points
        is_dup = false
        for up in unique_points
            if norm(p - up) < tol
                is_dup = true
                break
            end
        end
        if !is_dup
            push!(unique_points, p)
        end
    end

    return unique_points
end

"""
    order_segments_to_contour(segments; tol=1e-4)

Order line segments into a continuous contour.

# Arguments
- `segments`: Vector of (p1, p2) tuples representing line segments

# Returns
- Vector of [x, z] points forming the ordered contour
"""
function order_segments_to_contour(segments; tol=1e-4)
    if isempty(segments)
        return Vector{Float64}[]
    end

    # Build adjacency - find which segments connect
    n = length(segments)
    used = falses(n)

    # Start with first segment
    contour = [segments[1][1], segments[1][2]]
    used[1] = true
    n_used = 1

    # Iteratively find and add connecting segments
    max_iter = n * 2
    iter = 0
    while n_used < n && iter < max_iter
        iter += 1
        found = false

        current_end = contour[end]
        current_start = contour[1]

        for i in 1:n
            if used[i]
                continue
            end

            p1, p2 = segments[i]

            # Check if segment connects to end of contour
            if norm(p1 - current_end) < tol
                push!(contour, p2)
                used[i] = true
                n_used += 1
                found = true
                break
            elseif norm(p2 - current_end) < tol
                push!(contour, p1)
                used[i] = true
                n_used += 1
                found = true
                break
            # Check if segment connects to start of contour
            elseif norm(p1 - current_start) < tol
                pushfirst!(contour, p2)
                used[i] = true
                n_used += 1
                found = true
                break
            elseif norm(p2 - current_start) < tol
                pushfirst!(contour, p1)
                used[i] = true
                n_used += 1
                found = true
                break
            end
        end

        if !found
            # No connecting segment found, might have disjoint contours
            # Try to start a new contour segment (for complex geometries)
            for i in 1:n
                if !used[i]
                    push!(contour, segments[i][1])
                    push!(contour, segments[i][2])
                    used[i] = true
                    n_used += 1
                    break
                end
            end
        end
    end

    return contour
end

"""
    contour_to_airfoil(contour::Vector{Vector{Float64}})

Convert a 2D contour to normalized airfoil coordinates.

Assumes contour is in [x, z] format where x is chordwise, z is thickness direction.

# Returns
- (x, y): Normalized airfoil coordinates with chord = 1, LE at origin
"""
function contour_to_airfoil(contour::Vector{Vector{Float64}})
    if isempty(contour)
        return Float64[], Float64[]
    end

    # Extract x (chordwise) and y (thickness) coordinates
    x = [p[1] for p in contour]
    y = [p[2] for p in contour]

    # Find leading edge (minimum x)
    le_idx = argmin(x)
    x_le, y_le = x[le_idx], y[le_idx]

    # Find trailing edge (maximum x)
    te_idx = argmax(x)
    x_te, y_te = x[te_idx], y[te_idx]

    # Chord length
    chord = sqrt((x_te - x_le)^2 + (y_te - y_le)^2)

    if chord < 1e-10
        @warn "Degenerate airfoil with near-zero chord"
        return Float64[], Float64[]
    end

    # Rotation angle to align chord with x-axis
    angle = atan(y_te - y_le, x_te - x_le)
    cos_a, sin_a = cos(angle), sin(angle)

    # Transform to normalized coordinates
    x_norm = similar(x)
    y_norm = similar(y)
    for i in eachindex(x)
        dx = x[i] - x_le
        dy = y[i] - y_le
        x_norm[i] = (cos_a * dx + sin_a * dy) / chord
        y_norm[i] = (-sin_a * dx + cos_a * dy) / chord
    end

    return x_norm, y_norm
end

"""
    reorder_airfoil_selig(x::Vector, y::Vector)

Reorder airfoil coordinates to Selig format (TE upper -> LE -> TE lower).

This is the standard format expected by most airfoil tools.
"""
function reorder_airfoil_selig(x::Vector{T}, y::Vector{T}) where T
    n = length(x)
    if n < 3
        return x, y
    end

    # Find LE and TE indices
    le_idx = argmin(x)
    te_idx = argmax(x)

    # Determine if we need to reverse the ordering
    # Check if contour goes clockwise or counterclockwise from TE
    # by looking at the y-value trend after TE

    # Split into upper and lower based on y-values relative to camber line
    # For simplicity, assume points with y > mean(y at same x) are upper

    # First, separate points
    upper_pts = Tuple{T, T}[]
    lower_pts = Tuple{T, T}[]

    # Find approximate camber at TE and LE
    y_te = y[te_idx]
    y_le = y[le_idx]

    for i in eachindex(x)
        # Simple heuristic: use sign of y to classify
        # (This works for symmetric and slightly cambered airfoils)
        if y[i] >= y_le
            push!(upper_pts, (x[i], y[i]))
        else
            push!(lower_pts, (x[i], y[i]))
        end
    end

    # If classification is bad (unbalanced), try different approach
    if length(upper_pts) < n/4 || length(lower_pts) < n/4
        # Split based on contour direction
        # LE is at minimum x, so split contour there
        before_le = 1:le_idx
        after_le = le_idx:n

        # Determine which half is upper based on y values
        y_before = mean(y[before_le])
        y_after = mean(y[after_le])

        if y_before > y_after
            upper_pts = [(x[i], y[i]) for i in before_le]
            lower_pts = [(x[i], y[i]) for i in after_le]
        else
            upper_pts = [(x[i], y[i]) for i in after_le]
            lower_pts = [(x[i], y[i]) for i in before_le]
        end
    end

    # Sort upper by decreasing x (TE to LE)
    sort!(upper_pts, by=p -> -p[1])
    # Sort lower by increasing x (LE to TE)
    sort!(lower_pts, by=p -> p[1])

    # Combine: TE upper -> LE -> TE lower
    x_selig = T[]
    y_selig = T[]

    for (px, py) in upper_pts
        push!(x_selig, px)
        push!(y_selig, py)
    end
    # Add lower points (skip first if it's duplicate of LE)
    for (i, (px, py)) in enumerate(lower_pts)
        if i == 1 && !isempty(upper_pts)
            ux, uy = upper_pts[end]
            if abs(px - ux) < 1e-6 && abs(py - uy) < 1e-6
                continue  # Skip duplicate LE
            end
        end
        push!(x_selig, px)
        push!(y_selig, py)
    end

    return x_selig, y_selig
end

"""
    slice_mesh_at_plane(vertices, faces, point, normal; tol=1e-6)

Slice a mesh with the plane through `point` with unit `normal`, returning the 3D
segments `(p1, p2)` where the surface crosses it. The plane may have any
orientation, so the cut can follow a curved or swept span.
"""
function slice_mesh_at_plane(vertices, faces, point, normal; tol=1e-6)
    n = normalize(normal)
    p0 = Float64.(point)
    segments = Tuple{Vector{Float64}, Vector{Float64}}[]
    for face in faces
        v1, v2, v3 = vertices[face[1]], vertices[face[2]], vertices[face[3]]
        pts = Vector{Float64}[]
        for (va, vb) in ((v1, v2), (v2, v3), (v3, v1))
            da = dot(va - p0, n)
            db = dot(vb - p0, n)
            if da * db < 0
                t = da / (da - db)
                push!(pts, va .+ t .* (vb .- va))
            elseif abs(da) < tol
                push!(pts, Float64.(va))
            end
        end
        unique_pts = Vector{Float64}[]
        for p in pts
            any(norm(p - q) < tol for q in unique_pts) || push!(unique_pts, p)
        end
        length(unique_pts) == 2 && push!(segments, (unique_pts[1], unique_pts[2]))
    end
    return segments
end

"""
    span_edge_curves(vertices; n_samples=60) -> (arc_length, leading, trailing)

Build leading-/trailing-edge polylines by splitting the vertices into `n_samples`
constant-`y` slabs and taking the most forward (min `x`) and most aft (max `x`)
point in each. Returns the cumulative leading-edge arc length and the `3 × N` edge
matrices, ordered by `y`.
"""
function span_edge_curves(vertices; n_samples=60)
    ys = [v[2] for v in vertices]
    y_min, y_max = extrema(ys)
    centers = range(y_min, y_max, n_samples)
    half = (y_max - y_min) / n_samples
    le = Vector{Float64}[]
    te = Vector{Float64}[]
    for yc in centers
        in_slab = [v for v in vertices if yc - half ≤ v[2] ≤ yc + half]
        isempty(in_slab) && continue
        push!(le, Float64.(argmin(v -> v[1], in_slab)))
        push!(te, Float64.(argmax(v -> v[1], in_slab)))
    end
    leading = reduce(hcat, le)
    trailing = reduce(hcat, te)
    arc_length = zeros(size(leading, 2))
    for i in 2:length(arc_length)
        arc_length[i] = arc_length[i-1] + norm(leading[:, i] - leading[:, i-1])
    end
    return arc_length, leading, trailing
end

"""
    plane_contour_to_airfoil(contour3d, LE_point, TE_point, span_dir)

Project a 3D slice contour into the local chord (`TE_point - LE_point`) /
thickness (`span_dir × chord`, +z up) plane and return normalized Selig airfoil
coordinates `(x, y)`. Returns `nothing` for a degenerate slice.
"""
function plane_contour_to_airfoil(contour3d, LE_point, TE_point, span_dir)
    chord = TE_point - LE_point
    chord .-= dot(chord, span_dir) .* span_dir
    chord_len = norm(chord)
    chord_len < 1e-9 && return nothing
    chord ./= chord_len
    up = cross(span_dir, chord)
    norm(up) < 1e-12 && return nothing
    up ./= norm(up)
    up[3] < 0 && (up = -up)
    projected = [[dot(p - LE_point, chord), dot(p - LE_point, up)] for p in contour3d]
    x, y = contour_to_airfoil(projected)
    isempty(x) && return nothing
    return reorder_airfoil_selig(x, y)
end

"""
    perpendicular_sections(vertices, faces, n_sections; n_samples=60)

Extract `n_sections` airfoil cross-sections following a curved or swept span.
Stations are spaced at equal leading-edge arc-length intervals (panel centers);
each is sliced with a plane perpendicular to the local span and projected into the
chord/thickness plane, keeping the airfoil undistorted near curved tips. Returns
`(; LE_point, TE_point, x_airfoil, y_airfoil)` NamedTuples with un-normalized 3D
edge points.
"""
function perpendicular_sections(vertices, faces, n_sections; n_samples=60)
    arc_length, leading, trailing = span_edge_curves(vertices; n_samples)
    keep = [1; [i for i in 2:length(arc_length)
                if arc_length[i] > arc_length[i-1] + 1e-12]]
    arc_length, leading, trailing = arc_length[keep], leading[:, keep], trailing[:, keep]
    total = arc_length[end]
    total < 1e-9 && error("Degenerate span: leading-edge length ≈ 0")
    le_interp = ntuple(i -> linear_interpolation(arc_length, leading[i, :],
                       extrapolation_bc=Line()), 3)
    te_interp = ntuple(i -> linear_interpolation(arc_length, trailing[i, :],
                       extrapolation_bc=Line()), 3)
    delta = total / n_sections
    arc_step = total * 1e-3
    sections = NamedTuple[]
    for k in 1:n_sections
        sk = delta / 2 + (k - 1) * delta
        LE_point = [le_interp[i](sk) for i in 1:3]
        TE_point = [te_interp[i](sk) for i in 1:3]
        sp = min(sk + arc_step, total)
        sm = max(sk - arc_step, 0.0)
        span_dir = [le_interp[i](sp) - le_interp[i](sm) for i in 1:3]
        norm(span_dir) < 1e-12 && continue
        span_dir ./= norm(span_dir)
        segments = slice_mesh_at_plane(vertices, faces, LE_point, span_dir)
        isempty(segments) && continue
        contour3d = order_segments_to_contour(segments)
        length(contour3d) < 5 && continue
        airfoil = plane_contour_to_airfoil(contour3d, LE_point, TE_point, span_dir)
        airfoil === nothing && continue
        push!(sections, (; LE_point, TE_point,
                         x_airfoil=airfoil[1], y_airfoil=airfoil[2]))
    end
    return sections
end

"""
    slice_obj_at_positions(obj_path::String, y_positions::Vector)

Slice an OBJ mesh at multiple spanwise positions.

# Arguments
- `obj_path`: Path to .obj file
- `y_positions`: Vector of y-coordinates for slicing

# Returns
- Vector of (x, y) tuples, one per slice, with normalized airfoil coordinates
"""
function slice_obj_at_positions(obj_path::String, y_positions::Vector{T}) where T
    # Load mesh
    vertices, faces = read_faces(obj_path)

    # Slice at each position
    airfoils = Tuple{Vector{Float64}, Vector{Float64}}[]

    for y_pos in y_positions
        segments = slice_mesh_at_y(vertices, faces, Float64(y_pos))

        if isempty(segments)
            @warn "No intersection at y = $y_pos"
            push!(airfoils, (Float64[], Float64[]))
            continue
        end

        contour = order_segments_to_contour(segments)
        x, y = contour_to_airfoil(contour)

        if !isempty(x)
            x, y = reorder_airfoil_selig(x, y)
        end

        push!(airfoils, (x, y))
    end

    return airfoils
end

"""
    slice_obj_wing(obj_path::String, n_slices::Int)

Slice an OBJ wing mesh into n_slices evenly spaced airfoil sections.

Slices are positioned at panel centers: with n_slices panels, the positions are
at γ = δ/2, 3δ/2, 5δ/2, ... where δ = 1/n_slices. This avoids degenerate
tip profiles.

# Arguments
- `obj_path`: Path to .obj file
- `n_slices`: Number of slices (panels)

# Returns
- Vector of (x, y) tuples with normalized airfoil coordinates
- Vector of y_positions where slices were taken
"""
function slice_obj_wing(obj_path::String, n_slices::Int)
    # Load mesh to find bounds
    vertices, faces = read_faces(obj_path)

    # Find y-extent of mesh
    y_coords = [v[2] for v in vertices]
    y_min, y_max = minimum(y_coords), maximum(y_coords)

    # Slice at panel centers: δ/2, 3δ/2, 5δ/2, ... (n_slices-1/2)δ
    # where δ = span / n_slices
    span = y_max - y_min
    delta = span / n_slices
    y_positions = [y_min + delta/2 + i*delta for i in 0:(n_slices-1)]

    airfoils = slice_obj_at_positions(obj_path, y_positions)

    return airfoils, y_positions
end

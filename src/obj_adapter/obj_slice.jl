"""
OBJ mesh slicing for extracting 2D airfoil cross-sections.

This module provides functions to slice 3D wing meshes at specific spanwise
positions to extract 2D airfoil profiles.
"""

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
    march_edges(vertices, faces; step) -> (; le, te, point, tangent, arclen)

March the leading edge outward from mid-span in both directions in steps of arc
length `step`. Each cut is a vertical spanwise plane (both the chordwise and vertical
components of the running LE tangent dropped from the normal) so a tip that curls
downward can't tilt the plane toward horizontal, where its min-chord "LE" pick would
jump across the wing. Marching stops when the leading edge stops advancing spanwise.
Cuts sample mesh *edges*, so the picks are robust to vertex density. Returns, ordered
along the span, the LE/TE points, each cut's plane origin and tangent, and the
cumulative LE arc length. Build the airfoil for a chosen station with `build_section`.
"""
function march_edges(vertices, faces; step)
    ys = [v[2] for v in vertices]
    y_min, y_max = extrema(ys)
    y_mid = (y_min + y_max) / 2

    # Vertical spanwise plane (z dropped): a tip curl mustn't tilt the cut horizontal.
    cut(point, tangent) = begin
        segments = slice_mesh_at_plane(vertices, faces, point,
                                       normalize([0.0, tangent[2], 0.0]))
        isempty(segments) && return nothing
        crossings = Vector{Float64}[]
        for (start_point, end_point) in segments
            push!(crossings, start_point)
            push!(crossings, end_point)
        end
        return (; le=argmin(pt -> pt[1], crossings), te=argmax(pt -> pt[1], crossings))
    end

    mid_cut = cut([0.0, y_mid, 0.0], [0.0, 1.0, 0.0])
    mid_cut === nothing && error("Could not slice mid-span; check mesh / rotation.")

    march(direction) = begin
        rows = NamedTuple[]
        prev_le = mid_cut.le
        tangent = [0.0, direction, 0.0]
        for _ in 1:10_000
            probe = prev_le .+ step .* tangent
            found = cut(probe, tangent)
            if found === nothing
                # Bisect the last step to land the tip station on the outermost cut
                # that still yields a valid (non-degenerate) airfoil section.
                lo, hi, tip = 0.0, step, nothing
                for _ in 1:24
                    mid = (lo + hi) / 2
                    probe_mid = prev_le .+ mid .* tangent
                    here = cut(probe_mid, tangent)
                    # Reject a near-degenerate tip slice whose min-chord "LE" has jumped
                    # chordwise (an artifact that would skew the leading-edge arc length).
                    valid = here !== nothing &&
                            abs(here.le[1] - prev_le[1]) < step &&
                            build_section(vertices, faces, here.le, here.te,
                                          probe_mid, tangent) !== nothing
                    valid ? (lo = mid; tip = (; here.le, here.te, point=probe_mid)) :
                            (hi = mid)
                end
                tip === nothing || push!(rows, (; tip.le, tip.te, tip.point, tangent))
                break
            end
            le_step = found.le .- prev_le
            (norm(le_step) < 1e-9 || le_step[2] * direction <= 0.0) && break
            push!(rows, (; found.le, found.te, point=probe, tangent))
            tangent = normalize(le_step)
            prev_le = found.le
        end
        return rows
    end

    center = (; mid_cut.le, mid_cut.te, point=[0.0, y_mid, 0.0], tangent=[0.0, 1.0, 0.0])
    rows = vcat(reverse(march(-1.0)), [center], march(1.0))
    arclen = zeros(length(rows))
    for i in 2:length(rows)
        arclen[i] = arclen[i-1] + norm(rows[i].le - rows[i-1].le)
    end
    return (; le=[r.le for r in rows], te=[r.te for r in rows],
            point=[r.point for r in rows], tangent=[r.tangent for r in rows], arclen)
end

"""
    densify_contour(contour, max_edge) -> Vector{Vector{Float64}}

Insert evenly spaced points along contour edges longer than `max_edge`. Coarse mesh
triangles otherwise leave large hops in the slice cloud, which force the shrink
wrap's auto-raised rolling ball far up and over-smooth the wrapped airfoil.
"""
function densify_contour(contour, max_edge)
    out = Vector{Float64}[]
    for i in 1:length(contour)-1
        p1, p2 = contour[i], contour[i+1]
        n = clamp(ceil(Int, norm(p2 .- p1) / max_edge), 1, 100)
        for k in 0:n-1
            push!(out, p1 .+ (p2 .- p1) .* (k / n))
        end
    end
    push!(out, contour[end])
    return out
end

"""
    build_section(vertices, faces, le, te, point, tangent) -> section or nothing

Slice the mesh at one marched station ([`airfoil_frame`](@ref) drops the chordwise
tilt) and project it, producing the full
`(; LE_point, TE_point, span_dir, contour3d, x_airfoil, y_airfoil)`, or `nothing`.
"""
function build_section(vertices, faces, le, te, point, tangent)
    frame = airfoil_frame(le, te, tangent)
    frame === nothing && return nothing
    x_af, y_af, z_af = frame
    segments = slice_mesh_at_plane(vertices, faces, point, y_af)
    isempty(segments) && return nothing
    contour = order_segments_to_contour(segments)
    length(contour) < 5 && return nothing
    contour = densify_contour(contour, 0.005 * norm(te .- le))
    af = plane_contour_to_airfoil(contour, le, x_af, z_af)
    af === nothing && return nothing
    return (; LE_point=le, TE_point=te, span_dir=y_af, contour3d=contour,
            x_airfoil=af[1], y_airfoil=af[2])
end

"""
    airfoil_frame(LE_point, TE_point, span_tangent) -> (x_af, y_af, z_af)

Local airfoil axes for a slice, built so the slice plane contains the chord. With
`x̂ = [1,0,0]` and the LE `span_tangent`: `ẑ = x̂ × span` (up), `ŷ = ẑ × x̂` (spanwise
slice normal — the span projected into the `y`-`z` plane), `x̂_c = ŷ × ẑ` (chord,
oriented LE → TE). Slice the mesh with `y_af`. Returns `nothing` if the span is
parallel to global-x.
"""
function airfoil_frame(LE_point, TE_point, span_tangent)
    xg = [1.0, 0.0, 0.0]
    z_af = cross(xg, span_tangent)
    norm(z_af) < 1e-9 && return nothing
    z_af = normalize(z_af)
    z_af[3] < 0 && (z_af = -z_af)
    y_af = normalize(cross(z_af, xg))
    x_af = normalize(cross(y_af, z_af))
    dot(x_af, TE_point .- LE_point) < 0 && (x_af = -x_af)
    return x_af, y_af, z_af
end

"""
    plane_contour_to_airfoil(contour3d, LE_point, x_af, z_af)

Project a 3D slice contour into the local airfoil frame (chord axis `x_af`, up axis
`z_af`) and return normalized Selig coordinates `(x, y)`, or `nothing` if degenerate.
"""
function plane_contour_to_airfoil(contour3d, LE_point, x_af, z_af)
    projected = [[dot(p .- LE_point, x_af), dot(p .- LE_point, z_af)] for p in contour3d]
    x, y = contour_to_airfoil(projected)
    isempty(x) && return nothing
    return reorder_airfoil_selig(x, y)
end

"""
    perpendicular_sections(vertices, faces, n_sections; n_bins=60, rotation=I)

Extract `n_sections` airfoil cross-sections following a curved or swept span. The
leading edge is marched into `n_bins` stations (`march_edges`); the airfoil is
built (`build_section`) at the marched station nearest each equal
**leading-edge arc-length** target. Each section is
`(; LE_point, TE_point, span_dir, contour3d, x_airfoil, y_airfoil)`.

Stations closed to a point at the tips are skipped ([`station_indices`](@ref)), and
`wingtip_distance` insets the outermost sections a further arc length.

The slicer assumes `x` = chordwise, `y` = spanwise, `z` = up. Pass a `3×3` rotation
matrix to reorient a mesh stored in another convention before slicing.
"""
function perpendicular_sections(vertices, faces, n_sections; n_bins=60, rotation=I,
                                wingtip_distance=0.0)
    rotation === I || (vertices = [rotation * v for v in vertices])
    ys = [v[2] for v in vertices]
    m = march_edges(vertices, faces; step=(maximum(ys) - minimum(ys)) / n_bins)
    out = NamedTuple[]
    for i in station_indices(m, n_sections; wingtip_distance)
        sec = build_section(vertices, faces, m.le[i], m.te[i], m.point[i], m.tangent[i])
        sec === nothing || push!(out, sec)
    end
    return out
end

"""
    station_indices(march, n; wingtip_distance=0.0, min_chord_frac=0.01) -> Vector{Int}

Indices of the [`march_edges`](@ref) stations nearest `n` targets spread over the
leading-edge arc length. Stations whose chord has closed to less than
`min_chord_frac` of the longest one are left out of that range first, so a wing
tapering to a point puts its outermost sections on the last stations that still
have an airfoil to slice rather than on the point itself. The remaining first and
last targets sit a further `wingtip_distance` (arc length) inboard.
"""
function station_indices(march, n; wingtip_distance=0.0, min_chord_frac=0.01)
    chords = [norm(te .- le) for (le, te) in zip(march.le, march.te)]
    usable = findall(≥(min_chord_frac * maximum(chords)), chords)
    arclen = march.arclen
    inner, outer = arclen[first(usable)], arclen[last(usable)]
    d = clamp(wingtip_distance, 0.0, (outer - inner) / 2)
    n == 1 && return [argmin(abs.(arclen .- (inner + outer) / 2))]
    return [argmin(abs.(arclen .- t)) for t in range(inner + d, outer - d, n)]
end

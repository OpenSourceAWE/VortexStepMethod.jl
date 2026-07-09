"""
    ShrinkWrap(; clearance=0.006, min_concave_radius=0.02, cell_size=0.001,
               n_points=120, curvature_weight=0.05)

Distance-field shrink wrap: the rolling-ball offset of a raw slice point cloud,
extracted as one closed airfoil contour. A grid distance field is thresholded at
the rolling-ball radius (bridging cloud gaps and crevices narrower than about
twice `min_concave_radius`), flood-filled, eroded back to `clearance`, and the
resulting level set is traced with marching squares, faired (each pass clamped so
the contour keeps `clearance`) and resampled. Because the contour is
parameterized by arclength rather than `x`, the leading edge comes out genuinely
round — the offset of the cloud nose — the blunt trailing edge is capped by an arc
of radius `clearance`, and a single-membrane cloud becomes a thin capsule (at least
one `cell_size` half-thickness).

# Fields
- `clearance`: offset the contour keeps outside every cloud point; floored at one
  `cell_size` (the grid cannot represent a tighter wrap). Also the radius every
  convex corner is rounded at.
- `min_concave_radius`: rolling-ball radius — concave features narrower than about
  twice this are bridged by a fillet of roughly this radius; convex geometry is
  unaffected. Auto-raised so the ball can neither fall through the cloud's largest
  point gap nor pinch off between points during erosion, so sparse clouds get a
  correspondingly looser, smoother wrap.
- `cell_size`: distance-field grid resolution as a chord fraction; sets the
  geometric fidelity of the wrap.
- `n_points`: output stations per surface; the contour has `2*n_points - 1` points,
  cosine-clustered in arclength at the leading and trailing edges.
- `curvature_weight`: extra sampling measure per radian of contour turning (chord
  fraction), concentrating output points into corners so XFoil's spline can follow
  them; `0` gives plain cosine-in-arclength sampling.
"""
@with_kw struct ShrinkWrap
    clearance::Float64 = 0.006
    min_concave_radius::Float64 = 0.02
    cell_size::Float64 = 0.001
    n_points::Int = 120
    curvature_weight::Float64 = 0.05
end

"""
    distance_parabolas!(d, f, n, v, z)

One pass of the Felzenszwalb–Huttenlocher distance transform: writes into `d` the
lower envelope of the parabolas `(q - p)^2 + f[p]` over `p`, for `q in 1:n`.
`v` and `z` are scratch (length `n` and `n + 1`).
"""
function distance_parabolas!(d, f, n, v, z)
    k = 1
    v[1] = 1
    z[1] = -Inf
    z[2] = Inf
    for q in 2:n
        s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * (q - v[k]))
        while s <= z[k]
            k -= 1
            s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * (q - v[k]))
        end
        k += 1
        v[k] = q
        z[k] = s
        z[k+1] = Inf
    end
    k = 1
    for q in 1:n
        while z[k+1] < q
            k += 1
        end
        d[q] = (q - v[k])^2 + f[v[k]]
    end
    return d
end

"""
    squared_distance_transform!(field) -> field

In-place 2D squared Euclidean distance transform in grid-index units. On input
`field` holds `0` (or a local squared offset) at seeds and a large finite value
elsewhere; on output every node holds its squared distance to the nearest seed.
"""
function squared_distance_transform!(field::Matrix{Float64})
    nx, ny = size(field)
    m = max(nx, ny)
    d, f = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
    v, z = Vector{Int}(undef, m), Vector{Float64}(undef, m + 1)
    for i in 1:nx
        f[1:ny] .= @view field[i, :]
        distance_parabolas!(d, f, ny, v, z)
        field[i, :] .= @view d[1:ny]
    end
    for j in 1:ny
        f[1:nx] .= @view field[:, j]
        distance_parabolas!(d, f, nx, v, z)
        field[:, j] .= @view d[1:nx]
    end
    return field
end

"""
    flood_outside(blocked) -> BitMatrix

Mark the nodes reachable from the grid border without entering `blocked`
(4-connected flood fill); unreached nodes are the solid plus its enclosed holes.
"""
function flood_outside(blocked::AbstractMatrix{Bool})
    nx, ny = size(blocked)
    outside = falses(nx, ny)
    stack = Int[]
    visit(i, j) = if !blocked[i, j] && !outside[i, j]
        outside[i, j] = true
        push!(stack, (j - 1) * nx + i)
    end
    for i in 1:nx
        visit(i, 1)
        visit(i, ny)
    end
    for j in 1:ny
        visit(1, j)
        visit(nx, j)
    end
    while !isempty(stack)
        idx = pop!(stack)
        i, j = (idx - 1) % nx + 1, (idx - 1) ÷ nx + 1
        i > 1 && visit(i - 1, j)
        i < nx && visit(i + 1, j)
        j > 1 && visit(i, j - 1)
        j < ny && visit(i, j + 1)
    end
    return outside
end

"""
    trace_level_set(field, level) -> Vector{Vector{NTuple{2,Float64}}}

Marching-squares contours of `field .== level`, each returned as a closed loop of
grid-frame vertices (unit = one cell, node 1 at 1.0). Saddle cells are resolved by
the cell-center average. Loops touching the grid border are dropped.
"""
function trace_level_set(field::Matrix{Float64}, level::Float64)
    nx, ny = size(field)
    inside(v) = v > level
    frac(a, b) = (level - a) / (b - a)
    points = Dict{NTuple{3,Int},NTuple{2,Float64}}()
    links = Dict{NTuple{3,Int},Vector{NTuple{3,Int}}}()
    connect(a, b) = begin
        push!(get!(links, a, NTuple{3,Int}[]), b)
        push!(get!(links, b, NTuple{3,Int}[]), a)
    end
    for j in 1:ny-1, i in 1:nx-1
        v00, v10 = field[i, j], field[i+1, j]
        v01, v11 = field[i, j+1], field[i+1, j+1]
        b00, b10, b01, b11 = inside(v00), inside(v10), inside(v01), inside(v11)
        b00 == b10 == b01 == b11 && continue
        bottom, top = (i, j, 0), (i, j + 1, 0)
        left, right = (i, j, 1), (i + 1, j, 1)
        b00 != b10 && (points[bottom] = (i + frac(v00, v10), Float64(j)))
        b01 != b11 && (points[top] = (i + frac(v01, v11), Float64(j + 1)))
        b00 != b01 && (points[left] = (Float64(i), j + frac(v00, v01)))
        b10 != b11 && (points[right] = (Float64(i + 1), j + frac(v10, v11)))
        crossed = [e for (e, c) in ((bottom, b00 != b10), (top, b01 != b11),
                                    (left, b00 != b01), (right, b10 != b11)) if c]
        if length(crossed) == 4
            if inside((v00 + v10 + v01 + v11) / 4) == b00
                connect(bottom, right)
                connect(top, left)
            else
                connect(bottom, left)
                connect(top, right)
            end
        else
            connect(crossed[1], crossed[2])
        end
    end
    visited = Set{NTuple{3,Int}}()
    loops = Vector{Vector{NTuple{2,Float64}}}()
    for start in keys(links)
        start in visited && continue
        loop = NTuple{2,Float64}[]
        prev, current = start, start
        closed = false
        while true
            push!(visited, current)
            push!(loop, points[current])
            nbrs = links[current]
            length(nbrs) == 2 || break
            prev, current = current, (nbrs[1] == prev ? nbrs[2] : nbrs[1])
            current == start && (closed = true; break)
        end
        closed && push!(loops, loop)
    end
    return loops
end

"""
    grid_sampler(field, x0, y0, cell) -> f(px, py)

Bilinear interpolant of `field` at physical points; node `(i, j)` sits at
`(x0 + (i-1)*cell, y0 + (j-1)*cell)`. Queries are clamped to the grid.
"""
function grid_sampler(field::Matrix{Float64}, x0, y0, cell)
    nx, ny = size(field)
    return function (px, py)
        gi = clamp((px - x0) / cell + 1, 1.0, nx - 1e-9)
        gj = clamp((py - y0) / cell + 1, 1.0, ny - 1e-9)
        i, j = floor(Int, gi), floor(Int, gj)
        ti, tj = gi - i, gj - j
        return (field[i, j] * (1 - ti) + field[i+1, j] * ti) * (1 - tj) +
               (field[i, j+1] * (1 - ti) + field[i+1, j+1] * ti) * tj
    end
end

"""
    largest_linking_gap(x, y) -> Float64

Largest edge of the Euclidean minimum spanning tree of the points (Prim, O(n²)) —
the longest hop needed to keep the cloud connected. The rolling ball must not fall
through it.
"""
function largest_linking_gap(x, y)
    n = length(x)
    n <= 1 && return 0.0
    in_tree = falses(n)
    in_tree[1] = true
    best = [(x[k] - x[1])^2 + (y[k] - y[1])^2 for k in 1:n]
    gap2 = 0.0
    for _ in 2:n
        j, closest = 0, Inf
        for k in 1:n
            if !in_tree[k] && best[k] < closest
                j, closest = k, best[k]
            end
        end
        gap2 = max(gap2, closest)
        in_tree[j] = true
        for k in 1:n
            if !in_tree[k]
                d2 = (x[k] - x[j])^2 + (y[k] - y[j])^2
                d2 < best[k] && (best[k] = d2)
            end
        end
    end
    return sqrt(gap2)
end

"""
    resample_arc(ax, ay, n, curvature_weight) -> (x, y)

Resample the polyline `(ax, ay)` at `n` stations cosine-clustered in a measure that
blends arclength with `curvature_weight` extra length per radian of turning, so both
ends attract points and sharp features are resolved by several panels regardless of
their size (like XFoil's curvature-attracted paneling). Endpoints are preserved.
"""
function resample_arc(ax, ay, n, curvature_weight)
    m = length(ax)
    turn = zeros(m)
    for k in 2:m-1
        a = atan(ay[k] - ay[k-1], ax[k] - ax[k-1])
        b = atan(ay[k+1] - ay[k], ax[k+1] - ax[k])
        turn[k] = abs(rem(b - a, 2pi, RoundNearest))
    end
    w = zeros(m)
    for k in 2:m
        w[k] = w[k-1] + hypot(ax[k] - ax[k-1], ay[k] - ay[k-1]) +
               curvature_weight * (turn[k-1] + turn[k]) / 2
    end
    x_out, y_out = zeros(n), zeros(n)
    for (k, t) in enumerate(w[end] .* (1 .- cos.(range(0.0, pi, n))) ./ 2)
        idx = clamp(searchsortedlast(w, t), 1, m - 1)
        f = clamp((t - w[idx]) / max(w[idx+1] - w[idx], eps()), 0.0, 1.0)
        x_out[k] = (1 - f) * ax[idx] + f * ax[idx+1]
        y_out[k] = (1 - f) * ay[idx] + f * ay[idx+1]
    end
    return x_out, y_out
end

"""
    shrink_wrap(x, y, method::ShrinkWrap) -> (x, y)

Wrap the point cloud `(x, y)` into a clean closed airfoil in Selig order (TE upper →
LE → TE lower), following [`ShrinkWrap`](@ref): distance field on a `cell_size`
grid, closing with the rolling ball (`min_concave_radius`), offset outward by
`clearance`, traced as a single
closed contour and resampled to cosine panels in a curvature-weighted arclength
measure. The first and last point coincide at the
trailing edge (the TE cap is part of the contour). The output stays in the
normalized frame of the input cloud (chord slightly longer than 1, nose apex near
`x = -clearance`) and is ready to write as a `.dat` or fit with
[`LeastSquaresFit`](@ref).
"""
function shrink_wrap(x, y, method::ShrinkWrap)
    xn, yn, _ = normalize_airfoil(collect(float.(x)), collect(float.(y)))
    cell = method.cell_size
    gap = max(method.clearance, cell)
    linking = largest_linking_gap(xn, yn)
    ball = max(method.min_concave_radius, gap + 3cell,
               1.5 * (linking^2 / 4 + gap^2) / (2gap))
    pad = ball + 3cell
    x0, y0 = -pad, minimum(yn) - pad
    nx = ceil(Int, (1.0 + pad - x0) / cell) + 1
    ny = ceil(Int, (maximum(yn) + pad - y0) / cell) + 1

    far = Float64(nx^2 + ny^2)
    field = fill(far, nx, ny)
    for (px, py) in zip(xn, yn)
        gi, gj = (px - x0) / cell + 1, (py - y0) / cell + 1
        for i in (floor(Int, gi), floor(Int, gi) + 1),
            j in (floor(Int, gj), floor(Int, gj) + 1)

            (1 <= i <= nx && 1 <= j <= ny) || continue
            field[i, j] = min(field[i, j], (gi - i)^2 + (gj - j)^2)
        end
    end
    squared_distance_transform!(field)
    cloud_dist = sqrt.(field) .* cell

    outside = flood_outside(cloud_dist .<= ball)
    depth = [outside[i, j] ? 0.0 : far for i in 1:nx, j in 1:ny]
    squared_distance_transform!(depth)
    depth .= sqrt.(depth) .* cell

    loops = trace_level_set(depth, ball - gap)
    isempty(loops) && error("ShrinkWrap traced no contour; the cloud may be too" *
                            " sparse or thin for cell_size=$(cell).")
    shoelace(loop) = abs(sum(p[1] * loop[mod1(k + 1, length(loop))][2] -
                             loop[mod1(k + 1, length(loop))][1] * p[2]
                             for (k, p) in enumerate(loop))) / 2
    areas = shoelace.(loops)
    main = loops[argmax(areas)]
    if length(loops) > 1 && sort(areas)[end-1] > 0.01 * maximum(areas)
        @warn "ShrinkWrap discarded a comparable secondary contour; increase" *
              " min_concave_radius to bridge cloud gaps."
    end
    px = [x0 + (p[1] - 1) * cell for p in main]
    py = [y0 + (p[2] - 1) * cell for p in main]

    dist_at = grid_sampler(cloud_dist, x0, y0, cell)
    keep_clear(cx, cy) = begin
        d = dist_at(cx, cy)
        d >= gap && return (cx, cy)
        h = cell / 2
        gx = (dist_at(cx + h, cy) - dist_at(cx - h, cy)) / cell
        gy = (dist_at(cx, cy + h) - dist_at(cx, cy - h)) / cell
        len = hypot(gx, gy)
        len < 1e-9 && return (cx, cy)
        return (cx + (gap - d) * gx / len, cy + (gap - d) * gy / len)
    end
    m = length(px)
    for pass in 0:20
        ox, oy = copy(px), copy(py)
        for k in 1:m
            if pass > 0
                a, b = mod1(k - 1, m), mod1(k + 1, m)
                px[k] = ox[k] + 0.4 * (ox[a] - 2ox[k] + ox[b])
                py[k] = oy[k] + 0.4 * (oy[a] - 2oy[k] + oy[b])
            end
            px[k], py[k] = keep_clear(px[k], py[k])
        end
    end

    anchor = argmax(px)
    px, py = circshift(px, 1 - anchor), circshift(py, 1 - anchor)
    le = argmax((px .- px[1]) .^ 2 .+ (py .- py[1]) .^ 2)
    if sum(@view py[1:le]) / le < sum(@view py[le:end]) / (m - le + 1)
        reverse!(px)
        reverse!(py)
        px, py = circshift(px, 1), circshift(py, 1)
        le = argmax((px .- px[1]) .^ 2 .+ (py .- py[1]) .^ 2)
    end
    xu, yu = resample_arc(px[1:le], py[1:le], method.n_points,
                          method.curvature_weight)
    xl, yl = resample_arc(vcat(px[le:end], px[1]), vcat(py[le:end], py[1]),
                          method.n_points, method.curvature_weight)
    return vcat(xu, xl[2:end]), vcat(yu, yl[2:end])
end

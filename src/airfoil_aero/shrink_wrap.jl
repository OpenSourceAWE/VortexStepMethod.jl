"""
    ShrinkWrap(; clearance=0.006, min_concave_radius=0.02, min_clearance=0.001,
               n_points=120, curvature_weight=0.05)

Rolling-ball shrink wrap: the offset of a raw slice point cloud, traced exactly as
one closed airfoil contour. A disk of radius `min_concave_radius` is pivoted around
the outside of the cloud and the polygon through the points it touches is offset by
`clearance` ([`pivot_contour`](@ref)). Being traced rather than sampled on a grid,
the wrap has no resolution floor: the leading edge comes out genuinely round — the
offset of the cloud nose — the blunt trailing edge is capped by an arc of radius
`clearance`, and a single-membrane cloud becomes a capsule of exactly that
half-thickness.

# Fields
- `clearance`: offset the contour keeps outside every cloud point, and the radius
  every convex corner is rounded at; floored at `min_clearance` only when the cloud
  is an open single membrane rather than a closed loop.
- `min_concave_radius`: rolling-ball radius — concave features narrower than about
  twice this are bridged straight; convex geometry is unaffected. Raised where the
  cloud is sparse enough that the ball would otherwise fall through its largest point
  gap.
- `min_clearance`: floor on `clearance` for an open single-membrane cloud, so
  `clearance=0` still leaves it a capsule with two distinguishable sides instead of
  a bare curve; a closed loop uses its true `clearance`, so `clearance=0` hugs the
  cloud and keeps its sharp trailing edge sharp. Also accepted under its former name
  `cell_size`.
- `n_points`: output stations per surface; the contour has `2*n_points - 1` points,
  cosine-clustered in arclength at the leading and trailing edges.
- `curvature_weight`: extra sampling measure per radian of contour turning (chord
  fraction), concentrating output points into corners so XFoil's spline can follow
  them; `0` gives plain cosine-in-arclength sampling.
"""
struct ShrinkWrap
    clearance::Float64
    min_concave_radius::Float64
    min_clearance::Float64
    n_points::Int
    curvature_weight::Float64
end

function ShrinkWrap(; clearance=0.006, min_concave_radius=0.02, min_clearance=0.001,
                    cell_size=nothing, n_points=120, curvature_weight=0.05)
    return ShrinkWrap(clearance, min_concave_radius,
                      isnothing(cell_size) ? min_clearance : cell_size,
                      n_points, curvature_weight)
end

"""
    point_buckets(x, y, reach) -> Dict

Sparse uniform grid of the point indices, one bucket per `reach`-sized cell, so every
point within `reach` of a query lies in one of the nine buckets around it.
"""
function point_buckets(x, y, reach)
    buckets = Dict{NTuple{2,Int},Vector{Int}}()
    for k in eachindex(x)
        key = (floor(Int, x[k] / reach), floor(Int, y[k] / reach))
        push!(get!(buckets, key, Int[]), k)
    end
    return buckets
end

"""
    turn_measure(cross, dot) -> Float64

Pseudo-angle in `[0, 4)` of the direction with the given `cross` and `dot` against a
reference, rising monotonically with the counterclockwise angle it stands for. Ordering
directions by it is ordering them by angle, without the `atan` that dominates the
pivot's inner loop.
"""
function turn_measure(cross, dot)
    if cross >= 0
        return dot >= 0 ? cross / (cross + dot) : 1 + (-dot) / (cross - dot)
    end
    return dot <= 0 ? 2 + (-cross) / (-cross - dot) : 3 + dot / (dot - cross)
end

"""
    pivot_step(x, y, buckets, r, p, centre) -> (q, centre) or nothing

Rotate the empty disk of radius `r` counterclockwise about point `p`, starting from
`centre`, until its rim reaches a second point. Returns that point and the disk
centre touching both, or `nothing` when nothing lies within `2r` of `p`.
"""
function pivot_step(x, y, buckets, r, p, centre)
    px, py = x[p], y[p]
    ux, uy = centre[1] - px, centre[2] - py
    best_turn, best_q, best_centre = Inf, 0, centre
    i0, j0 = floor(Int, px / 2r), floor(Int, py / 2r)
    for i in i0-1:i0+1, j in j0-1:j0+1
        bucket = get(buckets, (i, j), nothing)
        bucket === nothing && continue
        for q in bucket
            q == p && continue
            dx, dy = x[q] - px, y[q] - py
            span = dx * dx + dy * dy
            (span < 1e-24 || span > 4r * r) && continue
            out = sqrt(max(r * r / span - 0.25, 0.0))
            mx, my = dx / 2, dy / 2
            ox, oy = -dy * out, dx * out
            for (vx, vy) in ((mx + ox, my + oy), (mx - ox, my - oy))
                (vx - ux)^2 + (vy - uy)^2 < 1e-24 && continue
                turn = turn_measure(ux * vy - uy * vx, ux * vx + uy * vy)
                (turn < 1e-9 || turn >= best_turn) && continue
                best_turn, best_q, best_centre = turn, q, (px + vx, py + vy)
            end
        end
    end
    best_q == 0 && return nothing
    return best_q, best_centre
end

"""
    push_arc!(px, py, cx, cy, radius, from, sweep)

Append the arc about `(cx, cy)` that starts at angle `from` and turns by the signed
angle `sweep`, endpoints included, stepped fine enough to stay within `1e-6` of the
true arc. A zero radius appends the centre once.
"""
function push_arc!(px, py, cx, cy, radius, from, sweep)
    steps = radius < 1e-12 ? 0 :
            clamp(ceil(Int, abs(sweep) / sqrt(8e-6 / radius)), 1, 400)
    for k in 0:steps
        angle = from + sweep * k / max(steps, 1)
        push!(px, cx + radius * cos(angle))
        push!(py, cy + radius * sin(angle))
    end
    return nothing
end

"""
    edge_normal(ax, ay, bx, by) -> Float64

Angle of the outward normal of the wrap edge running from `(ax, ay)` to `(bx, by)`,
which [`pivot_contour`](@ref) traces with the cloud on its left.
"""
edge_normal(ax, ay, bx, by) = atan(ax - bx, by - ay)

"""
    pivot_contour(x, y, r, clearance) -> (px, py)

Trace the alpha shape of the cloud `(x, y)` for alpha `r`, offset outward by
`clearance`, as a closed polyline. A disk of radius `r` is pivoted around the outside
of the cloud from its leftmost point; the points it touches are the polygon's
vertices in order, the ones inside a concavity narrower than the disk having been
skipped and bridged straight. The offset rounds each convex vertex with an arc of
radius `clearance` and chamfers each reflex one across the bisector, so it holds that
distance from every contact.
"""
function pivot_contour(x, y, r, clearance)
    start = argmin(x)
    buckets = point_buckets(x, y, 2r)
    contacts, centres = Int[], NTuple{2,Float64}[]
    p, centre = start, (x[start] - r, y[start])
    # A single-membrane stretch is touched once from each side, so the cycle closes on
    # the (point, centre) pair rather than on the point alone.
    for _ in 1:4*length(x)+100
        step = pivot_step(x, y, buckets, r, p, centre)
        step === nothing && break
        next, centre = step
        !isempty(contacts) && p == contacts[1] &&
            hypot(centre[1] - centres[1][1], centre[2] - centres[1][2]) < 1e-12 &&
            break
        push!(contacts, p)
        push!(centres, centre)
        p = next
    end
    px, py = Float64[], Float64[]
    n = length(contacts)
    n < 3 && return px, py
    tail = contacts[n]
    from = edge_normal(x[tail], y[tail], x[contacts[1]], y[contacts[1]])
    for k in 1:n
        i, j = contacts[k], contacts[mod1(k + 1, n)]
        to = edge_normal(x[i], y[i], x[j], y[j])
        turn = rem(to - from, 2pi, RoundNearest)
        if turn >= 0
            push_arc!(px, py, x[i], y[i], clearance, from, turn)
        else
            push!(px, x[i] + clearance * cos(from + turn / 2))
            push!(py, y[i] + clearance * sin(from + turn / 2))
        end
        from = to
    end
    return px, py
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
    smoothed_curvature(node_arclength, turn, band) -> Vector{Float64}

Turning per unit arclength at each node, averaged over a window of length `band`
centred on it. Averaging over a length leaves the measure independent of how densely
the trace happens to be stepped, which varies by a factor of several between the
`clearance` arcs and the straight edges between them.
"""
function smoothed_curvature(node_arclength, turn, band)
    m = length(turn)
    cumulative_turn = cumsum(turn)
    density = zeros(m)
    lo, hi = 1, 1
    for k in 1:m
        while lo < m && node_arclength[lo] < node_arclength[k] - band / 2
            lo += 1
        end
        while hi < m && node_arclength[hi+1] <= node_arclength[k] + band / 2
            hi += 1
        end
        span = node_arclength[hi] - node_arclength[lo]
        density[k] = span > 0 ?
            (cumulative_turn[hi] - cumulative_turn[lo] + turn[lo]) / span : 0.0
    end
    return density
end

"""
    enforce_min_spacing!(stations, floor_length) -> stations

Push the sampling `stations` apart until no gap is shorter than `floor_length`, taking
the room out of the gaps that have it to spare so the first and last station stay put.
`floor_length` is capped at the mean gap, which a uniform spacing already meets.
"""
function enforce_min_spacing!(stations, floor_length)
    n = length(stations)
    n < 3 && return stations
    gap = diff(stations)
    total = sum(gap)
    finish = stations[n]
    floor_length = min(floor_length, total / (n - 1))
    for _ in 1:20
        deficit = sum(max(0.0, floor_length - g) for g in gap)
        deficit <= 1e-15 * total && break
        spare = sum(max(0.0, g - floor_length) for g in gap)
        spare <= 0 && break
        shrink = 1 - min(deficit / spare, 1.0)
        for k in eachindex(gap)
            gap[k] = gap[k] < floor_length ? floor_length :
                     floor_length + (gap[k] - floor_length) * shrink
        end
    end
    gap .*= total / sum(gap)
    for k in 2:n
        stations[k] = stations[k-1] + gap[k-1]
    end
    stations[n] = finish
    return stations
end

"""
    resample_arc(ax, ay, n, curvature_weight; min_panel_fraction=0.25) -> (x, y)

Resample the polyline `(ax, ay)` at `n` stations cosine-clustered in a measure that
blends arclength with `curvature_weight` extra length per radian of turning, so both
ends attract points and sharp features are resolved by several panels regardless of
their size (like XFoil's curvature-attracted paneling). The turning is taken as a
curvature density smoothed over a few panel lengths ([`smoothed_curvature`](@ref)), so
the measure grows smoothly and neighbouring panels stay comparable in length. No panel
is left shorter than `min_panel_fraction` of the mean
([`enforce_min_spacing!`](@ref)). Endpoints are preserved.
"""
function resample_arc(ax, ay, n, curvature_weight; min_panel_fraction=0.25)
    m = length(ax)
    turn = zeros(m)
    for k in 2:m-1
        a = atan(ay[k] - ay[k-1], ax[k] - ax[k-1])
        b = atan(ay[k+1] - ay[k], ax[k+1] - ax[k])
        turn[k] = abs(rem(b - a, 2pi, RoundNearest))
    end
    node_arclength = zeros(m)
    for k in 2:m
        node_arclength[k] = node_arclength[k-1] +
                            hypot(ax[k] - ax[k-1], ay[k] - ay[k-1])
    end
    curvature = smoothed_curvature(node_arclength, turn,
                                   3 * node_arclength[end] / max(n - 1, 1))
    w = zeros(m)
    for k in 2:m
        segment = node_arclength[k] - node_arclength[k-1]
        w[k] = w[k-1] + segment *
               (1 + curvature_weight * (curvature[k-1] + curvature[k]) / 2)
    end
    stations = w[end] .* (1 .- cos.(range(0.0, pi, n))) ./ 2
    enforce_min_spacing!(stations, min_panel_fraction * w[end] / max(n - 1, 1))
    x_out, y_out = zeros(n), zeros(n)
    for (k, t) in enumerate(stations)
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
LE → TE lower), following [`ShrinkWrap`](@ref): the rolling ball
(`min_concave_radius`) is pivoted around the cloud, its contact side offset outward
by `clearance` ([`pivot_contour`](@ref)), and the resulting arcs are resampled to
cosine panels in a curvature-weighted arclength measure. The first and last point
coincide at the
trailing edge (the TE cap is part of the contour). A closed loop (first and last
cloud points coincident) keeps its true `clearance`, so `clearance=0` hugs the input
and leaves a sharp trailing edge sharp; an open single-membrane cloud is floored at
`min_clearance`. The output stays in the
normalized frame of the input cloud (chord slightly longer than 1, nose apex near
`x = -clearance`) and is ready to write as a `.dat` or fit with
[`LeastSquaresFit`](@ref).
"""
function shrink_wrap(x, y, method::ShrinkWrap)
    xn, yn, _ = normalize_airfoil(collect(float.(x)), collect(float.(y)))
    closed = hypot(xn[end] - xn[1], yn[end] - yn[1]) < 0.02
    gap = closed ? method.clearance : max(method.clearance, method.min_clearance)
    linking = largest_linking_gap(xn, yn)
    ball = max(method.min_concave_radius, 1.01 * linking / 2)
    px, py = pivot_contour(xn, yn, ball, gap)
    length(px) < 3 && error("ShrinkWrap traced no contour; the cloud may be too" *
                            " sparse or thin to pivot a ball of $(ball) around.")
    m = length(px)

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

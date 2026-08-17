"""
    ShrinkWrap(; clearance=0.006, min_concave_radius=0.02, min_clearance=0.001,
               n_points=120, curvature_weight=0.05)

Rolling-ball shrink wrap: the offset of a raw slice point cloud, traced exactly as
one closed airfoil contour. A disk of radius `min_concave_radius` is pivoted around
the outside of the cloud ([`pivot_contour`](@ref)); the wrap is what its contact
side sweeps, so it consists of arcs of radius `clearance` about the points the disk
touches, joined by arcs of the disk radius across the gaps it cannot enter. Being
built from arcs rather than sampled on a grid, the wrap has no resolution floor:
the leading edge comes out genuinely round — the offset of the cloud nose — the
blunt trailing edge is capped by an arc of radius `clearance`, and a
single-membrane cloud becomes a capsule of exactly that half-thickness.

# Fields
- `clearance`: offset the contour keeps outside every cloud point, and the radius
  every convex corner is rounded at; floored at `min_clearance`.
- `min_concave_radius`: rolling-ball radius — concave features narrower than about
  twice this are bridged by a fillet of roughly this radius; convex geometry is
  unaffected. Auto-raised so the ball can neither fall through the cloud's largest
  point gap nor dip below `clearance` between neighbouring points, so sparse clouds
  get a correspondingly looser, smoother wrap.
- `min_clearance`: floor on `clearance`, so `clearance=0` still leaves a single
  membrane a capsule with two distinguishable sides instead of a bare curve. Also
  accepted under its former name `cell_size`, which set the resolution of the
  distance field this used to be traced on.
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
angle `sweep`, stepped fine enough to stay within `1e-6` of the true arc. The closing
point is left off so consecutive arcs concatenate without duplicates.
"""
function push_arc!(px, py, cx, cy, radius, from, sweep)
    steps = radius < 1e-12 ? 1 :
            clamp(ceil(Int, abs(sweep) / sqrt(8e-6 / radius)), 1, 400)
    for k in 0:steps-1
        angle = from + sweep * k / steps
        push!(px, cx + radius * cos(angle))
        push!(py, cy + radius * sin(angle))
    end
    return nothing
end

"""
    pivot_contour(x, y, r, clearance) -> (px, py)

Trace the rolling-ball closing of the cloud `(x, y)` offset outward by `clearance`,
as a dense closed polyline. A disk of radius `r` is pivoted around the outside of the
cloud from its leftmost point; each point the disk touches contributes an arc of
radius `clearance` about itself, and each pivot an arc of radius `r - clearance`
about the empty disk's centre, spanning the gap the disk could not enter. Those arcs
are the exact wrap boundary, so nothing here has a resolution to choose.
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
    isempty(contacts) && return px, py
    for (k, i) in enumerate(contacts)
        incoming = centres[mod1(k - 1, length(centres))]
        outgoing = centres[k]
        from = atan(incoming[2] - y[i], incoming[1] - x[i])
        to = atan(outgoing[2] - y[i], outgoing[1] - x[i])
        push_arc!(px, py, x[i], y[i], clearance, from, mod(to - from, 2pi))
        j = contacts[mod1(k + 1, length(contacts))]
        from = atan(y[i] - outgoing[2], x[i] - outgoing[1])
        to = atan(y[j] - outgoing[2], x[j] - outgoing[1])
        push_arc!(px, py, outgoing[1], outgoing[2], r - clearance, from,
                  rem(to - from, 2pi, RoundNearest))
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
centred on it. Averaging over a length rather than over a node count is what makes it
independent of how the traced arcs happen to be stepped: [`push_arc!`](@ref) steps each
arc to a chord tolerance, so its segments scale as `sqrt(radius)` and the trace's node
spacing alternates by a factor of several between the `clearance` arcs and the wider
pivot arcs. Diffusing per node would carry that alternation into the measure and leave
the resampled panel lengths oscillating, which XFoil's viscous solver cannot follow.
A sharp corner still refines a band of panels gradually, the way XFoil's `PANGEN`
bunches panels on a *smoothed* curvature.
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
Cosine clustering thins as the square of the station count, so at the counts a wrap
emits it puts panels a small fraction of a percent of chord at the leading edge, which
XFoil's viscous solver cannot march through however smoothly they are graded.
"""
function enforce_min_spacing!(stations, floor_length)
    n = length(stations)
    n < 3 && return stations
    gap = diff(stations)
    total = sum(gap)
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
trailing edge (the TE cap is part of the contour). The output stays in the
normalized frame of the input cloud (chord slightly longer than 1, nose apex near
`x = -clearance`) and is ready to write as a `.dat` or fit with
[`LeastSquaresFit`](@ref).
"""
function shrink_wrap(x, y, method::ShrinkWrap)
    xn, yn, _ = normalize_airfoil(collect(float.(x)), collect(float.(y)))
    gap = max(method.clearance, method.min_clearance)
    linking = largest_linking_gap(xn, yn)
    ball = max(method.min_concave_radius, gap, 1.01 * linking / 2,
               1.5 * (linking^2 / 4 + gap^2) / (2gap))
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

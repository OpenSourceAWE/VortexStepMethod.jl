"""
    ShrinkWrap(; clearance=0.006, smoothing=1.0, n_points=120)

Discrete shrink-wrap that turns a raw slice point cloud into a clean closed airfoil
contour. Two height curves (upper and lower) share a cosine `x/c` grid and are each
Laplacian-smoothed and then clipped back outside the cloud, so the result encloses
every point with a `clearance` gap while fairing over interior noise. Where the cloud
collapses to a single membrane the surfaces hold `clearance` apart, giving a thin
cambered airfoil. The leading- and trailing-edge grid endpoints are never smoothed, so
they stay as sharp as the geometry closes to.

# Fields
- `clearance`: gap each surface keeps outside every point, and the half-thickness where
  the cloud collapses to a membrane (`2*clearance` total).
- `smoothing`: higher smooths harder — rounder corners, more fairing over interior
  points; lower hugs the cloud more tightly with sharper corners.
- `n_points`: grid stations per surface; the contour has `2*n_points - 1` points,
  cosine-clustered at the leading and trailing edges like `Xfoil.pane`.
"""
@with_kw struct ShrinkWrap
    clearance::Float64 = 0.006
    smoothing::Float64 = 1.0
    n_points::Int = 120
end

"""
    shrink_wrap(x, y, method::ShrinkWrap) -> (x, y)

Wrap the point cloud `(x, y)` into a clean closed airfoil in Selig order (TE upper → LE
→ TE lower) on a cosine `x/c` grid, following [`ShrinkWrap`](@ref). The output encloses
every input point with `method.clearance` clearance and is ready to write as a `.dat`
or fit with [`LeastSquaresFit`](@ref).

Each surface starts at the cloud envelope (the extreme points ± clearance, with each
point raising both bracketing stations so linear interpolation encloses it exactly),
then alternates one Laplacian smoothing pass with a clip back to that envelope. The
number of passes scales with `smoothing` and grid resolution so the smoothing length is
independent of `n_points`.
"""
function shrink_wrap(x, y, method::ShrinkWrap)
    xn, yn, _ = normalize_airfoil(collect(float.(x)), collect(float.(y)))
    n = method.n_points
    stations = (1 .- cos.(range(0.0, pi, n))) ./ 2
    lambda = 0.4

    fill_gaps!(v) = begin
        finite = findall(isfinite, v)
        isempty(finite) && return fill!(v, 0.0)
        v[1:finite[1]-1] .= v[finite[1]]
        v[finite[end]+1:end] .= v[finite[end]]
        for k in 2:length(finite)
            a, b = finite[k-1], finite[k]
            for j in a+1:b-1
                s = (stations[j] - stations[a]) / (stations[b] - stations[a])
                v[j] = (1 - s) * v[a] + s * v[b]
            end
        end
        return v
    end

    envelope(upper) = begin
        env = fill(upper ? -Inf : Inf, n)
        for (xi, yi) in zip(xn, yn)
            j = clamp(searchsortedlast(stations, clamp(xi, 0.0, 1.0)), 1, n - 1)
            target = upper ? yi + method.clearance : yi - method.clearance
            for k in (j, j + 1)
                env[k] = upper ? max(env[k], target) : min(env[k], target)
            end
        end
        return fill_gaps!(env)
    end

    smooth!(v) = begin
        old = copy(v)
        for j in 2:n-1
            v[j] = old[j] + lambda * (old[j-1] - 2old[j] + old[j+1])
        end
        return v
    end

    env_up, env_lo = envelope(true), envelope(false)
    yu, yl = copy(env_up), copy(env_lo)
    for _ in 1:max(1, round(Int, method.smoothing * n^2 / 25))
        yu .= max.(smooth!(yu), env_up)
        yl .= min.(smooth!(yl), env_lo)
    end
    return vcat(reverse(stations), stations[2:end]), vcat(reverse(yu), yl[2:end])
end

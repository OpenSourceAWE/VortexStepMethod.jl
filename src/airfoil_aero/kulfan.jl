"""
Kulfan CST (Class Shape Transformation) airfoil parameterization.

This module provides functions to:
- Fit Kulfan parameters from airfoil coordinates
- Generate airfoil coordinates from Kulfan parameters
"""

"""
    KulfanParameters

Kulfan CST parameters for an airfoil.

# Fields
- `upper_weights::Vector{Float64}`: 8 weights for upper surface
- `lower_weights::Vector{Float64}`: 8 weights for lower surface
- `leading_edge_weight::Float64`: Leading edge modification weight
- `TE_thickness::Float64`: Trailing edge thickness
"""
struct KulfanParameters
    upper_weights::Vector{Float64}
    lower_weights::Vector{Float64}
    leading_edge_weight::Float64
    TE_thickness::Float64
end

"""
    KulfanFitMethod

Abstract supertype for Kulfan-parameter fitting strategies. Concrete methods
[`LeastSquaresFit`](@ref) and [`EnvelopeFit`](@ref) are passed to
[`fit_kulfan_parameters`](@ref), which dispatches on them.
"""
abstract type KulfanFitMethod end

"""
    LeastSquaresFit(; n_weights=8)

Least-squares Kulfan fit: solves a single system that fits both surfaces
*through* the points, matching AeroSandbox's `get_kulfan_parameters` (the
parameterization NeuralFoil was trained on). Sensitive to noisy interior points;
use [`EnvelopeFit`](@ref) for those.
"""
@with_kw struct LeastSquaresFit <: KulfanFitMethod
    n_weights::Int = 8
end

"""
    EnvelopeFit(; min_distance=0.0, tightness=1.0, distance_penalty=1e4,
                n_grid=80, n_weights=8, ...)

Enclosing Kulfan fit: an outer envelope that keeps every point at least `min_distance`
inside the surfaces (a hard constraint, so the fit always encloses the geometry) while
hugging them as `tightness` asks. Robust to noisy interior points that pull
[`LeastSquaresFit`](@ref) inward. See [`fit_kulfan_parameters`](@ref) for the algorithm.

# Fields
- `min_distance`: clearance each surface keeps outside every point
- `distance_penalty`: weight of the hard clearance constraint
- `tightness`: weight of the hug objective — pulls the surfaces together (minimising
  thickness), each held out by the clearance constraint; traded off against a fixed
  arc-length (smoothing) term
- `tightness_dist`: chordwise distribution scaling `tightness`, sampled at equal `x/c`
  and linearly interpolated. Lower it over a zone to fair smoothly instead of hugging —
  e.g. over a kite's LE tube, where the clearance constraint pins the upper surface to
  the tube crown so the lower surface floats out and fairs. Default `ones(n_grid)`.
- `n_weights`: weights per surface (NeuralFoil uses 8)
- `n_grid`: chordwise stations for the distributions and the perimeter
- `n_iter`, `tol`: iteration controls; `regularization`: Tikhonov weight toward the LSQ fit
"""
@with_kw struct EnvelopeFit <: KulfanFitMethod
    min_distance::Float64 = 0.0
    n_weights::Int = 8
    tightness::Float64 = 1.0
    distance_penalty::Float64 = 1e4
    n_grid::Int = 80
    tightness_dist::Vector{Float64} = ones(n_grid)
    n_iter::Int = 30
    regularization::Float64 = 1e-6
    tol::Float64 = 1e-10
end

"""
    bernstein_basis(x::AbstractVector, n::Int)

Compute Bernstein polynomial basis matrix.

Returns matrix of shape (length(x), n+1) where each column is a Bernstein polynomial
B_i(x) = binomial(n,i) * x^i * (1-x)^(n-i)
"""
function bernstein_basis(x::AbstractVector{T}, n::Int) where T
    m = length(x)
    B = zeros(T, m, n + 1)
    for i in 0:n
        coeff = binomial(n, i)
        @. B[:, i+1] = coeff * x^i * (1 - x)^(n - i)
    end
    return B
end

"""
    class_function(x; N1=0.5, N2=1.0)

CST class function for airfoils: C(x) = x^N1 * (1-x)^N2

Standard values N1=0.5, N2=1.0 give round leading edge and pointed trailing edge.
"""
function class_function(x::AbstractVector{T}; N1=0.5, N2=1.0) where T
    # Avoid numerical issues at endpoints
    x_safe = clamp.(x, eps(T), 1 - eps(T))
    return x_safe.^N1 .* (1 .- x_safe).^N2
end

"""
    normalize_airfoil(x::Vector, y::Vector)

Normalize airfoil coordinates to unit chord with LE at origin.

Returns normalized (x, y) and transformation parameters.
"""
function normalize_airfoil(x::Vector{T}, y::Vector{T}) where T
    # Find leading edge (minimum x)
    le_idx = argmin(x)
    x_le, y_le = x[le_idx], y[le_idx]

    # Find trailing edge (maximum x)
    te_idx = argmax(x)
    x_te, y_te = x[te_idx], y[te_idx]

    # Chord length and angle
    chord = sqrt((x_te - x_le)^2 + (y_te - y_le)^2)
    angle = atan(y_te - y_le, x_te - x_le)

    # Transform: translate to LE at origin, rotate, scale
    x_norm = similar(x)
    y_norm = similar(y)
    cos_a, sin_a = cos(angle), sin(angle)

    for i in eachindex(x)
        dx = x[i] - x_le
        dy = y[i] - y_le
        x_rot = cos_a * dx + sin_a * dy
        y_rot = -sin_a * dx + cos_a * dy
        x_norm[i] = x_rot / chord
        y_norm[i] = y_rot / chord
    end

    return x_norm, y_norm, (x_le=x_le, y_le=y_le, chord=chord, angle=angle)
end

"""
    leading_edge_basis(x::AbstractVector, n_weights::Int)

Leading-edge-modification basis used by AeroSandbox/NeuralFoil:
`x * (1 - x)^(n_weights + 0.5)`. Added identically to both surfaces.
"""
function leading_edge_basis(x::AbstractVector{T}, n_weights::Int) where T
    return x .* max.(1 .- x, zero(T)).^(n_weights + 0.5)
end

"""
    fit_kulfan_parameters(x::Vector, y::Vector, method::KulfanFitMethod)
    fit_kulfan_parameters(x::Vector, y::Vector; n_weights=8)

Fit Kulfan CST parameters to airfoil coordinates. The `method` selects the
strategy: [`LeastSquaresFit`](@ref) fits *through* the points (default), while
[`EnvelopeFit`](@ref) wraps tightly *around* them and is robust to noisy interior
points. The keyword form is shorthand for `LeastSquaresFit(; n_weights)`.

Input coordinates are assumed to be in Selig order (TE upper -> LE -> TE lower).

# Returns
- `KulfanParameters`: Fitted parameters
"""
function fit_kulfan_parameters(x::Vector{T}, y::Vector{T}; n_weights::Int=8) where T
    return fit_kulfan_parameters(x, y, LeastSquaresFit(; n_weights))
end

"""
    fit_kulfan_parameters(x, y, method::LeastSquaresFit)

Least-squares fit matching AeroSandbox's `get_kulfan_parameters`: both surfaces
share a single least-squares system with a shared leading-edge weight and a
trailing-edge thickness.
"""
function fit_kulfan_parameters(x::Vector{T}, y::Vector{T},
                               method::LeastSquaresFit) where T
    n_weights = method.n_weights
    x_norm, y_norm, _ = normalize_airfoil(x, y)
    xv = clamp.(x_norm, zero(T), one(T))

    le_idx = argmin(xv)
    is_upper = (1:length(xv)) .<= le_idx

    CS = class_function(xv) .* bernstein_basis(xv, n_weights - 1)
    le_col = leading_edge_basis(xv, n_weights)
    te_col = ifelse.(is_upper, xv ./ 2, .-xv ./ 2)

    A = hcat((.!is_upper) .* CS, is_upper .* CS, le_col, te_col)
    coeffs = A \ y_norm
    TE_thickness = coeffs[end]

    if TE_thickness < 0
        A = hcat((.!is_upper) .* CS, is_upper .* CS, le_col)
        coeffs = A \ y_norm
        TE_thickness = zero(T)
    end

    lower_weights = coeffs[1:n_weights]
    upper_weights = coeffs[n_weights+1:2n_weights]
    leading_edge_weight = coeffs[2n_weights+1]

    return KulfanParameters(upper_weights, lower_weights, leading_edge_weight,
                            TE_thickness)
end

"""
    fit_kulfan_parameters(x, y, method::EnvelopeFit)

Fit Kulfan CST parameters to an outer envelope of the point cloud `(x, y)`: each
surface hugs the points but a hard constraint keeps it at least `min_distance`
outside every one of them, so the fit encloses the geometry even with noisy interior
points (unlike [`LeastSquaresFit`](@ref), which fits *through* them). `min_distance`
also sets the minimum thickness where the surfaces nearly coincide, preventing
collapse to zero. LE and TE are pinned at the extreme points by the CST class function.

Solved as "minimise `tightness`·Σ(upper − lower)² + (arc length) subject to
`upper(xᵢ) ≥ yᵢ + min_distance` and `lower(xᵢ) ≤ yᵢ - min_distance`" by iteratively
reweighted active-set least-squares. Each iteration stacks four blocks: the thickness
objective scaled by `tightness`·`tightness_dist(xᵢ)` (a two-surface quantity, so it
targets 0 and never references the point `y` values — interior points cannot erode the
clearance); a penalty row per violated clearance constraint (`distance_penalty`);
arc-length rows over `n_grid` chordwise segments reweighted by inverse length; and a
Tikhonov term toward the least-squares seed (`regularization`). Only the ratio of
`tightness` to the unit arc-length weight matters, so the arc-length weight is fixed at 1.

# Returns
- `KulfanParameters`: fitted parameters
"""
function fit_kulfan_parameters(x::Vector{T}, y::Vector{T},
                               method::EnvelopeFit) where T
    n_weights = method.n_weights
    d = T(method.min_distance)
    x_norm, y_norm, _ = normalize_airfoil(x, y)
    xv = clamp.(x_norm, zero(T), one(T))
    m = length(xv)
    n_params = 2n_weights + 2

    surface_design(xs) = begin
        cs = class_function(xs) .* bernstein_basis(xs, n_weights - 1)
        leb = leading_edge_basis(xs, n_weights)
        half = xs ./ 2
        zb = zeros(T, length(xs), n_weights)
        upper = hcat(zb, cs, leb, half)
        lower = hcat(cs, zb, leb, .-half)
        return upper, lower
    end

    upper_pts, lower_pts = surface_design(xv)
    A = vcat(upper_pts, .-lower_pts)
    b = vcat(y_norm .+ d, d .- y_norm)

    sample_dist(dist, xs) = begin
        itp = linear_interpolation(range(zero(T), one(T), length(dist)),
                                   T.(dist); extrapolation_bc=Flat())
        [itp(x) for x in xs]
    end
    tight_dist = sample_dist(method.tightness_dist, xv)
    thickness_rows = sqrt.(T(method.tightness) .* tight_dist) .* (upper_pts .- lower_pts)
    thickness_rhs = zeros(T, m)

    grid = (one(T) .- cos.(range(zero(T), T(pi), method.n_grid))) ./ 2
    upper_grid, lower_grid = surface_design(grid)
    diff_upper = upper_grid[2:end, :] .- upper_grid[1:end-1, :]
    diff_lower = lower_grid[2:end, :] .- lower_grid[1:end-1, :]
    perimeter_rows = vcat(diff_upper, diff_lower)
    dx2 = repeat((grid[2:end] .- grid[1:end-1]) .^ 2, 2)

    reg_rows = sqrt(T(method.regularization)) .* Matrix{T}(I, n_params, n_params)
    p0 = fit_kulfan_parameters(x, y; n_weights)
    theta0 = vcat(p0.lower_weights, p0.upper_weights, p0.leading_edge_weight,
                  p0.TE_thickness)
    # A near-zero-thickness point set (e.g. an open single-membrane slice) makes the
    # least-squares seed blow up; fall back to a thin flat seed so the envelope
    # constraints (min_distance) shape a clean thin airfoil instead.
    (all(isfinite, theta0) && maximum(abs, theta0) ≤ 10) ||
        (theta0 = vcat(zeros(T, 2n_weights + 1), 2d))
    reg_rhs = sqrt(T(method.regularization)) .* theta0

    theta = copy(theta0)
    for _ in 1:method.n_iter
        distance_weight = sqrt.(T(method.distance_penalty) .* ((A * theta .- b) .< 0))
        segment_length = sqrt.(dx2 .+ (perimeter_rows * theta) .^ 2)
        arclength_weight = sqrt.(one(T) ./ segment_length)
        M = vcat(thickness_rows, distance_weight .* A,
                 arclength_weight .* perimeter_rows, reg_rows)
        rhs = vcat(thickness_rhs, distance_weight .* b,
                   zeros(T, length(arclength_weight)), reg_rhs)
        theta_new = M \ rhs
        converged = maximum(abs.(theta_new .- theta)) < T(method.tol)
        theta = theta_new
        converged && break
    end

    lower_weights = theta[1:n_weights]
    upper_weights = theta[n_weights+1:2n_weights]
    leading_edge_weight = theta[2n_weights+1]
    TE_thickness = max(theta[2n_weights+2], zero(T))

    return KulfanParameters(upper_weights, lower_weights, leading_edge_weight,
                            TE_thickness)
end

"""
    kulfan_to_coordinates(params::KulfanParameters; n_points=100)

Generate airfoil coordinates from Kulfan parameters.

Returns (x, y) coordinates in Selig format (TE upper -> LE -> TE lower).
"""
function kulfan_to_coordinates(params::KulfanParameters; n_points::Int=100)
    # Cosine spacing for better LE resolution
    theta = range(0, pi, n_points)
    x = @. (1 - cos(theta)) / 2

    n_upper = length(params.upper_weights)
    n_lower = length(params.lower_weights)

    C = class_function(x)
    B_upper = bernstein_basis(x, n_upper - 1)
    B_lower = bernstein_basis(x, n_lower - 1)

    # Shape functions
    S_upper = B_upper * params.upper_weights
    S_lower = B_lower * params.lower_weights

    # Leading-edge modification (shared, added to both surfaces)
    le = params.leading_edge_weight .* leading_edge_basis(x, n_upper)

    # y coordinates
    y_upper = C .* S_upper .+ le .+ x .* (params.TE_thickness / 2)
    y_lower = C .* S_lower .+ le .- x .* (params.TE_thickness / 2)

    # Combine in Selig format
    x_coords = vcat(reverse(x), x[2:end])
    y_coords = vcat(reverse(y_upper), y_lower[2:end])

    return x_coords, y_coords
end

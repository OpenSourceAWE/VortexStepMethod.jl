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
    EnvelopeFit(; min_distance=0.0, n_weights=8, tightness=1.0,
                perimeter_weight=1.0, penalty=1e4, n_grid=80, n_iter=30,
                regularization=1e-6, tol=1e-10)

Enclosing Kulfan fit: wraps the airfoil tightly *around* the point cloud,
keeping every point at least `min_distance` inside the surfaces. Robust to noisy
interior points (e.g. the internal structure of a ram-air kite slice) that pull
[`LeastSquaresFit`](@ref) inward. See [`fit_kulfan_parameters`](@ref) for the
method and algorithm.

Two knobs trade off the fit:
- `tightness` pulls the surfaces onto the outermost points (minimises enclosed
  area); raise it to follow concave detail more closely.
- `perimeter_weight` minimises the surfaces' arc length (a shrink-wrap); raise it
  for a smoother, shorter perimeter at the cost of some tightness.

# Fields
- `min_distance`: required clearance between each surface and the points
- `n_weights`: weights per surface (NeuralFoil uses 8)
- `tightness`: weight of the enclosed-area objective (hug the points)
- `perimeter_weight`: weight of the arc-length objective (smooth, short perimeter)
- `penalty`: weight of each active clearance-constraint row
- `n_grid`: number of chordwise stations used to measure perimeter
- `n_iter`, `tol`: active-set / reweighting iteration controls
- `regularization`: Tikhonov weight toward the least-squares fit
"""
@with_kw struct EnvelopeFit <: KulfanFitMethod
    min_distance::Float64 = 0.0
    n_weights::Int = 8
    tightness::Float64 = 1.0
    perimeter_weight::Float64 = 1.0
    penalty::Float64 = 1e4
    n_grid::Int = 80
    n_iter::Int = 30
    regularization::Float64 = 1e-6
    tol::Float64 = 1e-10
end

"""
    fit_clearance(method::KulfanFitMethod) -> Float64

Chordwise/normal clearance a fit method wraps around the points: `min_distance`
for [`EnvelopeFit`](@ref), `0` otherwise. Used to inset the stored raw points to
the same frame as the fitted airfoil.
"""
fit_clearance(method::EnvelopeFit) = method.min_distance
fit_clearance(::KulfanFitMethod) = 0.0

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
    inset_airfoil(x::Vector, y::Vector, clearance) -> (x_inset, y_inset)

Normalize coordinates to unit chord, then uniformly scale them into the box
`[clearance, 1 - clearance]`. This leaves `clearance` of chordwise room at each end
so an airfoil fitted on `[0, 1]` can wrap around the leading and trailing edges
(and close the trailing edge) instead of being pinned to the extreme points. With
`clearance == 0` this is plain normalization.
"""
function inset_airfoil(x::Vector{T}, y::Vector{T}, clearance::Real) where T
    x_norm, y_norm, _ = normalize_airfoil(x, y)
    scale = one(T) - 2 * T(clearance)
    return T(clearance) .+ scale .* x_norm, scale .* y_norm
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

Fit Kulfan CST parameters so the airfoil *encloses* the point cloud `(x, y)` as
tightly as possible while keeping every point at least `method.min_distance`
inside the surfaces.

Unlike [`LeastSquaresFit`](@ref), which fits *through* the points and is pulled
around by noisy interior points (e.g. the internal structure of a ram-air kite
slice), this fits a tight outer envelope: the upper surface stays above every
point and the lower surface below it, so only the outermost points are active and
interior points are ignored. For flat/thin sections where the surfaces nearly
coincide, `min_distance` sets the minimum clearance between each surface and the
points, preventing the envelope from collapsing to zero thickness. The points are
first inset into `[min_distance, 1 - min_distance]` (see [`inset_airfoil`](@ref))
so the fitted `[0, 1]` airfoil has chordwise room to wrap around the leading and
trailing edges and to close the trailing edge, rather than being pinned to the
extreme points.

It approximates the program "minimise `tightness`·(enclosed area) +
`perimeter_weight`·(surface arc length) subject to
`upper(xᵢ) ≥ yᵢ + min_distance` and `lower(xᵢ) ≤ yᵢ - min_distance`" by an
iteratively reweighted active-set least-squares. Each iteration solves one
stacked least-squares (QR via `\\`) of four blocks: a squared-thickness objective
(`tightness`), a penalty row per currently-violated clearance constraint
(`penalty`), arc-length rows over an `n_grid`-station chordwise grid whose
segments are reweighted by inverse length each iteration so the quadratic descends
the true perimeter (`perimeter_weight`), and a Tikhonov term toward the
least-squares fit (`regularization`). Because the objective is area and perimeter
rather than distance to each point, interior points exert no inward pull, the
perimeter term keeps the surfaces smooth without a non-representable uniform
offset, and the method degrades gracefully when `min_distance` is too large for
the section to satisfy.

# Returns
- `KulfanParameters`: fitted parameters
"""
function fit_kulfan_parameters(x::Vector{T}, y::Vector{T},
                               method::EnvelopeFit) where T
    n_weights = method.n_weights
    d = T(method.min_distance)
    x_inset, y_norm = inset_airfoil(x, y, d)
    xv = clamp.(x_inset, zero(T), one(T))
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

    thickness_rows = sqrt(T(method.tightness)) .* (upper_pts .- lower_pts)
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
    reg_rhs = sqrt(T(method.regularization)) .* theta0

    theta = copy(theta0)
    for _ in 1:method.n_iter
        clearance_weight = sqrt.(T(method.penalty) .* ((A * theta .- b) .< 0))
        segment_length = sqrt.(dx2 .+ (perimeter_rows * theta) .^ 2)
        perimeter_weight = sqrt.(T(method.perimeter_weight) ./ segment_length)
        M = vcat(thickness_rows, clearance_weight .* A,
                 perimeter_weight .* perimeter_rows, reg_rows)
        rhs = vcat(thickness_rhs, clearance_weight .* b,
                   zeros(T, length(perimeter_weight)), reg_rhs)
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

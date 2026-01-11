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
    split_airfoil_surfaces(x::Vector, y::Vector)

Split airfoil coordinates into upper and lower surfaces.

Assumes coordinates start from TE upper, go around LE, end at TE lower (Selig format),
or a closed contour. Returns (x_upper, y_upper, x_lower, y_lower) with x sorted 0->1.
"""
function split_airfoil_surfaces(x::Vector{T}, y::Vector{T}) where T
    n = length(x)

    # Find leading edge (minimum x)
    le_idx = argmin(x)

    # Split into upper (indices 1:le_idx) and lower (indices le_idx:end)
    # Upper surface goes from TE to LE (reverse to get LE to TE)
    if le_idx > 1
        upper_indices = le_idx:-1:1
    else
        upper_indices = [1]
    end

    # Lower surface goes from LE to TE
    if le_idx < n
        lower_indices = le_idx:n
    else
        lower_indices = [n]
    end

    x_upper = x[upper_indices]
    y_upper = y[upper_indices]
    x_lower = x[lower_indices]
    y_lower = y[lower_indices]

    # Ensure sorted by x (LE to TE)
    upper_order = sortperm(x_upper)
    lower_order = sortperm(x_lower)

    return (x_upper[upper_order], y_upper[upper_order],
            x_lower[lower_order], y_lower[lower_order])
end

"""
    fit_kulfan_parameters(x::Vector, y::Vector; n_weights=8, outlier_threshold=0.05)

Fit Kulfan CST parameters to airfoil coordinates using least squares.

Uses a two-pass approach: first fit to identify outliers (interior mesh points),
then refit with outliers removed.

# Arguments
- `x::Vector`: x-coordinates (normalized 0-1)
- `y::Vector`: y-coordinates
- `n_weights::Int=8`: Number of weights per surface (NeuralFoil uses 8)
- `outlier_threshold::Real=0.05`: Distance threshold for outlier removal

# Returns
- `KulfanParameters`: Fitted parameters
"""
function fit_kulfan_parameters(x::Vector{T}, y::Vector{T}; n_weights::Int=8,
                               outlier_threshold::Real=0.05) where T
    # Normalize coordinates
    x_norm, y_norm, transform = normalize_airfoil(x, y)

    # First pass: fit to all points
    params_initial = _fit_kulfan_single_pass(x_norm, y_norm, n_weights)

    # Compute residuals and filter outliers
    x_clean, y_clean = _filter_outliers(x_norm, y_norm, params_initial,
                                        outlier_threshold)

    # Second pass: fit to cleaned points
    if length(x_clean) < length(x_norm) * 0.5
        # Too many points filtered - use original
        @warn "Outlier filtering removed too many points, using original"
        return params_initial
    end

    return _fit_kulfan_single_pass(x_clean, y_clean, n_weights)
end

"""
    _fit_kulfan_single_pass(x_norm, y_norm, n_weights)

Internal: single pass Kulfan fitting.
"""
function _fit_kulfan_single_pass(x_norm::Vector{T}, y_norm::Vector{T},
                                  n_weights::Int) where T
    # Split into upper and lower surfaces
    x_upper, y_upper, x_lower, y_lower = split_airfoil_surfaces(x_norm, y_norm)

    # Fit each surface separately
    upper_weights, le_weight_upper, te_upper = fit_surface(
        x_upper, y_upper, n_weights, true
    )
    lower_weights, le_weight_lower, te_lower = fit_surface(
        x_lower, y_lower, n_weights, false
    )

    # Average LE weight (should be similar for both surfaces)
    leading_edge_weight = (le_weight_upper + le_weight_lower) / 2

    # TE thickness is difference between upper and lower at TE
    TE_thickness = te_upper - te_lower

    return KulfanParameters(upper_weights, lower_weights, leading_edge_weight,
                           max(0.0, TE_thickness))
end

"""
    filter_airfoil_outliers(x, y, params, threshold)

Filter points that are more than `threshold` away from the fitted Kulfan curve.
Returns filtered (x, y) coordinates.
"""
function filter_airfoil_outliers(x::Vector{T}, y::Vector{T}, params::KulfanParameters,
                                 threshold::Real) where T
    # Generate dense Kulfan curve for distance computation
    x_fit, y_fit = kulfan_to_coordinates(params; n_points=200)

    # For each input point, find minimum distance to fitted curve
    keep = Bool[]
    for i in eachindex(x)
        min_dist = Inf
        for j in eachindex(x_fit)
            dist = sqrt((x[i] - x_fit[j])^2 + (y[i] - y_fit[j])^2)
            min_dist = min(min_dist, dist)
        end
        push!(keep, min_dist < threshold)
    end

    return x[keep], y[keep]
end

# Internal alias
const _filter_outliers = filter_airfoil_outliers

"""
    fit_surface(x::Vector, y::Vector, n_weights::Int, is_upper::Bool)

Fit Kulfan parameters to a single surface using least squares.

The CST equation is: y(x) = C(x) * S(x) + x * TE_offset
where C(x) = x^0.5 * (1-x) is the class function
and S(x) = sum(w_i * B_i(x)) is the shape function

We solve for weights using least squares: A * [w; TE_offset] = y
"""
function fit_surface(x::Vector{T}, y::Vector{T}, n_weights::Int,
                     is_upper::Bool) where T
    n = length(x)

    # Filter out points too close to endpoints (causes numerical issues)
    valid = (x .> 1e-6) .& (x .< 1 - 1e-6)
    x_fit = x[valid]
    y_fit = y[valid]
    n_fit = length(x_fit)

    if n_fit < n_weights + 1
        @warn "Not enough points for fitting, using fewer weights"
        n_weights = max(1, n_fit - 1)
    end

    # Build design matrix
    # y = C(x) * (B * w) + x * TE_offset
    # y = [C(x) .* B, x] * [w; TE_offset]
    C = class_function(x_fit)
    B = bernstein_basis(x_fit, n_weights - 1)  # n_weights basis functions

    # Design matrix: [C .* B, x]
    A = hcat(C .* B, x_fit)

    # Solve least squares
    params = A \ y_fit

    weights = params[1:n_weights]
    te_offset = params[end]

    # Leading edge weight is related to the curvature at LE
    # For CST, the LE modification is typically the first weight
    le_weight = weights[1]

    return weights, le_weight, te_offset
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

    # y coordinates
    y_upper = C .* S_upper .+ x .* (params.TE_thickness / 2)
    y_lower = C .* S_lower .- x .* (params.TE_thickness / 2)

    # Combine in Selig format
    x_coords = vcat(reverse(x), x[2:end])
    y_coords = vcat(reverse(y_upper), y_lower[2:end])

    return x_coords, y_coords
end

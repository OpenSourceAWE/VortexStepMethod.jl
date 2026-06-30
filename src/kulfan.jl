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
    leading_edge_basis(x::AbstractVector, n_weights::Int)

Leading-edge-modification basis used by AeroSandbox/NeuralFoil:
`x * (1 - x)^(n_weights + 0.5)`. Added identically to both surfaces.
"""
function leading_edge_basis(x::AbstractVector{T}, n_weights::Int) where T
    return x .* max.(1 .- x, zero(T)).^(n_weights + 0.5)
end

"""
    fit_kulfan_parameters(x::Vector, y::Vector; n_weights=8)

Fit Kulfan CST parameters to airfoil coordinates using a single least-squares
system, matching AeroSandbox's `get_kulfan_parameters` (the parameterization
NeuralFoil was trained on).

Upper and lower surfaces are fit together with a shared leading-edge weight and
a trailing-edge thickness. Input coordinates are assumed to be in Selig order
(TE upper -> LE -> TE lower).

# Arguments
- `x::Vector`: x-coordinates
- `y::Vector`: y-coordinates
- `n_weights::Int=8`: Number of weights per surface (NeuralFoil uses 8)

# Returns
- `KulfanParameters`: Fitted parameters
"""
function fit_kulfan_parameters(x::Vector{T}, y::Vector{T}; n_weights::Int=8) where T
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

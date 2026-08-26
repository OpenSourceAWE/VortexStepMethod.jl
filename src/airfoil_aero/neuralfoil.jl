"""
NeuralFoil - Pure Julia implementation of the NeuralFoil neural network for airfoil
aerodynamics prediction.

Based on: https://github.com/peterdsharpe/NeuralFoil

This module loads the pre-trained neural network weights and provides functions to
evaluate airfoil aerodynamics from Kulfan CST parameters.
"""

using NPZ

# Path to NeuralFoil weights (relative to this package)
const NEURALFOIL_WEIGHTS_DIR = joinpath(@__DIR__, "data")

# Global cache for loaded weights
const _NN_CACHE = Dict{String, Any}()

"""
    NeuralFoilModel

Holds the neural network weights and configuration.
"""
struct NeuralFoilModel
    weights::Vector{Matrix{Float32}}
    biases::Vector{Vector{Float32}}
    layer_indices::Vector{Int}
    mean_inputs::Vector{Float32}
    inv_cov_inputs::Matrix{Float64}
    n_inputs::Int
end

"""
    swish(x)

SiLU/Swish activation function: x * sigmoid(x) = x / (1 + exp(-x))
"""
function swish(x::T) where T
    return x / (one(T) + exp(-x))
end

function swish(x::AbstractArray)
    return x ./ (one(eltype(x)) .+ exp.(-x))
end

"""
    sigmoid(x)

Sigmoid activation function.
"""
function sigmoid(x::T) where T
    x_clipped = clamp(x, T(-500), T(500))
    return one(T) / (one(T) + exp(-x_clipped))
end

"""
    load_neuralfoil_model(model_size::String="xlarge"; weights_dir=nothing)

Load NeuralFoil neural network weights from .npz files.

# Arguments
- `model_size`: One of "xxsmall", "xsmall", "small", "medium", "large",
                "xlarge", "xxlarge", "xxxlarge"
- `weights_dir`: Directory containing the .npz files (defaults to package data dir)

# Returns
- `NeuralFoilModel`: Loaded model ready for evaluation
"""
function load_neuralfoil_model(model_size::String="xlarge"; weights_dir=nothing)
    # Check cache
    cache_key = "$(model_size)_$(weights_dir)"
    if haskey(_NN_CACHE, cache_key)
        return _NN_CACHE[cache_key]
    end

    if weights_dir === nothing
        weights_dir = NEURALFOIL_WEIGHTS_DIR
    end

    # Load neural network weights
    nn_path = joinpath(weights_dir, "nn-$(model_size).npz")
    if !isfile(nn_path)
        error("NeuralFoil weights not found at $nn_path. " *
              "Copy weights from NeuralFoil package.")
    end

    nn_data = npzread(nn_path)

    # Parse layer indices from keys like "net.0.weight", "net.2.weight", etc.
    layer_indices = Int[]
    for key in keys(nn_data)
        m = match(r"net\.(\d+)\.weight", key)
        if m !== nothing
            push!(layer_indices, parse(Int, m.captures[1]))
        end
    end
    sort!(unique!(layer_indices))

    # Load weights and biases in order
    weights = Matrix{Float32}[]
    biases = Vector{Float32}[]
    for idx in layer_indices
        push!(weights, Float32.(nn_data["net.$(idx).weight"]))
        push!(biases, Float32.(vec(nn_data["net.$(idx).bias"])))
    end

    # Load input distribution statistics
    dist_path = joinpath(weights_dir, "scaled_input_distribution.npz")
    if !isfile(dist_path)
        error("Input distribution file not found at $dist_path")
    end

    dist_data = npzread(dist_path)
    mean_inputs = Float32.(vec(dist_data["mean_inputs_scaled"]))
    inv_cov_inputs = Float64.(dist_data["inv_cov_inputs_scaled"])

    model = NeuralFoilModel(weights, biases, layer_indices, mean_inputs,
                           inv_cov_inputs, length(mean_inputs))

    _NN_CACHE[cache_key] = model
    return model
end

"""
    nn_forward(x::AbstractMatrix, model::NeuralFoilModel)

Forward pass through the neural network.

# Arguments
- `x`: Input matrix of shape (n_inputs, n_cases)
- `model`: Loaded NeuralFoil model

# Returns
- Output matrix of shape (n_outputs, n_cases)
"""
function nn_forward(x::AbstractMatrix{T}, model::NeuralFoilModel) where T
    n_layers = length(model.weights)

    for i in 1:n_layers
        x = model.weights[i] * x .+ model.biases[i]
        if i < n_layers  # No activation on last layer
            x = swish(x)
        end
    end

    return x
end

"""
    squared_mahalanobis_distance(x::AbstractMatrix, model::NeuralFoilModel)

Compute squared Mahalanobis distance from training distribution.

This is used to penalize predictions far from the training data.
"""
function squared_mahalanobis_distance(x::AbstractMatrix, model::NeuralFoilModel)
    n_cases = size(x, 2)
    result = zeros(n_cases)

    for i in 1:n_cases
        x_minus_mean = x[:, i] .- model.mean_inputs
        result[i] = dot(x_minus_mean, model.inv_cov_inputs * x_minus_mean)
    end

    return result
end

"""
    fill_case_input!(x, case, params::KulfanParameters, alpha_deg, Re, n_crit,
                     xtr_upper, xtr_lower)

Write one NeuralFoil input column: rows 1–18 the Kulfan shape, rows 19–25 the flow
condition (`alpha_deg` in degrees). The single point where the network's input layout
is defined, shared by every `prepare_inputs` method.
"""
function fill_case_input!(x::AbstractMatrix, case::Int, params::KulfanParameters,
                          alpha_deg, Re, n_crit, xtr_upper, xtr_lower)
    for i in 1:8
        x[i, case] = params.upper_weights[i]
        x[8 + i, case] = params.lower_weights[i]
    end
    x[17, case] = params.leading_edge_weight
    x[18, case] = params.TE_thickness * 50  # Scale factor from NeuralFoil
    x[19, case] = sind(2 * alpha_deg)
    x[20, case] = cosd(alpha_deg)
    x[21, case] = 1 - cosd(alpha_deg)^2
    x[22, case] = (log(Re) - 12.5) / 3.5
    x[23, case] = (n_crit - 9) / 4.5
    x[24, case] = xtr_upper
    x[25, case] = xtr_lower
    return nothing
end

"""
    prepare_inputs(params::KulfanParameters, alpha, Re;
                   n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)
    prepare_inputs(params::AbstractVector{KulfanParameters}, alpha, Re; kwargs...)

Prepare neural network inputs from Kulfan parameters and flow conditions. The vector
form takes one shape per case, so a whole wing's panels go through the network in a
single forward pass.

# Arguments
- `params`: Kulfan CST parameters
- `alpha`: Angle of attack in degrees (scalar or vector)
- `Re`: Reynolds number (scalar or vector)
- `n_crit`: Critical amplification factor (default 9.0)
- `xtr_upper`: Forced transition location on upper surface (0-1)
- `xtr_lower`: Forced transition location on lower surface (0-1)

# Returns
- Input matrix of shape (25, n_cases)
"""
function prepare_inputs(params::KulfanParameters, alpha, Re;
                        n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)
    alpha_vec = alpha isa Number ? [alpha] : collect(alpha)
    n_cases = length(alpha_vec)
    return prepare_inputs(fill(params, n_cases), alpha_vec, Re;
                          n_crit, xtr_upper, xtr_lower)
end

function prepare_inputs(params::AbstractVector{KulfanParameters}, alpha, Re;
                        n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)
    alpha_vec = alpha isa Number ? fill(alpha, length(params)) : collect(alpha)
    n_cases = length(alpha_vec)
    length(params) == n_cases || throw(ArgumentError(
        "prepare_inputs: $(length(params)) shapes for $n_cases angles of attack."))
    per_case(v) = v isa Number ? fill(v, n_cases) : collect(v)
    Re_vec, n_crit_vec = per_case(Re), per_case(n_crit)
    xtr_upper_vec, xtr_lower_vec = per_case(xtr_upper), per_case(xtr_lower)

    x = zeros(Float32, 25, n_cases)
    for i in 1:n_cases
        fill_case_input!(x, i, params[i], alpha_vec[i], Re_vec[i], n_crit_vec[i],
                         xtr_upper_vec[i], xtr_lower_vec[i])
    end
    return x
end

"""
    flip_inputs(x::AbstractMatrix)

Flip inputs for symmetry embedding (swap upper/lower, negate alpha).
"""
function flip_inputs(x::AbstractMatrix{T}) where T
    x_flip = copy(x)

    # Swap upper and lower weights with sign flip
    x_flip[1:8, :] .= -x[9:16, :]   # upper <- -lower
    x_flip[9:16, :] .= -x[1:8, :]   # lower <- -upper

    # Flip LE weight
    x_flip[17, :] .= -x[17, :]

    # Flip sin(2*alpha) (index 19)
    x_flip[19, :] .= -x[19, :]

    # Swap transition locations
    x_flip[24, :] .= x[25, :]  # xtr_upper <- xtr_lower
    x_flip[25, :] .= x[24, :]  # xtr_lower <- xtr_upper

    return x_flip
end

"""
    flip_outputs(y::AbstractMatrix)

Flip outputs back after evaluating with flipped inputs.
"""
function flip_outputs(y::AbstractMatrix{T}) where T
    y_flip = copy(y)
    N = (size(y, 1) - 6) ÷ 6

    # CL (row 2) and CM (row 4) flip sign
    y_flip[2, :] .= -y[2, :]
    y_flip[4, :] .= -y[4, :]

    # Transition locations swap (rows 5, 6)
    y_flip[5, :] .= y[6, :]
    y_flip[6, :] .= y[5, :]

    if N > 0
        # Swap upper/lower theta+H blocks (rows 7:6+2N ↔ 7+3N:6+5N)
        y_flip[7:(6+2N), :] .= y[(7+3N):(6+5N), :]
        y_flip[(7+3N):(6+5N), :] .= y[7:(6+2N), :]
        # Swap upper/lower ue/vinf blocks with sign flip (velocity mirrors)
        y_flip[(7+2N):(6+3N), :] .= -1 .* y[(7+5N):(6+6N), :]
        y_flip[(7+5N):(6+6N), :] .= -1 .* y[(7+2N):(6+3N), :]
    end

    return y_flip
end

"""
    compute_optimal_x_points(n) -> Vector{Float64}

NeuralFoil's boundary-layer station x/c: midpoints of a uniform [0,1] grid.
"""
function compute_optimal_x_points(n)
    s = range(0, 1, n + 1)
    return [(s[i] + s[i+1]) / 2 for i in 1:n]
end

"""
    NeuralFoilResult

Results from NeuralFoil aerodynamic analysis.
"""
struct NeuralFoilResult
    alpha::Vector{Float64}
    CL::Vector{Float64}
    CD::Vector{Float64}
    CM::Vector{Float64}
    analysis_confidence::Vector{Float64}
end

"""
    neuralfoil_aero(params::KulfanParameters, alpha, Re;
                    model_size="xlarge", weights_dir=nothing, kwargs...)

Compute aerodynamic coefficients using NeuralFoil.

# Arguments
- `params`: Kulfan CST parameters
- `alpha`: Angle of attack in degrees (scalar or vector)
- `Re`: Reynolds number
- `model_size`: Neural network size (default "xlarge")
- `weights_dir`: Directory containing weight files

# Keyword Arguments
- `n_crit`: Critical amplification factor (default 9.0)
- `xtr_upper`: Upper surface transition location (default 1.0)
- `xtr_lower`: Lower surface transition location (default 1.0)

# Returns
- `NeuralFoilResult`: CL, CD, CM and analysis confidence
"""
function neuralfoil_aero(params::KulfanParameters, alpha, Re;
                         model_size::String="xlarge", weights_dir=nothing,
                         n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)
    y = neuralfoil_fused_output(params, alpha, Re; model_size, weights_dir,
                                n_crit, xtr_upper, xtr_lower)
    cl, cd, cm, confidence = decode_coefficients(y)
    alpha_vec = alpha isa Number ? [Float64(alpha)] : Float64.(collect(alpha))
    return NeuralFoilResult(alpha_vec, cl, cd, cm, confidence)
end

"""
    decode_coefficients(y) -> (cl, cd, cm, confidence)

Turn a fused network output matrix into the integrated coefficients per case, undoing
NeuralFoil's output scaling. The single place that scaling is written down.
"""
function decode_coefficients(y::AbstractMatrix)
    return (Vector{Float64}(y[2, :] ./ 2),
            Vector{Float64}(clamp.(exp.((y[3, :] .- 2) .* 2), 0.0, 1.0)),
            Vector{Float64}(y[4, :] ./ 20),
            Vector{Float64}(sigmoid.(y[1, :])))
end

"""
    neuralfoil_fused_output(params, alpha, Re; model_size, weights_dir,
                            n_crit, xtr_upper, xtr_lower) -> Matrix

Forward pass with the symmetry embedding: evaluate the network, evaluate the
top/bottom-flipped case, flip its outputs back, and average. Returns the fused
output matrix (`n_outputs × n_cases`), with the Mahalanobis penalty already
applied to the confidence logit (row 1).
"""
function neuralfoil_fused_output(params, alpha, Re;
                                 model_size::String="xlarge", weights_dir=nothing,
                                 n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)
    x = prepare_inputs(params, alpha, Re; n_crit, xtr_upper, xtr_lower)
    return fused_output(x, load_neuralfoil_model(model_size; weights_dir))
end

"""
    fused_output(x, model) -> Matrix

Symmetry-fused forward pass over a prepared input matrix, see
[`neuralfoil_fused_output`](@ref). Takes the inputs already built so a caller that
assembles its own batch does not go back through [`prepare_inputs`](@ref).
"""
function fused_output(x::AbstractMatrix, model::NeuralFoilModel)
    y = nn_forward(x, model)
    y[1, :] .-= squared_mahalanobis_distance(x, model) ./ (2 * model.n_inputs)

    x_flip = flip_inputs(x)
    y_flip = nn_forward(x_flip, model)
    y_flip[1, :] .-= squared_mahalanobis_distance(x_flip, model) ./ (2 * model.n_inputs)

    return (y .+ flip_outputs(y_flip)) ./ 2
end

"""
    neuralfoil_section(params, alpha, Re; kwargs...) -> NamedTuple

Full NeuralFoil evaluation returning integrated coefficients and the surface
pressure distribution reconstructed from the predicted edge-velocity ratios
(`Cp = 1 - (ue/vinf)^2`) at NeuralFoil's `N` fixed station x/c.

Returns `(; alpha, cl, cd, cm, confidence, x, cp_upper, cp_lower, ue_upper,
ue_lower)`, with `x` of length `N` and the four matrices sized `N × n_alpha`. A
reconstruction should interpolate `ue` and square afterwards: it is linear in arc
length through a stagnation point, where `Cp` is quadratic.
"""
function neuralfoil_section(params::KulfanParameters, alpha, Re;
                            model_size::String="large", weights_dir=nothing,
                            n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)
    y = neuralfoil_fused_output(params, alpha, Re; model_size, weights_dir,
                                n_crit, xtr_upper, xtr_lower)
    alpha_vec = alpha isa Number ? [Float64(alpha)] : Float64.(collect(alpha))
    station_x, ue_upper, ue_lower = decode_surface_velocity(y)
    return (; alpha = alpha_vec,
            cl = Vector{Float64}(y[2, :] ./ 2),
            cd = Vector{Float64}(clamp.(exp.((y[3, :] .- 2) .* 2), 0.0, 1.0)),
            cm = Vector{Float64}(y[4, :] ./ 20),
            confidence = Vector{Float64}(sigmoid.(y[1, :])),
            x = station_x,
            cp_upper = Matrix{Float64}(1 .- ue_upper .^ 2),
            cp_lower = Matrix{Float64}(1 .- ue_lower .^ 2),
            ue_upper, ue_lower)
end

"""
    decode_surface_velocity(y) -> (station_x, ue_upper, ue_lower)

The edge-velocity ratios a fused output matrix carries, at NeuralFoil's own `N`
fixed stations. Each `ue` matrix is `N × n_cases`. Interpolate these rather than
`Cp = 1 - ue²`: `ue` is linear in arc length through a stagnation point, where
`Cp` is quadratic.
"""
function decode_surface_velocity(y::AbstractMatrix)
    N = (size(y, 1) - 6) ÷ 6
    return (compute_optimal_x_points(N),
            Matrix{Float64}(y[(7 + 2N):(6 + 3N), :]),
            Matrix{Float64}(y[(7 + 5N):(6 + 6N), :]))
end

"""
    neuralfoil_aero(x::Vector, y::Vector, alpha, Re; kwargs...)

Convenience function that fits Kulfan parameters from coordinates first.
"""
function neuralfoil_aero(x::Vector, y::Vector, alpha, Re; kwargs...)
    params = fit_kulfan_parameters(x, y)
    return neuralfoil_aero(params, alpha, Re; kwargs...)
end

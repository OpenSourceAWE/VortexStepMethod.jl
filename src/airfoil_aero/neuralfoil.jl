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
    layer_buffers(model, n_cases) -> Vector{Matrix{Float32}}

One activation matrix per network layer, `n_cases` wide: the storage a forward pass
writes through. [`NeuralFoilWorkspace`](@ref) holds two sets of them, one per symmetry.
"""
layer_buffers(model::NeuralFoilModel, n_cases::Int) =
    [zeros(Float32, size(w, 1), n_cases) for w in model.weights]

"""
    NeuralFoilWorkspace(model, n_cases)

Scratch a symmetry-fused forward pass runs inside, so a caller evaluating a batch of the
same width every solve allocates nothing at all. Sized once against `model` for at most
`n_cases` cases; [`fused_output!`](@ref) takes any batch up to that width, which is what
lets one workspace serve both a sampled polar batch and the narrower pressure batch.
"""
struct NeuralFoilWorkspace
    "Layer activations of the direct pass."
    direct::Vector{Matrix{Float32}}
    "Layer activations of the top/bottom-flipped pass."
    flipped::Vector{Matrix{Float32}}
    "The flipped input batch."
    inputs_flip::Matrix{Float32}
    "One input case less the training mean."
    centered::Vector{Float64}
    "NeuralFoil's boundary-layer station x/c, the axis a surface reconstruction reads."
    stations::Vector{Float64}
end

function NeuralFoilWorkspace(model::NeuralFoilModel, n_cases::Int)
    n_cases >= 1 || throw(ArgumentError(
        "A NeuralFoilWorkspace needs at least one case; got $n_cases."))
    n_outputs = size(last(model.weights), 1)
    return NeuralFoilWorkspace(layer_buffers(model, n_cases),
                               layer_buffers(model, n_cases),
                               zeros(Float32, model.n_inputs, n_cases),
                               zeros(model.n_inputs),
                               compute_optimal_x_points((n_outputs - 6) ÷ 6))
end

"How many cases a workspace was sized for."
case_capacity(work::NeuralFoilWorkspace) = size(work.inputs_flip, 2)

"""
    nn_forward(x::AbstractMatrix, model::NeuralFoilModel)

Forward pass through the neural network.

# Arguments
- `x`: Input matrix of shape (n_inputs, n_cases)
- `model`: Loaded NeuralFoil model

# Returns
- Output matrix of shape (n_outputs, n_cases)
"""
nn_forward(x::AbstractMatrix, model::NeuralFoilModel) =
    Matrix(nn_forward!(layer_buffers(model, size(x, 2)), x, model))

"""
    nn_forward!(layers, x, model, n_cases=size(x, 2)) -> AbstractMatrix

Forward pass writing each layer's activations into `layers`, returning a view of the last
one over the `n_cases` columns used. The matrix products go through BLAS in place, so a
pass over a batch the buffers were sized for allocates nothing.
"""
function nn_forward!(layers::Vector{Matrix{Float32}}, x::AbstractMatrix,
                     model::NeuralFoilModel, n_cases::Int=size(x, 2))
    n_layers = length(model.weights)
    out = view(layers[1], :, 1:n_cases)
    mul!(out, model.weights[1], view(x, :, 1:n_cases))
    add_bias!(out, model.biases[1], n_layers > 1)
    for i in 2:n_layers
        out = view(layers[i], :, 1:n_cases)
        mul!(out, model.weights[i], view(layers[i - 1], :, 1:n_cases))
        add_bias!(out, model.biases[i], i < n_layers)
    end
    return view(layers[n_layers], :, 1:n_cases)
end

"""
    add_bias!(out, bias, activate) -> out

Add a layer's bias to its activations in place, passing them through [`swish`](@ref)
unless this is the output layer, which carries no activation.
"""
function add_bias!(out::AbstractMatrix, bias::AbstractVector, activate::Bool)
    @inbounds for case in axes(out, 2), k in eachindex(bias)
        value = out[k, case] + bias[k]
        out[k, case] = activate ? swish(value) : value
    end
    return out
end

"""
    squared_mahalanobis_distance(x::AbstractMatrix, model::NeuralFoilModel)

Compute squared Mahalanobis distance from training distribution, one per case.

This is used to penalize predictions far from the training data.
"""
function squared_mahalanobis_distance(x::AbstractMatrix, model::NeuralFoilModel)
    centered = zeros(model.n_inputs)
    return [mahalanobis_case(x, case, model, centered) for case in axes(x, 2)]
end

"""
    mahalanobis_case(x, case, model, centered) -> Float64

One case's squared Mahalanobis distance, taking `centered` as the scratch the centred
input column is written into rather than allocating one per case.
"""
function mahalanobis_case(x::AbstractMatrix, case::Int, model::NeuralFoilModel,
                          centered::Vector{Float64})
    inv_cov = model.inv_cov_inputs
    @inbounds for i in eachindex(centered)
        centered[i] = x[i, case] - model.mean_inputs[i]
    end
    distance = 0.0
    @inbounds for j in eachindex(centered)
        weighted = 0.0
        for i in eachindex(centered)
            weighted += inv_cov[i, j] * centered[i]
        end
        distance += weighted * centered[j]
    end
    return distance
end

"""
    penalize_confidence!(y, x, model, centered) -> y

Subtract every case's Mahalanobis penalty from the confidence logit, row 1 of the network
output. This is what makes a shape far from the training distribution report a low
confidence, and it is applied to each symmetry before the two are fused.
"""
function penalize_confidence!(y::AbstractMatrix, x::AbstractMatrix,
                              model::NeuralFoilModel, centered::Vector{Float64})
    @inbounds for case in axes(y, 2)
        y[1, case] -= mahalanobis_case(x, case, model, centered) / (2 * model.n_inputs)
    end
    return y
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
flip_inputs(x::AbstractMatrix) = flip_inputs!(similar(x), x)

"""
    flip_inputs!(x_flip, x) -> x_flip

[`flip_inputs`](@ref) written into storage the caller owns: the upper and lower weights
swap with a sign change, the leading-edge weight and `sin(2·alpha)` negate, and the two
transition locations swap. Everything else carries over.
"""
function flip_inputs!(x_flip::AbstractMatrix, x::AbstractMatrix)
    x_flip .= x
    @inbounds for case in axes(x, 2)
        for i in 1:8
            x_flip[i, case] = -x[8 + i, case]
            x_flip[8 + i, case] = -x[i, case]
        end
        x_flip[17, case] = -x[17, case]
        x_flip[19, case] = -x[19, case]
        x_flip[24, case] = x[25, case]
        x_flip[25, case] = x[24, case]
    end
    return x_flip
end

"""
    flip_outputs(y::AbstractMatrix)

Flip outputs back after evaluating with flipped inputs.
"""
function flip_outputs(y::AbstractMatrix)
    n_stations = (size(y, 1) - 6) ÷ 6
    y_flip = similar(y)
    for row in axes(y, 1)
        source, sign = flipped_row(row, n_stations)
        @views y_flip[row, :] .= sign .* y[source, :]
    end
    return y_flip
end

"""
    flipped_row(row, n_stations) -> (source, sign)

Which row of a network output a top/bottom-flipped case reads `row` from, and with which
sign. The single place the mirror symmetry of the output layout is written down, shared
by [`flip_outputs`](@ref) and [`fuse_flipped!`](@ref).

`CL` and `CM` change sign, the two transition locations swap, and of the six
`n_stations`-long boundary-layer blocks the upper and lower halves swap — the momentum
thickness and shape factor as they are, the edge velocity mirrored, so it changes sign.
"""
function flipped_row(row::Int, n_stations::Int)
    (row == 2 || row == 4) && return (row, -1)
    row == 5 && return (6, 1)
    row == 6 && return (5, 1)
    (row <= 6 || n_stations == 0) && return (row, 1)
    block = (row - 7) ÷ n_stations
    block in (0, 1) && return (row + 3n_stations, 1)
    block == 2 && return (row + 3n_stations, -1)
    block in (3, 4) && return (row - 3n_stations, 1)
    block == 5 && return (row - 3n_stations, -1)
    return (row, 1)
end

"""
    fuse_flipped!(y, y_flip) -> y

Average a direct output with its flipped counterpart in place, undoing the flip on the
way ([`flipped_row`](@ref)). Reading `y_flip` while writing `y` is what lets the fusion
land in storage that is already there instead of a third matrix.
"""
function fuse_flipped!(y::AbstractMatrix, y_flip::AbstractMatrix)
    n_stations = (size(y, 1) - 6) ÷ 6
    @inbounds for row in axes(y, 1)
        source, sign = flipped_row(row, n_stations)
        for case in axes(y, 2)
            y[row, case] = (y[row, case] + sign * y_flip[source, case]) / 2
        end
    end
    return y
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
    n_cases = size(y, 2)
    return decode_coefficients!(zeros(n_cases), zeros(n_cases), zeros(n_cases),
                                zeros(n_cases), y)
end

"""
    decode_coefficients!(cl, cd, cm, confidence, y) -> (cl, cd, cm, confidence)

[`decode_coefficients`](@ref) into vectors the caller owns, so a live polar reads a
refresh out of the network without allocating four vectors for it.
"""
function decode_coefficients!(cl, cd, cm, confidence, y::AbstractMatrix)
    @inbounds for case in axes(y, 2)
        cl[case] = y[2, case] / 2
        cd[case] = clamp(exp((y[3, case] - 2) * 2), 0.0, 1.0)
        cm[case] = y[4, case] / 20
        confidence[case] = sigmoid(y[1, case])
    end
    return (cl, cd, cm, confidence)
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
fused_output(x::AbstractMatrix, model::NeuralFoilModel) =
    Matrix(fused_output!(NeuralFoilWorkspace(model, size(x, 2)), x, model))

"""
    fused_output!(work, x, model) -> AbstractMatrix

[`fused_output`](@ref) run entirely inside `work`, allocating nothing for a batch the
workspace has room for. The result is a view of the workspace's own storage and stays
valid until the next pass through it.
"""
function fused_output!(work::NeuralFoilWorkspace, x::AbstractMatrix,
                       model::NeuralFoilModel)
    n_cases = size(x, 2)
    n_cases <= case_capacity(work) || throw(ArgumentError(
        "A NeuralFoilWorkspace sized for $(case_capacity(work)) cases was handed " *
        "$n_cases."))
    x_flip = view(work.inputs_flip, :, 1:n_cases)
    flip_inputs!(x_flip, x)
    y = nn_forward!(work.direct, x, model, n_cases)
    penalize_confidence!(y, x, model, work.centered)
    y_flip = nn_forward!(work.flipped, x_flip, model, n_cases)
    penalize_confidence!(y_flip, x_flip, model, work.centered)
    return fuse_flipped!(y, y_flip)
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
    n_stations, upper, lower = surface_velocity_rows(y)
    return (compute_optimal_x_points(n_stations),
            Matrix{Float64}(y[upper, :]), Matrix{Float64}(y[lower, :]))
end

"""
    decode_surface_velocity!(ue_upper, ue_lower, y) -> (ue_upper, ue_lower)

[`decode_surface_velocity`](@ref) into matrices the caller owns. The station axis is
fixed by the network and is carried by [`NeuralFoilWorkspace`](@ref) instead, so a live
pressure refresh reads a pass out without allocating anything for it.
"""
function decode_surface_velocity!(ue_upper::AbstractMatrix, ue_lower::AbstractMatrix,
                                  y::AbstractMatrix)
    _, upper, lower = surface_velocity_rows(y)
    @inbounds for case in axes(y, 2)
        for (k, row) in enumerate(upper)
            ue_upper[k, case] = y[row, case]
        end
        for (k, row) in enumerate(lower)
            ue_lower[k, case] = y[row, case]
        end
    end
    return (ue_upper, ue_lower)
end

"""
    surface_velocity_rows(y) -> (n_stations, upper, lower)

The number of boundary-layer stations a fused output matrix carries and the row ranges
its upper and lower edge-velocity ratios sit in. Viewing those rows is how a caller that
already holds the output reads the velocities without copying them out
([`decode_surface_velocity`](@ref) is the copying form).
"""
function surface_velocity_rows(y::AbstractMatrix)
    n_stations = (size(y, 1) - 6) ÷ 6
    return (n_stations, (7 + 2n_stations):(6 + 3n_stations),
            (7 + 5n_stations):(6 + 6n_stations))
end

"""
    neuralfoil_aero(x::Vector, y::Vector, alpha, Re; kwargs...)

Convenience function that fits Kulfan parameters from coordinates first.
"""
function neuralfoil_aero(x::Vector, y::Vector, alpha, Re; kwargs...)
    params = fit_kulfan_parameters(x, y)
    return neuralfoil_aero(params, alpha, Re; kwargs...)
end

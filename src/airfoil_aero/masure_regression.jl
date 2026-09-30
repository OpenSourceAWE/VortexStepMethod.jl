"Airfoil parameters the masure regression takes, in its input order (alpha follows)."
const MASURE_PARAMETERS = ("t", "eta", "kappa", "delta", "lambda", "phi")

"Reynolds numbers a masure regression model exists for, with their file suffix."
const MASURE_REYNOLDS = Dict(1.0e6 => "1e6", 5.0e6 => "5e6", 2.0e7 => "2e7")

"""
    ExtraTreesForest

The regression trees predicting one output, nodes of all trees concatenated. A leaf
has `left == 0`.
"""
struct ExtraTreesForest
    "first node of each tree"
    roots::Vector{Int32}
    "child taken when the feature is at or below the threshold"
    left::Vector{Int32}
    "child taken when the feature is above the threshold"
    right::Vector{Int32}
    "input feature a node splits on"
    feature::Vector{Int32}
    "split threshold, in scaled input units"
    threshold::Vector{Float64}
    "prediction of a leaf"
    value::Vector{Float64}
end

"""
    MasureModel

A masure regression model for one Reynolds number: the input standardisation and one
[`ExtraTreesForest`](@ref) per output, in the order CD, CL, CM.
"""
struct MasureModel
    input_mean::Vector{Float64}
    input_scale::Vector{Float64}
    forests::Vector{ExtraTreesForest}
end

const MASURE_CACHE = Dict{String, MasureModel}()

"Index array `key` of an exported model, shifted from 0-based to 1-based."
one_based(data::AbstractDict, key::String) = data[key] .+ Int32(1)

function ExtraTreesForest(data::AbstractDict, prefix::String)
    return ExtraTreesForest(one_based(data, "$(prefix)_roots"),
                            one_based(data, "$(prefix)_left"),
                            one_based(data, "$(prefix)_right"),
                            one_based(data, "$(prefix)_feature"),
                            data["$(prefix)_threshold"], data["$(prefix)_value"])
end

"""
    load_masure_model(Re, ml_models_dir) -> MasureModel

Load the masure regression model for Reynolds number `Re` (1e6, 5e6 or 2e7) from
`ml_models_dir/ET_re<Re>.npz`, as written by `scripts/export_masure_models.py`.
"""
function load_masure_model(Re::Real, ml_models_dir::AbstractString)
    haskey(MASURE_REYNOLDS, Re) || error("No masure regression model for Re = $Re; " *
        "available: $(join(sort!(collect(keys(MASURE_REYNOLDS))), ", ")).")
    path = abspath(joinpath(ml_models_dir, "ET_re$(MASURE_REYNOLDS[Re]).npz"))
    haskey(MASURE_CACHE, path) && return MASURE_CACHE[path]
    isfile(path) || error("Masure regression model not found at $path. Download " *
        "the models from https://doi.org/10.5281/zenodo.16925758 and convert them " *
        "with scripts/export_masure_models.py.")
    data = npzread(path)
    return MASURE_CACHE[path] = MasureModel(data["input_mean"], data["input_scale"],
        [ExtraTreesForest(data, "output$k") for k in 0:2])
end

"""
    forest_predict(forest::ExtraTreesForest, x) -> Float64

Mean of the leaf values the scaled input `x` reaches in each tree.
"""
function forest_predict(forest::ExtraTreesForest, x::AbstractVector{Float32})
    total = 0.0
    for root in forest.roots
        node = root
        while forest.left[node] != 0
            node = x[forest.feature[node]] <= forest.threshold[node] ?
                forest.left[node] : forest.right[node]
        end
        total += forest.value[node]
    end
    return total / length(forest.roots)
end

"""
    masure_aero(model::MasureModel, params, alpha) -> (cl, cd, cm)

Lift, drag and moment coefficients of the LEI airfoil described by `params` (keys
[`MASURE_PARAMETERS`](@ref), `delta` in degrees) at each angle of attack in `alpha`
[deg].
"""
function masure_aero(model::MasureModel, params::AbstractDict, alpha::AbstractVector)
    x = [scale_input(model, params[name], j) for (j, name) in enumerate(MASURE_PARAMETERS)]
    push!(x, 0.0f0)
    coefficients = zeros(length(alpha), length(model.forests))
    for (i, angle) in enumerate(alpha)
        x[end] = scale_input(model, angle, length(x))
        for (k, forest) in enumerate(model.forests)
            coefficients[i, k] = forest_predict(forest, x)
        end
    end
    return coefficients[:, 2], coefficients[:, 1], coefficients[:, 3]
end

"""
    scale_input(model, value, j) -> Float32

Input `j` standardised as the model's scaler does, rounded to `Float32` as sklearn trees
compare it.
"""
scale_input(model::MasureModel, value::Real, j::Int) =
    Float32((Float64(value) - model.input_mean[j]) / model.input_scale[j])

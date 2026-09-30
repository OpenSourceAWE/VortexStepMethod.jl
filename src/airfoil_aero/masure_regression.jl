"""
Masure regression - pure Julia evaluation of the Extra-Trees models that predict the
polars of a parametric leading-edge-inflatable (LEI) airfoil.

The trained scikit-learn models are https://doi.org/10.5281/zenodo.16925758;
`scripts/export_masure_models.py` converts them to the `.npz` files read here.
"""

"Airfoil parameters the masure regression takes, in its input order (alpha follows)."
const MASURE_PARAMETERS = ("t", "eta", "kappa", "delta", "lambda", "phi")

"Reynolds numbers a masure regression model exists for, with their file suffix."
const MASURE_REYNOLDS = Dict(1.0e6 => "1e6", 5.0e6 => "5e6", 2.0e7 => "2e7")

const _MASURE_CACHE = Dict{Tuple{Float64,String}, Any}()

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

function ExtraTreesForest(data::AbstractDict, prefix::String)
    array(name) = data["$(prefix)_$(name)"]
    return ExtraTreesForest(array("roots") .+ Int32(1), array("left") .+ Int32(1),
                            array("right") .+ Int32(1), array("feature") .+ Int32(1),
                            array("threshold"), array("value"))
end

"""
    load_masure_model(Re, ml_models_dir) -> MasureModel

Load the masure regression model for Reynolds number `Re` (1e6, 5e6 or 2e7) from
`ml_models_dir/ET_re<Re>.npz`, as written by `scripts/export_masure_models.py`.
"""
function load_masure_model(Re::Real, ml_models_dir::AbstractString)
    haskey(MASURE_REYNOLDS, Re) || error("No masure regression model for Re = $Re; " *
        "available: $(join(sort!(collect(keys(MASURE_REYNOLDS))), ", ")).")
    path = joinpath(ml_models_dir, "ET_re$(MASURE_REYNOLDS[Re]).npz")
    get!(_MASURE_CACHE, (Float64(Re), abspath(path))) do
        isfile(path) || error("Masure regression model not found at $path. Download " *
            "the models from https://doi.org/10.5281/zenodo.16925758 and convert them " *
            "with scripts/export_masure_models.py.")
        data = npzread(path)
        MasureModel(data["input_mean"], data["input_scale"],
                    [ExtraTreesForest(data, "output$k") for k in 0:2])
    end::MasureModel
end

"""
    predict(forest::ExtraTreesForest, x) -> Float64

Mean of the leaf values the scaled input `x` reaches in each tree.
"""
function predict(forest::ExtraTreesForest, x::AbstractVector{Float32})
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
    raw = [Float64(params[name]) for name in MASURE_PARAMETERS]
    x = zeros(Float32, length(raw) + 1)
    coefficients = zeros(length(alpha), length(model.forests))
    for (i, angle) in enumerate(alpha)
        for (j, value) in enumerate((raw..., Float64(angle)))
            # sklearn trees compare the scaled input rounded to Float32.
            x[j] = Float32((value - model.input_mean[j]) / model.input_scale[j])
        end
        for (k, forest) in enumerate(model.forests)
            coefficients[i, k] = predict(forest, x)
        end
    end
    return coefficients[:, 2], coefficients[:, 1], coefficients[:, 3]
end

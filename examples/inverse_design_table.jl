# Sample the kite designs of inverse_design_setup.jl, keep one per cell of polar-property
# space, check that a nearest-match lookup moves smoothly, and write the table the
# interactive page reads. Run with threads: julia --project=examples -t auto <this file>
using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using CodecZlib
using JSON
using Random
using Sobol
include(joinpath(@__DIR__, "inverse_design_setup.jl"))

TABLE_FILE = joinpath(dirname(@__DIR__), "docs", "src", "assets",
                      "inverse_design_kites.json.gz")
N_POOL = 2^16
# Cells per polar property; the table keeps one design per occupied cell.
CELLS = 10
PROPERTY_KEYS = (cl_max="clMax", stall_angle="stallAngle", lift_to_drag="maxLiftToDrag",
                 cd_zero_lift="cdZeroLift")

"Design vectors spread evenly over the parameter ranges, and their polars."
function sample_pool(model, n)
    sequence = SobolSeq(length(DESIGN_PARAMETERS))
    designs = [from_unit(next!(sequence)) for _ in 1:n]
    polars = Vector{Any}(undef, n)
    Threads.@threads for i in 1:n
        polars[i] = design_polar(model, designs[i])
    end
    return designs, polars
end

"Property vectors scaled to [0, 1] over the range the designs span."
function scaled_properties(properties)
    values = reduce(hcat, [collect(p) for p in properties])
    lower, upper = minimum(values; dims=2), maximum(values; dims=2)
    return (values .- lower) ./ (upper .- lower), vec(lower), vec(upper)
end

"""
Indices of one design per occupied cell of the scaled property space: the one nearest the
V3's parameters, so neighbouring cells hold kites changed from the V3 alike.
"""
function select_per_cell(scaled, designs, x_v3; cells=CELLS)
    best = Dict{NTuple{4, Int}, Int}()
    change(i) = norm(to_unit(designs[i]) - to_unit(x_v3))
    for i in axes(scaled, 2)
        cell = Tuple(min.(floor.(Int, scaled[:, i] .* cells), cells - 1))
        if !haskey(best, cell) || change(i) < change(best[cell])
            best[cell] = i
        end
    end
    return sort!(collect(values(best)))
end

"Index in `subset` of the design nearest `target` in scaled property space."
nearest_design(scaled, subset, target) =
    argmin(i -> sum(k -> abs2(scaled[k, i] - target[k]), 1:4), subset)

"""
Sweeps of one slider at a time across its range in 100 steps, from settings drawn from the
designs, picking the design nearest in scaled property space. Returns the RMS parameter
change of each step over the parameter ranges, and how far each pick's moved property lies
from the slider.
"""
function slider_sweeps(scaled, designs, subset; n_starts=100, rng=MersenneTwister(1))
    jumps, misses = Float64[], Float64[]
    for _ in 1:n_starts, slider in 1:4
        settings = scaled[:, rand(rng, subset)]
        previous = nothing
        for value in range(0, 1; length=100)
            settings[slider] = value
            pick = nearest_design(scaled, subset, settings)
            push!(misses, abs(scaled[slider, pick] - value))
            isnothing(previous) || push!(jumps, norm(to_unit(designs[pick]) -
                to_unit(designs[previous])) / sqrt(length(DESIGN_PARAMETERS)))
            previous = pick
        end
    end
    return jumps, misses
end

"Values rounded to `digits`, NaN written as JSON null."
rounded(values, digits) = [isnan(v) ? nothing : round(v; digits) for v in values]

"One design of the table: its parameters, properties, polar and section edges."
function table_entry(x, polar, properties, sections)
    edges = [(le=rounded(sections.le[i], 3), te=rounded(sections.te[i], 3))
             for i in reverse(eachindex(sections.le))]
    return Dict(
        "parameters" => Dict(p.key => round(v; sigdigits=4)
                             for (p, v) in zip(DESIGN_PARAMETERS, x)),
        "properties" => Dict(PROPERTY_KEYS[k] => round(properties[k]; sigdigits=4)
                             for k in keys(PROPERTY_KEYS)),
        "cl" => rounded(polar.cl, 4), "cd" => rounded(polar.cd, 5),
        "sections" => [Dict("le" => e.le, "te" => e.te) for e in edges])
end

"Write the table of the designs in `subset`, with the V3 as reference, as gzipped JSON."
function write_table(path, designs, polars, properties, subset, x_v3, v3_polar)
    reference = table_entry(x_v3, v3_polar, polar_properties(v3_polar), V3_RIBS)
    reference["name"] = "TU Delft V3"
    table = Dict(
        "source" => "VortexStepMethod.jl inverse design: $(length(subset)) LEI kites of " *
            "$(round(V3_AREA; digits=1)) m² projected area, VSM with masure-regression " *
            "profiles at Re = 5e6",
        "alpha" => ALPHA,
        "parameters" => [Dict("key" => p.key, "label" => p.label, "unit" => p.unit)
                         for p in DESIGN_PARAMETERS],
        "designs" => [table_entry(designs[i], polars[i], properties[i],
                                  design_sections(designs[i])) for i in subset],
        "reference" => reference)
    mkpath(dirname(path))
    open(GzipCompressorStream, path, "w") do io
        JSON.print(io, table)
    end
end

model = load_masure_model(REYNOLDS, MASURE_MODELS_DIR)
x_v3, _ = fit_v3_shape()
v3_polar = kite_polar(masure_wing(model, V3_RIBS))
if isfile(POOL_FILE)
    designs, polars = read_pool(POOL_FILE)
else
    designs, polars = sample_pool(model, N_POOL)
    write_pool(POOL_FILE, designs, polars)
end
properties = polar_properties.(polars)
valid = findall(eachindex(polars)) do i
    all(isfinite, properties[i]) && all(>(0), filter(!isnan, polars[i].cd))
end
designs, polars, properties = designs[valid], polars[valid], properties[valid]
scaled, lower, upper = scaled_properties(properties)
println("$(length(valid)) of $N_POOL designs show all four properties and positive drag.")
println("property ranges: ", NamedTuple{keys(PROPERTY_KEYS)}(
    tuple(zip(round.(lower; sigdigits=3), round.(upper; sigdigits=3))...)))
for cells in (0, 6, 8, 10, 12, 14, 16)
    candidates = cells == 0 ? collect(axes(scaled, 2)) :
                 select_per_cell(scaled, designs, x_v3; cells)
    jumps, misses = slider_sweeps(scaled, designs, candidates)
    println(rpad(cells == 0 ? "whole pool" : "$(cells)^4 cells", 12),
            lpad(length(candidates), 6), " designs: jumps above 0.1 per sweep ",
            round(99 * mean(jumps .> 0.1); digits=2),
            ", median and 95th percentile slider miss ", round(median(misses); digits=3),
            " ", round(quantile(misses, 0.95); digits=3))
end
subset = select_per_cell(scaled, designs, x_v3)
write_table(TABLE_FILE, designs, polars, properties, subset, x_v3, v3_polar)
println("wrote $TABLE_FILE, $(round(filesize(TABLE_FILE) / 2^20; digits=2)) MiB")

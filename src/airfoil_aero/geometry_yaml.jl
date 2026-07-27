"""
    write_geometry_yaml(path, section_rows, airfoil_rows)

Write a geometry YAML via [`write_yaml`](@ref), one line per section/airfoil row.
`section_rows` are `[airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]`; `airfoil_rows`
are `[airfoil_id, type, info_dict]` where `info_dict` holds `dat_file`,
`csv_file_path`, and optionally `raw_dat_file`, `cp_file`, `cf_file`.
"""
function write_geometry_yaml(path::String, section_rows, airfoil_rows)
    data = Dict(
        "wing_sections" => Dict(
            "headers" => ["airfoil_id", "LE_x", "LE_y", "LE_z",
                          "TE_x", "TE_y", "TE_z"],
            "data" => section_rows),
        "wing_airfoils" => Dict(
            "headers" => ["airfoil_id", "type", "info_dict"],
            "data" => [[id, type, info] for (id, type, info) in airfoil_rows]))
    return write_yaml(path, data)
end

yaml_scalar(x::Bool) = string(x)
yaml_scalar(x::Integer) = string(x)
yaml_scalar(x::Real) = string(round(Float64(x); digits=3))
yaml_scalar(x::AbstractString) = "\"$(replace(x, '\\' => '/'))\""
yaml_scalar(::Nothing) = "null"
yaml_scalar(x) = string(x)

yaml_sorted(x::AbstractDict) = sort!(collect(x); by = p -> string(first(p)))
yaml_isscalar(x) = x isa Union{Real,AbstractString,Bool,Nothing}
yaml_scalardict(x) = x isa AbstractDict && all(yaml_isscalar, values(x))
yaml_flowable(x) = yaml_isscalar(x) || yaml_scalardict(x) ||
    (x isa AbstractVector && all(e -> yaml_isscalar(e) || yaml_scalardict(e), x))

yaml_flow(x) =
    yaml_isscalar(x) ? yaml_scalar(x) :
    x isa AbstractDict ?
        "{" * join(("$k: $(yaml_flow(v))" for (k, v) in yaml_sorted(x)), ", ") * "}" :
    x isa AbstractVector ? "[" * join(map(yaml_flow, x), ", ") * "]" :
    yaml_scalar(x)

function yaml_emit(io, x, indent)
    pad = "  "^indent
    if x isa AbstractDict
        for (k, v) in yaml_sorted(x)
            yaml_flowable(v) ? println(io, pad, k, ": ", yaml_flow(v)) :
                (println(io, pad, k, ":"); yaml_emit(io, v, indent + 1))
        end
    elseif x isa AbstractVector
        for el in x
            yaml_flowable(el) ? println(io, pad, "- ", yaml_flow(el)) :
                (println(io, pad, "-"); yaml_emit(io, el, indent + 1))
        end
    else
        println(io, pad, yaml_scalar(x))
    end
end

"""
    write_yaml(path, data)

Write `data` (nested `Dict`s, vectors, scalars) to `path` as YAML with every float
rounded to millimetre precision (3 decimals) and every leaf list on one line. The
single writer for all generated geometry YAMLs, so their formatting stays consistent.
Mapping keys are emitted in sorted order.
"""
function write_yaml(path::String, data)
    open(path, "w") do io
        yaml_emit(io, data, 0)
    end
    return path
end

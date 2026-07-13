"""
    CpPolar(data::CpData) -> CpPolar

Build the per-slice `cp(alpha, delta)` interpolants (linear, `Line()`
extrapolation) from raw [`CpData`](@ref).
"""
function CpPolar(data::CpData)
    build(grid) = [linear_interpolation((data.alpha_range, data.delta_range),
                                        grid[i, :, :]; extrapolation_bc=Line())
                   for i in 1:data.n_chord]
    return CpPolar(data, build(data.cp_upper), build(data.cp_lower))
end

"""
    CpData(polar::CpPolar) -> CpData

The raw samples underlying a [`CpPolar`](@ref) (strips the interpolants).
"""
CpData(polar::CpPolar) = polar.data

"""
    read_cp_data(path) -> Union{Nothing, CpData}

Load a section surface-pressure table from a single `cp.csv`; returns `nothing` if
`path` does not exist. Header: `alpha, delta, up@<x/c>…, lo@<x/c>…`. `alpha`/`delta`
are in degrees (converted to radians); `chord_x`/`n_chord` come from the `@<x/c>`
column headers (upper and lower must share the same slices). Each row is one
`(alpha, delta)` point; the rows must form a full `alpha × delta` grid.
"""
function read_cp_data(path::AbstractString)
    isfile(path) || return nothing
    lines = readlines(String(path))
    header = strip.(split(lines[1], ','))
    xof(prefix) = [parse(Float64, String(h[(length(prefix) + 1):end]))
                   for h in header if startswith(h, prefix)]
    chord_x = xof("up@")
    chord_x == xof("lo@") ||
        throw(ArgumentError("cp.csv upper/lower chord slices differ: $path"))
    n_chord = length(chord_x)
    n_chord > 0 || throw(ArgumentError("cp.csv has no up@/lo@ columns: $path"))

    rows = [strip.(split(l, ',')) for l in lines[2:end] if !isempty(strip(l))]
    alphas = [deg2rad(parse(Float64, String(r[1]))) for r in rows]
    deltas = [deg2rad(parse(Float64, String(r[2]))) for r in rows]
    alpha_range = sort(unique(alphas))
    delta_range = sort(unique(deltas))
    length(rows) == length(alpha_range) * length(delta_range) ||
        throw(ArgumentError("cp.csv rows do not form a full alpha × delta grid: $path"))

    cp_upper = fill(NaN, n_chord, length(alpha_range), length(delta_range))
    cp_lower = fill(NaN, n_chord, length(alpha_range), length(delta_range))
    for (k, r) in enumerate(rows)
        length(r) == 2 + 2n_chord ||
            throw(ArgumentError("cp.csv row $(k) has wrong column count: $path"))
        ia = searchsortedfirst(alpha_range, alphas[k])
        jd = searchsortedfirst(delta_range, deltas[k])
        vals = parse.(Float64, String.(r[3:end]))
        cp_upper[:, ia, jd] .= vals[1:n_chord]
        cp_lower[:, ia, jd] .= vals[(n_chord + 1):2n_chord]
    end
    return CpData(n_chord, chord_x, alpha_range, delta_range, cp_upper, cp_lower)
end

"""
    write_cp_data(path, data::CpData) -> path

Write a [`CpData`](@ref) to `path` as a `cp.csv` in the [`read_cp_data`](@ref)
format: header `alpha, delta, up@<x/c>…, lo@<x/c>…`, one row per `(alpha, delta)`
with `alpha`/`delta` in degrees.
"""
function write_cp_data(path::AbstractString, data::CpData)
    xhdr(prefix) = [string(prefix, round(x; digits=5)) for x in data.chord_x]
    header = ["alpha"; "delta"; xhdr("up@"); xhdr("lo@")]
    open(String(path), "w") do io
        println(io, join(header, ","))
        for (jd, delta) in enumerate(data.delta_range),
            (ia, alpha) in enumerate(data.alpha_range)
            row = [string(round(rad2deg(alpha); digits=6)),
                   string(round(rad2deg(delta); digits=6))]
            append!(row, string.(data.cp_upper[:, ia, jd]))
            append!(row, string.(data.cp_lower[:, ia, jd]))
            println(io, join(row, ","))
        end
    end
    return path
end

"""
    cp_distribution(polar::CpPolar, alpha, delta) -> (cp_upper, cp_lower)

Upper and lower surface `Cp` at every chord slice for `alpha` and `delta` (radians).
"""
function cp_distribution(polar::CpPolar, alpha, delta)
    return ([itp(alpha, delta) for itp in polar.cp_upper_interp],
            [itp(alpha, delta) for itp in polar.cp_lower_interp])
end

"""
    delta_cp(polar::CpPolar, alpha, delta) -> Vector{Float64}

Chordwise load `ΔCp = Cp_lower - Cp_upper` at every chord slice.
"""
function delta_cp(polar::CpPolar, alpha, delta)
    upper, lower = cp_distribution(polar, alpha, delta)
    return lower .- upper
end

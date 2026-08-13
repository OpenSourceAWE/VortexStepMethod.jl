"""
    SectionAero

Per-2D-section surface aerodynamics at native airfoil resolution: the closed contour
`(x, y)` per trailing-edge deflection, plus surface pressure `cp` and skin friction
`cf` per contour node over an `(alpha, delta)` grid (both radians). Replaces the old
chord-slice Cp table. Use [`section_surface`](@ref) to get the interpolated contour +
`cp`/`cf` at any `(alpha, delta)`.

Fields:
- `alpha_range`, `delta_range`: grid axes (radians).
- `x`, `y`: contour node coordinates, `n_node × n_delta`.
- `cp`, `cf`: surface pressure and skin friction, `n_node × n_alpha × n_delta`.
- `cp_interp`, `cf_interp`, `x_interp`, `y_interp`: per-node interpolants.
"""
struct SectionAero{A,B,C,D}
    alpha_range::Vector{Float64}
    delta_range::Vector{Float64}
    x::Matrix{Float64}
    y::Matrix{Float64}
    cp::Array{Float64,3}
    cf::Array{Float64,3}
    cp_interp::Vector{A}
    cf_interp::Vector{B}
    x_interp::Vector{C}
    y_interp::Vector{D}
end

"""
    SectionAero(alpha_range, delta_range, x, y, cp, cf) -> SectionAero

Build the per-node interpolants (linear, `Line()` extrapolation) from the raw grid. A
single-element `delta_range` is padded to two identical slices so every interpolant
shares one 2D `(alpha, delta)` form.
"""
function SectionAero(alpha_range, delta_range, x, y, cp, cf)
    alpha_range = collect(float.(alpha_range))
    delta_range = collect(float.(delta_range))
    if length(delta_range) == 1
        # linear_interpolation needs >= 2 delta nodes; duplicate the single slice.
        # The second delta coordinate is arbitrary: both slices are equal, so the
        # interpolant is constant in delta and its exact spacing never matters.
        delta_range = [delta_range[1], delta_range[1] + 1]
        x = hcat(x, x)
        y = hcat(y, y)
        cp = cat(cp, cp; dims=3)
        cf = cat(cf, cf; dims=3)
    end
    n_node = size(cp, 1)
    cp_interp = [linear_interpolation((alpha_range, delta_range), cp[i, :, :];
                                      extrapolation_bc=Line()) for i in 1:n_node]
    cf_interp = [linear_interpolation((alpha_range, delta_range), cf[i, :, :];
                                      extrapolation_bc=Line()) for i in 1:n_node]
    x_interp = [linear_interpolation(delta_range, x[i, :]; extrapolation_bc=Line())
                for i in 1:n_node]
    y_interp = [linear_interpolation(delta_range, y[i, :]; extrapolation_bc=Line())
                for i in 1:n_node]
    return SectionAero(alpha_range, delta_range, x, y, cp, cf,
                       cp_interp, cf_interp, x_interp, y_interp)
end

"""
    section_surface(aero::SectionAero, alpha, delta) -> (x, y, cp, cf)

Interpolated closed contour `(x, y)` and per-node surface pressure `cp` and skin
friction `cf` at `(alpha, delta)` (radians).
"""
function section_surface(aero::SectionAero, alpha, delta)
    x = [itp(delta) for itp in aero.x_interp]
    y = [itp(delta) for itp in aero.y_interp]
    cp = [itp(alpha, delta) for itp in aero.cp_interp]
    cf = [itp(alpha, delta) for itp in aero.cf_interp]
    return x, y, cp, cf
end

"""
    delta_suffix(delta) -> String

Filename tag for a non-zero trailing-edge deflection `delta` [rad], e.g. `d5`, `dm3`
(−3°), `d2p5` (2.5°). Uses millidegree precision (negatives as `m`, the decimal point
as `p`) so sub-degree deflections get distinct `.dat` files instead of colliding.
Shared by [`write_section_aero`](@ref VortexStepMethod.AirfoilAero.write_section_aero)
and [`read_section_aero`](@ref).
"""
function delta_suffix(delta)
    deg = round(rad2deg(delta); digits=3)
    s = isinteger(deg) ? string(Int(deg)) : string(deg)
    return "d" * replace(s, "-" => "m", "." => "p")
end

"""
    read_dat(path) -> (x, y)

Read Selig `.dat` airfoil coordinates (two whitespace-separated columns), skipping the
name/header line and any non-numeric lines.
"""
function read_dat(path::AbstractString)
    x = Float64[]
    y = Float64[]
    for ln in eachline(String(path))
        p = split(strip(ln))
        length(p) >= 2 || continue
        xv, yv = tryparse(Float64, p[1]), tryparse(Float64, p[2])
        (xv === nothing || yv === nothing) && continue
        push!(x, xv)
        push!(y, yv)
    end
    return x, y
end

"""
    read_node_table(path) -> (alpha, delta, values)

Read a per-node aero table into radian `alpha`/`delta` vectors (one entry per row) and a
`nrow × n_node` value matrix. The file suffix picks the format: `.arrow` (columns
`alpha`, `delta` in degrees and a per-row list column `values`), anything else CSV
(header `alpha, delta, n0, n1, …`, angles in degrees).
"""
function read_node_table(path::AbstractString)
    if endswith(String(path), ".arrow")
        # bytes, not the mmap Arrow.Table(path) takes: Windows locks a mapped file
        table = Arrow.Table(read(String(path)))
        values = Matrix{Float64}(undef, length(table.alpha), length(first(table.values)))
        for k in axes(values, 1)
            @inbounds values[k, :] .= table.values[k]
        end
        return deg2rad.(table.alpha), deg2rad.(table.delta), values
    end
    lines = [l for l in readlines(String(path)) if !isempty(strip(l))]
    n_row = length(lines) - 1
    n_col = count(==(','), lines[1]) + 1
    alpha = Vector{Float64}(undef, n_row)
    delta = Vector{Float64}(undef, n_row)
    values = Matrix{Float64}(undef, n_row, n_col - 2)
    for k in 1:n_row
        fields = split(lines[k + 1], ',')
        alpha[k] = deg2rad(parse(Float64, fields[1]))
        delta[k] = deg2rad(parse(Float64, fields[2]))
        @inbounds for j in 3:n_col
            values[k, j - 2] = parse(Float64, fields[j])
        end
    end
    return alpha, delta, values
end

"""
    write_node_rows(path, alpha, delta, values) -> path

Write a per-node aero table from radian `alpha`/`delta` vectors (one entry per row)
and a `nrow × n_node` value matrix. The file suffix picks the format, matching
[`read_node_table`](@ref): `.arrow` for the binary form, anything else CSV. Angles
are written in degrees either way.
"""
function write_node_rows(path::AbstractString, alpha, delta, values)
    if endswith(String(path), ".arrow")
        Arrow.write(String(path),
            (alpha=rad2deg.(alpha), delta=rad2deg.(delta),
             values=[values[k, :] for k in axes(values, 1)]))
        return path
    end
    open(String(path), "w") do io
        println(io, "alpha,delta," * join(("n$(j - 1)" for j in axes(values, 2)), ","))
        for k in axes(values, 1)
            row = [rad2deg(alpha[k]), rad2deg(delta[k])]
            append!(row, @view values[k, :])
            println(io, join(row, ","))
        end
    end
    return path
end

"""
    convert_node_table(src, dst) -> dst

Rewrite a per-node aero table in the format `dst`'s suffix names. Use to move an
existing dataset between CSV and Arrow without re-running the (slow) airfoil solver
that produced it.
"""
convert_node_table(src::AbstractString, dst::AbstractString) =
    write_node_rows(dst, read_node_table(src)...)

"""
    read_section_aero(dat_file, cp_file, cf_file) -> Union{Nothing, SectionAero}

Assemble a [`SectionAero`](@ref) from the human-readable files: the airfoil contour
(`dat_file`, plus `{stem}_{delta_suffix(δ)}.dat` per non-zero deflection) and the per-node `Cp`
and `cf` tables in the `.dat` node order, CSV or Arrow as their suffix says (see
[`read_node_table`](@ref)). Returns `nothing` if any file is missing. This is the single
loader for both provided and generated aero.
"""
function read_section_aero(dat_file::AbstractString, cp_file::AbstractString,
                           cf_file::AbstractString)
    (isfile(dat_file) && isfile(cp_file) && isfile(cf_file)) || return nothing
    alpha, delta, cp_rows = read_node_table(cp_file)
    _, _, cf_rows = read_node_table(cf_file)
    alpha_range = sort(unique(alpha))
    delta_range = sort(unique(delta))
    n_node = size(cp_rows, 2)
    stem = replace(String(dat_file), r"\.dat$" => "")
    x = fill(NaN, n_node, length(delta_range))
    y = fill(NaN, n_node, length(delta_range))
    for (jd, d) in enumerate(delta_range)
        xd, yd = read_dat(iszero(d) ? String(dat_file) :
                          "$(stem)_$(delta_suffix(d)).dat")
        x[:, jd] .= xd
        y[:, jd] .= yd
    end
    cp = fill(NaN, n_node, length(alpha_range), length(delta_range))
    cf = fill(NaN, n_node, length(alpha_range), length(delta_range))
    for r in eachindex(alpha)
        ia = searchsortedfirst(alpha_range, alpha[r])
        jd = searchsortedfirst(delta_range, delta[r])
        cp[:, ia, jd] .= cp_rows[r, :]
        cf[:, ia, jd] .= cf_rows[r, :]
    end
    return SectionAero(alpha_range, delta_range, x, y, cp, cf)
end

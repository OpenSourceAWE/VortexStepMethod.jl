"""
    KulfanParameters

Kulfan CST parameters for an airfoil: the weights of the class-shape transformation
each surface is built from, the shared leading-edge modification weight and the
trailing-edge thickness. Lives here rather than in `AirfoilAero` because a
[`Panel`](@ref) carries the shape it is currently flying (`live_shape`), which a
plot or a traction pattern reads without knowing how the shape was produced. Fit one
with `AirfoilAero.fit_kulfan_parameters`, deform one with `AirfoilAero.deform_kulfan`
and turn one into coordinates with `AirfoilAero.kulfan_to_coordinates`.

Mutable, so a live polar source can rewrite one shape every solve instead of building a
new one, and the panel pointing at it follows without being told
(`AirfoilAero.deform_kulfan!`). That also makes a panel's `live_shape` the very object it
was sampled from rather than a copy that compares equal to it.

# Fields
- `upper_weights::Vector{Float64}`: weights for upper surface
- `lower_weights::Vector{Float64}`: weights for lower surface
- `leading_edge_weight::Float64`: Leading edge modification weight
- `TE_thickness::Float64`: Trailing edge thickness
"""
mutable struct KulfanParameters
    upper_weights::Vector{Float64}
    lower_weights::Vector{Float64}
    leading_edge_weight::Float64
    TE_thickness::Float64
end

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
    read_dat_coordinates(path) -> (x, y)

Read Selig `.dat` airfoil coordinates (two whitespace-separated columns), skipping the
name/header line, comments, and any row that is not a finite pair. The single `.dat`
reader; [`write_dat`](@ref VortexStepMethod.AirfoilAero.write_dat) is its writer.
"""
function read_dat_coordinates(path::AbstractString)
    x = Float64[]
    y = Float64[]
    for line in eachline(String(path))
        fields = split(strip(line))
        length(fields) >= 2 || continue
        xp = tryparse(Float64, fields[1])
        yp = tryparse(Float64, fields[2])
        (xp === nothing || yp === nothing) && continue
        (isfinite(xp) && isfinite(yp)) || continue
        push!(x, xp)
        push!(y, yp)
    end
    return x, y
end

"""
    csv_fields(values) -> String

Join `values` into one comma-separated CSV line at 16 significant digits.
"""
csv_fields(values) = join((@sprintf("%.16g", v) for v in values), ",")

"""
    read_node_table(path) -> (alpha, delta, values, columns)

Read an aero table into a radian `alpha` vector (one entry per row), a radian `delta`
vector or `nothing` when the table has no `delta` column, the `nrow × ncol` matrix of
the value columns and their names. The file suffix picks the format: `.arrow` (angle
columns, a per-row list column `values`, the names in the `columns` metadata entry) or
anything else CSV (a header line naming every column, found case-insensitively in any
order). Angles are stored in degrees.
"""
function read_node_table(path::AbstractString)
    names, data = endswith(String(path), ".arrow") ? read_arrow_columns(path) :
                                                     read_csv_columns(path)
    lowercase_names = lowercase.(names)
    alpha_col = findfirst(==("alpha"), lowercase_names)
    alpha_col === nothing && error("Table $path has no alpha column")
    delta_col = findfirst(==("delta"), lowercase_names)
    value_cols = setdiff(eachindex(names), (alpha_col, delta_col))
    delta = delta_col === nothing ? nothing : deg2rad.(data[:, delta_col])
    return deg2rad.(data[:, alpha_col]), delta, view(data, :, value_cols), names[value_cols]
end

"""
    read_arrow_columns(path) -> (names, data)

Read an Arrow aero table into its column names and one `Float64` matrix: the angle
columns, then the list column `values` spread under the names its `columns` metadata
entry holds, `n0, n1, …` where it has none.
"""
function read_arrow_columns(path::AbstractString)
    # bytes, not the mmap Arrow.Table(path) takes: Windows locks a mapped file
    table = Arrow.Table(read(String(path)))
    angles = filter(!=(:values), Arrow.names(table))
    n_value = length(first(table.values))
    metadata = Arrow.getmetadata(table)
    columns = metadata === nothing || !haskey(metadata, "columns") ?
              ["n$(j - 1)" for j in 1:n_value] : split(metadata["columns"], ',')
    data = Matrix{Float64}(undef, length(table.values), length(angles) + n_value)
    for (j, angle) in enumerate(angles)
        data[:, j] .= getproperty(table, angle)
    end
    for k in axes(data, 1)
        @inbounds data[k, (length(angles) + 1):end] .= table.values[k]
    end
    return [String.(angles); String.(columns)], data
end

"""
    read_csv_columns(path) -> (names, data)

Read a CSV table with one header line into its column names and one `Float64` matrix.
Throws on a row whose field count differs from the header's.
"""
function read_csv_columns(path::AbstractString)
    lines = [strip(l) for l in readlines(String(path)) if !isempty(strip(l))]
    isempty(lines) && error("Table $path is empty")
    names = String.(strip.(split(lines[1], ',')))
    data = Matrix{Float64}(undef, length(lines) - 1, length(names))
    for k in axes(data, 1)
        fields = split(lines[k + 1], ',')
        length(fields) == length(names) ||
            error("Row $k of $path has $(length(fields)) fields, header $(length(names))")
        @inbounds for j in eachindex(fields)
            data[k, j] = parse(Float64, fields[j])
        end
    end
    return names, data
end

"""
    write_node_rows(path, alpha, delta, values; columns) -> path

Write an aero table from a radian `alpha` vector (one entry per row), a radian `delta`
vector or `nothing` for a table without a `delta` column, and a `nrow × ncol` value
matrix whose columns are named `columns` (default `n0, n1, …`). The file suffix picks
the format, matching [`read_node_table`](@ref): `.arrow`, or anything else CSV at 16
significant digits. Angles are written in degrees at 16 significant digits either way.
"""
function write_node_rows(path::AbstractString, alpha, delta, values;
                         columns=["n$(j - 1)" for j in axes(values, 2)])
    degrees(angle) = round.(rad2deg.(angle); sigdigits=16)
    angles = delta === nothing ? (alpha=degrees(alpha),) :
                                 (alpha=degrees(alpha), delta=degrees(delta))
    if endswith(String(path), ".arrow")
        Arrow.write(String(path),
            merge(angles, (values=[values[k, :] for k in axes(values, 1)],));
            metadata=["columns" => join(columns, ",")])
        return path
    end
    open(String(path), "w") do io
        println(io, join((keys(angles)..., columns...), ","))
        for k in axes(values, 1)
            println(io, csv_fields(angle[k] for angle in angles), ",",
                    csv_fields(@view values[k, :]))
        end
    end
    return path
end

"""
    convert_node_table(src, dst) -> dst

Rewrite an aero table, a per-node `Cp`/`cf` table or a polar, in the format `dst`'s
suffix names, keeping its column names. Use to move an existing dataset between CSV and
Arrow without re-running the (slow) airfoil solver that produced it.
"""
function convert_node_table(src::AbstractString, dst::AbstractString)
    alpha, delta, values, columns = read_node_table(src)
    return write_node_rows(dst, alpha, delta, values; columns)
end

"""
    read_section_aero(dat_file, cp_file, cf_file) -> Union{Nothing, SectionAero}

Assemble a [`SectionAero`](@ref) from the human-readable files: the airfoil contour
(`dat_file`, plus `{stem}_{delta_suffix(δ)}.dat` per non-zero deflection) and the
per-node `Cp` and `cf` tables in the `.dat` node order, CSV or Arrow as their suffix
says (see [`read_node_table`](@ref)). Returns `nothing` if one of the three named files
is missing, and `nothing` with a warning if a contour does not hold exactly the tables'
nodes. This is the single loader for both provided and generated aero.
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
        contour = iszero(d) ? String(dat_file) : "$(stem)_$(delta_suffix(d)).dat"
        xd, yd = isfile(contour) ? read_dat_coordinates(contour) : (Float64[], Float64[])
        if length(xd) != n_node
            @warn "$contour holds $(length(xd)) finite coordinates, not the $n_node " *
                  "nodes of its Cp table; this airfoil gets no surface aero."
            return nothing
        end
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

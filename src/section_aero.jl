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
single-element `delta_range` is padded to two constant slices so every interpolant
shares one 2D `(alpha, delta)` form.
"""
function SectionAero(alpha_range, delta_range, x, y, cp, cf)
    alpha_range = collect(float.(alpha_range))
    delta_range = collect(float.(delta_range))
    if length(delta_range) == 1
        delta_range = [delta_range[1], delta_range[1] + 1.0]
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

Read a per-node aero CSV (header `alpha, delta, n0, n1, …`; angles in degrees) into
radian `alpha`/`delta` vectors (one entry per row) and a `nrow × n_node` value matrix.
"""
function read_node_table(path::AbstractString)
    lines = [l for l in readlines(String(path)) if !isempty(strip(l))]
    rows = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
    alpha = deg2rad.([r[1] for r in rows])
    delta = deg2rad.([r[2] for r in rows])
    values = reduce(vcat, (permutedims(r[3:end]) for r in rows))
    return alpha, delta, values
end

"""
    read_section_aero(dat_file, cp_file, cf_file) -> Union{Nothing, SectionAero}

Assemble a [`SectionAero`](@ref) from the human-readable files: the airfoil contour
(`dat_file`, plus `{stem}_d{deg}.dat` per non-zero deflection) and the per-node `Cp`
and `cf` tables (`alpha, delta, n0…` in the `.dat` node order). Returns `nothing` if any
file is missing. This is the single loader for both provided and generated aero.
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
                          "$(stem)_d$(round(Int, rad2deg(d))).dat")
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

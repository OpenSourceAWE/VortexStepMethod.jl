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
    write_section_aero(path, aero::SectionAero) -> path

Write a [`SectionAero`](@ref) to `path` as an `.npz` (arrays `alpha_range`,
`delta_range`, `x`, `y`, `cp`, `cf`).
"""
function write_section_aero(path::AbstractString, aero::SectionAero)
    npzwrite(String(path), Dict("alpha_range" => aero.alpha_range,
        "delta_range" => aero.delta_range, "x" => aero.x, "y" => aero.y,
        "cp" => aero.cp, "cf" => aero.cf))
    return path
end

"""
    read_section_aero(path) -> Union{Nothing, SectionAero}

Load a [`SectionAero`](@ref) from an `.npz`; `nothing` if `path` does not exist.
"""
function read_section_aero(path::AbstractString)
    isfile(path) || return nothing
    d = npzread(String(path))
    return SectionAero(d["alpha_range"], d["delta_range"],
                       d["x"], d["y"], d["cp"], d["cf"])
end

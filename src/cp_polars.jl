"""
    CpPolar

Per-panel surface-pressure table for a 2D airfoil section, covering angle of
attack `alpha` and trailing-edge deflection `delta` (both in radians).

Stations are equally spaced in arc length: `n_panels` along the upper surface
followed by `n_panels` along the lower surface, each ordered from leading to
trailing edge. Station `i` tracks the same material point for every `delta`,
because arc-length fraction is preserved under a rigid trailing-edge rotation.

Fields:
- `n_panels`: stations per surface (total stations is `2 * n_panels`).
- `alpha_range`, `delta_range`: grid axes (radians).
- `station_frac`: arc-length fraction of each station (0 at LE, 1 at TE).
- `station_xy`: undeflected `(x, y)` of each station, `2*n_panels × 2`.
- `cp_grid`: Cp samples, `2*n_panels × length(alpha_range) × length(delta_range)`.
- `cp_interp`: one `cp(alpha, delta)` interpolant per station.
"""
struct CpPolar{I}
    n_panels::Int
    alpha_range::Vector{Float64}
    delta_range::Vector{Float64}
    station_frac::Vector{Float64}
    station_xy::Matrix{Float64}
    cp_grid::Array{Float64,3}
    cp_interp::Vector{I}
end

function CpPolar(n_panels, alpha_range, delta_range, station_frac, station_xy, cp_grid)
    interps = [linear_interpolation((alpha_range, delta_range), cp_grid[i, :, :];
                                    extrapolation_bc=Line()) for i in axes(cp_grid, 1)]
    return CpPolar(n_panels, collect(float.(alpha_range)), collect(float.(delta_range)),
                   station_frac, station_xy, cp_grid, interps)
end

"""
    arclength_fraction(x, y) -> Vector{Float64}

Cumulative arc length along an ordered polyline, normalized to `[0, 1]`.
"""
function arclength_fraction(x, y)
    s = zeros(length(x))
    for i in 2:length(x)
        s[i] = s[i-1] + hypot(x[i] - x[i-1], y[i] - y[i-1])
    end
    total = s[end]
    total > 0 && (s ./= total)
    return s
end

"""
    surface_arc(x, y, le_idx, upper) -> (s, idx)

Arc-length fraction `s` (LE→TE) and the source indices `idx` for one surface of
a Selig-ordered airfoil. `upper=true` walks `le_idx:-1:1`, `upper=false` walks
`le_idx:length(x)`.
"""
function surface_arc(x, y, le_idx, upper)
    idx = upper ? collect(le_idx:-1:1) : collect(le_idx:length(x))
    s = arclength_fraction(@view(x[idx]), @view(y[idx]))
    return s, idx
end

"""
    resample_surface(s, values, fracs) -> Vector{Float64}

Linearly interpolate `values` (given at arc-length fractions `s`) onto target
fractions `fracs`, dropping zero-length segments so the breakpoints are strictly
increasing.
"""
function resample_surface(s, values, fracs)
    keep = [true; diff(s) .> 0]
    itp = linear_interpolation(s[keep], values[keep]; extrapolation_bc=Line())
    return itp.(fracs)
end

"""
    generate_cp_polar(x, y; alpha_range, delta_range, n_panels, reynolds_number,
                      mach=0.0, crease_frac=0.9, npan=140, remove_nan=true) -> CpPolar

Build a per-panel surface-pressure table `cp[i](alpha, delta)` from Selig-ordered
airfoil coordinates `(x, y)`. For each `delta` the trailing edge is deflected with
[`turn_trailing_edge!`](@ref) and XFoil is solved over `alpha_range`; the Cp
distribution from `Xfoil.cpdump` is resampled onto `n_panels` equal arc-length
stations per surface. Non-converged points stay `NaN` and, when `remove_nan`,
are filled per station with [`interpolate_matrix_nans!`](@ref).
"""
function generate_cp_polar(x, y; alpha_range, delta_range, n_panels, reynolds_number,
        mach=0.0, crease_frac=0.9, npan=140, remove_nan=true)
    n_station = 2 * n_panels
    n_alpha = length(alpha_range)
    n_delta = length(delta_range)
    station_frac = [(i - 0.5) / n_panels for i in 1:n_panels]
    cp_grid = fill(NaN, n_station, n_alpha, n_delta)

    xb0 = collect(float.(x))
    yb0 = collect(float.(y))
    normalize_foil!(xb0, yb0)
    station_xy = baseline_station_xy(xb0, yb0, station_frac)

    for (jd, delta) in enumerate(delta_range)
        xb = copy(xb0)
        yb = copy(yb0)
        lower, upper = get_lower_upper(xb, yb, crease_frac)
        turn_trailing_edge!(delta, xb, yb, lower, upper, crease_frac)
        Xfoil.set_coordinates(xb, yb)
        xp, yp = Xfoil.pane(npan=npan)
        sweep_alpha_cp!(cp_grid, jd, xp, yp, alpha_range, reynolds_number, mach,
                        n_panels, station_frac)
    end

    if remove_nan
        fill_station_nans!(cp_grid)
    end
    return CpPolar(n_panels, alpha_range, delta_range, station_frac, station_xy, cp_grid)
end

"""
    generate_cp_polar(params::KulfanParameters; n_points=120, kwargs...) -> CpPolar

Generate a [`CpPolar`](@ref) from Kulfan/CST parameters, sampling coordinates with
[`kulfan_to_coordinates`](@ref) before delegating to the coordinate method.
"""
function generate_cp_polar(params::KulfanParameters; n_points::Int=120, kwargs...)
    x, y = kulfan_to_coordinates(params; n_points=n_points)
    return generate_cp_polar(x, y; kwargs...)
end

"""
    baseline_station_xy(x, y, station_frac) -> Matrix{Float64}

Undeflected `(x, y)` of each station: upper LE→TE first, then lower LE→TE.
"""
function baseline_station_xy(x, y, station_frac)
    le = argmin(x)
    xy = zeros(2 * length(station_frac), 2)
    n = length(station_frac)
    for (k, upper) in enumerate((true, false))
        s, idx = surface_arc(x, y, le, upper)
        rows = (k - 1) * n .+ (1:n)
        xy[rows, 1] = resample_surface(s, x[idx], station_frac)
        xy[rows, 2] = resample_surface(s, y[idx], station_frac)
    end
    return xy
end

"""
    sweep_alpha_cp!(cp_grid, jd, xp, yp, alpha_range, reynolds_number, mach,
                    n_panels, station_frac)

Solve XFoil over `alpha_range` for one deflected geometry `(xp, yp)` and store the
resampled Cp into `cp_grid[:, :, jd]`. Negative and positive angles are swept
outward from zero with a reinit at each side, mirroring [`run_solve_alpha`](@ref).
"""
function sweep_alpha_cp!(cp_grid, jd, xp, yp, alpha_range, reynolds_number, mach,
        n_panels, station_frac)
    le = argmin(xp)
    s_up, idx_up = surface_arc(xp, yp, le, true)
    s_lo, idx_lo = surface_arc(xp, yp, le, false)
    neg = sort(findall(<(0.0), alpha_range); rev=true)
    pos = sort(findall(>=(0.0), alpha_range))
    for order in (neg, pos)
        reinit = true
        for ia in order
            _, _, _, _, converged = Xfoil.solve_alpha(rad2deg(alpha_range[ia]),
                reynolds_number; iter=50, reinit=reinit, mach=mach, xtrip=(0.05, 0.05))
            reinit = false
            converged || continue
            _, cp = Xfoil.cpdump()
            length(cp) == length(xp) || continue
            cp_grid[1:n_panels, ia, jd] = resample_surface(s_up, cp[idx_up], station_frac)
            cp_grid[n_panels+1:end, ia, jd] =
                resample_surface(s_lo, cp[idx_lo], station_frac)
        end
    end
    return nothing
end

"""
    fill_station_nans!(cp_grid)

Fill `NaN`s in each station's `(alpha, delta)` slice with
[`interpolate_matrix_nans!`](@ref). Stations that never converged are skipped
with a warning rather than looping forever.
"""
function fill_station_nans!(cp_grid)
    for i in axes(cp_grid, 1)
        slice = cp_grid[i, :, :]
        if all(isnan, slice)
            @warn "Station $i has no converged Cp; leaving NaN."
            continue
        end
        interpolate_matrix_nans!(slice; prn=false)
        cp_grid[i, :, :] = slice
    end
    return nothing
end

"""
    cp_distribution(polar::CpPolar, alpha, delta) -> Vector{Float64}

Surface-pressure coefficients at all `2*n_panels` stations for `alpha` and `delta`
(radians): upper LE→TE first, then lower LE→TE.
"""
function cp_distribution(polar::CpPolar, alpha, delta)
    return [itp(alpha, delta) for itp in polar.cp_interp]
end

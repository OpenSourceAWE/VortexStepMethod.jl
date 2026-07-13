"""
    AbstractAirfoilSolver

Supertype for a 2D-airfoil solver backend. Concrete subtypes ([`XFoilSolver`](@ref),
[`NeuralFoilSolver`](@ref)) select, by dispatch, how a deformed section is analysed.
Add a backend by subtyping this and defining [`analyze_section`](@ref) (and
optionally [`analyze_sweep`](@ref)); the deformation and Cp resampling are shared.
"""
abstract type AbstractAirfoilSolver end

"""
    SectionSolution

Result of a single 2D analysis: integrated coefficients plus the surface pressure
as `(x/c, Cp)` per surface (each single-valued in x, leading edge to trailing
edge). `confidence` is `1.0`/`NaN` for XFoil (converged/not) or NeuralFoil's
`analysis_confidence`. The generator resamples the `(x, Cp)` pairs onto the shared
chord slices.
"""
struct SectionSolution
    alpha::Float64
    cl::Float64
    cd::Float64
    cm::Float64
    confidence::Float64
    x_upper::Vector{Float64}
    cp_upper::Vector{Float64}
    x_lower::Vector{Float64}
    cp_lower::Vector{Float64}
end

"""
    write_polar_csv(filepath, sols::Vector{SectionSolution})

Write a solver sweep (from any [`AbstractAirfoilSolver`](@ref)) to a `POLAR_VECTORS`
CSV (`alpha, Cd, Cs, Cl, Cm`; alpha in degrees). Non-converged angles (`NaN`) are
skipped, so this works for both NeuralFoil and XFoil sweeps.
"""
function write_polar_csv(filepath::String, sols::Vector{SectionSolution})
    open(filepath, "w") do io
        println(io, "alpha,Cd,Cs,Cl,Cm")
        for s in sols
            isnan(s.cl) && continue
            row = (rad2deg(s.alpha), s.cd, 0.0, s.cl, s.cm)
            println(io, join((@sprintf("%.4f", v) for v in row), ","))
        end
    end
    return filepath
end

"""
    write_polar_matrix_csv(filepath, alpha_range, delta_range, cl, cd, cm)

Write an `(alpha × delta)` sweep to a long-format `POLAR_MATRICES` CSV with columns
`alpha, delta, Cl, Cd, Cm` (both angles in degrees), one row per grid point. The
`delta` column is what marks the file as a matrix polar to the loader. `alpha_range`
and `delta_range` are in radians; `cl`/`cd`/`cm` are `length(alpha) × length(delta)`.
"""
function write_polar_matrix_csv(filepath::String, alpha_range, delta_range,
        cl::AbstractMatrix, cd::AbstractMatrix, cm::AbstractMatrix)
    open(filepath, "w") do io
        println(io, "alpha,delta,Cl,Cd,Cm")
        for (j, d) in enumerate(delta_range), (i, a) in enumerate(alpha_range)
            row = (rad2deg(a), rad2deg(d), cl[i, j], cd[i, j], cm[i, j])
            println(io, join((@sprintf("%.4f", v) for v in row), ","))
        end
    end
    return filepath
end

"""
    DeformedSection

A section deformed by a trailing-edge deflection: its coordinates `(x, y)` and the
[`LeastSquaresFit`](@ref) `kulfan` parameters fitted to them. XFoil consumes `(x, y)`
directly, NeuralFoil consumes `kulfan` — both describe the same shape.
"""
struct DeformedSection
    kulfan::KulfanParameters
    x::Vector{Float64}
    y::Vector{Float64}
end

"""
    deform_section(x, y, delta; crease_frac=0.9, thickness_frac=1.0,
                   flip_thickness_neg=true,
                   wrap_method=ShrinkWrap(clearance=0.0))
        -> DeformedSection

Deform the airfoil coordinates `(x, y)` by trailing-edge deflection `delta` (radians)
about the crease at `crease_frac` along the chord (pivoting through thickness at
`thickness_frac`, 1.0 = top surface), then re-[`shrink_wrap`](@ref) the deflected shape
into clean cosine panels and fit [`LeastSquaresFit`](@ref) Kulfan parameters to it.
XFoil consumes the coordinates directly, NeuralFoil the Kulfan parameters.

`flip_thickness_neg` folds a soft membrane about its lower surface for negative `delta`.
The re-wrap uses zero clearance (it hugs the deflected shape at grid resolution);
the rolling-ball wrap bridges the crease with a `min_concave_radius` fillet instead
of the overlapping panels that XFoil's own repaneling can hit
there. `delta == 0` returns the base coordinates unchanged.
"""
function deform_section(x, y, delta; crease_frac=0.9, thickness_frac=1.0,
                        flip_thickness_neg=true,
                        wrap_method::ShrinkWrap=ShrinkWrap(clearance=0.0))
    xd, yd = collect(float.(x)), collect(float.(y))
    if !iszero(delta)
        pivot = flip_thickness_neg && delta < 0 ? 1 - thickness_frac : thickness_frac
        lower, upper = get_lower_upper(xd, yd, crease_frac)
        turn_trailing_edge!(delta, xd, yd, lower, upper, crease_frac; thickness_frac=pivot)
        xd, yd = shrink_wrap(xd, yd, wrap_method)
    end
    kulfan = fit_kulfan_parameters(xd, yd, LeastSquaresFit())
    return DeformedSection(kulfan, xd, yd)
end

"""
    split_surfaces(x, cp) -> (x_upper, cp_upper, x_lower, cp_lower)

Split a Selig-ordered `(x, cp)` node set at the leading edge (`argmin(x)`) into
upper and lower surfaces, each ordered leading edge → trailing edge (x increasing).
"""
function split_surfaces(x, cp)
    le = argmin(x)
    upper = le:-1:1
    lower = le:length(x)
    return x[upper], cp[upper], x[lower], cp[lower]
end

"""
    resample_x(xs, values, targets) -> Vector{Float64}

Linearly interpolate `values` (given at `xs`) onto `targets`, sorting by `x` and
dropping zero-length segments so the breakpoints are strictly increasing.
"""
function resample_x(xs, values, targets)
    p = sortperm(xs)
    xs, values = xs[p], values[p]
    keep = [true; diff(xs) .> 0]
    itp = linear_interpolation(xs[keep], values[keep]; extrapolation_bc=Line())
    return itp.(targets)
end

"""
    analyze_sweep(solver, def, alpha_range, Re) -> Vector{SectionSolution}

Analyse a deformed section over `alpha_range` (radians). Generic fallback maps
[`analyze_section`](@ref); backends override for efficiency (XFoil reinit sweep,
NeuralFoil vectorization).
"""
function analyze_sweep(solver::AbstractAirfoilSolver, def::DeformedSection, alpha_range, Re)
    return [analyze_section(solver, def, a, Re) for a in alpha_range]
end

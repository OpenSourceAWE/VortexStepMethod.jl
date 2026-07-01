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
    DeformedSection

A section deformed by a trailing-edge deflection: the `EnvelopeFit` `kulfan`
parameters and their repaned coordinates `(x, y)`. NeuralFoil consumes `kulfan`,
XFoil consumes `(x, y)` — both describe the same shape.
"""
struct DeformedSection
    kulfan::KulfanParameters
    x::Vector{Float64}
    y::Vector{Float64}
end

"""
    deform_section(base, delta; crease_frac=0.9, thickness_frac=1.0,
                   flip_thickness_neg=true, fit_method=EnvelopeFit(),
                   n_points=120) -> DeformedSection

Deform a base airfoil (Kulfan parameters) by trailing-edge deflection `delta`
(radians) about the crease axis (`crease_frac` along chord, `thickness_frac`
through thickness — 1.0 = top surface), then refit with `fit_method` and repane.
The refit/repane produces a clean, even panel distribution for both backends.

With `flip_thickness_neg` a negative `delta` mirrors the pivot through thickness
(`thickness_frac -> 1 - thickness_frac`): a soft membrane wing that folds about its
top surface for positive deflection rotates about its lower surface for negative.
"""
function deform_section(base::KulfanParameters, delta;
                        crease_frac=0.9, thickness_frac=1.0, flip_thickness_neg=true,
                        fit_method::KulfanFitMethod=EnvelopeFit(), n_points::Int=120)
    if flip_thickness_neg && delta < 0
        thickness_frac = 1 - thickness_frac
    end
    x, y = kulfan_to_coordinates(base; n_points)
    lower, upper = get_lower_upper(x, y, crease_frac)
    turn_trailing_edge!(delta, x, y, lower, upper, crease_frac; thickness_frac)
    kulfan = fit_kulfan_parameters(x, y, fit_method)
    xr, yr = kulfan_to_coordinates(kulfan; n_points)
    return DeformedSection(kulfan, xr, yr)
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

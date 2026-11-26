# Forward model shared by inverse_design.jl and inverse_design_table.jl: a KiteShape with
# masure-regression sections, its VSM polar, the polar's properties and the inverse.
using LinearAlgebra
using DelimitedFiles
using Statistics
using Optim
using VortexStepMethod
using VortexStepMethod.AirfoilAero: load_masure_model, masure_aero

# Directory holding ET_re5e6.npz: the Re = 5e6 model of
# https://doi.org/10.5281/zenodo.16925758 converted with scripts/export_masure_models.py.
MASURE_MODELS_DIR = get(ENV, "MASURE_MODELS_DIR",
                        joinpath(dirname(@__DIR__), "data", "masure"))
REYNOLDS = 5e6
VA = 20.0
# Section angles of attack [°] the profile polars are tabulated at.
SECTION_ALPHA = collect(-10.0:1.0:25.0)
# Kite angles of attack [°] of every 3D polar, measured to the centre chord.
ALPHA = collect(-8.0:1.0:22.0)
N_SECTIONS = 19
# Sampled designs and their polars, written by inverse_design_table.jl.
POOL_FILE = joinpath(dirname(@__DIR__), "output", "inverse_design_pool.csv")
# Profile parameters held at the centre of the regression's training data.
FIXED_PROFILE = Dict("eta" => 0.23, "delta" => -3.8, "phi" => 0.65, "lambda" => 0.25)
# Design parameters: the KiteShape fields after projected_area, then the centre profile.
DESIGN_PARAMETERS = [
    (key="aspectRatio", label="Aspect ratio", unit="-", lower=2.5, upper=6.0),
    (key="arcAngle", label="Arc angle at the tips", unit="deg", lower=60.0, upper=90.0),
    (key="arcExponent", label="Arc distribution", unit="-", lower=0.8, upper=2.5),
    (key="taper", label="Tip chord / centre chord", unit="-", lower=0.15, upper=0.6),
    (key="sweep", label="Leading-edge sweep", unit="deg", lower=0.0, upper=30.0),
    (key="twist", label="Tip washout", unit="deg", lower=-4.0, upper=8.0),
    (key="tubeDiameter", label="Centre tube diameter / chord", unit="-",
     lower=0.05, upper=0.12),
    (key="camber", label="Centre camber / chord", unit="-", lower=0.04, upper=0.14),
]
LOWER = [p.lower for p in DESIGN_PARAMETERS]
UPPER = [p.upper for p in DESIGN_PARAMETERS]

"Columns of a CSV file with a header row, by name."
function read_columns(path)
    data, header = readdlm(path, ','; header=true)
    return Dict(name => Float64.(data[:, j]) for (j, name) in enumerate(vec(header)))
end

"Linear interpolation at `x` of `values` given at ascending `points`, extended linearly."
function interpolate_linear(points, values, x)
    i = clamp(searchsortedlast(points, x), 1, length(points) - 1)
    weight = (x - points[i]) / (points[i + 1] - points[i])
    return (1 - weight) * values[i] + weight * values[i + 1]
end

"""
The TU Delft V3 ribs and their profiles, pitched so the centre chord is level, the centre
leading edge at the origin.
"""
V3_RIBS = let columns = read_columns(joinpath(dirname(@__DIR__), "data", "TUDELFT_V3_KITE",
        "rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_19ribs.csv"))
    points(edge) = MVec3.(columns["$(edge)_x"], columns["$(edge)_y"], columns["$(edge)_z"])
    le, te = points("LE"), points("TE")
    centre = (length(le) + 1) ÷ 2
    chord = te[centre] - le[centre]
    pitch = atan(-chord[3], chord[1])
    rotation = [cos(pitch) 0 -sin(pitch); 0 1 0; sin(pitch) 0 cos(pitch)]
    level(point) = MVec3(rotation * (point - le[centre]))
    (le=level.(le), te=level.(te), tube=columns["d_tube"], camber=columns["camber"],
     pitch=rad2deg(pitch), centre)
end

"Arc-length fraction of each V3 rib along its leading edge in the y-z plane, 1 at +y."
V3_SPAN_FRACTION = let points = [MVec3(0, le[2], le[3]) for le in V3_RIBS.le]
    length_along = cumsum([0.0; norm.(diff(points))])
    (length_along[V3_RIBS.centre] .- length_along) ./ length_along[V3_RIBS.centre]
end

"Tube diameter and camber over their centre values along the V3 span, at fraction `s`."
function v3_profile_grading(s::Real)
    grading(values) = interpolate_linear(reverse(V3_SPAN_FRACTION), reverse(values), s) /
                      values[V3_RIBS.centre]
    return grading(V3_RIBS.tube), grading(V3_RIBS.camber)
end

"The masure-regression polar of a profile with tube diameter and camber over the chord."
function profile_polar(model, tube, camber)
    params = merge(FIXED_PROFILE, Dict("t" => tube, "kappa" => camber))
    cl, cd, cm = masure_aero(model, params, SECTION_ALPHA)
    return (deg2rad.(SECTION_ALPHA), cl, cd, cm)
end

"""
A VSM wing of `sections` (`le`, `te`, profile `tube` and `camber`) with masure-regression
profiles, two panels per section; mirrored profiles shared.
"""
function masure_wing(model, sections)
    wing = Wing(2 * (length(sections.le) - 1); spanwise_distribution=SPLIT_PROVIDED)
    polars = Dict{NTuple{2, Float64}, Any}()
    for i in eachindex(sections.le)
        profile = round.((sections.tube[i], sections.camber[i]); digits=6)
        polar = get!(() -> profile_polar(model, profile...), polars, profile)
        add_section!(wing, sections.le[i], sections.te[i], POLAR_VECTORS, polar)
    end
    refine!(wing)
    return wing
end

V3_AREA = let wing = Wing(2 * (length(V3_RIBS.le) - 1);
                         spanwise_distribution=SPLIT_PROVIDED)
    foreach((le, te) -> add_section!(wing, le, te, INVISCID), V3_RIBS.le, V3_RIBS.te)
    refine!(wing)
    calculate_projected_area(wing)
end

"The KiteShape of design vector `x` (in `DESIGN_PARAMETERS` order), at the V3's area."
design_shape(x) = KiteShape(; projected_area=V3_AREA, aspect_ratio=x[1], arc_angle=x[2],
    arc_exponent=x[3], taper=x[4], sweep=x[5], twist=x[6])

"Section edges and profile tube diameter and camber of design `x`."
function design_sections(x)
    sections = kite_sections(design_shape(x); n_sections=N_SECTIONS)
    grading = v3_profile_grading.(sections.s)
    return (; sections.le, sections.te, tube=x[7] .* first.(grading),
            camber=x[8] .* last.(grading))
end

design_wing(model, x) = masure_wing(model, design_sections(x))

"""
CL and CD of `wing` at each kite angle of attack in `ALPHA`, each solved from zero
circulation; NaN where VSM does not converge.
"""
function kite_polar(wing)
    body = BodyAerodynamics([wing])
    solver = Solver(wing.n_panels, wing.n_unrefined_sections; aerodynamic_model_type=VSM)
    cl = fill(NaN, length(ALPHA))
    cd = fill(NaN, length(ALPHA))
    for (i, alpha) in enumerate(deg2rad.(ALPHA))
        set_va!(body, apparent_wind(alpha, 0.0, VA))
        cfx, _, cfz = solve!(solver, body, nothing).force_coeffs
        solver.lr.converged || continue
        cl[i] = cfz * cos(alpha) - cfx * sin(alpha)
        cd[i] = cfx * cos(alpha) + cfz * sin(alpha)
    end
    return (; cl, cd)
end

design_polar(model, x) = kite_polar(design_wing(model, x))

"""
Maximum CL, the angle of attack it occurs at [°] (refined by a parabola through its
neighbours), maximum CL/CD and CD where CL crosses zero. NaN where the polar does not show
them inside `ALPHA`.
"""
function polar_properties(polar)
    cl, cd = polar
    peak = argmax(replace(cl, NaN => -Inf))
    stall_angle, cl_max = NaN, cl[peak]
    if 1 < peak < length(ALPHA)
        below, above = cl[peak - 1], cl[peak + 1]
        curvature = below - 2cl_max + above
        shift = (below - above) / (2curvature)
        stall_angle = ALPHA[peak] + shift * (ALPHA[2] - ALPHA[1])
        cl_max -= (below - above) * shift / 4
    end
    lift_to_drag = maximum(filter(!isnan, cl ./ cd); init=-Inf)
    crossing = findfirst(i -> cl[i] <= 0 < cl[i + 1], 1:length(ALPHA) - 1)
    cd_zero_lift = isnothing(crossing) ? NaN :
        interpolate_linear(cl[crossing:crossing + 1], cd[crossing:crossing + 1], 0.0)
    return (; cl_max, stall_angle, lift_to_drag, cd_zero_lift)
end

"Root-mean-square CL and CD differences between two polars."
function polar_misfit(polar, target)
    finite = .!isnan.(polar.cl) .& .!isnan.(target.cl)
    any(finite) || return (cl=Inf, cd=Inf)
    return (cl=sqrt(mean(abs2, polar.cl[finite] - target.cl[finite])),
            cd=sqrt(mean(abs2, polar.cd[finite] - target.cd[finite])))
end

"CL and CD errors weighed equal at these sizes in the inverse."
POLAR_SCALE = (cl=0.01, cd=0.001)

polar_cost(misfit) = (misfit.cl / POLAR_SCALE.cl)^2 + (misfit.cd / POLAR_SCALE.cd)^2

"Design vector from the unit cube the optimiser searches."
from_unit(u) = LOWER .+ clamp.(u, 0, 1) .* (UPPER .- LOWER)
to_unit(x) = (x .- LOWER) ./ (UPPER .- LOWER)

"`cost` of the design at `u` in the unit cube, plus its distance outside the cube."
function unit_cost(model, cost, u)
    x = from_unit(u)
    return cost(x, design_polar(model, x)) + 1e3 * sum(abs2, u .- clamp.(u, 0, 1))
end

"""
The design whose polar best matches `cost(design, polar)`, by Nelder-Mead from `start`
over the parameter ranges.
"""
function invert(model, cost, start; iterations=600)
    result = optimize(u -> unit_cost(model, cost, u), to_unit(start),
                      NelderMead(initial_simplex=Optim.AffineSimplexer(0.0, 0.08)),
                      Optim.Options(iterations=iterations))
    return from_unit(Optim.minimizer(result))
end

"Squared distance of the V3's rib edges from the generator's for shape parameters `x`."
function v3_shape_misfit(x)
    x[3] > 0.3 && 0.02 < x[4] <= 1 || return 1e3
    sections = kite_sections(design_shape(x); s=V3_SPAN_FRACTION)
    return sum(abs2, reduce(vcat, sections.le .- V3_RIBS.le)) +
           sum(abs2, reduce(vcat, sections.te .- V3_RIBS.te))
end

"""
The design vector of the generator closest to the V3 ribs, at the V3's centre profile, and
the RMS distance [m] of the rib edges from it.
"""
function fit_v3_shape()
    start = [3.5, 80.0, 1.5, 0.3, 15.0, 0.0]
    result = optimize(v3_shape_misfit, start, NelderMead(), Optim.Options(iterations=4000))
    rms_distance = sqrt(Optim.minimum(result) / (2 * length(V3_RIBS.le)))
    return [Optim.minimizer(result); V3_RIBS.tube[V3_RIBS.centre];
            V3_RIBS.camber[V3_RIBS.centre]], rms_distance
end

"Write each design vector followed by its CL and CD as one CSV row."
function write_pool(path, designs, polars)
    mkpath(dirname(path))
    rows = [[x; polar.cl; polar.cd]' for (x, polar) in zip(designs, polars)]
    writedlm(path, reduce(vcat, rows), ',')
end

"Design vectors and polars written by `write_pool`."
function read_pool(path)
    data = readdlm(path, ',')
    n_parameters, n_alpha = length(DESIGN_PARAMETERS), length(ALPHA)
    designs = [data[i, 1:n_parameters] for i in axes(data, 1)]
    polars = [(cl=data[i, n_parameters .+ (1:n_alpha)],
               cd=data[i, n_parameters + n_alpha .+ (1:n_alpha)]) for i in axes(data, 1)]
    return designs, polars
end

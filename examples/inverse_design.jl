# Inverse design of an LEI kite: from a target polar back to a geometry. Needs the sampled
# designs inverse_design_table.jl writes, and the masure models inverse_design_setup.jl
# names.
using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using CairoMakie
using VortexStepMethod.AirfoilAero: lei_poly_coeffs
include(joinpath(@__DIR__, "inverse_design_setup.jl"))

FIGURE_DIR = joinpath(dirname(@__DIR__), "docs", "src", "inverse_design")
LITERATURE_DIR = joinpath(dirname(@__DIR__), "data", "TUDELFT_V3_KITE",
                          "literature_results")
# Polars closer than this to each other count as the same polar.
TOLERANCE = (cl=0.01, cd=0.001)
# Sampled kites the inverse is started from to see how unique its answer is.
N_STARTS = 12
# Change of each design parameter that shows which polar property it moves.
STEPS = [2.0, 12.0, 1.0, 0.3, 16.0, 8.0, 0.04, 0.06]
COLORS = (reference=:black, first=:dodgerblue3, second=:orangered3)
PROPERTY_LABELS = ("max CL", "stall [°]", "max CL/CD", "CD zero lift")

model = load_masure_model(REYNOLDS, MASURE_MODELS_DIR)
designs, polars = read_pool(POOL_FILE)
x_v3, fit_distance = fit_v3_shape()
v3_polar = kite_polar(masure_wing(model, V3_RIBS))
v3_properties = polar_properties(v3_polar)

"Cost of a design whose polar should match `target`."
polar_target(target) = (x, polar) -> polar_cost(polar_misfit(polar, target))

"The sampled design with the lowest `cost`, where the inverse starts."
best_sample(cost) = designs[argmin(i -> cost(designs[i], polars[i]), eachindex(designs))]

"Design vectors side by side, one row per parameter, `columns` as `name => vector`."
function parameter_table(columns...; label=p -> p.key, width=10)
    rounded(j) = [string(round(x[j]; sigdigits=3)) for x in last.(columns)]
    header = rpad("", 34) * join(lpad.(first.(columns), width))
    rows = [rpad(label(p), 34) * join(lpad.(rounded(j), width))
            for (j, p) in enumerate(DESIGN_PARAMETERS)]
    return join([header; rows], "\n")
end

print_parameters(columns...) = println(parameter_table(columns...;
                                                       label=p -> "$(p.label) [$(p.unit)]"))

function print_misfit(name, polar, target)
    misfit = polar_misfit(polar, target)
    println("$name: RMS ΔCL = $(round(misfit.cl; sigdigits=2)), ",
            "RMS ΔCD = $(round(misfit.cd; sigdigits=2))")
end

"Front (y-z) and top (y-x, downstream down) views of kites drawn as their section edges."
function draw_kites!(front, top, kites, colors, labels)
    for (kite, color, label) in zip(kites, colors, labels)
        for (axis, (i, j), flip) in ((front, (2, 3), 1), (top, (2, 1), -1))
            lines!(axis, getindex.(kite.le, i), flip .* getindex.(kite.le, j); color, label)
            lines!(axis, getindex.(kite.te, i), flip .* getindex.(kite.te, j); color)
            for (le, te) in zip(kite.le, kite.te)
                lines!(axis, [le[i], te[i]], flip .* [le[j], te[j]]; color, linewidth=0.7)
            end
        end
    end
    front.aspect = top.aspect = DataAspect()
    hidedecorations!(front); hidedecorations!(top)
    hidespines!(front); hidespines!(top)
end

"CL, CD and CL/CD against angle of attack, one axis per cell."
function draw_polars!(cells, polars, colors, labels, linestyles)
    quantities = (("CL", p -> p.cl), ("CD", p -> p.cd), ("CL / CD", p -> p.cl ./ p.cd))
    for (cell, (name, quantity)) in zip(cells, quantities)
        axis = Axis(cell; xlabel="angle of attack [°]", ylabel=name)
        for (polar, color, label, linestyle) in zip(polars, colors, labels, linestyles)
            lines!(axis, ALPHA, quantity(polar); color, label, linestyle, linewidth=2)
        end
    end
end

"One slide: the kites' front and top views beside their polars and a table."
function kite_figure(title, kites, kite_polars, colors, labels, note;
                     linestyles=fill(:solid, length(kites)))
    figure = Figure(size=(1600, 900), fontsize=20)
    Label(figure[0, 1:3], title; fontsize=28, font=:bold, halign=:left)
    front, top = Axis(figure[1, 1]; title="front"), Axis(figure[2, 1]; title="top")
    draw_kites!(front, top, kites, colors, labels)
    axislegend(front; position=:cb, framevisible=false)
    draw_polars!([figure[1, 2], figure[2, 2], figure[1, 3]], kite_polars, colors, labels,
                 linestyles)
    Label(figure[2, 3], note; justification=:left, halign=:left, tellwidth=false,
          fontsize=17, font="DejaVu Sans Mono")
    colsize!(figure.layout, 1, Relative(0.42))
    rowsize!(figure.layout, 1, Auto(1))
    rowsize!(figure.layout, 2, Auto(1))
    return figure
end

mkpath(FIGURE_DIR)

# 0. Which parameters move which property: each stepped by STEPS around the V3 fit.
println(rpad("change across the step", 34), join(rpad.(PROPERTY_LABELS, 14)))
for (j, p) in enumerate(DESIGN_PARAMETERS)
    shifted(sign) = [x_v3[1:j - 1]; x_v3[j] + sign * STEPS[j] / 2; x_v3[j + 1:end]]
    change = map(-, polar_properties(design_polar(model, shifted(1))),
                 polar_properties(design_polar(model, shifted(-1))))
    println(rpad("$(p.label) ± $(STEPS[j] / 2)", 34),
            join(rpad.(round.(collect(change); sigdigits=2), 14)))
end

# 1. Round trip: hand only a polar to the inverse, first the V3 fit's, then the V3's own.
println("Generator fitted to the V3 ribs: RMS distance $(round(fit_distance; digits=3)) m")
fit_polar = design_polar(model, x_v3)
x_self = invert(model, polar_target(fit_polar), best_sample(polar_target(fit_polar)))
x_round = invert(model, polar_target(v3_polar), best_sample(polar_target(v3_polar)))
round_polar = design_polar(model, x_round)
print_parameters("V3 fit" => x_v3, "from fit" => x_self, "from V3" => x_round)
print_misfit("polar of the V3 fit vs the V3", fit_polar, v3_polar)
print_misfit("round trip from the fit vs the fit", design_polar(model, x_self), fit_polar)
print_misfit("round trip vs the V3", round_polar, v3_polar)
save(joinpath(FIGURE_DIR, "round_trip.png"), kite_figure(
    "Round trip: the inverse sees only the V3's polar", [V3_RIBS, design_sections(x_round)],
    [v3_polar, round_polar], [COLORS.reference, COLORS.first],
    ["TU Delft V3", "found from its polar"],
    parameter_table("V3 fit" => x_v3, "found" => x_round); linestyles=[:solid, :dash]))

# 2. How unique: the inverse run again from the sampled kites closest to the V3's polar,
# keeping the kites it returns within TOLERANCE of that polar.
closest = sortperm(polars; by=polar -> polar_cost(polar_misfit(polar, v3_polar)))
starts = designs[closest[1:N_STARTS]]
found = Vector{Vector{Float64}}(undef, N_STARTS)
Threads.@threads for i in 1:N_STARTS
    found[i] = invert(model, polar_target(v3_polar), starts[i])
end
found_polars = [design_polar(model, x) for x in found]
family = findall(found_polars) do polar
    misfit = polar_misfit(polar, v3_polar)
    misfit.cl < TOLERANCE.cl && misfit.cd < TOLERANCE.cd
end
println("$(length(family)) of $N_STARTS inverse runs fly within $(TOLERANCE.cl) in CL ",
        "and $(TOLERANCE.cd) in CD of the V3")
for (j, p) in enumerate(DESIGN_PARAMETERS)
    values = getindex.(found[family], j)
    println(rpad(p.label, 34), round(minimum(values); sigdigits=3), " to ",
            round(maximum(values); sigdigits=3), " of ", p.lower, " to ", p.upper)
end
pairs = [(i, j) for i in family for j in family if i < j]
first_kite, second_kite = argmax(((i, j),) -> norm(to_unit(found[i]) - to_unit(found[j])),
                                 pairs)
print_parameters("first" => found[first_kite], "second" => found[second_kite])
print_misfit("the two kites against each other", found_polars[first_kite],
             found_polars[second_kite])
save(joinpath(FIGURE_DIR, "uniqueness.png"), kite_figure(
    "Two kites, one polar: both within $(TOLERANCE.cl) in CL and $(TOLERANCE.cd) in CD " *
    "of the V3", design_sections.([found[first_kite], found[second_kite]]),
    [found_polars[first_kite], found_polars[second_kite]], [COLORS.first, COLORS.second],
    ["first", "second"],
    parameter_table("first" => found[first_kite], "second" => found[second_kite]);
    linestyles=[:solid, :dash]))

# 3. A design target the V3 does not meet: more lift, stalling later, at no cost in CL/CD
# or zero-lift drag. One-sided misses are charged; past the target the change from the V3
# is charged lightly, so the inverse returns the smallest change that meets it.
target = (cl_max=1.1 * v3_properties.cl_max, stall_angle=v3_properties.stall_angle + 2,
          lift_to_drag=v3_properties.lift_to_drag, cd_zero_lift=v3_properties.cd_zero_lift)
function target_cost(x, polar)
    p = polar_properties(polar)
    all(isfinite, p) || return 1e6
    miss = ((target.cl_max - p.cl_max) / 0.01, (target.stall_angle - p.stall_angle) / 0.2,
            (target.lift_to_drag - p.lift_to_drag) / 0.1,
            (p.cd_zero_lift - target.cd_zero_lift) / 0.001)
    return sum(abs2, max.(miss, 0)) + sum(abs2, to_unit(x) - to_unit(x_v3))
end
x_design = invert(model, target_cost, best_sample(target_cost))
designed_polar = design_polar(model, x_design)
design_properties = polar_properties(designed_polar)
print_parameters("V3 fit" => x_v3, "design" => x_design)

"Polar properties side by side, `columns` as `name => properties`."
property_table(columns...) = join(
    [rpad("", 16) * join(lpad.(first.(columns), 9));
     [rpad(label, 16) * join(lpad.([round(c[k]; sigdigits=3) for c in last.(columns)], 9))
      for (k, label) in zip(keys(target), PROPERTY_LABELS)]], "\n")
println(property_table("V3" => v3_properties, "target" => target,
                       "design" => design_properties))

# 4. A second profile model: the design and the V3 with Breukels' LEI regression sections.
function breukels_polar(sections)
    wing = Wing(2 * (length(sections.le) - 1); spanwise_distribution=SPLIT_PROVIDED)
    for i in eachindex(sections.le)
        add_section!(wing, sections.le[i], sections.te[i], POLY,
                     lei_poly_coeffs(sections.tube[i], sections.camber[i]))
    end
    refine!(wing)
    return kite_polar(wing)
end
println("with Breukels' sections\n", property_table(
    "V3" => polar_properties(breukels_polar(V3_RIBS)),
    "design" => polar_properties(breukels_polar(design_sections(x_design)))))
save(joinpath(FIGURE_DIR, "design.png"), kite_figure(
    "Design target: 10% more lift, stall 2° later, same CL/CD",
    [V3_RIBS, design_sections(x_design)], [v3_polar, designed_polar],
    [COLORS.reference, COLORS.first], ["TU Delft V3", "design"],
    property_table("V3" => v3_properties, "design" => design_properties)))

# 5. An independent check of the forward model: the V3 against RANS CFD and the wind
# tunnel, with profiles from the regression's Re = 1e6 model, near the tests' Reynolds
# numbers. Their angle of attack is measured in the CAD frame, the centre chord pitched up.
v3_low_re = kite_polar(masure_wing(load_masure_model(1e6, MASURE_MODELS_DIR), V3_RIBS))
function literature_polar(file)
    columns = read_columns(joinpath(LITERATURE_DIR, file))
    return (alpha=columns["alpha"] .+ V3_RIBS.pitch, cl=columns["CL"], cd=columns["CD"])
end
measured = [
    ("RANS CFD, Re = 1e6", literature_polar(
        "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv")),
    ("wind tunnel, Re = 5e5", literature_polar(
        "windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv"))]
figure = Figure(size=(1600, 600), fontsize=20)
Label(figure[0, 1:3], "Check: VSM with regression profiles against CFD and wind " *
      "tunnel, TU Delft V3"; fontsize=26, font=:bold, halign=:left)
for (k, (name, quantity)) in enumerate((("CL", p -> p.cl), ("CD", p -> p.cd),
                                       ("CL / CD", p -> p.cl ./ p.cd)))
    axis = Axis(figure[1, k]; xlabel="angle of attack [°]", ylabel=name)
    lines!(axis, ALPHA, quantity(v3_low_re); color=:black, linewidth=2,
           label="VSM, Re = 1e6")
    for ((label, polar), marker) in zip(measured, (:circle, :utriangle))
        scatter!(axis, polar.alpha, quantity(polar); marker, markersize=12, label)
    end
    k == 1 && axislegend(axis; position=:lt, framevisible=false)
end
save(joinpath(FIGURE_DIR, "check.png"), figure)
low_re_stall = polar_properties(v3_low_re).stall_angle
for (label, polar) in measured
    keep = ALPHA[1] .<= polar.alpha .<= low_re_stall
    error_cl = interpolate_linear.(Ref(ALPHA), Ref(v3_low_re.cl), polar.alpha[keep]) -
               polar.cl[keep]
    error_cd = interpolate_linear.(Ref(ALPHA), Ref(v3_low_re.cd), polar.alpha[keep]) -
               polar.cd[keep]
    println("VSM vs $label below stall ($(count(keep)) points): ",
            "mean ΔCL = ", round(mean(error_cl); sigdigits=2),
            ", RMS ΔCL = ", round(sqrt(mean(abs2, error_cl)); sigdigits=2),
            ", mean ΔCD = ", round(mean(error_cd); sigdigits=2),
            ", RMS ΔCD = ", round(sqrt(mean(abs2, error_cd)); sigdigits=2))
end
nothing

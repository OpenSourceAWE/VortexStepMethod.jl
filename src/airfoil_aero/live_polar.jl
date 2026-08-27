"""
    LivePolarSettings(; offsets=deg2rad.(-12.0:3.0:12.0), model_size="xlarge",
                      weights_dir=nothing, n_crit=9.0)

How a live polar is sampled. Every solve, each panel's deformed shape is evaluated at
its reference angle of attack plus each of `offsets`, and those values become the
panel's `SAMPLED` polar directly — there is no fit in between.

Sampling rather than fitting is what lets a panel hold a stall. A polynomial over the
same window averages the knee into a slope and past the peak returns a lift slope of
the wrong sign, which is a divergence the solve cannot recover from; sampled values
are the polar at their own angles and reproduce whatever shape lies between them.

The offsets set both what the polar resolves and how far the solve may move before it
reads a flat end (see [`polar_drift`](@ref)): their spacing is the resolution a stall
knee is caught at, their reach is the room the solve has. The default spacing matches
the 3° grid the offline polar tables are generated on, over a range wide enough that a
tip crossing into stall stays inside it.
"""
@with_kw struct LivePolarSettings
    "Angles of attack [rad] off the reference angle, ascending and straddling zero."
    offsets::Vector{Float64} = deg2rad.(collect(-12.0:3.0:12.0))
    "NeuralFoil network size."
    model_size::String = "xlarge"
    "Directory holding the network weights; `nothing` takes the packaged ones."
    weights_dir::Union{Nothing, String} = nothing
    "Critical amplification factor of the transition model."
    n_crit::Float64 = 9.0
end

"""
    LivePolars(base; settings=LivePolarSettings(), n_stations=60)

Live polar source for a wing whose panels each carry one undeformed airfoil in `base`.
Holds the fixed CST basis, the per-panel reference angles and the network input scratch.
A refresh is dominated by the forward pass itself; the deformation is a matvec against
the constant basis, about half a microsecond a panel. Drive it with
[`refresh_live_polars!`](@ref).
"""
mutable struct LivePolars
    "Sampling configuration."
    settings::LivePolarSettings
    "The fixed CST basis every panel's deflection is projected onto."
    basis::KulfanBasis
    "Undeformed Kulfan parameters per panel."
    base::Vector{KulfanParameters}
    "Deformed Kulfan parameters per panel, rewritten every refresh."
    deformed::Vector{KulfanParameters}
    "Reference angle [rad] the last refresh sampled about, per panel."
    alpha_ref::Vector{Float64}
    "One panel's sample angles [rad], rebuilt per panel inside a refresh."
    knots::Vector{Float64}
    "`25 × (n_panels · n_samples)` network input scratch."
    inputs::Matrix{Float32}
    "`25 × n_panels` network input scratch for the surface-pressure pass."
    pressure_inputs::Matrix{Float32}
    "Lowest confidence NeuralFoil reported over each panel's samples, last refresh."
    confidence::Vector{Float64}
end

function LivePolars(base::AbstractVector{KulfanParameters};
                    settings::LivePolarSettings=LivePolarSettings(),
                    n_stations::Int=60)
    offsets = settings.offsets
    length(offsets) >= 2 || throw(ArgumentError(
        "LivePolars needs at least two sample offsets; got $(length(offsets))."))
    issorted(offsets) || throw(ArgumentError(
        "LivePolars needs ascending sample offsets."))
    first(offsets) <= 0 <= last(offsets) || throw(ArgumentError(
        "LivePolars sample offsets must straddle zero, so the reference angle lies " *
        "inside the sampled range; got $(rad2deg.(extrema(offsets))) deg."))
    n_panels = length(base)
    return LivePolars(settings, KulfanBasis(; n_stations,
                                            n_weights=length(base[1].upper_weights)),
                      collect(base), collect(base), zeros(n_panels),
                      zeros(length(offsets)),
                      zeros(Float32, 25, n_panels * length(offsets)),
                      zeros(Float32, 25, n_panels), zeros(n_panels))
end

"""
    panel_kulfan_parameters(panels; delta=0.0) -> Vector{KulfanParameters}

Fit the undeformed Kulfan parameters of every panel from the surface contour its
`section_aero` carries, the starting point [`LivePolars`](@ref) deforms.

The contour is [`shrink_wrap`](@ref)ped first, the same route the offline polar
generator takes: fitting a raw slice directly returns weights that oscillate by an
order of magnitude more than the deformation ever will, and neighbouring panels then
disagree enough to cost the VSM solve its convergence. The wrap costs ~15 ms a panel,
which is why this is done once at build time and never inside the loop.
"""
function panel_kulfan_parameters(panels; delta=0.0)
    return map(panels) do panel
        panel.section_aero === nothing && throw(ArgumentError(
            "Live polars need a section_aero contour on every panel."))
        x, y, _, _ = section_surface(panel.section_aero, 0.0, delta)
        xw, yw = shrink_wrap(collect(float.(x)), collect(float.(y)), ShrinkWrap())
        return fit_kulfan_parameters(xw, yw, LeastSquaresFit())
    end
end

"""
    polar_drift(live, alpha) -> Float64

How far the largest panel angle of attack has drifted off the reference angle its polar
was sampled about, as a fraction of the reach the samples give it. Above `1` a panel is
being evaluated past the last sample, where the polar is held flat and says nothing
about how the panel is really behaving, so this is the diagnostic a caller's solve loop
reports or refreshes on. Asymmetric offsets are measured by their shorter side.
"""
function polar_drift(live::LivePolars, alpha::AbstractVector)
    length(alpha) == length(live.alpha_ref) || throw(ArgumentError(
        "polar_drift: $(length(alpha)) angles for $(length(live.alpha_ref)) panels."))
    offsets = live.settings.offsets
    reach = min(-first(offsets), last(offsets))
    return maximum(abs.(alpha .- live.alpha_ref)) / reach
end

"""
    deform_live_shapes!(live::LivePolars, deflection) -> Vector{KulfanParameters}

Deform every base airfoil by its camber increment and store the result in
`live.deformed`, which is returned. `nothing` leaves the base shapes in place.
The deformation is an analytic perturbation of one fixed weight vector, so it
never refits and never inherits the non-uniqueness of a fit.
"""
function deform_live_shapes!(live::LivePolars, deflection)
    for i in eachindex(live.base)
        live.deformed[i] = isnothing(deflection) ? live.base[i] :
            deform_kulfan(live.basis, live.base[i], deflection[i])
    end
    return live.deformed
end

"""
    apply_live_shapes!(live::LivePolars, panels; deflection=nothing)

Write each panel's deformed airfoil into `live_shape` without touching its polar,
so the shape a panel reports is the one its current deformation gives. Replaying a
log re-derives the drawn airfoil this way: the frame is never solved, so its polars
would say nothing, and the network pass they cost is the whole price of a refresh.
"""
function apply_live_shapes!(live::LivePolars, panels; deflection=nothing)
    length(panels) == length(live.base) || throw(ArgumentError(
        "apply_live_shapes!: $(length(panels)) panels for $(length(live.base)) " *
        "base airfoils."))
    deform_live_shapes!(live, deflection)
    for (panel, shape) in zip(panels, live.deformed)
        panel.live_shape = shape
    end
    return nothing
end

"""
    refresh_live_polars!(live, panels, alpha_ref, reynolds; deflection=nothing)
        -> Float64

Regenerate every panel's polar from its current shape and write it in as a `SAMPLED`
polar (see [`set_sampled_polar!`](@ref VortexStepMethod.set_sampled_polar!)). Per panel:
deform the base airfoil by `deflection` (a chord-normalized deflection on `live.basis.x`,
or `nothing` to keep the base shape), evaluate NeuralFoil at `alpha_ref .+ offsets`, and
hand those values straight to the panel.

The write is in place — same knot count, same vectors — so a refresh every solve costs
the forward pass and nothing else.

Each panel keeps the deformed shape it was evaluated at as its `live_shape`, so a plot
draws the airfoil the network actually saw.

`alpha_ref` [rad] and `reynolds` are per panel or scalar. Every panel and sample goes
through the network in one forward pass. Returns the lowest analysis confidence over
the batch, a value near zero meaning the deformed shape has left the region the network
was trained on.

Call it after the mesh has been rebuilt from the structure and before the solve: a mesh
rebuild re-seeds each panel's aero from its section, which would drop the polar.
"""
function refresh_live_polars!(live::LivePolars, panels, alpha_ref, reynolds;
                              deflection=nothing)
    n_panels = length(live.base)
    length(panels) == n_panels || throw(ArgumentError(
        "refresh_live_polars!: $(length(panels)) panels for $n_panels base airfoils."))
    per_panel(v) = v isa Number ? fill(float(v), n_panels) : collect(float.(v))
    alpha_vec, re_vec = per_panel(alpha_ref), per_panel(reynolds)
    live.alpha_ref .= alpha_vec

    offsets = live.settings.offsets
    n_samples = length(offsets)
    deform_live_shapes!(live, deflection)
    for i in 1:n_panels
        for k in 1:n_samples
            fill_case_input!(live.inputs, (i - 1) * n_samples + k, live.deformed[i],
                             rad2deg(alpha_vec[i] + offsets[k]), re_vec[i],
                             live.settings.n_crit, 1.0, 1.0)
        end
    end

    model = load_neuralfoil_model(live.settings.model_size;
                                  weights_dir=live.settings.weights_dir)
    cl, cd, cm, confidence = decode_coefficients(fused_output(live.inputs, model))

    for i in 1:n_panels
        samples = ((i - 1) * n_samples + 1):(i * n_samples)
        live.knots .= alpha_vec[i] .+ offsets
        live.confidence[i] = minimum(view(confidence, samples))
        set_sampled_polar!(panels[i], live.knots, view(cl, samples),
                           view(cd, samples), view(cm, samples);
                           shape=live.deformed[i])
    end
    return minimum(confidence)
end

"""
    live_xfoil_solver(live::LivePolars) -> XFoilSolver

XFoil settings matching the ones the live polars were evaluated at: the same critical
amplification factor, and free transition on both surfaces, which is what
[`refresh_live_polars!`](@ref) feeds the network. A comparison run at other settings
measures the settings rather than the network.
"""
live_xfoil_solver(live::LivePolars) =
    XFoilSolver(; ncrit=live.settings.n_crit, xtrip=(1.0, 1.0))

"""
    compare_live_polar(live::LivePolars, panel, panel_idx, alpha, reynolds;
                       solver=live_xfoil_solver(live)) -> NamedTuple

One viscous solve of panel `panel_idx`'s current deformed shape against the live polar
that panel is flying, at the same angle and Reynolds. Returns
`(; panel_idx, alpha, requested, reynolds, confidence, cl, cd, cm)`, each coefficient a
`(live, reference)` pair.

XFoil is marched out from zero in `ramp_step` increments rather than jumped to the
angle, and the comparison is made at the nearest angle it converged at — `alpha` is
that angle and `requested` the one asked for. Both coefficients are read there, so the
pair stays like for like even when the march stops short.

Run it when `confidence` is low. A confidence is the network's opinion of its own
inputs — a shape far from what it was trained on scores badly whether or not the answer
is wrong — so it says to go and check, not what the check will find. `cl` here is read
off the panel's `SAMPLED` polar rather than from a fresh network call, so what is
compared is what the solver actually flew.

A reference that will not converge anywhere on the ramp comes back `NaN` rather than
throwing, so a sweep over panels does not stop at the first one that fails.
"""
function compare_live_polar(live::LivePolars, panel, panel_idx, alpha, reynolds;
                            solver=live_xfoil_solver(live), ramp_step=deg2rad(1.0))
    shape = live.deformed[panel_idx]
    x, y = kulfan_to_coordinates(shape)
    section = DeformedSection(shape, collect(float.(x)), collect(float.(y)))
    angles = xfoil_ramp(alpha, ramp_step)
    solutions = analyze_sweep(solver, section, angles, reynolds)
    usable = [k for k in eachindex(angles) if !isnan(solutions[k].cl)]
    k = isempty(usable) ? length(angles) :
        usable[argmin(abs.(angles[usable] .- alpha))]
    at, reference = angles[k], solutions[k]
    return (; panel_idx, alpha=at, requested=alpha, reynolds,
            confidence=live.confidence[panel_idx],
            cl=(calculate_cl(panel, at), reference.cl),
            cd=(calculate_cd(panel, at), reference.cd),
            cm=(calculate_cm(panel, at), reference.cm))
end

"""
    xfoil_ramp(alpha, step) -> Vector{Float64}

Angles [rad] from zero out to `alpha` in `step` increments, `alpha` itself last. The
viscous march refuses a blunt inflated section jumped straight to its angle, and
`analyze_sweep` reinitialises at zero and walks outward, so handing it the whole ramp
is what gets the reference solved at all.
"""
function xfoil_ramp(alpha, step)
    step = abs(step)
    (step <= 0 || !isfinite(alpha)) && return [float(alpha)]
    angles = collect(range(0.0, float(alpha); step=alpha < 0 ? -step : step))
    isempty(angles) && return [float(alpha)]
    last(angles) ≈ alpha || push!(angles, float(alpha))
    return angles
end

"""
    refresh_live_pressure!(cp, live, contour_x, contour_y, leading_edge, alpha,
                           reynolds) -> Vector{Vector{Float64}}

Fill `cp[i]` with the surface pressure of panel `i`'s current deformed shape at
`alpha[i]`, reconstructed on the contour nodes `(contour_x[i], contour_y[i])` (see
[`contour_pressure`](@ref)). `cp` is written in place and returned.

This is the pressure half of a live polar. [`refresh_live_polars!`](@ref) makes the
panel forces follow the deformed shape; without this the pattern that spreads those
forces over the structure would still come from the undeformed section, so a
deformation would change how hard a panel pulls but not where it pulls.

One batched forward pass over all panels, at the converged angle of attack rather
than at the sampled ones — a panel's `Cp` is wanted at exactly one angle, and
evaluating there is both cheaper than storing the samples and exact. Call it after
the solve has converged, with the contours the traction pattern is indexed on —
their deformed `y`, since the reconstruction runs in arc length and a deformed nose
is a different distance around. Deform the shapes first, which
[`refresh_live_polars!`](@ref) already did for this solve.
"""
function refresh_live_pressure!(cp, live::LivePolars, contour_x, contour_y,
                                leading_edge, alpha, reynolds)
    n_panels = length(live.base)
    length(cp) == length(contour_x) == length(contour_y) ==
        length(leading_edge) == n_panels ||
        throw(ArgumentError("refresh_live_pressure!: $(length(cp)) pressure and " *
            "$(length(contour_x)) contour entries for $n_panels panels."))
    per_panel(v) = v isa Number ? fill(float(v), n_panels) : collect(float.(v))
    alpha_vec, re_vec = per_panel(alpha), per_panel(reynolds)
    for i in 1:n_panels
        fill_case_input!(live.pressure_inputs, i, live.deformed[i],
                         rad2deg(alpha_vec[i]), re_vec[i], live.settings.n_crit,
                         1.0, 1.0)
    end
    model = load_neuralfoil_model(live.settings.model_size;
                                  weights_dir=live.settings.weights_dir)
    station_x, ue_upper, ue_lower = decode_surface_velocity(
        fused_output(live.pressure_inputs, model))
    for i in 1:n_panels
        arc = contour_arc(contour_x[i], contour_y[i], leading_edge[i])
        cp[i] .= contour_pressure(station_x, view(ue_upper, :, i),
                                  view(ue_lower, :, i), contour_x[i], arc,
                                  leading_edge[i])
    end
    return cp
end

"""
    live_surface_friction!(cf, contour_x, reynolds) -> Vector{Vector{Float64}}

Fill `cf[i]` with the skin friction of panel `i`'s contour nodes at `reynolds[i]`,
by the same flat-plate closure ([`flat_plate_cf`](@ref)) the offline tables carry —
NeuralFoil does not predict skin friction. It depends on chord fraction and Reynolds
only, not on the shape or the angle of attack, so the live value differs from the
tabulated one purely by being at the panel's own flight Reynolds instead of the one
the tables were generated at.
"""
function live_surface_friction!(cf, contour_x, reynolds)
    length(cf) == length(contour_x) || throw(ArgumentError(
        "live_surface_friction!: $(length(cf)) friction and " *
        "$(length(contour_x)) contour entries."))
    re_vec = reynolds isa Number ? fill(float(reynolds), length(cf)) :
        collect(float.(reynolds))
    for i in eachindex(cf)
        cf[i] .= (flat_plate_cf(clamp(x, 0.0, 1.0), re_vec[i])
                  for x in contour_x[i])
    end
    return cf
end

"""
    contour_shape_matrix(contour_x, n_weights) -> Matrix{Float64}

The CST shape matrix `C(x)·B(x)` at chord fractions `contour_x`, mapping a change in
Kulfan weights to the normal offset it produces there. Build it once per contour: the
chord fractions of a panel's contour nodes never move, only the weights do.
"""
function contour_shape_matrix(contour_x, n_weights)
    x = clamp.(collect(float.(contour_x)), 0.0, 1.0)
    return class_function(x) .* bernstein_basis(x, n_weights - 1)
end

"""
    live_shape_offset!(offset, live::LivePolars, shape) -> Vector{Vector{Float64}}

Fill `offset[i]` with the normal offset, over chord, that panel `i`'s current
deformation adds to its contour, given the panel's `shape` matrix from
[`contour_shape_matrix`](@ref). Written in place and returned.

A deflection deforms both surfaces by the same camber increment, so one offset serves
the upper and lower halves of a contour alike. Added to a reference contour it gives
the surface the network was actually evaluated on, which is what a traction pattern
has to be draped over for its normals and segment areas to mean anything.
"""
function live_shape_offset!(offset, live::LivePolars, shape)
    length(offset) == length(shape) == length(live.base) || throw(ArgumentError(
        "live_shape_offset!: $(length(offset)) offset and $(length(shape)) shape " *
        "entries for $(length(live.base)) panels."))
    for i in eachindex(live.base)
        mul!(offset[i], shape[i],
             live.deformed[i].upper_weights .- live.base[i].upper_weights)
    end
    return offset
end

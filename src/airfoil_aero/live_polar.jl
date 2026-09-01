"""
    LivePolarSettings(; offsets=deg2rad.(-12.0:3.0:12.0), model_size="xlarge",
                      weights_dir=nothing, n_crit=9.0)

How a live polar is sampled. Every solve, each panel's deformed shape is evaluated at
its reference angle of attack plus each of `offsets`, and those values become the
panel's polar table directly — there is no fit in between.

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
Holds the fixed CST basis, the per-panel reference angles, the network and every buffer a
refresh needs, so a refresh writes through storage that is already there and allocates
nothing at all. A refresh is dominated by the forward pass itself; the deformation is a
matvec against the constant basis, about half a microsecond a panel. Drive it with
[`refresh_live_polars!`](@ref).
"""
mutable struct LivePolars
    "Sampling configuration."
    settings::LivePolarSettings
    "The fixed CST basis every panel's deflection is projected onto."
    basis::KulfanBasis
    "Undeformed Kulfan parameters per panel."
    base::Vector{KulfanParameters}
    "Deformed Kulfan parameters per panel, rewritten in place every refresh."
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
    "The network the settings name, held so a refresh never goes to the model cache."
    model::NeuralFoilModel
    "Forward-pass scratch, sized for the sampled batch and so wide enough for either pass."
    work::NeuralFoilWorkspace
    "Per-panel Reynolds number of the current pass."
    reynolds::Vector{Float64}
    "Per-panel angle [rad] the surface-pressure pass was evaluated at."
    alpha_at::Vector{Float64}
    "Lift coefficient of every sample of the last refresh."
    sample_cl::Vector{Float64}
    "Drag coefficient of every sample of the last refresh."
    sample_cd::Vector{Float64}
    "Moment coefficient of every sample of the last refresh."
    sample_cm::Vector{Float64}
    "Analysis confidence of every sample of the last refresh."
    sample_confidence::Vector{Float64}
    "Upper-surface edge velocity per station and panel, last surface-pressure pass."
    ue_upper::Matrix{Float64}
    "Lower-surface edge velocity per station and panel, last surface-pressure pass."
    ue_lower::Matrix{Float64}
    "Signed arc length of the contour one panel's pressure is being spread over."
    arc::Vector{Float64}
    "Buffers the pressure reconstruction walks from panel to panel through."
    pressure_scratch::ContourPressureScratch
    "The representable part of one panel's deflection."
    residual::Vector{Float64}
    "One panel's weight change off its base shape."
    weight_delta::Vector{Float64}
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
    n_panels, n_samples = length(base), length(offsets)
    n_weights = length(base[1].upper_weights)
    n_cases = n_panels * n_samples
    model = load_neuralfoil_model(settings.model_size; weights_dir=settings.weights_dir)
    work = NeuralFoilWorkspace(model, n_cases)
    n_network_stations = length(work.stations)
    deformed = [KulfanParameters(copy(p.upper_weights), copy(p.lower_weights),
                                 p.leading_edge_weight, p.TE_thickness) for p in base]
    return LivePolars(settings, KulfanBasis(; n_stations, n_weights),
                      collect(base), deformed, zeros(n_panels), zeros(n_samples),
                      zeros(Float32, 25, n_cases), zeros(Float32, 25, n_panels),
                      zeros(n_panels), model, work, zeros(n_panels), zeros(n_panels),
                      zeros(n_cases), zeros(n_cases), zeros(n_cases), zeros(n_cases),
                      zeros(n_network_stations, n_panels),
                      zeros(n_network_stations, n_panels),
                      Float64[], ContourPressureScratch(), zeros(n_stations),
                      zeros(n_weights))
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
    return maximum(i -> abs(alpha[i] - live.alpha_ref[i]),
                   eachindex(alpha, live.alpha_ref)) / reach
end

"""
    deform_live_shapes!(live::LivePolars, deflection) -> Vector{KulfanParameters}

Deform every base airfoil by its camber increment and store the result in
`live.deformed`, which is returned. `nothing` writes the base shapes back. The
deformation is an analytic perturbation of one fixed weight vector, so it never refits
and never inherits the non-uniqueness of a fit.

Each entry keeps the object it already was, its weights overwritten, so a panel holding
one from an earlier refresh follows the current shape and no frame allocates a new one.
"""
function deform_live_shapes!(live::LivePolars, deflection)
    for i in eachindex(live.base)
        shape, base = live.deformed[i], live.base[i]
        if isnothing(deflection)
            shape.upper_weights .= base.upper_weights
            shape.lower_weights .= base.lower_weights
            shape.leading_edge_weight = base.leading_edge_weight
            shape.TE_thickness = base.TE_thickness
        else
            deform_kulfan!(shape, live.basis, base, deflection[i], live.residual)
        end
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

Regenerate every panel's polar table from its current shape and write it in
(see [`set_polar!`](@ref VortexStepMethod.set_polar!)). Per panel:
deform the base airfoil by `deflection` (a chord-normalized deflection on `live.basis.x`,
or `nothing` to keep the base shape), evaluate NeuralFoil at `alpha_ref .+ offsets`, and
hand those values straight to the panel.

The write is in place — same knot count, same vectors — and so is everything around it:
the shapes, the network inputs, both symmetries of the forward pass and the decoded
coefficients all land in storage `live` already holds, so a refresh every solve costs the
forward pass and nothing else.

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
    live.alpha_ref .= alpha_ref
    live.reynolds .= reynolds

    offsets = live.settings.offsets
    n_samples = length(offsets)
    deform_live_shapes!(live, deflection)
    for i in 1:n_panels
        for k in 1:n_samples
            fill_case_input!(live.inputs, (i - 1) * n_samples + k, live.deformed[i],
                             rad2deg(live.alpha_ref[i] + offsets[k]), live.reynolds[i],
                             live.settings.n_crit, 1.0, 1.0)
        end
    end

    decode_coefficients!(live.sample_cl, live.sample_cd, live.sample_cm,
                         live.sample_confidence,
                         fused_output!(live.work, live.inputs, live.model))

    for i in 1:n_panels
        samples = ((i - 1) * n_samples + 1):(i * n_samples)
        live.knots .= live.alpha_ref[i] .+ offsets
        live.confidence[i] = minimum(view(live.sample_confidence, samples))
        set_polar!(panels[i], live.knots, view(live.sample_cl, samples),
                   view(live.sample_cd, samples), view(live.sample_cm, samples);
                   shape=live.deformed[i])
    end
    return minimum(live.sample_confidence)
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
off the panel's polar table rather than from a fresh network call, so what is
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

One batched forward pass over all panels — through the same workspace the polar refresh
uses, whose inputs it has already consumed — at the converged angle of attack rather
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
    live.alpha_at .= alpha
    live.reynolds .= reynolds
    for i in 1:n_panels
        fill_case_input!(live.pressure_inputs, i, live.deformed[i],
                         rad2deg(live.alpha_at[i]), live.reynolds[i],
                         live.settings.n_crit, 1.0, 1.0)
    end
    decode_surface_velocity!(live.ue_upper, live.ue_lower,
                             fused_output!(live.work, live.pressure_inputs, live.model))
    for i in 1:n_panels
        arc = contour_arc!(live.arc, contour_x[i], contour_y[i], leading_edge[i])
        contour_pressure!(cp[i], live.pressure_scratch, live.work.stations,
                          view(live.ue_upper, :, i), view(live.ue_lower, :, i),
                          contour_x[i], arc, leading_edge[i])
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
    for i in eachindex(cf)
        re = reynolds isa Number ? float(reynolds) : float(reynolds[i])
        nodes = contour_x[i]
        @inbounds for k in eachindex(cf[i], nodes)
            cf[i][k] = flat_plate_cf(clamp(nodes[k], 0.0, 1.0), re)
        end
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
        live.weight_delta .= live.deformed[i].upper_weights .- live.base[i].upper_weights
        mul!(offset[i], shape[i], live.weight_delta)
    end
    return offset
end

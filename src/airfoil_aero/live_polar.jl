"""
    LivePolarSettings(; order=2, half_window=deg2rad(4), n_samples=5,
                      model_size="xlarge", weights_dir=nothing, n_crit=9.0)

How a live polar is sampled and fitted. Every solve, each panel's deformed shape is
evaluated at `n_samples` angles of attack spanning `α_ref ± half_window` and a
polynomial of `order` is least-squares fitted through them, becoming the panel's
`TAYLOR` polar.

Order and window belong together. Order 2 over `±4°` tracks a full-polar solve to
within a few tenths of a percent and bends with the stall knee where a straight line
cuts across it; higher orders carry large coefficients that diverge as soon as the
solve steps outside the window, which is why the window is a correctness limit and not
an accuracy knob — see [`polar_drift`](@ref).

Fitting is by least squares over the window, never by finite differences at a tabulated
source's own knot spacing: a piecewise-linear polar has zero curvature inside a segment
and a spike at every knot, so a difference at the knot spacing returns pure
discretization artefact.
"""
@with_kw struct LivePolarSettings
    "Polynomial order in `α - α_ref`; 2 unless the window is narrowed to match."
    order::Int = 2
    "Half width [rad] of the fit window, also the drift the solve may not leave."
    half_window::Float64 = deg2rad(4.0)
    "Angles of attack sampled per panel per refresh."
    n_samples::Int = 5
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
Holds the fixed CST basis, the per-panel expansion points and the network input scratch.
A refresh is dominated by the forward pass itself; the deformation is a matvec against
the constant basis, about half a microsecond a panel. Drive it with
[`refresh_live_polars!`](@ref).
"""
mutable struct LivePolars
    "Sampling and fitting configuration."
    settings::LivePolarSettings
    "The fixed CST basis every panel's deflection is projected onto."
    basis::KulfanBasis
    "Undeformed Kulfan parameters per panel."
    base::Vector{KulfanParameters}
    "Deformed Kulfan parameters per panel, rewritten every refresh."
    deformed::Vector{KulfanParameters}
    "Expansion point [rad] the last fit was built about, per panel."
    alpha_ref::Vector{Float64}
    "Sample offsets [rad] off the expansion point, shared by every panel."
    offsets::Vector{Float64}
    "`(order + 1) × n_samples` least-squares solver for the shared offset grid."
    fit::Matrix{Float64}
    "`25 × (n_panels · n_samples)` network input scratch."
    inputs::Matrix{Float32}
end

function LivePolars(base::AbstractVector{KulfanParameters};
                    settings::LivePolarSettings=LivePolarSettings(),
                    n_stations::Int=60)
    settings.n_samples > settings.order || throw(ArgumentError(
        "LivePolars needs more samples than the fit order; got " *
        "$(settings.n_samples) and $(settings.order)."))
    offsets = collect(range(-settings.half_window, settings.half_window,
                            settings.n_samples))
    vandermonde = [d^(k - 1) for d in offsets, k in 1:(settings.order + 1)]
    n_panels = length(base)
    return LivePolars(settings, KulfanBasis(; n_stations,
                                            n_weights=length(base[1].upper_weights)),
                      collect(base), collect(base), zeros(n_panels), offsets,
                      pinv(vandermonde),
                      zeros(Float32, 25, n_panels * settings.n_samples))
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

How far the largest panel angle of attack has drifted off the expansion point its
polar was fitted about, as a fraction of the fit window. Above `1` a panel is being
evaluated outside the window and the fit must be rebuilt before its answer is used, so
this is the guard [`refresh_live_polars!`](@ref) leaves to the caller's solve loop.
"""
function polar_drift(live::LivePolars, alpha::AbstractVector)
    length(alpha) == length(live.alpha_ref) || throw(ArgumentError(
        "polar_drift: $(length(alpha)) angles for $(length(live.alpha_ref)) panels."))
    return maximum(abs.(alpha .- live.alpha_ref)) / live.settings.half_window
end

"""
    refresh_live_polars!(live, panels, alpha_ref, reynolds; deflection=nothing)
        -> Float64

Refit every panel's polar from its current shape and write it in as a `TAYLOR` polar
(see [`set_taylor_polar!`](@ref VortexStepMethod.set_taylor_polar!)). Per panel:
deform the base airfoil by `deflection` (a chord-normalized deflection on
`live.basis.x`, or `nothing` to keep the base shape), evaluate NeuralFoil at
`alpha_ref ± half_window`, and least-squares fit the polynomial.

Each panel keeps the deformed shape it was evaluated at as its `live_shape`, so a plot
draws the airfoil the network actually saw.

`alpha_ref` [rad] and `reynolds` are per panel or scalar. Every panel and sample goes
through the network in one forward pass. Returns the lowest analysis confidence over
the batch, a value near zero meaning the deformed shape has left the region the network
was trained on.

Call it after the mesh has been rebuilt from the structure and before the solve: a mesh
rebuild re-seeds each panel's aero from its section, which would drop the fit.
"""
function refresh_live_polars!(live::LivePolars, panels, alpha_ref, reynolds;
                              deflection=nothing)
    n_panels = length(live.base)
    length(panels) == n_panels || throw(ArgumentError(
        "refresh_live_polars!: $(length(panels)) panels for $n_panels base airfoils."))
    per_panel(v) = v isa Number ? fill(float(v), n_panels) : collect(float.(v))
    alpha_vec, re_vec = per_panel(alpha_ref), per_panel(reynolds)
    live.alpha_ref .= alpha_vec

    n_samples = live.settings.n_samples
    for i in 1:n_panels
        live.deformed[i] = isnothing(deflection) ? live.base[i] :
            deform_kulfan(live.basis, live.base[i], deflection[i])
        for k in 1:n_samples
            fill_case_input!(live.inputs, (i - 1) * n_samples + k, live.deformed[i],
                             rad2deg(alpha_vec[i] + live.offsets[k]), re_vec[i],
                             live.settings.n_crit, 1.0, 1.0)
        end
    end

    model = load_neuralfoil_model(live.settings.model_size;
                                  weights_dir=live.settings.weights_dir)
    cl, cd, cm, confidence = decode_coefficients(fused_output(live.inputs, model))

    for i in 1:n_panels
        samples = ((i - 1) * n_samples + 1):(i * n_samples)
        set_taylor_polar!(panels[i], alpha_vec[i], live.fit * cl[samples],
                          live.fit * cd[samples], live.fit * cm[samples];
                          window=live.settings.half_window, shape=live.deformed[i])
    end
    return minimum(confidence)
end

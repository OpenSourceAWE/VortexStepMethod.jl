"""Length [`smooth_norm`](@ref) adds in quadrature to keep a norm positive at zero."""
const SMOOTH_FLOOR = 1e-12

"""
    smooth_norm(vec, floor=SMOOTH_FLOOR)

Norm floored in quadrature: positive and differentiable at the origin, and
branch-free, which is what makes the functions below traceable.
"""
@inline smooth_norm(vec, floor=SMOOTH_FLOOR) = sqrt(sum(abs2, vec) + floor^2)
@inline smooth_norm(x::Number, floor=SMOOTH_FLOOR) = sqrt(x^2 + floor^2)

"""
    panel_span_vector(le_1, te_1, le_2, te_2)

Quarter-chord vector from section 2 to section 1: the panel's spanwise axis, its
length the span width.
"""
@inline panel_span_vector(le_1, te_1, le_2, te_2) =
    (0.75 .* le_1 .+ 0.25 .* te_1) .- (0.75 .* le_2 .+ 0.25 .* te_2)

@inline panel_chord(le_1, te_1, le_2, te_2) =
    0.5 * (smooth_norm(te_1 .- le_1) + smooth_norm(te_2 .- le_2))

"""
    panel_chord_weight(width_prev, width_own, width_next)

Section 1's share of the edge blend setting a panel's chord direction, weighted by
panel spacing so a panel leans towards a narrower neighbour. `nothing` for a
missing neighbour at a tip; a lone panel gets `0.5`.
"""
@inline function panel_chord_weight(width_prev, width_own, width_next)
    share = if isnothing(width_prev) && isnothing(width_next)
        one(width_own) / 2
    elseif isnothing(width_prev)
        width_own / (width_own + width_next)
    elseif isnothing(width_next)
        width_prev / (width_prev + width_own)
    else
        0.25 * (width_prev / (width_prev + width_own) +
                width_own / (width_own + width_next) + 1)
    end
    return 1 - share
end

"""
    panel_axes(le_1, te_1, le_2, te_2, chord_weight=0.5, orient=1)

Airfoil frame and size of the panel between two sections, as
`(; x_airf, y_airf, z_airf, chord, width)`. `x_airf` is chordwise, `y_airf` runs
along the quarter-chord line the bound vortex sits on — not square to `x_airf` on
a swept panel — and `z_airf` is `x_airf × y_airf` normalised. `chord_weight`
([`panel_chord_weight`](@ref)) is section 1's share of the chord-direction blend;
`orient` is `±1` and flips `y_airf` and `z_airf` together.
"""
@inline function panel_axes(le_1, te_1, le_2, te_2, chord_weight=0.5, orient=1)
    lean = chord_weight - 0.5
    chord_vec = (0.5 .* (te_1 .+ te_2) .- 0.5 .* (le_1 .+ le_2)) .+
        lean .* ((te_1 .- te_2) .- (le_1 .- le_2))
    span_vec = panel_span_vector(le_1, te_1, le_2, te_2)
    width = smooth_norm(span_vec)
    x_airf = chord_vec ./ smooth_norm(chord_vec)
    y_airf = orient .* (span_vec ./ width)
    z_cross = cross(x_airf, span_vec)
    z_airf = orient .* (z_cross ./ smooth_norm(z_cross))
    return (; x_airf, y_airf, z_airf,
            chord=panel_chord(le_1, te_1, le_2, te_2), width)
end

"""
    section_pitch_rate(delta_va, z_airf, chord)
    section_pitch_rate(velocity_leading, velocity_trailing, z_airf, chord)

Rate a section rotates about its own spanwise axis, positive nose-up. `delta_va`
is the trailing minus leading edge apparent wind; apparent wind is
`wind - velocity`, so that is the leading minus trailing edge velocity, hence the
reversed order in the four-argument form. Chordwise wind variation enters here
too. Builds a `pitch_rate_dist` for [`set_va!`](@ref) on a deforming structure,
where no single body rate describes every section.
"""
@inline section_pitch_rate(delta_va, z_airf, chord) =
    ifelse(chord > 0, dot(delta_va, z_airf) / smooth_norm(chord), zero(chord))
@inline section_pitch_rate(velocity_leading, velocity_trailing, z_airf, chord) =
    section_pitch_rate(velocity_leading .- velocity_trailing, z_airf, chord)

"""
    flow_curvature_cm(pitch_rate, chord, v_rel)

Thin-airfoil quarter-chord moment increment of a section pitching about its own
spanwise axis: `Δcm = -(π/4)·q̂` with `q̂ = q·c/(2·v_rel)`, `q` positive nose-up.
Pivot-independent, and lift needs no matching correction because the inflow is
already sampled at three-quarter chord.
"""
@inline flow_curvature_cm(pitch_rate, chord, v_rel) =
    ifelse(v_rel > 0, -0.25π * pitch_rate * chord / (2 * smooth_norm(v_rel)),
           zero(pitch_rate * chord))

"""
    effective_alpha(alpha, deficiency)

Angle the polars are read at: the geometric inflow angle less an unsteady lag. The
geometric angle still turns the force, so a lag shifts the coefficients only.
"""
@inline effective_alpha(alpha, deficiency) = alpha - deficiency

"""
    panel_inflow(axes, va_1, va_2, v_ind, dva_1=nothing, dva_2=nothing,
                 deficiency=0)

Flow a panel sees, as `(; v_eff, alpha, alpha_eff, v_span, pitch_rate)`, where
`v_span` is the effective velocity across the span. `axes` is a
[`panel_axes`](@ref) result. `dva_1`/`dva_2` are the sections' trailing minus
leading edge apparent wind, giving the [`section_pitch_rate`](@ref); `nothing`
leaves it zero. `deficiency` feeds [`effective_alpha`](@ref).
"""
@inline function panel_inflow(axes, va_1, va_2, v_ind, dva_1=nothing,
                              dva_2=nothing, deficiency=0)
    (; x_airf, y_airf, z_airf, chord) = axes
    v_eff = 0.5 .* (va_1 .+ va_2) .+ v_ind
    alpha = atan(dot(v_eff, z_airf), dot(v_eff, x_airf))
    v_span = cross(v_eff, y_airf)
    pitch_rate = isnothing(dva_1) ? zero(chord) :
        section_pitch_rate(0.5 .* (dva_1 .+ dva_2), z_airf, chord)
    return (; v_eff, alpha, alpha_eff=effective_alpha(alpha, deficiency),
            v_span, pitch_rate)
end

"""
    dynamic_pressure(rho_1, rho_2, v_span)

Panel dynamic pressure from its two section densities and the
[`panel_inflow`](@ref) span-wise velocity, as that vector or its magnitude.
"""
@inline dynamic_pressure(rho_1, rho_2, v_span::AbstractVector) =
    0.25 * (rho_1 + rho_2) * dot(v_span, v_span)
@inline dynamic_pressure(rho_1, rho_2, v_rel::Number) =
    0.25 * (rho_1 + rho_2) * v_rel^2

"""
    panel_force_directions(axes, alpha_dir, spanwise)

Lift and drag unit vectors, `(; dir_lift, dir_drag)`. `alpha_dir` is the angle that
turns the force, `spanwise` the wing's spanwise direction.
"""
@inline function panel_force_directions(axes, alpha_dir, spanwise)
    (; x_airf, y_airf, z_airf) = axes
    dir_inflow = cos(alpha_dir) .* x_airf .+ sin(alpha_dir) .* z_airf
    lift_cross = cross(dir_inflow, y_airf)
    dir_lift = lift_cross ./ smooth_norm(lift_cross)
    drag_cross = cross(spanwise, dir_lift)
    return (; dir_lift, dir_drag=drag_cross ./ smooth_norm(drag_cross))
end

"""Panel pitching moment per unit span about `y_airf`, positive nose-up."""
@inline panel_moment(cm, q_dyn, chord) = cm * q_dyn * chord^2

"""
    panel_couple_force(cm, q_dyn, chord, width, scale=1)

[`panel_moment`](@ref) as a force couple: the magnitude two opposed forces one
chord apart need to produce it. A particle model places the moment this way, along
the panel normal at its leading and trailing edge.
"""
@inline panel_couple_force(cm, q_dyn, chord, width, scale=1) =
    scale * width * cm * q_dyn * chord

"""
    panel_loads(axes, dirs, q_dyn, cl, cd, cm, scale=1)

Panel load from its polar coefficients, [`panel_axes`](@ref) and
[`panel_force_directions`](@ref), as `(; lift, drag, moment, force,
pitching_moment)`. The first three are per unit span; `force` and
`pitching_moment` are the whole panel's, `scale` included.
"""
@inline function panel_loads(axes, dirs, q_dyn, cl, cd, cm, scale=1)
    (; chord, width) = axes
    lift = cl * q_dyn * chord
    drag = cd * q_dyn * chord
    moment = panel_moment(cm, q_dyn, chord)
    force = (scale * width) .* (lift .* dirs.dir_lift .+ drag .* dirs.dir_drag)
    return (; lift, drag, moment, force, pitching_moment=scale * width * moment)
end

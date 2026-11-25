"""
    KiteShape(; projected_area=19.75, aspect_ratio=3.5, arc_angle=80, arc_exponent=1.6,
                taper=0.3, sweep=18, twist=0)

Planform and arc of a leading-edge-inflatable kite, in the terms a kite designer uses.
[`kite_sections`](@ref) turns it into section leading and trailing edges.

The span is walked by `s`, the arc length from the centre over the half arc length,
1 at the +y tip and -1 at the -y tip. The arc's tangent leans
`arc_angle * |s|^arc_exponent` from the horizontal, the chord is
`c0 * (taper + (1 - taper) * sqrt(1 - s^2))`, the leading edge moves downstream by
`tan(sweep) * b/2 * s^2` and the section turns nose down by `twist * |s|` about its
leading edge. Span and chord are scaled so the untwisted
kite has the projected area and aspect ratio asked for.
"""
Base.@kwdef struct KiteShape
    "area projected on the horizontal plane [m²]"
    projected_area::Float64 = 19.75
    "projected span squared over projected area [-]"
    aspect_ratio::Float64 = 3.5
    "angle of the arc's tangent from the horizontal at the tips [°]"
    arc_angle::Float64 = 80.0
    "how the arc is spread over the span: 1 bends evenly, more keeps the centre flat [-]"
    arc_exponent::Float64 = 1.6
    "tip chord over centre chord [-]"
    taper::Float64 = 0.3
    "angle from the centre leading edge to the tip leading edge, seen from above [°]"
    sweep::Float64 = 18.0
    "nose-down rotation of the tip sections against the centre [°]"
    twist::Float64 = 0.0
end

"Angle [rad] of the arc's tangent from the horizontal at span fraction `s`."
tangent_angle(shape::KiteShape, s::Real) =
    deg2rad(shape.arc_angle) * abs(s)^shape.arc_exponent

"Arc position (y, z) at span fraction `s` for a half arc length of 1, by Simpson's rule."
function arc_position(shape::KiteShape, s::Real; n_intervals=64)
    h = abs(s) / n_intervals
    y = z = 0.0
    for k in 0:n_intervals
        weight = k == 0 || k == n_intervals ? 1 : (isodd(k) ? 4 : 2)
        y += weight * cos(tangent_angle(shape, k * h))
        z += weight * sin(tangent_angle(shape, k * h))
    end
    return sign(s) * y * h / 3, -z * h / 3
end

"Chord at span fraction `s` over the centre chord."
chord_fraction(shape::KiteShape, s::Real) =
    shape.taper + (1 - shape.taper) * sqrt(max(1 - s^2, 0.0))

"""
    kite_sections(shape::KiteShape; n_sections=19, s=range(1, -1; length=n_sections))
        -> (; le, te, s)

Leading and trailing edges [m] of the sections at span fractions `s`, by default
`n_sections` evenly spaced in arc length from the +y tip to the -y tip, as `add_section!`
takes them. x points downstream, z up, the centre leading edge at the origin.
"""
function kite_sections(shape::KiteShape; n_sections::Int=19,
                       s=range(1.0, -1.0; length=n_sections))
    half_arc_span, _ = arc_position(shape, 1.0)
    fine_s = range(0.0, 1.0; length=201)
    projected_chord = [chord_fraction(shape, u) * cos(tangent_angle(shape, u))
                       for u in fine_s]
    chord_integral = step(fine_s) *
        (sum(projected_chord) - (projected_chord[1] + projected_chord[end]) / 2)
    span = sqrt(shape.aspect_ratio * shape.projected_area)
    half_arc_length = span / (2 * half_arc_span)
    centre_chord = shape.projected_area / (2 * half_arc_length * chord_integral)
    le = MVec3[]
    te = MVec3[]
    for fraction in s
        y, z = arc_position(shape, fraction)
        leading_edge = MVec3(tan(deg2rad(shape.sweep)) * span / 2 * fraction^2,
                             half_arc_length * y, half_arc_length * z)
        theta = tangent_angle(shape, fraction)
        normal = MVec3(0.0, sign(fraction) * sin(theta), cos(theta))
        twist = deg2rad(shape.twist) * abs(fraction)
        chord = centre_chord * chord_fraction(shape, fraction)
        push!(le, leading_edge)
        push!(te, leading_edge +
                  chord * (cos(twist) * MVec3(1, 0, 0) + sin(twist) * normal))
    end
    return (; le, te, s=collect(s))
end

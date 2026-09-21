
const KINEMATIC_VISCOSITY = 1.460e-5 # [m²/s] for the atmosphere at sea level.
                                     # see: https://en.wikipedia.org/wiki/Reynolds_number

"""
    turn_trailing_edge!(angle, x, y, lower_turn, upper_turn, crease_frac; thickness_frac=nothing)

Deflect airfoil trailing edge by rotating coordinates behind the crease line.
Positive angle deflects downward.

The pivot's through-thickness location is set by `thickness_frac` (0 = bottom
surface, 0.5 = mid, 1 = top): `y_pivot = lower_turn + thickness_frac*(upper_turn -
lower_turn)`. When `thickness_frac` is given, points are only rotated (no crease
cleanup) — intended for a following Kulfan refit that wraps the crease robustly.
When `thickness_frac === nothing` (legacy), the pivot is the top surface for
downward deflection and the bottom for upward, and folded-over crease points are
removed and the band averaged so the coordinates stay XFoil-paneable.
"""
function turn_trailing_edge!(angle, x, y, lower_turn, upper_turn, crease_frac;
                             thickness_frac=nothing)
    legacy = thickness_frac === nothing
    turn_distance = upper_turn - lower_turn
    flap_sign = angle > 0 ? 1 : -1
    y_turn = legacy ? (angle > 0 ? upper_turn : lower_turn) :
                      lower_turn + thickness_frac * turn_distance
    smooth_idx = Int[]
    rm_idx = Int[]
    for i in eachindex(x)
        if legacy
            if crease_frac - turn_distance < x[i] < crease_frac + turn_distance &&
               flap_sign * y[i] > 0
                push!(smooth_idx, i)
            elseif flap_sign * y[i] < 0 && crease_frac > x[i] > crease_frac - turn_distance
                push!(rm_idx, i)
            end
        end
        if x[i] > crease_frac
            x_rel = x[i] - crease_frac
            y_rel = y[i] - y_turn
            x[i] = crease_frac + x_rel * cos(angle) + y_rel * sin(angle)
            y[i] = y_turn - x_rel * sin(angle) + y_rel * cos(angle)
            if legacy && angle > 0 && x[i] < crease_frac - turn_distance/2 && y[i] > lower_turn
                push!(rm_idx, i)
            elseif legacy && angle < 0 && x[i] < crease_frac - turn_distance/2 &&
                   y[i] < upper_turn
                push!(rm_idx, i)
            end
        end
    end

    legacy || return nothing
    lower_i, upper_i = minimum(smooth_idx), maximum(smooth_idx)
    for i in smooth_idx
        window = min(i - lower_i + 1, upper_i - i + 1)
        x[i] = mean(x[i-window:i+window])
    end
    deleteat!(x, rm_idx)
    deleteat!(y, rm_idx)
    nothing
end


"""
    get_lower_upper(x, y, crease_frac) -> (lower, upper)

Heights of the lower and upper surfaces of the closed contour `(x, y)` at the chord
station `x = crease_frac`: the lowest and highest crossing of that vertical line with
the contour's segments. Throws an `ArgumentError` when the contour crosses it fewer
than twice.
"""
function get_lower_upper(x, y, crease_frac)
    lower, upper = Inf, -Inf
    n_crossings = 0
    for i in firstindex(x):lastindex(x)-1
        min(x[i], x[i+1]) <= crease_frac < max(x[i], x[i+1]) || continue
        t = (crease_frac - x[i]) / (x[i+1] - x[i])
        y_crossing = y[i] + t * (y[i+1] - y[i])
        lower, upper = min(lower, y_crossing), max(upper, y_crossing)
        n_crossings += 1
    end
    n_crossings >= 2 || throw(ArgumentError("the contour crosses x = $crease_frac " *
        "$n_crossings times; a closed one crosses it twice"))
    return lower, upper
end


const SPEED_OF_SOUND = 343 # [m/s] at 20 °C, see: https://en.wikipedia.org/wiki/Speed_of_sound
const KINEMATIC_VISCOSITY = 1.460e-5 # [m²/s] for the atmosphere at sea level.
                                     # see: https://en.wikipedia.org/wiki/Reynolds_number

"""
    normalize_foil!(x, y)

Scale airfoil coordinates to unit chord length, with x ∈ [0,1].
"""
function normalize_foil!(x, y)
    x_min = minimum(x)
    x_max = maximum(x)
    for i in eachindex(x)
        x[i] = (x[i] - x_min) / (x_max - x_min)
        y[i] = (y[i] - x_min) / (x_max - x_min)
    end
end

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
    solve_alpha!(cls, cds, cms, alpha_range, alpha_idxs, delta, re, x_, y_, 
                lower, upper, kite_speed, speed_of_sound, crease_frac)

Calculate aerodynamic coefficients for specified angles of attack using XFoil.
Results are stored in pre-allocated cls, cds, cms vectors.
"""
function solve_alpha!(cls, cds, cms, alpha_range, alpha_idxs, delta, re, x_, y_, lower, upper, kite_speed, speed_of_sound, crease_frac)
    x = deepcopy(x_)
    y = deepcopy(y_)
    turn_trailing_edge!(delta, x, y, lower, upper, crease_frac)
    Xfoil.set_coordinates(x, y)
    Xfoil.pane(npan=140)
    reinit = true
    for (alpha, alpha_idx) in zip(alpha_range, alpha_idxs)
        # Solve for the given angle of attack
        cl, cd, _, cm, converged = Xfoil.solve_alpha(rad2deg(alpha), re; iter=50, reinit=reinit, mach=kite_speed/speed_of_sound, xtrip=(0.05, 0.05))
        reinit = false
        if converged
            cls[alpha_idx] = cl
            cds[alpha_idx] = cd
            cms[alpha_idx] = cm
        end
    end
    return nothing
end

"""
    run_solve_alpha(alpha_range, delta, re, x_, y_, lower, upper, kite_speed, 
                   speed_of_sound, crease_frac) -> (cls, cds, cms)

Run XFoil analysis for full alpha range, handling positive and negative angles separately.
"""
function run_solve_alpha(alpha_range, delta, re, x_, y_, lower, upper, kite_speed, speed_of_sound, crease_frac)
    @info "solving alpha with trailing edge angle: $(rad2deg(delta)) degrees"
    cls = Float64[NaN for _ in alpha_range]
    cds = Float64[NaN for _ in alpha_range]
    cms = Float64[NaN for _ in alpha_range]
    neg_idxs = sort(findall(alpha_range .< 0.0), rev=true)
    neg_alpha_range = alpha_range[neg_idxs]
    pos_idxs = sort(findall(alpha_range .>= 0.0))
    pos_alpha_range = alpha_range[pos_idxs]
    
    solve_alpha!(cls, cds, cms, neg_alpha_range, neg_idxs, delta, 
                        re, x_, y_, lower, upper, kite_speed, speed_of_sound, crease_frac)
    solve_alpha!(cls, cds, cms, pos_alpha_range, pos_idxs, delta, 
                        re, x_, y_, lower, upper, kite_speed, speed_of_sound, crease_frac)
    return cls, cds, cms
end

"""
    get_lower_upper(x, y, crease_frac) -> (lower, upper)

Find y-coordinates where upper/lower surfaces intersect the hinge line.
"""
function get_lower_upper(x, y, crease_frac)
    lower_trailing_edge = 0.0
    upper_trailing_edge = 0.0
    min_lower_distance = Inf
    min_upper_distance = Inf
    for (xi, yi) in zip(x, y)
        if yi < 0
            lower_distance = abs(xi - crease_frac)
            if lower_distance < min_lower_distance
                min_lower_distance = lower_distance
                lower_trailing_edge = yi
            end
        else
            upper_distance = abs(xi - crease_frac)
            if upper_distance < min_upper_distance
                min_upper_distance = upper_distance
                upper_trailing_edge = yi
            end
        end
    end
    return lower_trailing_edge, upper_trailing_edge
end

"""
    create_polars(; dat_path, cl_polar_path, cd_polar_path, cm_polar_path, wind_vel, area, 
                 width, crease_frac, alpha_range, delta_range, remove_nan=true)

Generate complete aerodynamic coefficient tables for an airfoil with trailing edge deflection.

This function:
1. Reads airfoil coordinates from a .dat file
2. Normalizes the coordinates to unit chord
3. Computes Reynolds number based on normal flight conditions
4. Generates coefficient matrices using XFoil for combinations of:
   - Angles of attack (alpha_range)
   - Trailing edge deflections (delta_range)
5. Optionally interpolates any NaN values from failed computations
6. Saves results to CSV files

# Arguments
- `dat_path`: Path to airfoil coordinate file
- `cl_polar_path`: Output path for lift coefficient CSV
- `cd_polar_path`: Output path for drag coefficient CSV
- `cm_polar_path`: Output path for moment coefficient CSV
- `wind_vel`: Reference velocity (m/s)
- `area`: Wing area (m²)
- `width`: Wing span (m)
- `crease_frac`: Chordwise location of hinge line (0-1)
- `alpha_range`: Vector of angles of attack to analyze (radians)
- `delta_range`: Vector of trailing edge deflections to analyze (radians)
- `remove_nan`: Whether to interpolate NaN values in results (default: true)

# Outputs
Creates three CSV files containing lift, drag, and moment coefficient matrices.
Each file includes headers with angle labels and uses degrees for readability.
"""
function create_polars(; dat_path, cl_polar_path, cd_polar_path, cm_polar_path, wind_vel, area, 
    width, crease_frac, alpha_range, delta_range, remove_nan=true
)
    @info "Creating polars. This can take several minutes."
    tic()

    cl_matrix = zeros(length(alpha_range), length(delta_range))
    cd_matrix = zeros(length(alpha_range), length(delta_range))
    cm_matrix = zeros(length(alpha_range), length(delta_range))

    kite_speed = wind_vel
    chord_length = area / width
    local reynolds_number = kite_speed * chord_length / KINEMATIC_VISCOSITY # https://en.wikipedia.org/wiki/Reynolds_number

    # Read airfoil coordinates from a file.
    local x = Float64[]
    local y = Float64[]
    f = open(dat_path, "r")
    try
        for line in eachline(f)
            entries = split(chomp(line))
            try
                push!(x, parse(Float64, entries[1]))
                push!(y, parse(Float64, entries[2]))
            catch err
                if err isa ArgumentError
                    println(entries)
                else
                    rethrow(err)
                end
            end
        end
    finally
        close(f)
    end
    normalize_foil!(x, y)
    Xfoil.set_coordinates(x, y)
    x, y = Xfoil.pane(npan=140)
    lower, upper = get_lower_upper(x, y, crease_frac)

    for j in eachindex(delta_range)
        cl_matrix[:, j], cd_matrix[:, j], cm_matrix[:, j] = run_solve_alpha(alpha_range, delta_range[j], 
                        reynolds_number, x, y, lower, upper, kite_speed, SPEED_OF_SOUND, crease_frac)
    end
    cl_matrix = Matrix{Float64}(cl_matrix)
    cd_matrix = Matrix{Float64}(cd_matrix)
    cm_matrix = Matrix{Float64}(cm_matrix)

    println("Generated lift matrix:")
    display(cl_matrix)
    println("Generated drag matrix:")
    display(cd_matrix)
    println("Generated moment matrix:")
    display(cm_matrix)
    
    @info "Relative trailing_edge height: $(upper - lower)"
    @info "Reynolds number for flying speed of $kite_speed is $reynolds_number"

    if remove_nan
        interpolate_matrix_nans!(cl_matrix)
        interpolate_matrix_nans!(cd_matrix)
        interpolate_matrix_nans!(cm_matrix)
    end
    
    write_aero_matrix(cl_polar_path, cl_matrix, alpha_range, delta_range, "C_l")
    write_aero_matrix(cd_polar_path, cd_matrix, alpha_range, delta_range, "C_d")
    write_aero_matrix(cm_polar_path, cm_matrix, alpha_range, delta_range, "C_m")
        
    toc()
end


const SPEED_OF_SOUND = 343 # [m/s] at 20 °C, see: https://en.wikipedia.org/wiki/Speed_of_sound
const KINEMATIC_VISCOSITY = 1.460e-5 # [m²/s] for the atmosphere at sea level.
                                     # see: https://en.wikipedia.org/wiki/Reynolds_number

function normalize_foil!(x, y)
    x_min = minimum(x)
    x_max = maximum(x)
    for i in eachindex(x)
        x[i] = (x[i] - x_min) / (x_max - x_min)
        y[i] = (y[i] - x_min) / (x_max - x_min)
    end
end

function turn_trailing_edge!(angle, x, y, lower_turn, upper_turn, crease_frac)
    turn_distance = upper_turn - lower_turn
    smooth_idx = []
    rm_idx = []

    sign = angle > 0 ? 1 : -1
    y_turn = angle > 0 ? upper_turn : lower_turn
    for i in eachindex(x)
        if crease_frac - turn_distance < x[i] < crease_frac + turn_distance && sign * y[i] > 0
            append!(smooth_idx, i)
        elseif sign * y[i] < 0 && crease_frac > x[i] > crease_frac - turn_distance
            append!(rm_idx, i)
        end
        if x[i] > crease_frac
            x_rel = x[i] - crease_frac
            y_rel = y[i] - y_turn
            x[i] = crease_frac + x_rel * cos(angle) + y_rel * sin(angle)
            y[i] = y_turn - x_rel * sin(angle) + y_rel * cos(angle)
            if angle > 0 && x[i] < crease_frac - turn_distance/2 && y[i] > lower_turn
                append!(rm_idx, i)
            elseif angle < 0 && x[i] < crease_frac - turn_distance/2 && y[i] < upper_turn
                append!(rm_idx, i)
            end
        end
    end

    #TODO: lower and upper is slightly off because of smoothing
    lower_i, upper_i = minimum(smooth_idx), maximum(smooth_idx)
    for i in smooth_idx
        window = min(i - lower_i + 1, upper_i - i + 1)
        x[i] = mean(x[i-window:i+window])
    end
    deleteat!(x, rm_idx)
    deleteat!(y, rm_idx)
    nothing
end

function solve_alpha!(cls, cds, cms, alpha_range, alpha_idxs, delta, re, x_, y_, lower, upper, kite_speed, speed_of_sound, crease_frac)
    x = deepcopy(x_)
    y = deepcopy(y_)
    turn_trailing_edge!(delta, x, y, lower, upper, crease_frac)
    Xfoil.set_coordinates(x, y)
    x, y = Xfoil.pane(npan=140)
    reinit = true
    for (alpha, alpha_idx) in zip(alpha_range, alpha_idxs)
        converged = false
        cl = 0.0
        cd = 0.0
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

function create_polars(; dat_path, polar_path, wind_vel, area, width, crease_frac, alpha_range, delta_range)
    @info "Creating polars. This can take several minutes."
    tic()

    cl_matrix = zeros(length(alpha_range), length(delta_range))
    cd_matrix = zeros(length(alpha_range), length(delta_range))
    cm_matrix = zeros(length(alpha_range), length(delta_range))

    kite_speed = wind_vel
    chord_length = area / width
    local reynolds_number = kite_speed * chord_length / KINEMATIC_VISCOSITY # https://en.wikipedia.org/wiki/Reynolds_number

    # Read airfoil coordinates from a file.
    local x, y = open(dat_path, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            try
                push!(x, parse(Float64, entries[1]))
                push!(y, parse(Float64, entries[2]))
            catch ArgumentError
                println(entries)
            end
        end
        x, y
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
    
    serialize(polar_path, (alpha_range, delta_range, cl_matrix, cd_matrix, cm_matrix))
        
    toc()
end


"""
    write_aero_matrix(filepath::String, matrix::Matrix{Float64}, 
                      alpha_range::Vector{Float64}, delta_range::Vector{Float64};
                      label::String="C_l")

Write an aerodynamic coefficient matrix to CSV with angle labels.
The first row contains flap deflection angles, first column contains angles of attack.

# Arguments
- `filepath`: Path to output CSV file
- `matrix`: Matrix of aerodynamic coefficients
- `alpha_range`: Vector of angle of attack values in radians
- `delta_range`: Vector of flap deflection angles in radians
- `label`: Coefficient label for the header (default: "C_l")
"""
function write_aero_matrix(filepath::String, matrix::Matrix{Float64}, 
                         alpha_range::Vector{Float64}, delta_range::Vector{Float64};
                         label::String="C_l")
    open(filepath, "w") do io
        # Write header with delta values
        println(io, "$label/delta," * join(["δ=$(round(rad2deg(δ), digits=1))°" for δ in delta_range], ","))
        
        # Write data rows with alpha values and coefficients
        for i in eachindex(alpha_range)
            row = "α=$(round(rad2deg(alpha_range[i]), digits=1))°," * join(matrix[i,:], ",")
            println(io, row)
        end
    end
end

"""
    read_aero_matrix(filepath::String) -> (Matrix{Float64}, Vector{Float64}, Vector{Float64})

Read an aerodynamic coefficient matrix from CSV with angle labels.
Returns the coefficient matrix and corresponding angle ranges.

# Returns
- `matrix`: Matrix of aerodynamic coefficients
- `alpha_range`: Vector of angle of attack values in radians
- `delta_range`: Vector of flap deflection angles in radians
"""
function read_aero_matrix(filepath::String)
    lines = readlines(filepath)
    
    # Parse header to get delta values
    header = split(lines[1], ',')
    delta_values = map(header[2:end]) do δ_str
        # Extract number between "δ=" and "°"
        m = match(r"δ=(-?\d+\.?\d*)°", δ_str)
        deg2rad(parse(Float64, m.captures[1]))
    end
    
    # Initialize matrix
    n_rows = length(lines) - 1  # Subtract header
    n_cols = length(delta_values)
    matrix = zeros(n_rows, n_cols)
    alpha_values = zeros(n_rows)
    
    # Parse data rows
    for (i, line) in enumerate(lines[2:end])
        entries = split(line, ',')
        # Extract alpha value
        m = match(r"α=(-?\d+\.?\d*)°", entries[1])
        alpha_values[i] = deg2rad(parse(Float64, m.captures[1]))
        # Parse coefficient values
        matrix[i,:] .= parse.(Float64, entries[2:end])
    end
    
    return matrix, alpha_values, delta_values
end


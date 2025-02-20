using Distributed, Timers, Serialization, SharedArrays
using Xfoil
using Pkg
using ControlPlots

function distributed_create_polars(foil_file, v_wind, middle_length, tip_length)
    procs = addprocs()

    function normalize!(x, y)
        x_min = minimum(x)
        x_max = maximum(x)
        for i in eachindex(x)
            x[i] = (x[i] - x_min) / (x_max - x_min)
            y[i] = (y[i] - x_min) / (x_max - x_min)
        end
    end

    @eval @everywhere begin
        using Xfoil, Statistics, SharedArrays

        function turn_trailing_edge!(angle, x, y, lower_turn, upper_turn)
            theta = deg2rad(angle)
            x_turn = 0.75
            turn_distance = upper_turn - lower_turn
            smooth_idx = []
            rm_idx = []

            sign = theta > 0 ? 1 : -1
            y_turn = theta > 0 ? upper_turn : lower_turn
            for i in eachindex(x)
                if x_turn - turn_distance < x[i] < x_turn + turn_distance && sign * y[i] > 0
                    append!(smooth_idx, i)
                elseif sign * y[i] < 0 && x_turn > x[i] > x_turn - turn_distance
                    append!(rm_idx, i)
                end
                if x[i] > x_turn
                    x_rel = x[i] - x_turn
                    y_rel = y[i] - y_turn
                    x[i] = x_turn + x_rel * cos(theta) + y_rel * sin(theta)
                    y[i] = y_turn - x_rel * sin(theta) + y_rel * cos(theta)
                    if theta > 0 && x[i] < x_turn - turn_distance/2 && y[i] > lower_turn
                        append!(rm_idx, i)
                    elseif theta < 0 && x[i] < x_turn - turn_distance/2 && y[i] < upper_turn
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

        function calculate_constants(d_trailing_edge_angle_, x_, y_, cp, lower, upper)
            d_trailing_edge_angle = deg2rad(d_trailing_edge_angle_)
            x = deepcopy(x_)
            y = deepcopy(y_)
            c_te = 0.0
            if d_trailing_edge_angle > 0
                x_ref, y_ref = 0.75, upper
            else
                x_ref, y_ref = 0.75, lower
            end
            
            # straighten out the trailing_edge in order to find the trailing edge torque constant
            for i in eachindex(x)
                x_rel = x[i] - x_ref
                y_rel = y[i] - y_ref
                x[i] = x_ref + x_rel * cos(-d_trailing_edge_angle) + y_rel * sin(-d_trailing_edge_angle)
                y[i] = y_ref - x_rel * sin(-d_trailing_edge_angle) + y_rel * cos(-d_trailing_edge_angle)
            end
            
            x2 = []
            y2 = []
            for i in 1:(length(x)-1)
                if x[i] > x_ref && lower < y[i] < upper
                    push!(x2, x[i])
                    push!(y2, y[i])
                    dx = x[i+1] - x[i]
                    cp_avg = (cp[i] + cp[i+1]) / 2
                    c_te -= dx * cp_avg * (x[i] - x_ref) # clockwise torque around (x_ref, y_ref)
                end
            end
            return c_te
        end

        function solve_alpha!(cls, cds, c_tes, alphas, alpha_idxs, d_trailing_edge_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound)
            x = deepcopy(x_)
            y = deepcopy(y_)
            turn_trailing_edge!(d_trailing_edge_angle, x, y, lower, upper)
            Xfoil.set_coordinates(x, y)
            x, y = Xfoil.pane(npan=140)
            @show d_trailing_edge_angle
            reinit = true
            for (alpha, alpha_idx) in zip(alphas, alpha_idxs)
                converged = false
                cl = 0.0
                cd = 0.0
                # Solve for the given angle of attack
                cl, cd, _, cm, converged = Xfoil.solve_alpha(alpha, re; iter=50, reinit=reinit, mach=kite_speed/speed_of_sound, xtrip=(0.05, 0.05)) # TODO: use 5% point
                reinit = false
                if converged
                    _, cp = Xfoil.cpdump()
                    c_te = calculate_constants(d_trailing_edge_angle, x, y, cp, lower, upper)
                    cls[alpha_idx] = cl
                    cds[alpha_idx] = cd
                    c_tes[alpha_idx] = c_te
                end
            end
            return nothing
        end

        function run_solve_alpha(alphas, d_trailing_edge_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound)
            cls = Float64[NaN for _ in alphas]
            cds = Float64[NaN for _ in alphas]
            c_tes = Float64[NaN for _ in alphas]
            neg_idxs = sort(findall(alphas .< 0.0), rev=true)
            neg_alphas = alphas[neg_idxs]
            pos_idxs = sort(findall(alphas .>= 0.0))
            pos_alphas = alphas[pos_idxs]
            solve_alpha!(cls, cds, c_tes, neg_alphas, neg_idxs, d_trailing_edge_angle, 
                                re, x_, y_, lower, upper, kite_speed, speed_of_sound)
            solve_alpha!(cls, cds, c_tes, pos_alphas, pos_idxs, d_trailing_edge_angle, 
                                re, x_, y_, lower, upper, kite_speed, speed_of_sound)
            return cls, cds, c_tes
        end
    end

    function get_lower_upper(x, y)
        lower_trailing_edge = 0.0
        upper_trailing_edge = 0.0
        min_lower_distance = Inf
        min_upper_distance = Inf
        for (xi, yi) in zip(x, y)
            if yi < 0
                lower_distance = abs(xi - 0.75)
                if lower_distance < min_lower_distance
                    min_lower_distance = lower_distance
                    lower_trailing_edge = yi
                end
            else
                upper_distance = abs(xi - 0.75)
                if upper_distance < min_upper_distance
                    min_upper_distance = upper_distance
                    upper_trailing_edge = yi
                end
            end
        end
        return lower_trailing_edge, upper_trailing_edge
    end

    function remove_nothing(matrix)
        # Identify rows and columns containing `nothing`
        rows_to_keep = [all(!isnothing, matrix[i, :]) for i in 1:size(matrix)[1]]
        cols_to_keep = [all(!isnothing, matrix[:, j]) for j in 1:size(matrix)[2]]

        return matrix[rows_to_keep, cols_to_keep]
    end

    println("Creating polars")
    if !endswith(foil_file, ".dat")
        foil_file *= ".dat"
    end
    polar_file = foil_file[1:end-4] * "_polar.bin"

    # alphas = -90:1.0:90
    alphas = -1.0:1.0:1.0
    # d_trailing_edge_angles = -90:1.0:90
    d_trailing_edge_angles = -1.0:1.0:1.0
    cl_matrix = SharedArray{Float64}((length(alphas), length(d_trailing_edge_angles)), init = (a) -> fill!(a, NaN))
    cd_matrix = SharedArray{Float64}((length(alphas), length(d_trailing_edge_angles)), init = (a) -> fill!(a, NaN))
    c_te_matrix = SharedArray{Float64}((length(alphas), length(d_trailing_edge_angles)), init = (a) -> fill!(a, NaN))
    
    kite_speed = v_wind
    speed_of_sound = 343
    reynolds_number = kite_speed * (middle_length + tip_length)/2 / 1.460e-5

    # Read airfoil coordinates from a file.
    x, y = open(foil_file, "r") do f
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
    normalize!(x, y)
    Xfoil.set_coordinates(x, y)
    x, y = Xfoil.pane(npan=140)
    lower, upper = get_lower_upper(x, y)

    try
        @sync @distributed for j in eachindex(d_trailing_edge_angles)
            cl_matrix[:, j], cd_matrix[:, j], c_te_matrix[:, j] = run_solve_alpha(alphas, d_trailing_edge_angles[j], 
                            reynolds_number, x, y, lower, upper, kite_speed, speed_of_sound)
        end
    catch e
        println(e)
    finally
        println("Removing processes")
        rmprocs(procs)
    end

    println("cl_matrix")
    [println(cl_matrix[i, :]) for i in eachindex(alphas)]

    println("Relative trailing_edge height: ", upper - lower)
    println("Reynolds number for flying speed of $kite_speed is $reynolds_number")

    serialize(polar_file, (alphas, d_trailing_edge_angles, Matrix(cl_matrix), Matrix(cd_matrix), Matrix(c_te_matrix)))
    toc()
end

distributed_create_polars("data/centre_line_profile_with_billow", 9.81, 2.0, 1.0)
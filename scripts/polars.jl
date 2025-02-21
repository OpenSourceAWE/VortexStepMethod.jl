using Distributed, Timers, Serialization, SharedArrays
using Interpolations
using Xfoil
using Pkg
using ControlPlots
using Logging

@info "Creating polars. This can take several minutes."

procs = addprocs()

try
    function normalize!(x, y)
        x_min = minimum(x)
        x_max = maximum(x)
        for i in eachindex(x)
            x[i] = (x[i] - x_min) / (x_max - x_min)
            y[i] = (y[i] - x_min) / (x_max - x_min)
        end
    end

    @eval @everywhere using VortexStepMethod, Xfoil, Statistics, SharedArrays

    alphas = -5:1.0:30
    # alphas = -1.0:1.0:1.0
    d_trailing_edge_angles = -10:1.0:50
    # d_trailing_edge_angles = -1.0:1.0:1.0
    cl_matrix = SharedArray{Float64}((length(alphas), length(d_trailing_edge_angles)))
    cd_matrix = SharedArray{Float64}((length(alphas), length(d_trailing_edge_angles)))
    cm_matrix = SharedArray{Float64}((length(alphas), length(d_trailing_edge_angles)))

    @everywhere begin
        function turn_trailing_edge!(angle, x, y, lower_turn, upper_turn, x_turn)
            theta = deg2rad(angle)
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
        
        function solve_alpha!(cls, cds, cms, alphas, alpha_idxs, d_trailing_edge_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound, x_turn)
            x = deepcopy(x_)
            y = deepcopy(y_)
            turn_trailing_edge!(d_trailing_edge_angle, x, y, lower, upper, x_turn)
            Xfoil.set_coordinates(x, y)
            x, y = Xfoil.pane(npan=140)
            @info "d_trailing_edge_angle: $d_trailing_edge_angle"
            reinit = true
            for (alpha, alpha_idx) in zip(alphas, alpha_idxs)
                converged = false
                cl = 0.0
                cd = 0.0
                # Solve for the given angle of attack
                cl, cd, _, cm, converged = Xfoil.solve_alpha(alpha, re; iter=50, reinit=reinit, mach=kite_speed/speed_of_sound, xtrip=(0.05, 0.05))
                reinit = false
                if converged
                    cls[alpha_idx] = cl
                    cds[alpha_idx] = cd
                    cms[alpha_idx] = cm
                end
            end
            return nothing
        end
        
        function run_solve_alpha(alphas, d_trailing_edge_angle, re, x_, y_, lower, upper, kite_speed, speed_of_sound, x_turn)
            @info "solving alpha"
            cls = Float64[NaN for _ in alphas]
            cds = Float64[NaN for _ in alphas]
            cms = Float64[NaN for _ in alphas]
            neg_idxs = sort(findall(alphas .< 0.0), rev=true)
            neg_alphas = alphas[neg_idxs]
            pos_idxs = sort(findall(alphas .>= 0.0))
            pos_alphas = alphas[pos_idxs]
            solve_alpha!(cls, cds, cms, neg_alphas, neg_idxs, d_trailing_edge_angle, 
                                re, x_, y_, lower, upper, kite_speed, speed_of_sound, x_turn)
            solve_alpha!(cls, cds, cms, pos_alphas, pos_idxs, d_trailing_edge_angle, 
                                re, x_, y_, lower, upper, kite_speed, speed_of_sound, x_turn)
            return cls, cds, cms
        end
    end

    function get_lower_upper(x, y)
        lower_trailing_edge = 0.0
        upper_trailing_edge = 0.0
        min_lower_distance = Inf
        min_upper_distance = Inf
        for (xi, yi) in zip(x, y)
            if yi < 0
                lower_distance = abs(xi - x_turn)
                if lower_distance < min_lower_distance
                    min_lower_distance = lower_distance
                    lower_trailing_edge = yi
                end
            else
                upper_distance = abs(xi - x_turn)
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

    function replace_nan!(matrix)
        rows, cols = size(matrix)
        distance = max(rows, cols)
        for i in 1:rows
            for j in 1:cols
                if isnan(matrix[i, j])
                    neighbors = []
                    for d in 1:distance
                        found = false
                        if i-d >= 1 && !isnan(matrix[i-d, j]);
                            push!(neighbors, matrix[i-d, j])
                            found = true
                        end
                        if i+d <= rows && !isnan(matrix[i+d, j])
                            push!(neighbors, matrix[i+d, j])
                            found = true
                        end
                        if j-d >= 1 && !isnan(matrix[i, j-d])
                            push!(neighbors, matrix[i, j-d])
                            found = true
                        end
                        if j+d <= cols && !isnan(matrix[i, j+d])
                            push!(neighbors, matrix[i, j+d])
                            found = true
                        end
                        if found; break; end
                    end
                    if !isempty(neighbors)
                        matrix[i, j] = sum(neighbors) / length(neighbors)
                    end
                end
            end
        end
        return nothing
    end
    
    kite_speed = v_wind
    speed_of_sound = 343
    chord_length = area / width
    local reynolds_number = kite_speed * chord_length / 1.460e-5 # wikipedia

    # Read airfoil coordinates from a file.
    local x, y = open(foil_path, "r") do f
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

    @everywhere begin
        x = $x
        y = $y
        x_turn = $x_turn
        reynolds_number = $reynolds_number
    end

    @sync @distributed for j in eachindex(d_trailing_edge_angles)
        cl_matrix[:, j], cd_matrix[:, j], cm_matrix[:, j] = run_solve_alpha(alphas, d_trailing_edge_angles[j], 
                        reynolds_number, x, y, lower, upper, kite_speed, speed_of_sound, x_turn)
    end
    @info "cl_matrix" begin
        display(cl_matrix)
    end
    
    function replace_nan!(matrix)
        rows, cols = size(matrix)
        distance = max(rows, cols)
        for i in 1:rows
            for j in 1:cols
                if isnan(matrix[i, j])
                    neighbors = []
                    for d in 1:distance
                        found = false
                        if i-d >= 1 && !isnan(matrix[i-d, j]);
                            push!(neighbors, matrix[i-d, j])
                            found = true
                        end
                        if i+d <= rows && !isnan(matrix[i+d, j])
                            push!(neighbors, matrix[i+d, j])
                            found = true
                        end
                        if j-d >= 1 && !isnan(matrix[i, j-d])
                            push!(neighbors, matrix[i, j-d])
                            found = true
                        end
                        if j+d <= cols && !isnan(matrix[i, j+d])
                            push!(neighbors, matrix[i, j+d])
                            found = true
                        end
                        if found; break; end
                    end
                    if !isempty(neighbors)
                        matrix[i, j] = sum(neighbors) / length(neighbors)
                    end
                end
            end
        end
        return nothing
    end
    
    replace_nan!(cl_matrix)
    replace_nan!(cd_matrix)
    replace_nan!(cm_matrix)
    
    cl_interp_ = extrapolate(scale(interpolate(cl_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
    cd_interp_ = extrapolate(scale(interpolate(cd_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
    cm_interp_ = extrapolate(scale(interpolate(cm_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
    
    cl_interp = (α, β) -> cl_interp_(α, β)
    cd_interp = (α, β) -> cd_interp_(α, β)
    cm_interp = (α, β) -> cm_interp_(α, β)
    
    @info "Relative trailing_edge height: $upper - $lower"
    @info "Reynolds number for flying speed of $kite_speed is $reynolds_number"
    
    toc()
    serialize(polar_path, (cl_interp, cd_interp, cm_interp))

catch e
    @info "Removing processes"
    rmprocs(procs)
    throw(e)
finally
    @info "Removing processes"
    rmprocs(procs)
end    


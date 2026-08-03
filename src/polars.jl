
"""
    interpolate_matrix_nans!(matrix::Matrix{Float64}; prn=true)

Replace NaN values in a matrix by interpolating from nearest non-NaN neighbors.
Uses an expanding search radius until valid neighbors are found.

# Arguments
- `matrix`: Matrix containing NaN values to be interpolated
"""
function interpolate_matrix_nans!(matrix::Matrix{Float64}; prn=true)
    rows, cols = size(matrix)
    nans_found = 0
    while any(isnan, matrix)
        for i in 1:rows, j in 1:cols
            if isnan(matrix[i,j])
                # Search in expanding radius until we find valid neighbors
                radius = 1
                values = Float64[]
                weights = Float64[]
                
                while isempty(values) && radius < max(rows, cols)
                    # Check all points at current Manhattan distance
                    for di in -radius:radius, dj in -radius:radius
                        if abs(di) + abs(dj) == radius  # Points exactly at distance 'radius'
                            ni, nj = i + di, j + dj
                            if 1 ≤ ni ≤ rows && 1 ≤ nj ≤ cols && !isnan(matrix[ni,nj])
                                # Weight by inverse distance
                                dist = sqrt(di^2 + dj^2)
                                push!(values, matrix[ni,nj])
                                push!(weights, 1/dist)
                            end
                        end
                    end
                    radius += 1
                end
                
                if !isempty(values)
                    # Calculate weighted average of found values
                    matrix[i,j] = sum(values .* weights) / sum(weights)
                    nans_found += 1
                else
                    throw(ArgumentError("Could not remove NaN"))
                end
            end
        end
    end
    if nans_found > 0 && prn
        @info "Removed $nans_found NaNs from the matrix."
    end
    return matrix
end


"""
    read_aero_matrix(filepath::AbstractString) -> (Matrix{Float64}, Vector{Float64}, Vector{Float64})

Read an aerodynamic coefficient matrix from CSV with angle labels.
Returns the coefficient matrix and corresponding angle ranges.

# Returns
- `matrix`: Matrix of aerodynamic coefficients
- `alpha_range`: Vector of angle of attack values in radians
- `delta_range`: Vector of flap deflection angles in radians
"""
function read_aero_matrix(filepath::AbstractString)
    lines = readlines(String(filepath))
    
    # Parse header to get delta values
    header = split(lines[1], ',')
    delta_values = map(header[2:end]) do δ_str
        # Extract number between "δ=" and "°"
        m = match(r"δ=(-?\d+\.?\d*)°", δ_str)
        m === nothing && throw(ArgumentError("Invalid delta header entry: $δ_str"))
        δ_cap = m.captures[1]
        (δ_cap isa AbstractString) || throw(ArgumentError("Missing delta value in header entry: $δ_str"))
        deg2rad(parse(Float64, String(δ_cap)))
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
        m === nothing && throw(ArgumentError("Invalid alpha row entry: $(entries[1])"))
        α_cap = m.captures[1]
        (α_cap isa AbstractString) || throw(ArgumentError("Missing alpha value in row entry: $(entries[1])"))
        alpha_values[i] = deg2rad(parse(Float64, String(α_cap)))
        # Parse coefficient values
        matrix[i,:] .= parse.(Float64, entries[2:end])
    end
    
    return matrix, alpha_values, delta_values
end

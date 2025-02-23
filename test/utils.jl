using LinearAlgebra
using Logging
using Statistics
using VortexStepMethod: Wing, WingAerodynamics, Solver

include("thesis_oriol_cayon.jl")

"""
Create an array with cosine spacing, from min to max values, with n points
"""
function cosspace(min_val::Real, max_val::Real, n_points::Int64)
    mean_val = (max_val + min_val) / 2
    amp = (max_val - min_val) / 2
    return mean_val .+ amp .* cos.(range(π, 0, length=n_points))
end

"""
Generate 3D coordinates of a rectangular wing with twist and dihedral.

# Arguments
- `chord::Vector{Float64}`: Chord lengths of wing panels
- `span::Float64`: Total wing span
- `twist::Vector{Float64}`: Twist angles in radians
- `beta::Vector{Float64}`: Dihedral angles in radians
- `N::Int64`: Number of spanwise panels
- `dist::String`: Distribution type ("cos" or "lin")

# Returns
- `coord::Matrix{Float64}`: 2N×3 matrix of coordinates
"""
function generate_coordinates_rect_wing(chord, span, twist, beta, N, dist)
    coord = zeros(2 * N, 3)
    span_points = if dist == "cos"
        cosspace(-span/2, span/2, N)
    elseif dist == "lin"
        range(-span/2, span/2, length=N)
    else
        error("Unknown distribution type: $dist")
    end

    for i in 1:N
        coord[2i-1, :] = [
            -0 * chord[i] * cos(twist[i]),
            span_points[i],
            0 * chord[i] * sin(twist[i]) - abs(span_points[i] * sin(beta[i]))
        ]
        coord[2i, :] = [
            1 * chord[i] * cos(twist[i]),
            span_points[i],
            -1 * chord[i] * sin(twist[i]) - abs(span_points[i] * sin(beta[i]))
        ]
    end
    return coord
end

"""
Generate 3D coordinates of a curved wing.

# Arguments
- `chord::Float64`: Wing chord
- `span::Float64`: Wing span
- `theta::Float64`: Angular extent of curvature (radians)
- `R::Float64`: Radius of curvature
- `N::Int64`: Number of spanwise panels
- `dist::String`: Distribution type ("cos", "lin", or "cos2")

# Returns
- `coord::Matrix{Float64}`: 2N×3 matrix of coordinates
"""
function generate_coordinates_curved_wing(chord, span, theta, R, N, dist)
    coord = zeros(2 * N, 3)
    theta_points = if dist == "cos"
        cosspace(-theta, theta, N)
    elseif dist == "lin"
        range(-theta, theta, length=N)
    elseif dist == "cos2"
        theta1 = cosspace(-theta, -theta/N/10, div(N,2))
        theta2 = cosspace(theta/N/10, theta, div(N,2))
        vcat(theta1, theta2)
    else
        error("Unknown distribution type: $dist")
    end

    for i in 1:N
        coord[2i-1, :] = [0, R * sin(theta_points[i]), R * cos(theta_points[i])]
        coord[2i, :] = [chord, R * sin(theta_points[i]), R * cos(theta_points[i])]
    end
    return coord
end

"""
Generate 3D coordinates of an elliptical wing.

# Arguments
- `max_chord::Float64`: Maximum chord length
- `span::Float64`: Wing span
- `N::Int64`: Number of spanwise panels
- `dist::String`: Distribution type ("cos" or "lin")

# Returns
- `coord::Matrix{Float64}`: 2N×3 matrix of coordinates
"""
function generate_coordinates_el_wing(max_chord, span, N, dist)
    coord = zeros(2 * N, 3)
    start = span * 1e-5
    y_arr = if dist == "cos"
        cosspace(-span/2 + start, span/2 - start, N)
    elseif dist == "lin"
        range(-span/2 + start, span/2 - start, length=N)
    else
        error("Unknown distribution type: $dist")
    end

    c_arr = 2 .* sqrt.(1 .- (y_arr ./ (span/2)).^2) .* max_chord/2

    for i in 1:N
        coord[2i-1, :] = [-0.25 * c_arr[i], y_arr[i], 0]
        coord[2i, :] = [0.75 * c_arr[i], y_arr[i], 0]
    end
    return coord
end

"""
Compare all elements in two lists of dictionaries
"""
function assert_list_dict_equality(variable1, variable_expected; atol=1e-5)
    for (var1, var_expected) in zip(variable1, variable_expected)
        for (key1, key2) in zip(keys(var1), keys(var_expected))
            @debug "Comparing key $key1"
            @debug "variable1[$key1] = $(var1[key1])"
            @debug "variable_expected[$key1] = $(var_expected[key1])"
            @test isapprox(var1[key1], var_expected[key1], atol=atol)
        end
    end
end

"""
Compare all elements in two lists of lists of dictionaries
"""
function assert_list_list_dict_equality(variable1, variable_expected; atol=1e-5)
    for (i, (list1, list_expected)) in enumerate(zip(variable1, variable_expected))
        @info "Comparing list $i"
        @info "list1: $list1"
        @info "list_expected: $list_expected"
        
        for (j, (dict1, dict_expected)) in enumerate(zip(list1, list_expected))
            @info "Comparing dictionary $j"
            @info "dict1: $dict1"
            @info "dict_expected: $dict_expected"
            
            @test keys(dict1) == keys(dict_expected)
            
            for (key1, key2) in zip(keys(dict1), keys(dict_expected))
                @info "Comparing key $key1 with $key2"
                @info "dict1[$key1]: $(dict1[key1])"
                @info "dict_expected[$key1]: $(dict_expected[key1])"
                
                @test key1 == key2
                if dict1[key1] isa String
                    @test dict1[key1] == dict_expected[key1]
                else
                    @test isapprox(dict1[key1], dict_expected[key1], atol=atol)
                end
            end
        end
    end
end

"""
Print a matrix with formatted output
"""
function print_matrix(matrix::AbstractMatrix, name::String="Matrix")
    println("$name:")
    display(round.(matrix, digits=3))
    println()
end

"""
Flip created coordinates in pairs.

Takes a matrix of coordinates and flips the order of coordinate pairs while maintaining
the pairing relationship.

# Arguments
- `coord::Matrix{Float64}`: Matrix of coordinates (2N×3)

# Returns
- `Matrix{Float64}`: Flipped coordinate matrix with same shape as input
"""
function flip_created_coord_in_pairs(coord::Matrix{Float64})
    n_pairs = size(coord, 1) ÷ 2
    reshaped = reshape(coord, (2, n_pairs, 3))
    flipped = reverse(reshaped, dims=2)     
    return reshape(flipped, (2*n_pairs, 3))     
end

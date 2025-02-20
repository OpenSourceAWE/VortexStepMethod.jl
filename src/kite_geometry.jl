
function read_faces(filename)
    vertices = []
    faces = []
    
    open(filename) do file
        for line in eachline(file)
            if startswith(line, "v ") && !startswith(line, "vt") && !startswith(line, "vn")
                parts = split(line)
                x = parse(Float64, parts[2])
                y = parse(Float64, parts[3])
                z = parse(Float64, parts[4])
                push!(vertices, [x, y, z])
            elseif startswith(line, "f ")
                parts = split(line)
                # Handle both f v1 v2 v3 and f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 formats
                indices = map(p -> parse(Int, split(p, '/')[1]), parts[2:4])
                push!(faces, indices)
            end
        end
    end
    return vertices, faces
end

function find_circle_center_and_radius(vertices)
    r = zeros(2)
    v_min = zeros(3)
    v_tip = zeros(3)
    v_min .= Inf

    # find the vertex with smallest x in the middle of the kite
    for v in vertices
        if abs(v[2]) ≤ 0.1
            if v[1] < v_min[1]
                v_min .= v
            end
        end
    end
    @show v_min

    # Find vertex furthest in y, -z direction
    max_score = -Inf
    v_tip .= 0.0
    for v in vertices
        # Score each vertex based on y and -z components
        # Higher y and lower z gives higher score
        score = v[2] - v[3]  # y - z
        if score > max_score
            max_score = score
            v_tip .= v
        end
    end

    function r_diff!(du, u, p)
        z = u[1]
        r .= Inf
        r[1] = sqrt(v_min[2]^2 + (v_min[3] - z)^2)
        r[2] = sqrt(v_tip[2]^2 + (v_tip[3] - z)^2)
        du[1] = r[1] - r[2]
        return nothing
    end

    prob = NonlinearProblem(r_diff!, [v_min[3]-0.1], nothing)
    result = NonlinearSolve.solve(prob, NewtonRaphson(; autodiff=AutoFiniteDiff(; relstep = 1e-3, absstep = 1e-3)); abstol = 1e-2)
    r_diff!(zeros(1), result, nothing)
    z = result[1]

    gamma_tip = atan(-v_tip[2], (v_tip[3] - z))

    return z, r[1], gamma_tip
end

# Create interpolations for max and min x coordinates
function create_interpolations(vertices, circle_center_z, radius, gamma_tip)
    gamma_range = range(-abs(gamma_tip)+1e-6, abs(gamma_tip)-1e-6, 100)
    vz_centered = [v[3] - circle_center_z for v in vertices]
    
    max_xs = zeros(length(gamma_range))
    min_xs = zeros(length(gamma_range))
    areas  = zeros(length(gamma_range))
    
    for (j, gamma) in enumerate(gamma_range)
        max_xs[j] = -Inf
        min_xs[j] = Inf
        for (i, v) in enumerate(vertices)
            # Rotate y coordinate to check box containment
            rotated_y = v[2] * cos(gamma) - vz_centered[i] * sin(gamma)
            if gamma ≤ 0.0 && -0.5 ≤ rotated_y ≤ 0.0
                max_xs[j] = max(max_xs[j], v[1])
                min_xs[j] = min(min_xs[j], v[1])
            elseif gamma > 0.0 && 0.0 ≤ rotated_y ≤ 0.5
                max_xs[j] = max(max_xs[j], v[1])
                min_xs[j] = min(min_xs[j], v[1])
            end
        end
        area = abs(max_xs[j] - min_xs[j]) * gamma_range.step * radius
        last_area = j > 1 ? areas[j-1] : 0.0
        areas[j] = last_area + area
    end

    te_interp = linear_interpolation(gamma_range, max_xs)
    le_interp = linear_interpolation(gamma_range, min_xs)
    area_interp = linear_interpolation(gamma_range, areas)
    
    return (le_interp, te_interp, area_interp, gamma_range, max_xs, min_xs)
end

# Calculate center of mass for a triangular surface mesh
function calculate_com(vertices, faces)
    area_total = 0.0
    com = zeros(3)
    
    for face in faces
        v1 = vertices[face[1]]
        v2 = vertices[face[2]]
        v3 = vertices[face[3]]
        
        # Calculate triangle area and centroid
        normal = cross(v2 - v1, v3 - v1)
        area = norm(normal) / 2
        centroid = (v1 + v2 + v3) / 3
        
        area_total += area
        com += area * centroid
    end
    
    return com / area_total
end

# Calculate inertia tensor for a triangular surface mesh with given mass
function calculate_inertia_tensor(vertices, faces, mass, com)
    # Initialize inertia tensor
    I = zeros(3, 3)
    total_area = 0.0
    
    for face in faces
        v1 = vertices[face[1]] .- com
        v2 = vertices[face[2]] .- com
        v3 = vertices[face[3]] .- com
        
        # Calculate triangle area
        normal = cross(v2 - v1, v3 - v1)
        area = norm(normal) / 2
        total_area += area
        
        # Calculate contribution to inertia tensor
        for i in 1:3
            for j in 1:3
                # Vertices relative to center of mass
                points = [v1, v2, v3]
                
                # Calculate contribution to inertia tensor
                for p in points
                    if i == j
                        # Diagonal terms
                        I[i,i] += area * (sum(p.^2) - p[i]^2)
                    else
                        # Off-diagonal terms
                        I[i,j] -= area * (p[i] * p[j])
                    end
                end
            end
        end
    end
    
    # Scale by mass/total_area to get actual inertia tensor
    return (mass / total_area) * I / 3
end


"""
    KiteWing

Represents a curved wing that inherits from Wing with additional geometric properties.

# Fields
- All fields from Wing:
  - `n_panels::Int`: Number of panels in aerodynamic mesh
  - `spanwise_panel_distribution::String`: Panel distribution type
  - `spanwise_direction::Vector{Float64}`: Wing span direction vector
  - `sections::Vector{Section}`: List of wing sections
- Additional fields:
  - `center_of_mass::Vector{Float64}`: Center of mass coordinates
  - `circle_center_z::Vector{Float64}`: Center of circle coordinates
  - `radius::Float64`: Radius of curvature

# Distribution types
Same as Wing
"""
mutable struct KiteWing <: AbstractWing
    n_panels::Int
    spanwise_panel_distribution::String
    spanwise_direction::Vector{Float64}
    sections::Vector{Section}
    
    # Additional fields for KiteWing
    mass::Float64
    center_of_mass::Vector{Float64}
    circle_center_z::Float64
    gamma_tip::Float64
    inertia_tensor::Matrix{Float64}
    radius::Float64
    le_interp::Function
    te_interp::Function
    area_interp::Function

    function KiteWing(obj_path, dat_path=nothing; mass=1.0, n_panels=54, spanwise_panel_distribution="linear", spanwise_direction=[0.0, 1.0, 0.0])
        !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && @error "Spanwise direction has to be [0.0, 1.0, 0.0]"
        if !isfile(obj_path)
            error("OBJ file not found: $obj_path")
        end
        vertices, faces = read_faces(obj_path)
        center_of_mass = calculate_com(vertices, faces)
        inertia_tensor = calculate_inertia_tensor(vertices, faces, mass, center_of_mass)
        
        circle_center_z, radius, gamma_tip = find_circle_center_and_radius(vertices)
        le_interp_, te_interp_, area_interp_, gammas, max_xs, min_xs = create_interpolations(vertices, circle_center_z, radius, gamma_tip)
        le_interp = (γ) -> le_interp_(γ)
        te_interp = (γ) -> te_interp_(γ)
        area_interp = (γ) -> area_interp_(γ)

        sections = Section[]
        if !isnothing(dat_path)
            
        end
        new(
            n_panels, spanwise_panel_distribution, spanwise_direction, sections,
            mass, center_of_mass, circle_center_z, gamma_tip, inertia_tensor, radius,
            le_interp, te_interp, area_interp
        )
    end
end

"""
    add_section!(wing::KiteWing, LE_point::PosVector, 
                TE_point::PosVector, aero_input::Vector{Any})

Add a new section to the wing.
"""
function add_section!(wing::KiteWing, gamma, aero_input; α=0.0)
    LE_point = [0.0, 0.0, wing.circle_center_z] .+ [wing.le_interp(gamma), sin(gamma) * wing.radius, cos(gamma) * wing.radius]
    if !isapprox(α, 0.0)
        local_y_vec = [0.0, sin(-gamma), cos(gamma)] × [1.0, 0.0, 0.0]
        TE_point = LE_point .+ rotate_v_around_k([wing.te_interp(gamma) - wing.le_interp(gamma), 0.0, 0.0], local_y_vec, α)
    else
        TE_point = LE_point .+ [wing.te_interp(gamma) - wing.le_interp(gamma), 0.0, 0.0]
    end
    push!(wing.sections, Section(LE_point, TE_point, aero_input))
    nothing
end

function rotate_v_around_k(v, k, θ)
    k = normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ)  + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end
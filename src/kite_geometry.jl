"""
    read_faces(filename)

Read vertices and faces from an OBJ file.

# Arguments
- `filename::String`: Path to .obj file

# Returns
- Tuple of (vertices, faces) where:
  - vertices: Vector of 3D coordinates [x,y,z]
  - faces: Vector of triangle vertex indices
"""
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
                indices = map(p -> parse(Int64, split(p, '/')[1]), parts[2:4])
                push!(faces, indices)
            end
        end
    end
    return vertices, faces
end

"""
    find_circle_center_and_radius(vertices)

Find the center and radius of the kite's curvature circle.

# Arguments
- `vertices`: Vector of 3D point coordinates

# Returns
- Tuple of (z_center, radius, gamma_tip) where:
  - z_center: Z-coordinate of circle center
  - radius: Circle radius
  - gamma_tip: Angle of the kite tip from z-axis
"""
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

    # Find vertex furthest in -y, -z direction
    max_score = -Inf
    v_tip .= 0.0
    for v in vertices
        # Score each vertex based on -y and -z components
        # lower y and lower z gives higher score
        score = -v[2] - v[3]  # - y - z
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
    @assert gamma_tip > 0.0

    return z, r[1], gamma_tip
end

"""
    create_interpolations(vertices, circle_center_z, radius, gamma_tip)

Create interpolation functions for leading/trailing edges and area.

# Arguments
- `vertices`: Vector of 3D point coordinates
- `circle_center_z`: Z-coordinate of circle center
- `radius`: Circle radius
- `gamma_tip`: Maximum angular extent

# Returns
- Tuple of (le_interp, te_interp, area_interp) interpolation functions
- Where le_interp and te_interp are tuples themselves, containing the x, y and z interpolations
"""
function create_interpolations(vertices, circle_center_z, radius, gamma_tip)
    gamma_range = range(-gamma_tip+1e-6, gamma_tip-1e-6, 100)
    vz_centered = [v[3] - circle_center_z for v in vertices]
    
    trailing_edges = zeros(3, length(gamma_range))
    leading_edges = zeros(3, length(gamma_range))
    areas  = zeros(length(gamma_range))
    
    for (j, gamma) in enumerate(gamma_range)
        trailing_edges[1, j] = -Inf
        leading_edges[1, j] = Inf
        for (i, v) in enumerate(vertices)
            # Rotate y coordinate to check box containment
            rotated_y = v[2] * cos(gamma) - vz_centered[i] * sin(gamma)
            if gamma ≤ 0.0 && -0.5 ≤ rotated_y ≤ 0.0
                if v[1] > trailing_edges[1, j]
                    trailing_edges[:, j] .= v
                end
                if v[1] < leading_edges[1, j]
                    leading_edges[:, j] .= v
                end
            elseif gamma > 0.0 && 0.0 ≤ rotated_y ≤ 0.5
                if v[1] > trailing_edges[1, j]
                    trailing_edges[:, j] .= v
                end
                if v[1] < leading_edges[1, j]
                    leading_edges[:, j] .= v
                end
            end
        end
        area = norm(leading_edges[:, j] - trailing_edges[:, j]) * gamma_range.step * radius
        last_area = j > 1 ? areas[j-1] : 0.0
        areas[j] = last_area + area
    end

    le_interp = ntuple(i -> linear_interpolation(gamma_range, leading_edges[i, :],
                                           extrapolation_bc=Line()), 3)
    te_interp = ntuple(i -> linear_interpolation(gamma_range, trailing_edges[i, :],
                                           extrapolation_bc=Line()), 3)
    area_interp = linear_interpolation(gamma_range, areas, extrapolation_bc=Line())
    
    return (le_interp, te_interp, area_interp)
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

"""
    calculate_inertia_tensor(vertices, faces, mass, com)

Calculate the inertia tensor for a triangulated surface mesh, assuming a thin shell with uniform 
surface density.

# Arguments
- `vertices`: Vector of 3D point coordinates representing mesh vertices
- `faces`: Vector of triangle indices, each defining a face of the mesh
- `mass`: Total mass of the shell in kg
- `com`: Center of mass coordinates [x,y,z]

# Method
Uses the thin shell approximation where:
1. Mass is distributed uniformly over the surface area
2. Each triangle contributes to the inertia based on its area and position
3. For each triangle vertex p, contribution to diagonal terms is: area * (sum(p²) - p_i²)
4. For off-diagonal terms: area * (-`p_i` * `p_j`)
5. Final tensor is scaled by mass/(3*total_area) to get correct units

# Returns
- 3×3 matrix representing the inertia tensor in kg⋅m²
"""
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
    interpolate_matrix_nans!(matrix::Matrix{Float64})

Replace NaN values in a matrix by interpolating from nearest non-NaN neighbors.
Uses an expanding search radius until valid neighbors are found.

# Arguments
- `matrix`: Matrix containing NaN values to be interpolated
"""
function interpolate_matrix_nans!(matrix::Matrix{Float64})
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
    @info "Removed $nans_found NaNs from the matrix."
    return matrix
end


"""
    KiteWing

Represents a curved wing that inherits from Wing with additional geometric properties.

# Fields
- All fields from Wing:
  - `n_panels::Int64`: Number of panels in aerodynamic mesh
  - `spanwise_panel_distribution`::PanelDistribution: see: [PanelDistribution](@ref)
  - `spanwise_direction::MVec3`: Wing span direction vector
  - `sections::Vector{Section}`: List of wing sections, see: [Section](@ref)
  -  refined_sections::Vector{Section}
  - `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not
- Additional fields:
  - `center_of_mass::Vector{Float64}`: Center of mass coordinates
  - `circle_center_z::Vector{Float64}`: Center of circle coordinates
  - gamma_tip::Float64: Angle between the body frame z axis and the vector going from the kite circular shape center to the wing tip.
  - `inertia_tensor`::Matrix{Float64}: see: [`calculate_inertia_tensor`](@ref)
  - radius::Float64: Radius of curvature
  - le_interp::NTuple{3, Extrapolation}: see: [Extrapolation](https://juliamath.github.io/Interpolations.jl/stable/extrapolation/)
  - te_interp::NTuple{3, Extrapolation}
  - area_interp::Extrapolation

"""
mutable struct KiteWing <: AbstractWing
    n_panels::Int64
    spanwise_panel_distribution::PanelDistribution
    spanwise_direction::MVec3
    sections::Vector{Section}
    refined_sections::Vector{Section}
    remove_nan::Bool
    
    # Additional fields for KiteWing
    non_deformed_sections::Vector{Section}
    mass::Float64
    center_of_mass::Vector{Float64}
    circle_center_z::Float64
    gamma_tip::Float64
    inertia_tensor::Matrix{Float64}
    radius::Float64
    le_interp::NTuple{3, Extrapolation}
    te_interp::NTuple{3, Extrapolation}
    area_interp::Extrapolation
    alpha_dist::Vector{Float64}
    beta_dist::Vector{Float64}
end

"""
    KiteWing(obj_path, dat_path; alpha=0.0, crease_frac=0.75, wind_vel=10., mass=1.0, 
             n_panels=54, n_sections=n_panels+1, spanwise_panel_distribution=UNCHANGED, 
             spanwise_direction=[0.0, 1.0, 0.0], remove_nan::Bool=true)

Constructor for a [KiteWing](@ref) that allows to use an `.obj` and a `.dat` file as input.

# Parameters
- obj_path: Path to the `.obj` file used for creating the geometry
- dat_path: Path to the `.dat` file, a standard format for 2d foil geometry

# Keyword Parameters
- alpha=0.0: Angle of attack of each segment relative to the x axis [rad]
- crease_frac=0.75: The x coordinate around which the trailing edge rotates on a normalized 2d foil, 
                    used in the xfoil polar generation
- wind_vel=10.0: Apparent wind speed in m/s, used in the xfoil polar generation
- mass=1.0: Mass of the wing in kg, used for the inertia calculations 
- `n_panels`=54: Number of panels.
- `n_sections`=n_panels+1: Number of sections (there is a section on each side of each panel.)
- `spanwise_panel_distribution`=UNCHANGED: see: [PanelDistribution](@ref)
- `spanwise_direction`=[0.0, 1.0, 0.0]
- `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not
"""
function KiteWing(obj_path, dat_path; alpha=0.0, crease_frac=0.75, wind_vel=10., mass=1.0, 
                  n_panels=54, n_sections=n_panels+1, spanwise_panel_distribution=UNCHANGED, 
                  spanwise_direction=[0.0, 1.0, 0.0], remove_nan=true)

    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))

    # Load or create polars
    (!endswith(dat_path, ".dat")) && (dat_path *= ".dat")
    (!isfile(dat_path)) && error("DAT file not found: $dat_path")
    polar_path = dat_path[1:end-4] * "_polar.bin"

    (!endswith(obj_path, ".obj")) && (obj_path *= ".obj")
    (!isfile(obj_path)) && error("OBJ file not found: $obj_path")
    info_path = obj_path[1:end-4] * "_info.bin"

    if !ispath(polar_path) || !ispath(info_path)
        @info "Reading $obj_path"
        vertices, faces = read_faces(obj_path)
        center_of_mass = calculate_com(vertices, faces)
        !(abs(center_of_mass[2]) < 0.01) && @error "Center of mass $center_of_mass has to lie on x-axis."
        inertia_tensor = calculate_inertia_tensor(vertices, faces, mass, center_of_mass)

        circle_center_z, radius, gamma_tip = find_circle_center_and_radius(vertices)
        le_interp, te_interp, area_interp = create_interpolations(vertices, circle_center_z, radius, gamma_tip)
        @info "Writing $info_path"
        serialize(info_path, (center_of_mass, inertia_tensor, circle_center_z, radius, gamma_tip, 
            le_interp, te_interp, area_interp))

        width = 2gamma_tip * radius
        area = area_interp(gamma_tip)
        @eval Main begin
            foil_path, polar_path, v_wind, area, width, x_turn =
                $dat_path, $polar_path, $wind_vel, $gamma_tip, $width, $crease_frac
            include(joinpath(dirname(@__FILE__), "../scripts/polars.jl"))
        end
    end

    @info "Loading polars and kite info from $polar_path and $info_path"
    (alpha_range, beta_range, cl_matrix::Matrix, cd_matrix::Matrix, cm_matrix::Matrix) = deserialize(polar_path)
    if remove_nan
        interpolate_matrix_nans!(cl_matrix)
        interpolate_matrix_nans!(cd_matrix)
        interpolate_matrix_nans!(cm_matrix)
    end

    (center_of_mass, inertia_tensor, circle_center_z, radius, gamma_tip, 
        le_interp, te_interp, area_interp) = deserialize(info_path)

    # Create sections
    sections = Section[]
    for gamma in range(-gamma_tip, gamma_tip, n_sections)
        aero_data = (collect(alpha_range), collect(beta_range), cl_matrix, cd_matrix, cm_matrix)
        LE_point = [le_interp[i](gamma) for i in 1:3]
        TE_point = [te_interp[i](gamma) for i in 1:3]
        push!(sections, Section(LE_point, TE_point, POLAR_MATRICES, aero_data))
    end

    KiteWing(n_panels, spanwise_panel_distribution, spanwise_direction, sections, sections, remove_nan, sections,
        mass, center_of_mass, circle_center_z, gamma_tip, inertia_tensor, radius,
        le_interp, te_interp, area_interp, zeros(n_panels), zeros(n_panels))
end

"""
    deform!(wing::KiteWing, alphas::AbstractVector, betas::AbstractVector; width)

Deform wing by applying left and right alpha and beta.

# Arguments
- `wing`: KiteWing to deform
- `alphas`: [left, right] the angle between of the kite and the body x-axis in radians
- `betas`: [left, right] the deformation of the trailing edges
- `width`: Transition width in meters to smoothe out the transition from left to right deformation

# Effects
Updates wing.sections with deformed geometry
"""
function deform!(wing::KiteWing, alphas::AbstractVector, betas::AbstractVector; width)
    local_y = zeros(MVec3)
    chord = zeros(MVec3)

    for (i, gamma) in enumerate(range(-wing.gamma_tip, wing.gamma_tip, wing.n_panels))
        normalized_gamma = (gamma * wing.radius / width + 0.5)  # Maps [-0.5, 0.5] to [0, 1]
        wing.alpha_dist[i] = if normalized_gamma <= 0.0
            alphas[1]
        elseif normalized_gamma >= 1.0
            alphas[2]
        else
            alphas[1] * (1.0 - normalized_gamma) + alphas[2] * normalized_gamma
        end
        wing.beta_dist[i] = if normalized_gamma <= 0.0
            betas[1]
        elseif normalized_gamma >= 1.0
            betas[2]
        else
            betas[1] * (1.0 - normalized_gamma) + betas[2] * normalized_gamma
        end

        section1 = wing.non_deformed_sections[i]
        section2 = wing.non_deformed_sections[i+1]
        local_y .= normalize(section1.LE_point - section2.LE_point)
        chord .= section1.TE_point .- section1.LE_point
        wing.sections[i].TE_point .= section1.LE_point .+ rotate_v_around_k(chord, local_y, wing.alpha_dist[i])
    end
    return nothing
end

"""
    rotate_v_around_k(v, k, θ)

Rotate vector v around axis k by angle θ using Rodrigues' rotation formula.

# Arguments
- `v`: Vector to rotate
- `k`: Rotation axis (will be normalized)
- `θ`: Rotation angle in radians

# Returns
- Rotated vector
"""
function rotate_v_around_k(v, k, θ)
    k = normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ)  + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end

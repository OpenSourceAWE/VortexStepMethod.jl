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
    vertices = Vector{Float64}[]
    faces = Vector{Int64}[]

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

    function r_diff!(du, u, _)
        z_local = u[1]
        r .= Inf
        r[1] = sqrt(v_min[2]^2 + (v_min[3] - z_local)^2)
        r[2] = sqrt(v_tip[2]^2 + (v_tip[3] - z_local)^2)
        du[1] = r[1] - r[2]
        return nothing
    end

    prob = NonlinearProblem(r_diff!, [v_min[3]-0.1], nothing)
    result = NonlinearSolve.solve(prob, NewtonRaphson(; autodiff=AutoFiniteDiff(; relstep = 1e-3, absstep = 1e-3)); abstol = 1e-2)
    r_diff!(zeros(1), result, nothing)
    z_center = result[1]

    gamma_tip = atan(-v_tip[2], (v_tip[3] - z_center))
    @assert gamma_tip > 0.0

    return z_center, r[1], gamma_tip
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
function create_interpolations(vertices, circle_center_z, radius, gamma_tip, R=I(3); interp_steps=40)
    gamma_range = range(-gamma_tip+gamma_tip/interp_steps*2,
                        gamma_tip-gamma_tip/interp_steps*2, interp_steps)
    stepsize = gamma_range.step.hi
    vz_centered = [v[3] - circle_center_z for v in vertices]

    te_gammas = zeros(length(gamma_range))
    le_gammas = zeros(length(gamma_range))
    trailing_edges = zeros(3, length(gamma_range))
    leading_edges = zeros(3, length(gamma_range))
    areas  = zeros(length(gamma_range))

    n_slices = length(gamma_range)
    for (j, gamma) in enumerate(gamma_range)
        trailing_edges[1, j] = -Inf
        leading_edges[1, j] = Inf

        # Determine if this is a tip slice and get search parameters
        is_first_tip = (j == 1)
        is_last_tip = (j == n_slices)

        if is_first_tip || is_last_tip
            # Tip slices: use directional search within adjacent slice region
            gamma_search = is_first_tip ? gamma_range[1] : gamma_range[end]
            max_te_score = -Inf
            max_le_score = -Inf

            for (i, v) in enumerate(vertices)
                gamma_v = atan(-v[2], vz_centered[i])

                # Check if vertex is in the adjacent slice region
                in_range = if gamma_search ≤ 0
                    gamma_search - stepsize ≤ gamma_v ≤ gamma_search
                else
                    gamma_search ≤ gamma_v ≤ gamma_search + stepsize
                end

                if in_range
                    if is_first_tip
                        # TE: furthest in [X, Y, -Z] direction
                        te_score = v[1] + v[2] - v[3]
                        if te_score > max_te_score
                            trailing_edges[:, j] .= v
                            te_gammas[j] = gamma_v
                            max_te_score = te_score
                        end
                        # LE: furthest in [-X, Y, -Z] direction
                        le_score = -v[1] + v[2] - v[3]
                        if le_score > max_le_score
                            leading_edges[:, j] .= v
                            le_gammas[j] = gamma_v
                            max_le_score = le_score
                        end
                    else  # is_last_tip
                        # TE: furthest in [X, -Y, -Z] direction
                        te_score = v[1] - v[2] - v[3]
                        if te_score > max_te_score
                            trailing_edges[:, j] .= v
                            te_gammas[j] = gamma_v
                            max_te_score = te_score
                        end
                        # LE: furthest in [-X, -Y, -Z] direction
                        le_score = -v[1] - v[2] - v[3]
                        if le_score > max_le_score
                            leading_edges[:, j] .= v
                            le_gammas[j] = gamma_v
                            max_le_score = le_score
                        end
                    end
                end
            end
        else
            # Interior slices: use standard min/max x-coordinate search
            for (i, v) in enumerate(vertices)
                gamma_v = atan(-v[2], vz_centered[i])
                if gamma ≤ 0 && gamma - stepsize ≤ gamma_v ≤ gamma
                    if v[1] > trailing_edges[1, j]
                        trailing_edges[:, j] .= v
                        te_gammas[j] = gamma_v
                    end
                    if v[1] < leading_edges[1, j]
                        leading_edges[:, j] .= v
                        le_gammas[j] = gamma_v
                    end
                elseif gamma > 0 && gamma ≤ gamma_v ≤ gamma + stepsize
                    if v[1] > trailing_edges[1, j]
                        trailing_edges[:, j] .= v
                        te_gammas[j] = gamma_v
                    end
                    if v[1] < leading_edges[1, j]
                        leading_edges[:, j] .= v
                        le_gammas[j] = gamma_v
                    end
                end
            end
        end

        area = norm(leading_edges[:, j] - trailing_edges[:, j]) * stepsize * radius
        last_area = j > 1 ? areas[j-1] : 0.0
        areas[j] = last_area + area
    end

    for j in eachindex(gamma_range)
        leading_edges[:, j] .= R * leading_edges[:, j]
        trailing_edges[:, j] .= R * trailing_edges[:, j]
    end

    le_interp = ntuple(i -> linear_interpolation(le_gammas, leading_edges[i, :],
                                           extrapolation_bc=Line()), 3)
    te_interp = ntuple(i -> linear_interpolation(te_gammas, trailing_edges[i, :],
                                           extrapolation_bc=Line()), 3)
    area_interp = linear_interpolation(gamma_range, areas, extrapolation_bc=Line())

    return (le_interp, te_interp, area_interp)
end

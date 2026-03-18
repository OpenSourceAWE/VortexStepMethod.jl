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

"""
    refine_obj_wing!(wing::AbstractWing; recompute_mapping=true)

Refine OBJ wing by computing position deltas and applying them to refined sections.

This method enables deformation support for OBJ wings by:
1. Recalculating evenly-spaced gammas for unrefined sections
2. Computing what unrefined sections SHOULD be (from interpolations)
3. Computing deltas between current and interpolated positions
4. Creating refined sections from interpolations + interpolated deltas
5. Computing panel mapping

# Arguments
- `wing::AbstractWing`: OBJ wing with le_interp/te_interp
- `recompute_mapping::Bool=true`: Whether to recompute refined_panel_mapping

# Effects
Updates wing.refined_sections and wing.non_deformed_sections in-place.
"""
function refine_obj_wing!(wing::AbstractWing; recompute_mapping=true)
    n_unrefined = wing.n_unrefined_sections
    n_refined = wing.n_panels + 1

    # 1. Calculate evenly spaced gammas for unrefined sections
    unrefined_gammas = range(-wing.gamma_tip, wing.gamma_tip, n_unrefined)

    # 2. Recalculate what unrefined sections SHOULD be from interpolations
    interpolated_unrefined_le = Matrix{Float64}(undef, n_unrefined, 3)
    interpolated_unrefined_te = Matrix{Float64}(undef, n_unrefined, 3)
    for (i, gamma) in enumerate(unrefined_gammas)
        interpolated_unrefined_le[i, :] .= [wing.le_interp[j](gamma) for j in 1:3]
        interpolated_unrefined_te[i, :] .= [wing.te_interp[j](gamma) for j in 1:3]
    end

    # 3. Compute deltas: current - interpolated
    deltas_le = Matrix{Float64}(undef, n_unrefined, 3)
    deltas_te = Matrix{Float64}(undef, n_unrefined, 3)
    for i in 1:n_unrefined
        deltas_le[i, :] .= wing.unrefined_sections[i].LE_point -
                           interpolated_unrefined_le[i, :]
        deltas_te[i, :] .= wing.unrefined_sections[i].TE_point -
                           interpolated_unrefined_te[i, :]
    end

    # 4. Create refined sections with interpolated deltas
    refined_gammas = range(-wing.gamma_tip, wing.gamma_tip, n_refined)
    if isempty(wing.refined_sections)
        wing.refined_sections = [Section() for _ in 1:n_refined]
    end

    for (idx, gamma) in enumerate(refined_gammas)
        # Get base position from interpolation
        base_le = [wing.le_interp[i](gamma) for i in 1:3]
        base_te = [wing.te_interp[i](gamma) for i in 1:3]

        # Find surrounding unrefined sections for delta interpolation
        unrefined_idx = searchsortedlast(collect(unrefined_gammas), gamma)
        unrefined_idx = clamp(unrefined_idx, 1, n_unrefined - 1)

        # Linear interpolation weight
        gamma_left = unrefined_gammas[unrefined_idx]
        gamma_right = unrefined_gammas[unrefined_idx + 1]
        t = (gamma - gamma_left) / (gamma_right - gamma_left)

        # Interpolate deltas
        delta_le = (1 - t) * deltas_le[unrefined_idx, :] +
                   t * deltas_le[unrefined_idx + 1, :]
        delta_te = (1 - t) * deltas_te[unrefined_idx, :] +
                   t * deltas_te[unrefined_idx + 1, :]

        # Apply deltas to get final position
        final_le = base_le + delta_le
        final_te = base_te + delta_te

        # Update refined section
        aero_model = wing.unrefined_sections[1].aero_model
        aero_data = wing.unrefined_sections[1].aero_data
        VortexStepMethod.reinit!(wing.refined_sections[idx], final_le, final_te,
                aero_model, aero_data)
    end

    # 5. Compute panel mapping and update non_deformed_sections
    recompute_mapping && VortexStepMethod.compute_refined_panel_mapping!(wing)
    VortexStepMethod.update_non_deformed_sections!(wing)

    return nothing
end

"""
    center_to_com!(vertices, faces)

Calculate center of mass of a mesh and translate vertices so that COM is at origin.

# Arguments
- `vertices`: Vector of 3D point coordinates
- `faces`: Vector of vertex indices for each face (can be triangular or non-triangular)

# Returns
- Vector representing the original center of mass before translation

# Notes
- Non-triangular faces are automatically triangulated into triangles
- Assumes uniform surface density
"""
function center_to_com!(vertices, faces; prn=true)
    area_total = 0.0
    com = zeros(3)

    for face in faces
        if length(face) == 3
            # Triangle case
            v1 = vertices[face[1]]
            v2 = vertices[face[2]]
            v3 = vertices[face[3]]

            # Calculate triangle area and centroid
            normal = cross(v2 - v1, v3 - v1)
            area = norm(normal) / 2
            centroid = (v1 + v2 + v3) / 3

            area_total += area
            com -= area * centroid
        else
            throw(ArgumentError("Triangulate faces in a CAD program first"))
        end
    end

    com = com / area_total
    !(abs(com[2]) < 0.01) && throw(ArgumentError("Center of mass $com of .obj file has to lie on the xz-plane."))
    prn && @info "Centering vertices of .obj file to the center of mass: $com"
    com[2] = 0.0
    for v in vertices
        v .+= com
    end
    return com
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

function calc_inertia_y_rotation(I_b_tensor)
    # Function for nonlinear solver - off-diagonal element should be zero
    function eq!(F, theta, _)
        # Rotation matrix around y-axis
        R_y = [
            cos(theta[1])  0  sin(theta[1]);
            0              1  0;
            -sin(theta[1]) 0  cos(theta[1])
        ]
        # Transform inertia tensor
        I_rotated = R_y * I_b_tensor * R_y'
        # We want the off-diagonal xz elements to be zero
        F[1] = I_rotated[1,3]
    end

    theta0 = [0.0]
    prob = NonlinearProblem(eq!, theta0, nothing)
    sol = NonlinearSolve.solve(prob, NewtonRaphson())
    theta_opt = sol.u[1]

    R_b_p = [
        cos(theta_opt)  0  sin(theta_opt);
        0               1  0;
        -sin(theta_opt) 0  cos(theta_opt)
    ]
    # Calculate diagonalized inertia tensor
    I_diag = R_b_p * I_b_tensor * R_b_p'
    @assert isapprox(I_diag[1,3], 0.0, atol=1e-5)
    return I_diag, R_b_p
end


"""
    ObjWing(obj_path, dat_path; kwargs...)

Create a deformable wing model from 3D geometry (.obj) and airfoil data (.dat) files.

This constructor builds a complete aerodynamic model by:
1. Loading wing geometry from the .obj file
2. Creating aerodynamic polars from the airfoil .dat file (or loading existing)
3. Computing inertial properties and coordinate transformations
4. Setting up control surfaces and panel distribution

The resulting Wing supports deformation through unrefined_deform! and deform! functions.

# Arguments
- `obj_path`: Path to .obj file containing 3D wing geometry
- `dat_path`: Path to .dat file containing 2D airfoil profile

# Keyword Arguments
- `crease_frac=0.9`: Normalized trailing edge hinge location (0-1)
- `wind_vel=10.0`: Reference wind velocity for XFoil analysis (m/s)
- `mass=1.0`: Wing mass (kg)
- `n_panels=56`: Number of aerodynamic panels across wingspan
- `n_unrefined_sections`: Number of unrefined sections for deformation control (default: inferred from geometry)
- `align_to_principal=false`: Align body frame to principal axes of inertia
- `spanwise_distribution=UNCHANGED`: Panel distribution type (forced to UNCHANGED for ObjWing)
- `remove_nan=true`: Interpolate NaN values in aerodynamic data
- `use_prior_polar=false`: Reuse prior refined/panel polar mapping on geometry updates
- `alpha_range=deg2rad.(-5:1:20)`: Angle of attack range for polars (rad)
- `delta_range=deg2rad.(-5:1:20)`: Trailing edge deflection range for polars (rad)
- `prn=true`: Print informational messages

# Returns
A fully initialized `Wing` instance ready for aerodynamic simulation with deformation support.

# Example
```julia
# Create a deformable wing from geometry files
wing = ObjWing(
    "path/to/wing.obj",
    "path/to/airfoil.dat";
    mass=1.5,
    n_panels=40,
    n_unrefined_sections=4
)

# Apply deformation
unrefined_deform!(wing, deg2rad.([5, 10, 5, 0]), deg2rad.([-5, 0, -5, 0]))
```
"""
function ObjWing(
    obj_path, dat_path;
    crease_frac=0.9, wind_vel=10., mass=1.0,
    n_panels=56, n_unrefined_sections=nothing,
    spanwise_distribution=UNCHANGED,
    spanwise_direction=[0.0, 1.0, 0.0], remove_nan=true, use_prior_polar=false, align_to_principal=false,
    alpha_range=deg2rad.(-5:1:20), delta_range=deg2rad.(-5:1:20), prn=true,
    interp_steps=n_panels+1
)

    # Set default: evenly spaced unrefined sections including both tips
    if isnothing(n_unrefined_sections)
        # Default to having same number of unrefined sections as refined (no interpolation needed)
        n_unrefined_sections = n_panels + 1
    end

    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))

    # Load or create polars
    (!endswith(dat_path, ".dat")) && (dat_path *= ".dat")
    (!isfile(dat_path)) && error("DAT file not found: $dat_path")
    cl_polar_path = dat_path[1:end-4] * "_cl_polar.csv"
    cd_polar_path = dat_path[1:end-4] * "_cd_polar.csv"
    cm_polar_path = dat_path[1:end-4] * "_cm_polar.csv"

    (!endswith(obj_path, ".obj")) && (obj_path *= ".obj")
    (!isfile(obj_path)) && error("OBJ file not found: $obj_path")

    ! prn || @info "Reading $obj_path"
    vertices, faces = read_faces(obj_path)
    T_cad_body = center_to_com!(vertices, faces; prn)
    inertia_tensor = calculate_inertia_tensor(vertices, faces, mass, zeros(3))

    if align_to_principal
        inertia_tensor, R_cad_body = calc_inertia_y_rotation(inertia_tensor)
    else
        R_cad_body = Matrix{Float64}(I, 3, 3)
    end
    circle_center_z, radius, gamma_tip = find_circle_center_and_radius(vertices)
    le_interp, te_interp, area_interp = create_interpolations(vertices, circle_center_z, radius, gamma_tip, R_cad_body; interp_steps)

    ! prn || @info "Loading 2d polars from $cl_polar_path, $cd_polar_path and $cm_polar_path"
    try
        if !ispath(cl_polar_path) || !ispath(cd_polar_path) || !ispath(cm_polar_path)
            width = 2gamma_tip * radius
            area = area_interp(gamma_tip)
            create_polars(; dat_path, cl_polar_path, cd_polar_path, cm_polar_path, wind_vel,
                area, width, crease_frac, alpha_range, delta_range, remove_nan)
        end

        cl_matrix, _, _ = read_aero_matrix(cl_polar_path)
        cd_matrix, _, _ = read_aero_matrix(cd_polar_path)
        cm_matrix, alpha_range, delta_range = read_aero_matrix(cm_polar_path)

        if remove_nan
            any(isnan.(cl_matrix)) && interpolate_matrix_nans!(cl_matrix; prn)
            any(isnan.(cd_matrix)) && interpolate_matrix_nans!(cd_matrix; prn)
            any(isnan.(cm_matrix)) && interpolate_matrix_nans!(cm_matrix; prn)
        end

        # Create unrefined sections (evenly spaced including both tips)
        sections = Section[]
        aero_data = (collect(alpha_range), collect(delta_range), cl_matrix, cd_matrix, cm_matrix)
        for gamma in range(-gamma_tip, gamma_tip, n_unrefined_sections)
            LE_point = [le_interp[i](gamma) for i in 1:3]
            TE_point = [te_interp[i](gamma) for i in 1:3]
            push!(sections, Section(LE_point, TE_point, POLAR_MATRICES, aero_data))
        end

        panel_props = PanelProperties{n_panels}()
        cache = [PreallocationTools.LazyBufferCache()]

        wing = Wing(n_panels, Int16(n_unrefined_sections), spanwise_distribution, panel_props, MVec3(spanwise_direction),
            sections, Section[], remove_nan, use_prior_polar, 0.0,  # billowing_angle
            Int16[],  # refined_panel_mapping empty
            Section[], zeros(n_panels), zeros(n_panels),  # non_deformed, theta, delta
            mass, gamma_tip, inertia_tensor, T_cad_body, R_cad_body, radius,
            le_interp, te_interp, area_interp, cache)

        # Auto-refine for backward compatibility
        refine_obj_wing!(wing; recompute_mapping=true)
        reinit!(wing)

        wing

    catch e
        if e isa BoundsError
            @error "Delete $cl_polar_path, $cd_polar_path and $cm_polar_path and try again."
        end
        rethrow(e)
    end
end

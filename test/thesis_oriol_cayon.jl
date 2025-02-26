"""
Post-process results to get global forces and aerodynamic coefficients

Parameters
----------
Fmag : Lift, Drag and Moment magnitudes
aero_coeffs : alpha, cl, cd, cm
ringvec : List of dictionaries containing the vectors that define each ring
Uinf : Wind speed velocity vector
controlpoints : List of dictionaries with the variables needed to define each wing section
Atot : Planform area

Returns
-------
F_rel : Lift and drag forces relative to the local angle of attack
F_gl : Lift and drag forces relative to the wind direction
Ltot : Total lift
Dtot : Total drag
CL : Global CL
CD : Global CD
"""
function output_results(Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot)
    rho = 1.225
    alpha = aero_coeffs[:, 1]
    F_rel = []
    F_gl = []
    Fmag_gl = []
    SideF = []
    Ltot = 0.0
    Dtot = 0.0
    SFtot = 0.0
    
    for i in eachindex(alpha)
        r0 = ringvec[i]["r0"]
        # Relative wind speed direction
        dir_urel = cos(alpha[i]) * controlpoints[i]["tangential"] +
                   sin(alpha[i]) * controlpoints[i]["normal"]
        dir_urel = dir_urel / norm(dir_urel)
        
        # Lift direction relative to Urel
        dir_L = cross(dir_urel, r0)
        dir_L = dir_L / norm(dir_L)
        
        # Drag direction relative to Urel
        dir_D = cross([0.0, 1.0, 0.0], dir_L)
        dir_D = dir_D / norm(dir_D)
        
        # Lift and drag relative to Urel
        L_rel = dir_L * Fmag[i, 1]
        D_rel = dir_D * Fmag[i, 2]
        push!(F_rel, [L_rel, D_rel])
        
        # Lift direction relative to the wind speed
        dir_L_gl = cross(Uinf, [0.0, 1.0, 0.0])
        dir_L_gl = dir_L_gl / norm(dir_L_gl)
        
        # Lift and drag relative to the windspeed
        L_gl = vector_projection(L_rel, dir_L_gl) + vector_projection(D_rel, dir_L_gl)
        D_gl = vector_projection(L_rel, Uinf) + vector_projection(D_rel, Uinf)
        push!(F_gl, [L_gl, D_gl])
        
        push!(Fmag_gl, [
            dot(L_rel, dir_L_gl) + dot(D_rel, dir_L_gl),
            dot(L_rel, Uinf / norm(Uinf)) + dot(D_rel, Uinf / norm(Uinf))
        ])
        
        push!(SideF, dot(L_rel, [0.0, 1.0, 0.0]) + dot(D_rel, [0.0, 1.0, 0.0]))
    end

    # Calculate total aerodynamic forces
    for i in eachindex(Fmag_gl)
        Ltot += Fmag_gl[i][1] * norm(ringvec[i]["r0"])
        Dtot += Fmag_gl[i][2] * norm(ringvec[i]["r0"])
        SFtot += SideF[i] * norm(ringvec[i]["r0"])
    end

    v_a = norm(Uinf)
    CL = Ltot / (0.5 * v_a^2 * Atot * rho)
    CD = Dtot / (0.5 * v_a^2 * Atot * rho)
    CS = SFtot / (0.5 * v_a^2 * Atot * rho)

    return F_rel, F_gl, Ltot, Dtot, CL, CD, CS
end

"""
Find the projection of a vector into a direction

Parameters
----------
v : vector to be projected
u : direction

Returns
-------
proj : projection of the vector v onto u
"""
function vector_projection(v, u)
    unit_u = u / norm(u)
    proj = dot(v, unit_u) * unit_u
    return proj
end


"""
Create geometry structures necessary for solving the system of circualtion

Parameters
----------
coordinates : coordinates the nodes (each section is defined by two nodes,
                                     the first is the LE, so each section
                                     defined by a pair of coordinates)
Uinf : Wind speed vector
N : Number of sections
ring_geo :  - '3fil': Each horsehoe is defined by 3 filaments
            - '5fil': Each horseshoe is defined by 5 filaments
model : VSM: Vortex Step method/ LLT: Lifting Line Theory

Returns
-------
controlpoints :  List of dictionaries with the variables needed to define each wing section  
rings : List of list with the definition of each vortex filament
wingpanels : List with the points defining each wing pannel
ringvec : List of dictionaries containing the vectors that define each ring
coord_L : coordinates of the aerodynamic centers of each wing panel
"""
function create_geometry_general(coordinates, Uinf, N, ring_geo, model)
    filaments = []
    controlpoints = []
    rings = []
    wingpanels = []
    ringvec = []
    coord_L = []

    # Go through all wing panels
    for i in 1:N-1
        # Identify points defining the panel
        section = Dict(
            "p1" => coordinates[2i-1, :],
            "p2" => coordinates[2i+1, :],
            "p3" => coordinates[2i+2, :],
            "p4" => coordinates[2i, :]
        )
        push!(wingpanels, section)

        di = norm(
            coordinates[2i-1, :] * 0.75 + 
            coordinates[2i, :] * 0.25 - 
            (coordinates[2i+1, :] * 0.75 + coordinates[2i+2, :] * 0.25)
        )

        if i == 1
            diplus = norm(
            coordinates[2*(i+1)-1, :] * 0.75 + 
            coordinates[2*(i+1), :] * 0.25 - 
            (coordinates[2*(i+1)+1, :] * 0.75 + coordinates[2*(i+1)+2, :] * 0.25)
            )
            ncp = di / (di + diplus)
        elseif i == N-1
            dimin = norm(
            coordinates[2*(i-1)-1, :] * 0.75 + 
            coordinates[2*(i-1), :] * 0.25 - 
            (coordinates[2*(i-1)+1, :] * 0.75 + coordinates[2*(i-1)+2, :] * 0.25)
            )
            ncp = dimin / (dimin + di)
        else
            dimin = norm(
            coordinates[2*(i-1)-1, :] * 0.75 + 
            coordinates[2*(i-1), :] * 0.25 - 
            (coordinates[2*(i-1)+1, :] * 0.75 + coordinates[2*(i-1)+2, :] * 0.25)
            )
            diplus = norm(
            coordinates[2*(i+1)-1, :] * 0.75 + 
            coordinates[2*(i+1), :] * 0.25 - 
            (coordinates[2*(i+1)+1, :] * 0.75 + coordinates[2*(i+1)+2, :] * 0.25)
            )
            ncp = 0.25 * (dimin/(dimin + di) + di/(di + diplus) + 1)
        end

        ncp = 1 - ncp
        chord = norm((section["p2"] + section["p1"])/2 - (section["p3"] + section["p4"])/2)
        
        LLpoint = (section["p2"] * (1-ncp) + section["p1"] * ncp) * 3/4 + 
                  (section["p3"] * (1-ncp) + section["p4"] * ncp) * 1/4
        VSMpoint = (section["p2"] * (1-ncp) + section["p1"] * ncp) * 1/4 + 
                   (section["p3"] * (1-ncp) + section["p4"] * ncp) * 3/4
        push!(coord_L, LLpoint)

        # Define bound vortex filament
        bound = Dict(
            "id" => "bound",
            "x1" => section["p1"] * 3/4 + section["p4"] * 1/4,
            "x2" => section["p2"] * 3/4 + section["p3"] * 1/4,
            "Gamma" => 0
        )
        push!(filaments, bound)

        x_airf = cross(VSMpoint - LLpoint, section["p2"] - section["p1"])
        x_airf = x_airf / norm(x_airf)
        y_airf = VSMpoint - LLpoint
        y_airf = y_airf / norm(y_airf)
        z_airf = bound["x2"] - bound["x1"]
        z_airf = z_airf / norm(z_airf)
        airf_coord = hcat(x_airf, y_airf, z_airf)

        normal = x_airf
        tangential = y_airf
        
        if model === VSM
            cp = Dict(
                "coordinates" => VSMpoint,
                "chord" => chord,
                "normal" => normal,
                "tangential" => tangential,
                "airf_coord" => airf_coord,
                "coordinates_aoa" => LLpoint
            )
            push!(controlpoints, cp)
        elseif model === LLT
            cp = Dict(
                "coordinates" => LLpoint,
                "chord" => chord,
                "normal" => normal,
                "tangential" => tangential,
                "airf_coord" => airf_coord
            )
            push!(controlpoints, cp)
        end

        temp = Dict(
            "r0" => bound["x2"] - bound["x1"],
            "r1" => cp["coordinates"] - bound["x1"],
            "r2" => cp["coordinates"] - bound["x2"],
            "r3" => cp["coordinates"] - (bound["x2"] + bound["x1"])/2
        )
        push!(ringvec, temp)

        temp = Uinf / norm(Uinf)
        if ring_geo == "3fil"
            # create trailing filaments, at x1 of bound filament
            temp1 = Dict("dir" => temp, "id" => "trailing_inf1", "x1" => bound["x1"], "Gamma" => 0)
            push!(filaments, temp1)

            # create trailing filaments, at x2 of bound filament 
            temp1 = Dict("x1" => bound["x2"], "dir" => temp, "id" => "trailing_inf2", "Gamma" => 0)
            push!(filaments, temp1)
        elseif ring_geo == "5fil"
            temp1 = Dict(
                "x1" => section["p4"],
                "x2" => bound["x1"],
                "Gamma" => 0,
                "id" => "trailing1"
            )
            push!(filaments, temp1)

            temp1 = Dict(
                "dir" => temp,
                "id" => "trailing_inf1",
                "x1" => section["p4"],
                "Gamma" => 0
            )
            push!(filaments, temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = Dict(
                "x2" => section["p3"],
                "x1" => bound["x2"],
                "Gamma" => 0,
                "id" => "trailing1"
            )
            push!(filaments, temp1)

            temp1 = Dict(
                "x1" => section["p3"],
                "dir" => temp,
                "id" => "trailing_inf2",
                "Gamma" => 0
            )
            push!(filaments, temp1)
        end

        push!(rings, filaments)
        filaments = []
    end

    coord_L = hcat(coord_L...)
    return controlpoints, rings, wingpanels, ringvec, coord_L
end


"""
Solve the VSM or LLM by finding the distribution of Gamma

Parameters
----------
ringvec : List of dictionaries containing the vectors that define each ring
controlpoints : List of dictionaries with the variables needed to define each wing section 
rings : List of list with the definition of each vortex filament
Uinf : Wind speed velocity vector
data_airf : 2D airfoil data with alpha, Cl, Cd, Cm 
recalc_alpha : True if you want to recalculate the induced angle of attack at 1/4 of the chord (VSM)
Gamma0 : Initial Guess of Gamma
model : VSM: Vortex Step method/ LLT: Lifting Line Theory

Returns
-------
MatrixU, MatrixV, MatrixW : Induction matrices
"""
function thesis_induction_matrix_creation(ringvec, controlpoints, rings, Uinf, Gamma0, data_airf, conv_crit, model)
    nocore = false  # To shut down core corrections input true
    
    # Initialization of parameters
    N = length(rings)
    MatrixU = zeros(N, N)
    MatrixV = zeros(N, N)  
    MatrixW = zeros(N, N)
    U_2D = zeros(N, 3)

    coord_cp = [controlpoints[icp]["coordinates"] for icp in 1:N]
    chord = [controlpoints[icp]["chord"] for icp in 1:N] 
    airf_coord = [controlpoints[icp]["airf_coord"] for icp in 1:N]

    for icp in 1:N
        if model === VSM
            # Velocity induced by infinite bound vortex with Gamma = 1
            U_2D[icp,:] = velocity_induced_bound_2D(ringvec[icp])
        end

        for jring in 1:N
            rings[jring] = update_Gamma_single_ring(rings[jring], 1, 1)
            
            # Calculate velocity induced by ring at control point
            velocity_induced = velocity_induced_single_ring_semiinfinite(
                rings[jring], coord_cp[icp], model, norm(Uinf)
            )
            
            # If CORE corrections are deactivated
            if nocore
                velocity_induced = velocity_induced_single_ring_semiinfinite_nocore(
                    rings[jring], coord_cp[icp], model
                )
            end

            # AIC Matrix 
            MatrixU[icp,jring] = velocity_induced[1]
            MatrixV[icp,jring] = velocity_induced[2]
            MatrixW[icp,jring] = velocity_induced[3]

            # Different from old thesis code as this was considered wrong
            if icp == jring
                if model === VSM
                    MatrixU[icp,jring] -= U_2D[icp,1]
                    MatrixV[icp,jring] -= U_2D[icp,2]
                    MatrixW[icp,jring] -= U_2D[icp,3]
                end
            end
        end
    end

    return MatrixU, MatrixV, MatrixW
end


"""
Calculate the velocity induced by a trailing vortex filament in a point in space

Vortex core correction from:
    Rick Damiani et al. "A vortex step method for nonlinear airfoil polar data as implemented in
KiteAeroDyn".

Parameters
----------
XV1 : Point A of the vortex filament (vector)
XV2 : Point B of the vortex filament (vector)
XVP : Controlpoint (vector)
gamma : Strength of the vortex (scalar)
Uinf : Inflow velocity modulus (scalar)

Returns
-------
vel_ind : induced velocity by the trailing fil. (vector)
"""
function velocity_3D_trailing_vortex(XV1, XV2, XVP, gamma, Uinf)
    r0 = XV2 - XV1  # Vortex filament
    r1 = XVP - XV1  # Controlpoint to one end of the vortex filament
    r2 = XVP - XV2  # Controlpoint to one end of the vortex filament

    alpha0 = 1.25643  # Oseen parameter
    nu = 1.48e-5  # Kinematic viscosity of air
    r_perp = dot(r1, r0) * r0 / (norm(r0)^2)  # Vector from XV1 to XVP perpendicular to the core radius
    epsilon = sqrt(4 * alpha0 * nu * norm(r_perp) / Uinf)  # Cut-off radius

    # Cross products used for later computations
    r1Xr0 = cross(r1, r0)
    r2Xr0 = cross(r2, r0)

    if norm(r1Xr0) / norm(r0) > epsilon  # Perpendicular distance from XVP to vortex filament (r0)
        r1Xr2 = cross(r1, r2)
        vel_ind = gamma / (4π) * r1Xr2 / (norm(r1Xr2)^2) * 
                 dot(r0, r1/norm(r1) - r2/norm(r2))
    else
        # The control point is placed on the edge of the radius core
        # proj stands for the vectors respect to the new controlpoint
        r1_proj = dot(r1, r0) * r0 / (norm(r0)^2) + epsilon * r1Xr0 / norm(r1Xr0)
        r2_proj = dot(r2, r0) * r0 / (norm(r0)^2) + epsilon * r2Xr0 / norm(r2Xr0)
        r1Xr2_proj = cross(r1_proj, r2_proj)
        vel_ind_proj = gamma / (4π) * r1Xr2_proj / (norm(r1Xr2_proj)^2) * 
                      dot(r0, r1_proj/norm(r1_proj) - r2_proj/norm(r2_proj))
        vel_ind = norm(r1Xr0) / (norm(r0) * epsilon) * vel_ind_proj
    end
    
    return vel_ind
end


"""
Calculate the velocity induced by a semiinfinite trailing vortex filament in a point in space

Vortex core correction from:
    Rick Damiani et al. "A vortex step method for nonlinear airfoil polar data as implemented in
KiteAeroDyn".

Parameters
----------
XV1 : Point A of the vortex filament (vector)
Vf : Direction vector of the filament
XVP : Controlpoint (vector)
GAMMA : Strength of the vortex (scalar)
Uinf : Inflow velocity modulus (scalar)

Returns
-------
vel_ind : induced velocity by the trailing fil. (vector)
"""
function velocity_3D_trailing_vortex_semiinfinite(XV1, Vf, XVP, GAMMA, Uinf)
    r1 = XVP - XV1  # Vector from XV1 to XVP
    r1XVf = cross(r1, Vf)

    alpha0 = 1.25643  # Oseen parameter
    nu = 1.48e-5  # Kinematic viscosity of air
    r_perp = dot(r1, Vf) * Vf  # Vector from XV1 to XVP perpendicular to the core radius
    epsilon = sqrt(4 * alpha0 * nu * norm(r_perp) / Uinf)  # Cut-off radius

    if norm(r1XVf) / norm(Vf) > epsilon
        # determine scalar
        K = GAMMA / (4π * norm(r1XVf)^2) * (1 + dot(r1, Vf) / norm(r1))
        # determine the three velocity components
        vel_ind = K * r1XVf
    else
        r1_proj = dot(r1, Vf) * Vf + epsilon * (r1/norm(r1) - Vf) / norm(r1/norm(r1) - Vf)
        r1XVf_proj = cross(r1_proj, Vf)
        K = GAMMA / (4π * norm(r1XVf_proj)^2) * (1 + dot(r1_proj, Vf) / norm(r1_proj))
        # determine the three velocity components
        vel_ind = K * r1XVf_proj
    end
    return vel_ind
end

"""
Calculate induced velocity for 2D bound vortex
"""
function velocity_induced_bound_2D(ringvec)
    r0 = ringvec["r0"]
    r3 = ringvec["r3"]

    cross_prod = cross(r0, r3)
    ind_vel = cross_prod / (dot(cross_prod, cross_prod)) / (2π) * norm(r0)

    return ind_vel
end


"""
Calculates the velocity induced by a ring at a certain controlpoint

Parameters
----------
ring : List of dictionaries defining the filaments of a vortex ring
controlpoint : Dictionary defining a controlpoint
model : VSM: Vortex Step method/ LLT: Lifting Line Theory 
Uinf : Wind speed vector

Returns
-------
velind : Induced velocity
"""
function velocity_induced_single_ring_semiinfinite(ring, controlpoint, model, Uinf)
    velind = [0.0, 0.0, 0.0]
    for filament in ring
        GAMMA = filament["Gamma"]
        XV1 = filament["x1"]
        XVP = controlpoint

        if filament["id"] == "trailing_inf1"
            Vf = filament["dir"]
            tempvel = velocity_3D_trailing_vortex_semiinfinite(
                XV1, Vf, XVP, GAMMA, Uinf
            )
        elseif filament["id"] == "trailing_inf2"
            Vf = filament["dir"]
            tempvel = velocity_3D_trailing_vortex_semiinfinite(
                XV1, Vf, XVP, -GAMMA, Uinf
            )
        elseif filament["id"] == "bound"
            if model === VSM
                XV2 = filament["x2"]
                tempvel = velocity_3D_bound_vortex(XV1, XV2, XVP, GAMMA)
            else
                tempvel = [0.0, 0.0, 0.0]
            end
        else
            XV2 = filament["x2"]
            tempvel = velocity_3D_trailing_vortex(XV1, XV2, XVP, GAMMA, Uinf)
        end

        velind[1] += tempvel[1]
        velind[2] += tempvel[2]
        velind[3] += tempvel[3]
    end

    return velind
end

"""
Calculate the velocity induced by a bound vortex filament in a point in space

Vortex core correction from:
    Rick Damiani et al. "A vortex step method for nonlinear airfoil polar data as implemented in
KiteAeroDyn".

Parameters
----------
XV1 : Point A of Bound vortex (vector)
XV2 : Point B of Bound vortex (vector) 
XVP : Control point (vector)
gamma : Strength of the vortex (scalar)

Returns
-------
vel_ind : Induced velocity (vector)
"""
function velocity_3D_bound_vortex(XV1, XV2, XVP, gamma)
    r0 = XV2 - XV1  # Vortex filament
    r1 = XVP - XV1  # Controlpoint to one end of the vortex filament
    r2 = XVP - XV2  # Controlpoint to one end of the vortex filament

    # Cross products used for later computations
    r1Xr0 = cross(r1, r0)
    r2Xr0 = cross(r2, r0)

    epsilon = 0.05 * norm(r0)  # Cut-off radius

    if norm(r1Xr0) / norm(r0) > epsilon  # Perpendicular distance from XVP to vortex filament (r0)
        r1Xr2 = cross(r1, r2)
        vel_ind = gamma / (4π) * r1Xr2 / (norm(r1Xr2)^2) * 
                 dot(r0, r1/norm(r1) - r2/norm(r2))
    else
        # The control point is placed on the edge of the radius core
        # proj stands for the vectors respect to the new controlpoint
        r1_proj = dot(r1, r0) * r0 / (norm(r0)^2) + epsilon * r1Xr0 / norm(r1Xr0)
        r2_proj = dot(r2, r0) * r0 / (norm(r0)^2) + epsilon * r2Xr0 / norm(r2Xr0)
        r1Xr2_proj = cross(r1_proj, r2_proj)
        vel_ind_proj = gamma / (4π) * r1Xr2_proj / (norm(r1Xr2_proj)^2) * 
                      dot(r0, r1_proj/norm(r1_proj) - r2_proj/norm(r2_proj))
        vel_ind = norm(r1Xr0) / (norm(r0) * epsilon) * vel_ind_proj
    end
    
    return vel_ind
end

"""
Update Gamma of all the filaments in a horshoe ring
"""
function update_Gamma_single_ring(ring, GammaNew, WeightNew)
    # Runs through all filaments
    for filament in ring
        filament["Gamma"] = filament["Gamma"] * (1 - WeightNew) + WeightNew * GammaNew
    end
    return ring
end

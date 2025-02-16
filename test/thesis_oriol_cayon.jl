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

    Umag = norm(Uinf)
    CL = Ltot / (0.5 * Umag^2 * Atot * rho)
    CD = Dtot / (0.5 * Umag^2 * Atot * rho)
    CS = SFtot / (0.5 * Umag^2 * Atot * rho)

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

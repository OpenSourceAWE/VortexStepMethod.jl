"""
    lei_poly_coeffs(tube_diameter, camber) -> (cl_coeffs, cd_coeffs, cm_coeffs)

Breukels leading-edge-inflatable α-polynomial coefficients (α in **degrees**) for a
section with normalized `tube_diameter` and `camber`. `cl_coeffs` is a cubic (4
coefficients, ascending order), `cd_coeffs`/`cm_coeffs` are quadratic (3). Feed the
result to a core `POLY` section as `aero_data = (cl_coeffs, cd_coeffs, cm_coeffs)`.
"""
function lei_poly_coeffs(tube_diameter, camber)
    t = tube_diameter
    k = camber
    C = Dict(
        20 => -0.008011, 21 => -0.000336, 22 => 0.000992,
        23 => 0.013936, 24 => -0.003838, 25 => -0.000161,
        26 => 0.001243, 27 => -0.009288, 28 => -0.002124,
        29 => 0.012267, 30 => -0.002398, 31 => -0.000274,
        32 => 0.0, 33 => 0.0, 34 => 0.0,
        35 => -3.371000, 36 => 0.858039, 37 => 0.141600,
        38 => 7.201140, 39 => -0.676007, 40 => 0.806629,
        41 => 0.170454, 42 => -0.390563, 43 => 0.101966
    )
    S = Dict{Int64,Float64}()
    S[9] = C[20]*t^2 + C[21]*t + C[22]
    S[10] = C[23]*t^2 + C[24]*t + C[25]
    S[11] = C[26]*t^2 + C[27]*t + C[28]
    S[12] = C[29]*t^2 + C[30]*t + C[31]
    S[13] = C[32]*t^2 + C[33]*t + C[34]
    S[14] = C[35]*t^2 + C[36]*t + C[37]
    S[15] = C[38]*t^2 + C[39]*t + C[40]
    S[16] = C[41]*t^2 + C[42]*t + C[43]

    cl_coeffs = [
        S[9]*k + S[10],
        S[11]*k + S[12],
        S[13]*k + S[14],
        S[15]*k + S[16],
    ]
    cd_coeffs = [
        ((0.546094*t + 0.022247)*k^2 +
         (-0.071462*t - 0.006527)*k +
         (0.002733*t + 0.000686)),
        0.0,
        ((0.123685*t + 0.143755)*k +
         (0.495159*t^2 - 0.105362*t + 0.033468)),
    ]
    cm_coeffs = [
        ((-0.284793*t - 0.026199)*k +
         (-0.024060*t + 0.000559)),
        0.0,
        ((-1.787703*t + 0.352443)*k +
         (-0.839323*t + 0.137932)),
    ]
    return cl_coeffs, cd_coeffs, cm_coeffs
end

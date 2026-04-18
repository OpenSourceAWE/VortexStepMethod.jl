"""
    parse_literature_column(col)

Parse a column of mixed-type data (Real or String) into Float64,
returning NaN for unparseable values.
"""
function parse_literature_column(col)
    [v isa Real ? Float64(v) :
        (y = tryparse(Float64, strip(string(v)));
         isnothing(y) ? NaN : y) for v in col]
end

"""
    extract_literature_polar_data(raw_data, path;
        angle_type="angle_of_attack")

Extract polar data from literature CSV data (as returned by
`readdlm`). Returns a named tuple with `polar_data` (array of
[angle, cl, cd, cs] vectors) and `cmx`, `cmy`, `cmz` vectors.

Supports both tuple format `(table, header)` from
`readdlm(...; header=true)` and matrix format where the first
row contains headers.
"""
function extract_literature_polar_data(raw_data, path;
        angle_type="angle_of_attack")
    table, header = if raw_data isa Tuple
        raw_table, raw_header = raw_data
        raw_table, lowercase.(strip.(string.(vec(raw_header))))
    else
        raw_data[2:end, :],
            lowercase.(strip.(string.(raw_data[1, :])))
    end

    # Find angle column based on angle_type
    contains_token(token, value) = occursin(token, lowercase(strip(string(value))))
    angle_idx = if angle_type == "side_slip"
        findfirst(
            x -> contains_token("beta", x) ||
                 contains_token("side_slip", x), header)
    else
        findfirst(
            x -> contains_token("alpha", x) ||
                 lowercase(strip(string(x))) == "aoa", header)
    end
    cl_idx = findfirst(x -> contains_token("cl", x), header)
    cd_idx = findfirst(x -> contains_token("cd", x), header)
    cs_idx = findfirst(x -> contains_token("cs", x), header)
    cmx_idx = findfirst(x -> contains_token("cmx", x), header)
    cmy_idx = findfirst(x -> contains_token("cmy", x), header)
    cmz_idx = findfirst(x -> contains_token("cmz", x), header)

    (isnothing(angle_idx) || isnothing(cl_idx) ||
        isnothing(cd_idx)) &&
        throw(ArgumentError(
            "Literature CSV must contain alpha/aoa (or " *
            "beta for side_slip), cl and cd columns: $path"))

    n_rows = size(table, 1)
    get_col(idx) = parse_literature_column(table[:, idx])
    nan_col() = fill(NaN, n_rows)

    return (
        polar_data=[
            get_col(angle_idx),
            get_col(cl_idx),
            get_col(cd_idx),
            cs_idx === nothing ? zeros(n_rows) : get_col(cs_idx)
        ],
        cmx=cmx_idx === nothing ? nan_col() : get_col(cmx_idx),
        cmy=cmy_idx === nothing ? nan_col() : get_col(cmy_idx),
        cmz=cmz_idx === nothing ? nan_col() : get_col(cmz_idx)
    )
end

"""
    generate_polar_data(solver, body_aero, angle_range;
        angle_type="angle_of_attack", angle_of_attack=0.0,
        side_slip=0.0, v_a=10.0)

Sweep over `angle_range` (degrees), solving at each angle. Returns
a named tuple `(polar_data, cmx, cmy, cmz, rey)` where `polar_data`
is `[angle, cl, cd, cs, gamma_dist, cl_dist, cd_dist, cs_dist,
reynolds]`.
"""
function generate_polar_data(
    solver,
    body_aero::BodyAerodynamics,
    angle_range;
    angle_type="angle_of_attack",
    angle_of_attack=0.0,
    side_slip=0.0,
    v_a=10.0
)
    n_panels = length(body_aero.panels)
    n_angles = length(angle_range)

    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    cmx = fill(NaN, n_angles)
    cmy = fill(NaN, n_angles)
    cmz = fill(NaN, n_angles)
    gamma_distribution = zeros(n_angles, n_panels)
    cl_distribution = zeros(n_angles, n_panels)
    cd_distribution = zeros(n_angles, n_panels)
    cs_distribution = zeros(n_angles, n_panels)
    reynolds_number = zeros(n_angles)

    for (i, angle_i) in enumerate(angle_range)
        if angle_type == "angle_of_attack"
            α = deg2rad(angle_i)
            β = side_slip
        elseif angle_type == "side_slip"
            α = angle_of_attack
            β = deg2rad(angle_i)
        else
            throw(ArgumentError(
                "angle_type must be 'angle_of_attack' " *
                "or 'side_slip'"))
        end

        set_va!(body_aero,
            [cos(α) * cos(β), sin(β), sin(α)] * v_a)

        results = solve(solver, body_aero,
            gamma_distribution[i, :])

        cl[i] = results["cl"]
        cd[i] = results["cd"]
        cs[i] = results["cs"]
        cmx[i] = get(results, "cmx", NaN)
        cmy[i] = get(results, "cmy", NaN)
        cmz[i] = get(results, "cmz", NaN)
        gamma_distribution[i, :] =
            results["gamma_distribution"]
        cl_distribution[i, :] = results["cl_distribution"]
        cd_distribution[i, :] = results["cd_distribution"]
        cs_distribution[i, :] = results["cs_distribution"]
        reynolds_number[i] = results["Rey"]
    end

    polar_data = [
        angle_range, cl, cd, cs,
        gamma_distribution, cl_distribution,
        cd_distribution, cs_distribution,
        reynolds_number
    ]

    return (polar_data=polar_data, cmx=cmx, cmy=cmy, cmz=cmz,
            rey=reynolds_number[1])
end

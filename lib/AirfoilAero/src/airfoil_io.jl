"""
    read_dat_coordinates(path) -> (x, y)

Read airfoil coordinates from a Selig-format `.dat` file, skipping header and
comment lines.
"""
function read_dat_coordinates(path::String)
    x = Float64[]
    y = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || !(isdigit(s[1]) || s[1] == '-' || s[1] == '.')) && continue
        parts = split(s)
        length(parts) >= 2 || continue
        xp = tryparse(Float64, parts[1])
        yp = tryparse(Float64, parts[2])
        (xp === nothing || yp === nothing) && continue
        push!(x, xp)
        push!(y, yp)
    end
    return x, y
end

"""
    write_dat(filepath, name, x, y) -> filepath

Write airfoil coordinates to a Selig-format `.dat` file.
"""
function write_dat(filepath::String, name::String, x::Vector, y::Vector)
    open(filepath, "w") do io
        println(io, name)
        for i in eachindex(x)
            println(io, @sprintf("%.8f %.8f", x[i], y[i]))
        end
    end
    return filepath
end

"""
    write_polar_csv(filepath, result::NeuralFoilResult)

Write a NeuralFoil result to a `POLAR_VECTORS` CSV (`alpha, Cd, Cs, Cl, Cm`).
"""
function write_polar_csv(filepath::String, result::NeuralFoilResult)
    n = length(result.alpha)
    open(filepath, "w") do io
        println(io, "alpha,Cd,Cs,Cl,Cm")
        for i in 1:n
            row = (result.alpha[i], result.CD[i], 0.0, result.CL[i], result.CM[i])
            println(io, join((@sprintf("%.16g", v) for v in row), ","))
        end
    end
    return filepath
end

"""
    generate_polar_from_coordinates(x, y, output_path; Re, alpha_range=-180:1:180, kwargs...)

Fit Kulfan parameters to `(x, y)`, evaluate NeuralFoil over `alpha_range` (degrees),
and write the polar CSV. Extra `kwargs` pass to [`neuralfoil_aero`](@ref). Returns
the `NeuralFoilResult`.
"""
function generate_polar_from_coordinates(x::Vector, y::Vector, output_path::String;
                                         Re::Real, alpha_range=-180:1:180, kwargs...)
    alphas = collect(Float64, alpha_range)
    params = fit_kulfan_parameters(x, y)
    result = neuralfoil_aero(params, alphas, Float64(Re); kwargs...)
    write_polar_csv(output_path, result)
    return result
end

"""
    generate_polar_from_dat(dat_path, output_path; Re, kwargs...)

Read a `.dat` airfoil and generate its polar CSV via
[`generate_polar_from_coordinates`](@ref). This is AirfoilAero's `.dat` → CSV entry
point.
"""
function generate_polar_from_dat(dat_path::String, output_path::String;
                                 Re::Real, kwargs...)
    isfile(dat_path) || error("DAT file not found: $dat_path")
    x, y = read_dat_coordinates(dat_path)
    isempty(x) && error("No valid coordinates found in $dat_path")
    return generate_polar_from_coordinates(x, y, output_path; Re, kwargs...)
end

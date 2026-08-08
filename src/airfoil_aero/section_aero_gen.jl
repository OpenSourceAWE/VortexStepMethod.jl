"""
    generate_airfoil_aero(solver, base; alpha_range, delta_range, reynolds_number,
                          crease_frac=0.9, remove_nan=true) -> (SectionAero, sols)

Run **one** solver sweep of a base airfoil (Kulfan) over the `(alpha, delta)` grid
(radians) and return both the [`SectionAero`](@ref) (contour + `Cp` + `cf` per node) and
the raw `sols::Vector{Vector{SectionSolution}}` (one inner vector per delta). The `sols`
also carry `cl/cd/cm`, so a caller can write the polar from the same sweep — this is how
`obj_to_yaml` avoids a second sweep. Non-converged points stay `NaN` and, when
`remove_nan`, are filled per node with `interpolate_matrix_nans!`.
"""
function generate_airfoil_aero(solver::AbstractAirfoilSolver, base::KulfanParameters;
        alpha_range, delta_range, reynolds_number, crease_frac=0.9, remove_nan=true)
    x0, y0 = kulfan_to_coordinates(base)
    n_alpha, n_delta = length(alpha_range), length(delta_range)
    sols = Vector{Vector{SectionSolution}}(undef, n_delta)
    for (jd, delta) in enumerate(delta_range)
        def = deform_section(x0, y0, delta; crease_frac)
        sols[jd] = analyze_sweep(solver, def, alpha_range, reynolds_number)
    end

    lens = [length(s.cp) for col in sols for s in col if !isempty(s.cp)]
    isempty(lens) && throw(ErrorException("generate_airfoil_aero: no converged angle"))
    # the modal node count, so a single degenerate wrap cannot drop every other column
    n_node = argmax(u -> count(==(u), lens), unique(lens))

    x = fill(NaN, n_node, n_delta)
    y = fill(NaN, n_node, n_delta)
    cp = fill(NaN, n_node, n_alpha, n_delta)
    cf = fill(NaN, n_node, n_alpha, n_delta)
    for jd in 1:n_delta, (ia, s) in enumerate(sols[jd])
        (isempty(s.cp) || length(s.cp) != n_node) && continue
        cp[:, ia, jd] .= s.cp
        cf[:, ia, jd] .= s.cf
        if any(isnan, view(x, :, jd))
            x[:, jd] .= s.x
            y[:, jd] .= s.y
        end
    end

    if remove_nan
        for i in 1:n_node
            fill_node_nans!(cp, i)
            fill_node_nans!(cf, i)
        end
    end
    aero = SectionAero(collect(float.(alpha_range)), collect(float.(delta_range)),
                       x, y, cp, cf)
    return aero, sols
end

"""
    generate_airfoil_aero(solver, x::Vector, y::Vector; kwargs...) -> (SectionAero, sols)

Convenience: [`shrink_wrap`](@ref) the coordinates and fit base Kulfan parameters
([`LeastSquaresFit`](@ref)) first.
"""
function generate_airfoil_aero(solver::AbstractAirfoilSolver, x::Vector, y::Vector;
        kwargs...)
    xw, yw = shrink_wrap(x, y, ShrinkWrap())
    return generate_airfoil_aero(solver, fit_kulfan_parameters(xw, yw, LeastSquaresFit());
                                 kwargs...)
end

"""
    generate_section_aero(solver, base_or_coords; kwargs...) -> SectionAero

Just the [`SectionAero`](@ref) from [`generate_airfoil_aero`](@ref) (drops the raw sols).
"""
generate_section_aero(solver::AbstractAirfoilSolver, base::KulfanParameters; kwargs...) =
    generate_airfoil_aero(solver, base; kwargs...)[1]

function generate_section_aero(solver::AbstractAirfoilSolver, x::Vector, y::Vector;
        kwargs...)
    xw, yw = shrink_wrap(x, y, ShrinkWrap())
    return generate_section_aero(solver, fit_kulfan_parameters(xw, yw, LeastSquaresFit());
                                 kwargs...)
end

"""
    fill_node_nans!(grid, i)

Fill `NaN`s in node `i`'s `(alpha, delta)` matrix `grid[i, :, :]` with
`interpolate_matrix_nans!`; an all-`NaN` node is left untouched.
"""
function fill_node_nans!(grid, i)
    m = grid[i, :, :]
    all(isnan, m) && return grid
    interpolate_matrix_nans!(m; prn=false)
    grid[i, :, :] = m
    return grid
end

"""
    write_node_table(path, aero, values) -> path

Write a per-node aero table (`Cp` or `cf`, shaped `n_node × n_alpha × n_delta`), one row
per `(alpha, delta)` with angles in degrees. The file suffix picks the format: `.arrow`
(columns `alpha`, `delta` and a per-row list column `values`, an order of magnitude
faster to load), anything else a human-readable CSV with header `alpha, delta, n0, n1, …`
(node columns in the contour node order). [`read_node_table`](@ref) reads both.
"""
function write_node_table(path::AbstractString, aero::SectionAero, values)
    if endswith(String(path), ".arrow")
        grid = [(jd, ia) for jd in eachindex(aero.delta_range)
                         for ia in eachindex(aero.alpha_range)]
        Arrow.write(String(path),
            (alpha=[rad2deg(aero.alpha_range[ia]) for (_, ia) in grid],
             delta=[rad2deg(aero.delta_range[jd]) for (jd, _) in grid],
             values=[values[:, ia, jd] for (jd, ia) in grid]))
        return path
    end
    open(String(path), "w") do io
        println(io, "alpha,delta," * join(("n$(k - 1)" for k in 1:size(values, 1)), ","))
        for jd in eachindex(aero.delta_range), ia in eachindex(aero.alpha_range)
            row = [rad2deg(aero.alpha_range[ia]), rad2deg(aero.delta_range[jd])]
            append!(row, values[:, ia, jd])
            println(io, join(row, ","))
        end
    end
    return path
end

"""
    write_section_aero(dat_prefix, aero::SectionAero; table_prefix=dat_prefix,
                       table_format=:csv) -> (dat, cp_table, cf_table)

Write a [`SectionAero`](@ref) as human-readable files. The airfoil contours share
`dat_prefix`: `{dat_prefix}.dat` (contour at `delta=0`) plus
`{dat_prefix}_{delta_suffix(δ)}.dat` per non-zero deflection. The per-node `Cp`/`cf`
tables share `table_prefix` (defaults to `dat_prefix`): `{table_prefix}_cp.{ext}` and
`{table_prefix}_cf.{ext}`, where `ext` is `table_format`, either `:csv` (default,
readable) or `:arrow` (binary, an order of magnitude faster to load). Pass a separate
`table_prefix` to keep the surface tables out of the airfoil-shape directory.
`read_section_aero` reads either format back, detecting it from the suffix. The single
writer for surface aero (submodule side); loading lives in the main package.
"""
function write_section_aero(dat_prefix::AbstractString, aero::SectionAero;
                            table_prefix::AbstractString=dat_prefix,
                            table_format::Symbol=:csv)
    table_format in (:csv, :arrow) ||
        throw(ArgumentError("table_format must be :csv or :arrow, got :$table_format"))
    mkpath(dirname(dat_prefix))
    mkpath(dirname(table_prefix))
    for (jd, d) in enumerate(aero.delta_range)
        iszero(d) && continue
        write_dat("$(dat_prefix)_$(delta_suffix(d)).dat", "section",
                  aero.x[:, jd], aero.y[:, jd])
    end
    # never write `$dat_prefix.dat` all-NaN (read_dat_coordinates would drop it to empty)
    finite(jd) = any(isfinite, view(aero.x, :, jd))
    jz = findfirst(iszero, aero.delta_range)
    jdat = jz !== nothing && finite(jz) ? jz : findfirst(finite, eachindex(aero.delta_range))
    write_dat("$dat_prefix.dat", "section", aero.x[:, jdat], aero.y[:, jdat])
    cp_path = "$(table_prefix)_cp.$table_format"
    cf_path = "$(table_prefix)_cf.$table_format"
    write_node_table(cp_path, aero, aero.cp)
    write_node_table(cf_path, aero, aero.cf)
    return "$dat_prefix.dat", cp_path, cf_path
end

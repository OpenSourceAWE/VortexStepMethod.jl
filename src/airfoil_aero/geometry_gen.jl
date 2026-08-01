"""
    generate_airfoils(airfoils, output_dir; Re, alpha_range=-180:1:180,
                      delta_range=nothing, aero_solver=NeuralFoilSolver(),
                      reuse_valid_airfoils=true, crease_frac=0.75, verbose=true)
        -> (airfoil_rows, ok)

Run the 2D solver over a set of already-shrink-wrapped airfoils and write the per
section files each geometry route references. Shared by the `.obj` and Surfplan
adapters: the front-end supplies the clean airfoils and their leading/trailing-edge
placement; this writes the surface pressure/friction tables, polars and airfoil
`.dat`s, and returns the `wing_airfoils` rows plus the ids that solved.

`airfoils` is a vector of `(; id, x_fit, y_fit, x_raw, y_raw)`: `x_fit`/`y_fit` is
the wrapped airfoil the solver analyses; `x_raw`/`y_raw` the raw points it enclosed.

Writes into `output_dir` (indexed by `id`): `airfoils/{id}.dat` (wrapped shape),
`airfoils/{id}_raw.dat` (raw points), `airfoils/{id}_cp.csv` / `_cf.csv` (per-node
surface pressure and skin friction), and `polars/{id}.csv` (`POLAR_VECTORS`, or a
`POLAR_MATRICES` grid when `delta_range` is set). `aero_solver` selects the backend
([`NeuralFoilSolver`](@ref) default, [`XFoilSolver`](@ref) opt-in).

With `reuse_valid_airfoils=true` an airfoil the solver cannot converge is skipped
(the caller maps its sections to the nearest solved id); otherwise it errors.
"""
function generate_airfoils(airfoils, output_dir::String;
        Re::Real, alpha_range=-180:1:180, delta_range=nothing,
        aero_solver::AbstractAirfoilSolver=NeuralFoilSolver(),
        reuse_valid_airfoils::Bool=true, crease_frac=0.75, verbose::Bool=true)
    mkpath(joinpath(output_dir, "airfoils"))
    mkpath(joinpath(output_dir, "polars"))
    airfoil_rows = Vector{Any}[]
    ok = Int[]
    for af in airfoils
        j = af.id
        raw_rel = joinpath("airfoils", "$(j)_raw.dat")
        csv_rel = joinpath("polars", "$j.csv")
        try
            alphas = deg2rad.(collect(Float64, alpha_range))
            deltas = isnothing(delta_range) ? [0.0] :
                deg2rad.(collect(Float64, delta_range))
            aero, sols = generate_airfoil_aero(aero_solver, af.x_fit, af.y_fit;
                alpha_range=alphas, delta_range=deltas,
                reynolds_number=Float64(Re), crease_frac)
            clvals = collect(sol.cl for sol in sols[1])
            all(isnan, clvals) && error("solver produced no converged points")
            if isnothing(delta_range)
                write_polar_csv(joinpath(output_dir, csv_rel), sols[1])
            else
                coeff(f) = [f(sols[jd][ia]) for ia in eachindex(alphas),
                            jd in eachindex(deltas)]
                write_polar_matrix_csv(joinpath(output_dir, csv_rel), alphas, deltas,
                    coeff(s -> s.cl), coeff(s -> s.cd), coeff(s -> s.cm))
            end
            paths = write_section_aero(joinpath(output_dir, "airfoils", "$j"), aero)
            dat_rel, cp_rel, cf_rel = (relpath(p, output_dir) for p in paths)
            write_dat(joinpath(output_dir, raw_rel), "section_$(j)_raw",
                      af.x_raw, af.y_raw)
            push!(airfoil_rows, Any[j, "polar_vectors",
                Dict("dat_file" => dat_rel, "raw_dat_file" => raw_rel,
                     "csv_file_path" => csv_rel, "cp_file" => cp_rel,
                     "cf_file" => cf_rel)])
            push!(ok, j)
            finite_cl = [v for v in clvals if !isnan(v)]
            cl_max = isempty(finite_cl) ? NaN : maximum(finite_cl)
            verbose && println("  Airfoil $j: CL_max=$(round(cl_max, digits=2))")
        catch e
            reuse_valid_airfoils ||
                error("Airfoil $j polar generation failed: $(sprint(showerror, e))")
            @warn "Airfoil $j did not solve, reusing the nearest valid airfoil: " *
                sprint(showerror, e)
        end
    end
    return airfoil_rows, ok
end

using LinearAlgebra
using Test
using VortexStepMethod
using VortexStepMethod: calculate_AIC_matrices!, quarter_chord_AIC!, relaxation_factor_limit
using YAML

const PYTHON_REFERENCE = YAML.load_file(
    joinpath(@__DIR__, "..", "data", "python_reference", "reference.yaml"))

"""Rows of a nested list as a matrix."""
as_matrix(rows) = reduce(vcat, permutedims(Float64.(row)) for row in rows)

"""`POLAR_VECTORS` aero data of a reference section's `[alpha, cl, cd, cm]` table."""
function polar_vectors(table)
    polar = as_matrix(table)
    return (polar[:, 1], polar[:, 2], polar[:, 3], polar[:, 4])
end

"""Body and solver of a reference `case`, with the Python call's inflow and options."""
function reference_setup(case)
    n_panels = case["n_panels"]
    wing = Wing(n_panels; spanwise_distribution=LINEAR)
    for section in case["sections"]
        add_section!(wing, Float64.(section["LE"]), Float64.(section["TE"]), POLAR_VECTORS,
            polar_vectors(section["polar"]))
    end
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    va = case["va"]
    if first(va) isa AbstractVector
        va_low, va_high = Float64.(va[1]), Float64.(va[2])
        va_dist = [va_low[k] + (va_high[k] - va_low[k]) * (i - 1) / (n_panels - 1)
                   for i in 1:n_panels, k in 1:3]
        set_va!(body_aero, va_dist)
    else
        rates = case["body_rates"]
        omega = isnothing(rates) ? zeros(3) : Float64.(rates)
        set_va!(body_aero, Float64.(va), omega; reference_point=Float64.(case["pivot"]))
    end
    model = get(case["options"], "aerodynamic_model_type", "VSM") == "LLT" ? LLT : VSM
    solver = Solver(n_panels, length(case["sections"]);
        aerodynamic_model_type=model, rtol=1e-10, max_iterations=400_000,
        reference_point=Float64.(case["reference_point"]),
        is_with_viscous_drag_correction=get(case["options"],
            "is_with_viscous_drag_correction", false))
    return body_aero, solver, model
end

"""Julia's `(i, j, k)` influence array in the `(k, i, j)` nesting of the reference."""
reference_order(matrices) = [matrices[k][i][j] for i in eachindex(matrices[1]),
    j in eachindex(matrices[1]), k in eachindex(matrices)]

@testset "The default solve agrees with the Python VSM v2.3.0" begin
    @testset "Solver defaults are the Python ones" begin
        defaults = PYTHON_REFERENCE["defaults"]
        solver = Solver(4, 2)
        @test solver.max_iterations == defaults["max_iterations"]
        @test solver.core_radius_fraction == defaults["core_radius_fraction"]
        @test solver.correct_aoa == defaults["is_aoa_corrected"]
        @test solver.is_with_attached_trailed_force ==
            defaults["is_with_attached_trailed_vortex_force"]
        @test isnothing(solver.relaxation_factor)
    end

    @testset "$(case["name"])" for case in PYTHON_REFERENCE["cases"]
        body_aero, solver, model = reference_setup(case)
        panels = body_aero.panels
        sol = solve!(solver, body_aero; log=false)
        @test solver.lr.converged

        force_scale = norm(case["force"])
        moment_scale = norm(case["moment"])
        @testset "panels have the perpendicular chord and an orthonormal frame" begin
            @test [panel.chord for panel in panels] ≈ Float64.(case["chord"]) rtol = 1e-12
            @test [panel.width for panel in panels] ≈ Float64.(case["width"]) rtol = 1e-9
            for (i, panel) in enumerate(panels)
                @test panel.x_airf ≈ Float64.(case["chordwise"][i]) atol = 1e-12
                @test panel.y_airf ≈ Float64.(case["spanwise"][i]) atol = 1e-12
                @test panel.z_airf ≈ Float64.(case["normal"][i]) atol = 1e-12
            end
        end
        @testset "each panel feels the Python inflow" begin
            @test sol.va_vec_dist ≈ as_matrix(case["panel_va"]) rtol = 1e-12
        end
        if !isnothing(case["AIC"])
            @testset "influence matrices match" begin
                @test body_aero.AIC ≈ reference_order(case["AIC"]) rtol = 1e-10
                correction = zeros(size(body_aero.AIC))
                quarter_chord_AIC!(correction, body_aero, solver.core_radius_fraction,
                    solver.br.va_dist, solver.br.va_unit_dist)
                @test correction ≈ reference_order(case["AIC_quarter_chord"]) rtol = 1e-10
            end
        end
        @testset "circulation and angles match" begin
            @test sol.gamma_distribution ≈ Float64.(case["gamma"]) rtol = 1e-6
            @test solver.lr.alpha_dist ≈ Float64.(case["alpha"]) rtol = 1e-6
            @test sol.alpha_dist ≈ Float64.(case["alpha_at_ac"]) rtol = 1e-6
        end
        @testset "panel forces and moments match" begin
            forces = as_matrix(case["F_distribution"])'
            moments = as_matrix(case["M_distribution"])'
            @test maximum(abs, sol.f_body_3D .- forces) < 1e-6 * force_scale
            @test maximum(abs, sol.m_body_3D .- moments) < 1e-6 * moment_scale
            @test sol.cl_distribution ≈ Float64.(case["cl_distribution"]) rtol = 1e-6
            @test sol.cd_distribution ≈ Float64.(case["cd_distribution"]) rtol = 1e-5
        end
        @testset "totals and coefficients match" begin
            @test sol.force ≈ Float64.(case["force"]) atol = 1e-6 * force_scale
            @test sol.moment ≈ Float64.(case["moment"]) atol = 1e-6 * moment_scale
            @test [sol.lift, sol.drag, sol.side] ≈ Float64.(case["lift_drag_side"]) atol =
                1e-6 * force_scale
            @test [sol.cl, sol.cd, sol.cs] ≈ Float64.(case["force_coefficients"]) atol =
                1e-6
            @test sol.moment_coeffs ≈ Float64.(case["moment_coefficients"]) atol = 1e-5
            @test sol.q_ref ≈ case["q_ref"] rtol = 1e-12
            @test sol.projected_area ≈ case["projected_area"] rtol = 1e-12
            @test sol.area_all_panels ≈ case["area_all_panels"] rtol = 1e-12
        end
        if case["name"] != "v3"
            @testset "relaxation is 0.8 of the stability limit" begin
                @test 0.8 * relaxation_factor_limit(panels, model) ≈
                    case["relaxation_factor"] rtol = 1e-3
            end
        end
    end
end

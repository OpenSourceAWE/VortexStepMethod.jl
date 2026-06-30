using Test
using VortexStepMethod
using VortexStepMethod: KulfanParameters, fit_kulfan_parameters, kulfan_to_coordinates,
                       neuralfoil_aero

read_dat_coords(path) = begin
    x = Float64[]; y = Float64[]
    for ln in eachline(path)
        s = strip(ln)
        (isempty(s) || !(isdigit(s[1]) || s[1] == '-' || s[1] == '.')) && continue
        p = split(s); length(p) >= 2 || continue
        push!(x, parse(Float64, p[1])); push!(y, parse(Float64, p[2]))
    end
    x, y
end

@testset "Kulfan fit and NeuralFoil" begin
    @testset "Round-trip recovers known parameters" begin
        truth = KulfanParameters(fill(0.2, 8), fill(-0.2, 8), 0.0, 0.0)
        x, y = kulfan_to_coordinates(truth; n_points=120)
        fit = fit_kulfan_parameters(collect(x), collect(y))
        @test maximum(abs.(fit.upper_weights .- truth.upper_weights)) < 1e-5
        @test maximum(abs.(fit.lower_weights .- truth.lower_weights)) < 1e-5
        @test abs(fit.leading_edge_weight) < 1e-5
        @test abs(fit.TE_thickness) < 1e-10
    end

    dat = joinpath(@__DIR__, "..", "..", "data", "ram_air_kite", "ram_air_kite_foil.dat")
    xr, yr = read_dat_coords(dat)
    params = fit_kulfan_parameters(xr, yr)

    @testset "Fit matches aerosandbox get_kulfan_parameters" begin
        # Reference values from aerosandbox 4.2.9 get_kulfan_parameters(n_weights_per_side=8)
        ref_upper = [0.19725544843504217, 0.010915219977460575, 0.5941307473734625,
                     -0.27110843470789603, 0.5355058274193715, -0.17717585452753137,
                     0.23021043555024268, -0.00473923052893821]
        ref_lower = [-0.3653153052872525, -0.5523615615783967, 0.13925743643180583,
                     -0.5611284190191309, 0.11545728044627838, -0.3801887792186299,
                     -0.022461606251837227, -0.176490579153197]
        @test maximum(abs.(params.upper_weights .- ref_upper)) < 1e-9
        @test maximum(abs.(params.lower_weights .- ref_lower)) < 1e-9
        @test abs(params.leading_edge_weight - 0.9180525576877088) < 1e-9
        @test params.TE_thickness ≈ 0.0 atol = 1e-12
    end

    @testset "NeuralFoil matches Python neuralfoil (xlarge)" begin
        # Reference CL/CD/CM from neuralfoil 0.3.2 get_aero_from_kulfan_parameters, Re=5e5
        alphas = Float64.(-10:2:20)
        ref_CL = [-1.1223364415631283, -0.9446879205845478, -0.7397, -0.5262, -0.3103,
                  -0.091, 0.1345, 0.3571, 0.5788, 0.8821, 1.0837, 1.2544, 1.4,
                  1.4226, 1.3811, 1.2739]
        res = neuralfoil_aero(params, alphas, 5e5; model_size="xlarge")
        # Compare against full-precision endpoints + rounded interior to 1e-3
        @test res.CL[1] ≈ -1.1223364415631283 rtol = 1e-4
        @test res.CL[6] ≈ -0.091 atol = 1e-3
        @test res.CL[14] ≈ 1.4226 atol = 1e-3
        @test maximum(abs.(round.(res.CL; digits=4) .- ref_CL)) < 2e-3
    end
end

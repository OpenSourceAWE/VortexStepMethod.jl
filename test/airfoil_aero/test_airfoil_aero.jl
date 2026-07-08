using Test
using VortexStepMethod.AirfoilAero
import VortexStepMethod
using VortexStepMethod.AirfoilAero: KulfanParameters, LeastSquaresFit, ShrinkWrap,
                       shrink_wrap, fit_kulfan_parameters, kulfan_to_coordinates,
                       neuralfoil_aero, class_function, bernstein_basis,
                       leading_edge_basis, normalize_airfoil

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

    dat = joinpath(@__DIR__, "data", "test_airfoil.dat")
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

    @testset "Shrink-wrap encloses points with clearance" begin
        xn, yn, _ = normalize_airfoil(collect(float.(xr)), collect(float.(yr)))
        mid = (xn .> 0.05) .& (xn .< 0.95)

        interp(xs, ys, px) = begin
            px <= xs[1] && return ys[1]
            px >= xs[end] && return ys[end]
            i = searchsortedlast(xs, px)
            t = (px - xs[i]) / (xs[i+1] - xs[i])
            (1 - t) * ys[i] + t * ys[i+1]
        end

        function min_clearance(method)
            xw, yw = shrink_wrap(xr, yr, method)
            le = argmin(xw)
            xu, yu = reverse(xw[1:le]), reverse(yw[1:le])
            xl, yl = xw[le:end], yw[le:end]
            upper = [interp(xu, yu, px) for px in xn]
            lower = [interp(xl, yl, px) for px in xn]
            return minimum((upper .- yn)[mid]), minimum((yn .- lower)[mid])
        end

        env_upper, env_lower = min_clearance(ShrinkWrap(clearance=0.005))
        @test env_upper > 0.004
        @test env_lower > 0.004

        wide_upper, _ = min_clearance(ShrinkWrap(clearance=0.02))
        @test wide_upper > env_upper
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

@testset "Polar matrix generation (create_2d_polars)" begin
    dat = joinpath(@__DIR__, "data", "test_airfoil.dat")
    alpha_range = deg2rad.(-2:1:2)
    delta_range = deg2rad.(-1:1:1)

    @testset "$name backend" for (name, solver) in (
            ("NeuralFoil", NeuralFoilSolver(model_size="xlarge")),
            ("XFoil", XFoilSolver()))
        work = mktempdir()
        cl_path = joinpath(work, "cl.csv")
        cd_path = joinpath(work, "cd.csv")
        cm_path = joinpath(work, "cm.csv")
        create_2d_polars(; dat_path=dat, cl_polar_path=cl_path, cd_polar_path=cd_path,
            cm_polar_path=cm_path, wind_vel=15.0, area=20.0, width=8.0,
            crease_frac=0.75, alpha_range, delta_range, solver)

        cl, a, d = VortexStepMethod.read_aero_matrix(cl_path)
        @test size(cl) == (length(alpha_range), length(delta_range))
        @test a ≈ collect(alpha_range)
        @test d ≈ collect(delta_range)
        @test all(isfinite, cl)

        cd, _, _ = VortexStepMethod.read_aero_matrix(cd_path)
        @test all(x -> x > 0, cd)

        cm, _, _ = VortexStepMethod.read_aero_matrix(cm_path)
        @test size(cm) == (length(alpha_range), length(delta_range))
    end
end

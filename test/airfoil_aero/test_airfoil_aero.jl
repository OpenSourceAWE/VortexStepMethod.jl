using Test
using VortexStepMethod.AirfoilAero
import VortexStepMethod
using VortexStepMethod.AirfoilAero: KulfanParameters, LeastSquaresFit, ShrinkWrap,
                       shrink_wrap, fit_kulfan_parameters, kulfan_to_coordinates,
                       neuralfoil_aero, class_function, bernstein_basis,
                       leading_edge_basis, normalize_airfoil
using VortexStepMethod: SectionAero, section_surface, read_section_aero
using VortexStepMethod.AirfoilAero: write_section_aero

seg_dist(px, py, ax, ay, bx, by) = begin
    vx, vy = bx - ax, by - ay
    t = clamp(((px - ax) * vx + (py - ay) * vy) / max(vx^2 + vy^2, eps()), 0.0, 1.0)
    hypot(px - (ax + t * vx), py - (ay + t * vy))
end

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
        cloud_to_wrap(xw, yw) = minimum(
            minimum(seg_dist(xn[p], yn[p], xw[k], yw[k], xw[k+1], yw[k+1])
                    for k in 1:length(xw)-1) for p in eachindex(xn))

        xw, yw = shrink_wrap(xr, yr, ShrinkWrap(clearance=0.005))
        @test cloud_to_wrap(xw, yw) > 0.004
        @test (xw[1], yw[1]) == (xw[end], yw[end])
        @test minimum(xw) < -0.003

        xd, yd = shrink_wrap(xr, yr, ShrinkWrap(clearance=0.02))
        @test cloud_to_wrap(xd, yd) > 0.016

        xt, yt = shrink_wrap(xr, yr, ShrinkWrap(clearance=0.0))
        @test cloud_to_wrap(xt, yt) < 0.002
        # A closed loop at zero clearance keeps its trailing edge sharp instead of
        # stamping it into a round cap of `min_clearance` radius.
        @test maximum(xt) ≈ 1.0 atol = 1e-6
    end

    @testset "Shrink-wrap at zero clearance lands on a clean airfoil" begin
        clean = KulfanParameters(fill(0.2, 8), fill(-0.2, 8), 0.0, 0.0)
        xc, yc = kulfan_to_coordinates(clean; n_points=60)
        xn, yn, _ = normalize_airfoil(collect(float.(xc)), collect(float.(yc)))
        xw, yw = shrink_wrap(collect(xc), collect(yc), ShrinkWrap(clearance=0.0))
        @test maximum(minimum(seg_dist(xw[k], yw[k], xn[p], yn[p], xn[p+1], yn[p+1])
                              for p in 1:length(xn)-1) for k in eachindex(xw)) < 1e-9
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

    work = mktempdir()
    cl_path = joinpath(work, "cl.csv")
    cd_path = joinpath(work, "cd.csv")
    cm_path = joinpath(work, "cm.csv")
    create_2d_polars(; dat_path=dat, cl_polar_path=cl_path, cd_polar_path=cd_path,
        cm_polar_path=cm_path, wind_vel=15.0, area=20.0, width=8.0,
        crease_frac=0.75, alpha_range, delta_range,
        solver=NeuralFoilSolver(model_size="xlarge"))

    cl, a, d = VortexStepMethod.read_aero_matrix(cl_path)
    @test size(cl) == (length(alpha_range), length(delta_range))
    @test a ≈ collect(alpha_range)
    @test d ≈ collect(delta_range)
    @test all(isfinite, cl)
    @test maximum(abs.(cl[end, :] .- cl[1, :])) > 0.1
    @test maximum(abs.(cl[:, end] .- cl[:, 1])) > 0.02

    cd, _, _ = VortexStepMethod.read_aero_matrix(cd_path)
    @test all(x -> x > 0, cd)

    cm, _, _ = VortexStepMethod.read_aero_matrix(cm_path)
    @test size(cm) == (length(alpha_range), length(delta_range))
end

@testset "SectionAero table round-trip and interpolation" begin
    alpha_range = deg2rad.([-5.0, 0.0, 5.0, 10.0])
    delta_range = deg2rad.([-3.0, 0.0, 3.0])
    xc = [1.0, 0.5, 0.0, 0.5, 1.0]
    yc = [0.0, 0.06, 0.0, -0.04, 0.0]
    n_node = length(xc)
    x = repeat(xc, 1, length(delta_range))
    y = repeat(yc, 1, length(delta_range))
    val(i, ia, jd) = 100i + 10ia + jd
    cp = [float(val(i, ia, jd)) for i in 1:n_node,
          ia in eachindex(alpha_range), jd in eachindex(delta_range)]
    cf = cp ./ 1000
    aero = SectionAero(alpha_range, delta_range, x, y, cp, cf)

    prefix = joinpath(mktempdir(), "af")
    dat, cp_csv, cf_csv = write_section_aero(prefix, aero)
    back = read_section_aero(dat, cp_csv, cf_csv)
    xb, yb, cpb, cfb = section_surface(back, alpha_range[2], delta_range[1])
    @test isapprox(xb, xc; atol=1e-6)           # contour via .dat (8-decimal)
    @test isapprox(yb, yc; atol=1e-6)
    @test cpb ≈ cp[:, 2, 1]
    @test cfb ≈ cf[:, 2, 1]
    @test read_section_aero("no.dat", "no_cp.csv", "no_cf.csv") === nothing

    _, _, cpa, _ = section_surface(aero, alpha_range[3], delta_range[2])
    @test cpa ≈ cp[:, 3, 2]

    dat_a, cp_arrow, cf_arrow = write_section_aero(prefix, aero; table_format=:arrow)
    @test endswith(cp_arrow, ".arrow") && endswith(cf_arrow, ".arrow")
    arrow = read_section_aero(dat_a, cp_arrow, cf_arrow)
    @test arrow.alpha_range == back.alpha_range
    @test arrow.delta_range == back.delta_range
    @test isequal(arrow.cp, back.cp)
    @test isequal(arrow.cf, back.cf)
    @test_throws ArgumentError write_section_aero(prefix, aero; table_format=:parquet)
end

@testset "generate_section_aero builds a surface table" begin
    truth = KulfanParameters(fill(0.15, 8), fill(-0.15, 8), 0.1, 0.0)
    alpha_range = deg2rad.(-4.0:2.0:4.0)
    delta_range = deg2rad.([0.0, 5.0])
    aero = generate_section_aero(NeuralFoilSolver(model_size="medium"), truth;
        alpha_range, delta_range, reynolds_number=5e5)
    @test aero isa SectionAero
    x, y, cp, cf = section_surface(aero, alpha_range[2], delta_range[1])
    @test length(cp) == size(aero.cp, 1)
    @test all(isfinite, cp) && all(isfinite, cf)
    @test all(cf .>= 0)

    _, _, cp_lo_a, _ = section_surface(aero, alpha_range[1], delta_range[1])
    _, _, cp_hi_a, _ = section_surface(aero, alpha_range[end], delta_range[1])
    @test maximum(abs.(cp_lo_a .- cp_hi_a)) > 0.1

    xdat, ydat = read_dat_coords(joinpath(@__DIR__, "data", "test_airfoil.dat"))
    aero2 = generate_section_aero(NeuralFoilSolver(model_size="medium"), xdat, ydat;
        alpha_range, delta_range=deg2rad.([0.0, 5.0]), reynolds_number=5e5)
    @test aero2 isa SectionAero
    @test size(aero2.cp, 1) > 0
end

@testset "NeuralFoil physical invariants" begin
    sym = KulfanParameters(fill(0.2, 8), fill(-0.2, 8), 0.0, 0.0)

    sweep = neuralfoil_aero(sym, [-4.0, 0.0, 4.0], 5e5; model_size="xlarge")
    @test abs(sweep.CL[2]) < 0.02
    @test sweep.CL[1] < sweep.CL[2] < sweep.CL[3]
    @test sweep.CL[1] ≈ -sweep.CL[3] atol = 0.02

    xs, ys = kulfan_to_coordinates(sym; n_points=120)
    xs, ys = collect(xs), collect(ys)
    # geometric camber line; the weight sum is not a symmetry proxy once
    # shrink_wrap resamples the two surfaces at differing x-stations
    function max_camber(k)
        cx, cy = kulfan_to_coordinates(k; n_points=120)
        cy = collect(cy)
        return maximum(abs, reverse(cy[1:120]) .+ cy[120:end]) / 2
    end
    base = deform_section(xs, ys, 0.0)
    flap = deform_section(xs, ys, deg2rad(10.0); crease_frac=0.75)
    @test max_camber(base.kulfan) < 5e-3
    @test max_camber(flap.kulfan) > 0.02

    @test minimum(flap.x) ≈ 0 atol = 0.01
    @test maximum(flap.x) ≈ 1 atol = 0.01
    @test abs(flap.y[argmax(flap.x)]) < 0.01
    @test abs(flap.y[argmin(flap.x)]) < 0.01

    sec = neuralfoil_section(sym, [6.0], 5e5; model_size="xlarge")
    load = sec.cp_lower[:, 1] .- sec.cp_upper[:, 1]
    @test sum(load) > 0
    @test count(>(0), load) > length(load) ÷ 2
end

@testset "deform_section node count is delta-independent" begin
    # A base whose paneling differs from shrink_wrap's fixed output (2·n_points-1)
    # must still come back at that same length for delta == 0, or a delta sweep
    # mixes lengths and generate_airfoil_aero NaN-s the delta == 0 column.
    sym = KulfanParameters(fill(0.15, 8), fill(-0.15, 8), 0.0, 0.0)
    xb, yb = kulfan_to_coordinates(sym)   # default paneling, != shrink_wrap output
    base = deform_section(xb, yb, 0.0; crease_frac=0.8)
    @test length(base.x) != length(xb)    # delta == 0 is re-wrapped, not passed through
    for d in (deg2rad(-5.0), deg2rad(5.0), deg2rad(15.0))
        def = deform_section(xb, yb, d; crease_frac=0.8)
        @test length(def.x) == length(base.x)
        @test length(def.y) == length(base.y)
    end
end

@testset "generate_polar_from_coordinates POLAR_VECTORS sweep" begin
    x, y = read_dat_coords(joinpath(@__DIR__, "data", "test_airfoil.dat"))
    csv = joinpath(mktempdir(), "polar.csv")
    sols = generate_polar_from_coordinates(x, y, csv;
        Re=5e5, alpha_range=-4:2:4, solver=NeuralFoilSolver(model_size="medium"))
    @test sols isa AbstractVector
    @test isfile(csv)
    header = lowercase(readline(csv))
    @test occursin("alpha", header)
    @test !occursin("delta", header)
end

@testset "turn_trailing_edge! legacy crease cleanup" begin
    x, y = read_dat_coords(joinpath(@__DIR__, "data", "test_airfoil.dat"))
    crease_frac = 0.7
    for angle in (deg2rad(10.0), deg2rad(-10.0))
        xd, yd = collect(float.(x)), collect(float.(y))
        lower, upper = get_lower_upper(xd, yd, crease_frac)
        @test lower < upper
        n0 = length(xd)
        turn_trailing_edge!(angle, xd, yd, lower, upper, crease_frac)
        @test length(xd) == length(yd)
        @test length(xd) <= n0
        @test all(isfinite, xd) && all(isfinite, yd)
    end
end

@testset "load_neuralfoil_model missing weights errors" begin
    @test_throws ErrorException load_neuralfoil_model("nonexistent_size")

    weights_dir = joinpath(dirname(pathof(VortexStepMethod)), "airfoil_aero", "data")
    partial = mktempdir()
    cp(joinpath(weights_dir, "nn-medium.npz"), joinpath(partial, "nn-medium.npz"))
    @test_throws ErrorException load_neuralfoil_model("medium"; weights_dir=partial)
end

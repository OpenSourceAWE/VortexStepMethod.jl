using VortexStepMethod
using VortexStepMethod.AirfoilAero: read_dat_coordinates, fit_kulfan_parameters,
                                    LeastSquaresFit, DeformedSection, analyze_section,
                                    analyze_sweep, deform_section, neuralfoil_aero,
                                    generate_aero_matrices, NeuralFoilSolver, XFoilSolver
using Test

@testset "NeuralFoil vs XFoil agree on a flat wing (matched settings)" begin
    # NACA 0012 (squarely inside NeuralFoil's training distribution). NeuralFoil is
    # trained on XFoil, so with matched settings — here both tripped at 5% chord — the
    # two backends agree tightly, including across the trailing-edge deflection now that
    # the shrink wrap emits XFoil-friendly panels (no near-zero panels at the hinge, so
    # XFoil needs no internal repaneling). This is the suite's only XFoil path;
    # everything else uses NeuralFoil.
    dat = joinpath(@__DIR__, "data", "naca0012.dat")
    x, y = read_dat_coordinates(dat)

    Re = 1e6
    xtr = 0.05
    nf_solver = NeuralFoilSolver(model_size="xlarge", xtr_upper=xtr, xtr_lower=xtr)
    xf_solver = XFoilSolver(xtrip=(xtr, xtr))

    alpha_range = deg2rad.(-4:2:6)
    delta_range = deg2rad.(0:1:1)
    nf = generate_aero_matrices(nf_solver, x, y; alpha_range, delta_range, Re)
    xf = generate_aero_matrices(xf_solver, x, y; alpha_range, delta_range, Re)

    @test maximum(abs, nf[1][:, 1] .- xf[1][:, 1]) < 0.01   # Cl, base airfoil
    @test maximum(abs, nf[1] .- xf[1]) < 0.03               # Cl, incl. 1-deg deflection
    @test maximum(abs, nf[2] .- xf[2]) < 0.01               # Cd

    function flat_wing((cl, cd, cm))
        span, chord = 20.0, 1.0
        wing = Wing(4, spanwise_distribution=LINEAR)
        data = (collect(alpha_range), collect(delta_range), cl, cd, cm)
        add_section!(wing, [0.0, span/2, 0.0], [chord, span/2, 0.0], POLAR_MATRICES, data)
        add_section!(wing, [0.0, -span/2, 0.0], [chord, -span/2, 0.0], POLAR_MATRICES, data)
        refine!(wing)
        return wing
    end
    function lift(matrices)
        body = BodyAerodynamics([flat_wing(matrices)])
        aoa = deg2rad(4.0)
        set_va!(body, [cos(aoa), 0.0, sin(aoa)] * 15.0)
        return collect(solve!(Solver(body; aerodynamic_model_type=VSM), body).force_coeffs)
    end
    fc_nf, fc_xf = lift(nf), lift(xf)
    @test fc_nf[3] > 0 && fc_xf[3] > 0          # positive lift at positive alpha
    @test abs(fc_nf[3] - fc_xf[3]) < 0.03       # same lift from either backend's polars

    @testset "trip point and Reynolds move both backends and keep them matched" begin
        kf = fit_kulfan_parameters(collect(x), collect(y), LeastSquaresFit())
        base = DeformedSection(kf, collect(x), collect(y))
        nfcd(xt, re) = neuralfoil_aero(kf, [4.0], re; model_size="xlarge",
                                       xtr_upper=xt, xtr_lower=xt).CD[1]
        xfcd(xt, re) = analyze_section(XFoilSolver(xtrip=(xt, xt)), base, deg2rad(4.0), re).cd

        # forcing transition at 5% raises drag versus free transition, in both solvers,
        # and they stay matched at either transition setting
        @test nfcd(0.05, Re) - nfcd(1.0, Re) > 0.0015
        @test xfcd(0.05, Re) - xfcd(1.0, Re) > 0.0015
        @test isapprox(nfcd(0.05, Re), xfcd(0.05, Re); atol=0.002)
        @test isapprox(nfcd(1.0, Re), xfcd(1.0, Re); atol=0.002)

        # raising the Reynolds number lowers drag in both solvers, still matched
        @test nfcd(0.05, 1e6) - nfcd(0.05, 5e6) > 0.001
        @test xfcd(0.05, 1e6) - xfcd(0.05, 5e6) > 0.001
        @test isapprox(nfcd(0.05, 5e6), xfcd(0.05, 5e6); atol=0.002)
    end

    @testset "deflected section solves with and without XFoil repaneling" begin
        # The shrink-wrap smoothing keeps the hinge crease from collapsing panels to
        # near-zero length, so XFoil converges on the deflected shape either way.
        deflected = deform_section(x, y, deg2rad(3.0); crease_frac=0.75)
        alphas = deg2rad.(-2:2:4)
        xf_no = [s.cl for s in
                 analyze_sweep(XFoilSolver(xtrip=(xtr, xtr), repanel=false), deflected, alphas, Re)]
        xf_yes = [s.cl for s in
                  analyze_sweep(XFoilSolver(xtrip=(xtr, xtr), repanel=true), deflected, alphas, Re)]
        nf_cl = neuralfoil_aero(deflected.kulfan, rad2deg.(alphas), Re;
                                model_size="xlarge", xtr_upper=xtr, xtr_lower=xtr).CL
        @test all(isfinite, xf_no)
        @test all(isfinite, xf_yes)
        # the default (no repaneling) tracks NeuralFoil; XFoil's own repaneling
        # re-clusters the hinge and drifts, so it is not held to the same bound
        @test maximum(abs, nf_cl .- xf_no) < 0.03
    end
end

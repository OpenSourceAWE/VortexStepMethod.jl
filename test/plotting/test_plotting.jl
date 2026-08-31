module FakeMakieNoCurrentBackend end

module FakeMakieCurrentBackendThrows
    current_backend() = error("boom")
end

module FakeMakieReturnsCairoModule
    import CairoMakie
    current_backend() = CairoMakie
end

module FakeMakieReturnsOtherModule
    import Base
    current_backend() = Base
end
using CairoMakie
using MakieControlPlots
# MakieControlPlots activates GLMakie on load; force CairoMakie so headless
# figure saving uses the software backend.
CairoMakie.activate!()

using VortexStepMethod
using Test

const makie_ext = Base.get_extension(VortexStepMethod, :VortexStepMethodMakieExt)

# Helper to robustly delete files on platforms with occasional file locks
safe_rm(path) = begin
    if isfile(path)
        try
            rm(path; force=true)
        catch
            sleep(0.2)
            try
                rm(path; force=true)
            catch
                # last resort, ignore
            end
        end
    end
    nothing
end

global ram_wing = ram_air_matrix_wing(; n_panels=20, n_sections=4,
                          alpha_range=deg2rad.(-1:1.0:1),
                          delta_range=deg2rad.(-1:1.0:1))

function create_body_aero()
    n_panels = 20          # Number of panels
    span = 20.0            # Wing span [m]
    chord = 1.0            # Chord length [m]
    v_a = 20.0             # Magnitude of inflow velocity [m/s]
    alpha_deg = 30.0       # Angle of attack [degrees]
    alpha = deg2rad(alpha_deg)

    wing = Wing(n_panels, spanwise_distribution=LINEAR)

    add_section!(wing,
        [0.0, span/2, 0.0],
        [chord, span/2, 0.0],
        INVISCID)
    add_section!(wing,
        [0.0, -span/2, 0.0],
        [chord, -span/2, 0.0],
        INVISCID)

    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
    set_va!(body_aero, vel_app)
    body_aero
end

@testset "Plotting (Makie)" begin
    save_dir = tempdir()
    body_aero = create_body_aero()

    fig = plot_geometry(
        body_aero,
        "Rectangular_wing_geometry";
        data_type=".png",
        save_path=save_dir,
        is_save=true,
        is_show=false)
    @test fig isa Figure

    @test hasmethod(VortexStepMethod.show_plot, Tuple{Figure})
    @test_throws MethodError VortexStepMethod.show_plot(nothing)
    @test_nowarn VortexStepMethod.show_plot(fig)

    @test isfile(joinpath(save_dir,
                          "Rectangular_wing_geometry_angled_view.png"))
    safe_rm(joinpath(save_dir,
                     "Rectangular_wing_geometry_angled_view.png"))
    @test isfile(joinpath(save_dir,
                          "Rectangular_wing_geometry_front_view.png"))
    safe_rm(joinpath(save_dir,
                     "Rectangular_wing_geometry_front_view.png"))
    @test isfile(joinpath(save_dir,
                          "Rectangular_wing_geometry_side_view.png"))
    safe_rm(joinpath(save_dir,
                     "Rectangular_wing_geometry_side_view.png"))
    @test isfile(joinpath(save_dir,
                          "Rectangular_wing_geometry_top_view.png"))
    safe_rm(joinpath(save_dir,
                     "Rectangular_wing_geometry_top_view.png"))

    # Initialize the solvers
    vsm_solver = Solver(body_aero; aerodynamic_model_type=VSM)
    llt_solver = Solver(body_aero; aerodynamic_model_type=LLT)

    # Solve the VSM and LLT
    results_vsm = solve(vsm_solver, body_aero)
    results_llt = solve(llt_solver, body_aero)

    # Plot spanwise distributions
    y_coordinates = [panel.aero_center[2]
                     for panel in body_aero.panels]

    fig = plot_distribution(
        [y_coordinates, y_coordinates],
        [results_vsm, results_llt],
        ["VSM", "LLT"],
        title="Spanwise Distributions",
        is_show=false
    )
    @test fig isa Figure

    # Plot polar curves
    v_a = 20.0
    angle_range = range(0, 20, 20)
    fig = plot_polars(
        [llt_solver, vsm_solver],
        [body_aero, body_aero],
        ["VSM", "LLT"],
        angle_range=angle_range,
        angle_type="angle_of_attack",
        v_a=v_a,
        title="Rectangular Wing Polars",
        data_type=".png",
        save_path=save_dir,
        is_save=true,
        is_show=false
    )
    @test fig isa Figure
    @test isfile(joinpath(save_dir, "Rectangular_Wing_Polars.png"))
    safe_rm(joinpath(save_dir, "Rectangular_Wing_Polars.png"))

    # Plot polars with CL vs CD (cl_over_cd=false)
    fig = plot_polars(
        [llt_solver, vsm_solver],
        [body_aero, body_aero],
        ["VSM", "LLT"],
        angle_range=angle_range,
        angle_type="angle_of_attack",
        v_a=v_a,
        title="Polars CL vs CD",
        is_save=false,
        is_show=false,
        cl_over_cd=false
    )
    @test fig isa Figure

    # Plot combined analysis with cl_over_cd
    fig = plot_combined_analysis(
        vsm_solver, body_aero, results_vsm;
        angle_range=angle_range,
        angle_type="angle_of_attack",
        angle_of_attack=30.0,
        v_a=v_a,
        title="Combined Analysis",
        is_save=false,
        is_show=false,
        cl_over_cd=true
    )
    @test fig isa Figure

    # Plot combined analysis with cl_over_cd=false
    fig = plot_combined_analysis(
        vsm_solver, body_aero, results_vsm;
        angle_range=angle_range,
        angle_type="angle_of_attack",
        angle_of_attack=30.0,
        v_a=v_a,
        title="Combined CL vs CD",
        is_save=false,
        is_show=false,
        cl_over_cd=false
    )
    @test fig isa Figure

    # Test polar data plotting
    body_aero = BodyAerodynamics([ram_wing])
    fig = plot_polar_data(body_aero; is_show=false)
    @test fig isa Figure

    fig_rect = plot_polar_data(body_aero;
        alphas=collect(deg2rad.(-5:1.0:15)),
        delta_tes=collect(deg2rad.(-3:1.0:5)),
        is_show=false)
    @test fig_rect isa Figure

    # Edge cases: empty panels and distributed apparent wind
    body_aero_empty = create_body_aero()
    empty!(body_aero_empty.panels)
    @test_throws Exception plot_geometry(
        body_aero_empty,
        "Rectangular_wing_geometry_empty_panels";
        is_save=false,
        is_show=false,
    )

    body_aero_distributed = create_body_aero()
    n_panels = length(body_aero_distributed.panels)
    va_distribution = repeat([12.0 0.0 1.0], n_panels, 1)
    set_va!(body_aero_distributed, va_distribution)

    @test body_aero_distributed.has_distributed_va
    fig = plot_geometry(
        body_aero_distributed,
        "Rectangular_wing_geometry_distributed_va";
        is_save=false,
        is_show=false,
    )
    @test fig !== nothing

    literature_csv = joinpath(tempdir(), "polar_literature_aoa.csv")
    open(literature_csv, "w") do io
        write(io, "AOA,cl,cd,cs\n")
        write(io, "0.0,0.1,0.01,0.0\n")
        write(io, "5.0,0.5,0.02,0.01\n")
        write(io, "10.0,0.9,0.04,0.02\n")
    end

    try
        fig = plot_polars(
            Solver[],
            BodyAerodynamics[],
            ["Literature"],
            literature_path_list=[literature_csv],
            title="Literature AOA Header",
            is_save=false,
            is_show=false,
        )
        @test fig !== nothing
    finally
        safe_rm(literature_csv)
    end

    # Unit tests for shared extract_literature_polar_data
    using DelimitedFiles

    # Tuple parsing branch (e.g. readdlm(...; header=true) shape)
    tuple_table = [0.0 0.10 0.010; 5.0 0.20 0.020]
    tuple_header = [" AoA " "CL" "CD"]
    tuple_result = VortexStepMethod.extract_literature_polar_data(
        (tuple_table, tuple_header), "tuple.csv")
    @test tuple_result.polar_data[1] == tuple_table[:, 1]
    @test tuple_result.polar_data[2] == tuple_table[:, 2]
    @test tuple_result.polar_data[3] == tuple_table[:, 3]
    @test tuple_result.polar_data[4] == zeros(size(tuple_table, 1))
    @test all(isnan, tuple_result.cmx)
    @test all(isnan, tuple_result.cmy)
    @test all(isnan, tuple_result.cmz)

    # Matrix parsing branch (header in first row + explicit CS column)
    matrix_data = Any[
        "alpha" "cl" "cd" "cs";
        0.0 0.11 0.011 0.001;
        4.0 0.21 0.021 0.002
    ]
    matrix_result = VortexStepMethod.extract_literature_polar_data(
        matrix_data, "matrix.csv")
    @test Float64.(matrix_result.polar_data[1]) == [0.0, 4.0]
    @test Float64.(matrix_result.polar_data[2]) == [0.11, 0.21]
    @test Float64.(matrix_result.polar_data[3]) == [0.011, 0.021]
    @test Float64.(matrix_result.polar_data[4]) == [0.001, 0.002]

    # Missing required columns should throw a clear ArgumentError
    bad_data = Any[
        "aoa" "cl" "cs";
        0.0 0.1 0.0
    ]
    @test_throws ArgumentError VortexStepMethod.extract_literature_polar_data(
        bad_data, "bad.csv")

    # CM coefficient extraction from literature data
    cm_csv = tempname() * "_lit_cm.csv"
    open(cm_csv, "w") do io_cm
        write(io_cm,
            "alpha,cl,cd,cs,cmx,cmy,cmz\n" *
            "0.0,0.1,0.01,0.0,0.001,0.002,0.003\n" *
            "5.0,0.5,0.02,0.01,0.004,0.005,0.006\n")
    end
    cm_result = VortexStepMethod.extract_literature_polar_data(
        readdlm(cm_csv, ','), cm_csv)
    @test cm_result.polar_data[1] == [0.0, 5.0]
    @test Float64.(cm_result.cmx) == [0.001, 0.004]
    @test Float64.(cm_result.cmy) == [0.002, 0.005]
    @test Float64.(cm_result.cmz) == [0.003, 0.006]
    safe_rm(cm_csv)

    # angle_type="side_slip" literature loading
    beta_csv = tempname() * "_lit_beta.csv"
    open(beta_csv, "w") do io_beta
        write(io_beta,
            "alpha,beta,cl,cd,cs\n" *
            "7.4,0.0,0.7,0.06,0.0\n" *
            "7.4,5.0,0.68,0.07,0.01\n")
    end
    beta_result = VortexStepMethod.extract_literature_polar_data(
        readdlm(beta_csv, ','), beta_csv;
        angle_type="side_slip")
    @test beta_result.polar_data[1] == [0.0, 5.0]
    safe_rm(beta_csv)

    # Integration: literature CSV with AoA alias and no CS
    lit_no_cs_path = tempname() * "_lit_no_cs.csv"
    open(lit_no_cs_path, "w") do io_no_cs
        write(io_no_cs, "aoa,cl,cd\n0.0,0.10,0.010\n5.0,0.20,0.020\n")
    end
    fig_lit_no_cs = plot_polars(
        Any[],
        Any[],
        ["Literature no CS"];
        literature_path_list=[lit_no_cs_path],
        is_save=false,
        is_show=false
    )
    @test fig_lit_no_cs !== nothing
    safe_rm(lit_no_cs_path)

    # Integration: missing CD column should fail
    lit_bad_path = tempname() * "_lit_bad.csv"
    open(lit_bad_path, "w") do io_bad
        write(io_bad, "alpha,cl\n0.0,0.10\n5.0,0.20\n")
    end
    @test_throws ArgumentError plot_polars(
        Any[],
        Any[],
        ["Literature bad"];
        literature_path_list=[lit_bad_path],
        is_save=false,
        is_show=false
    )
    safe_rm(lit_bad_path)

    # Test show_moments=true with literature data
    cm_lit_path = tempname() * "_lit_moments.csv"
    open(cm_lit_path, "w") do io_cm_lit
        write(io_cm_lit,
            "alpha,cl,cd,cs,cmx,cmy,cmz\n" *
            "0.0,0.1,0.01,0.0,0.001,0.002,0.003\n" *
            "5.0,0.5,0.02,0.01,0.004,0.005,0.006\n")
    end
    fig_moments = plot_polars(
        Any[],
        Any[],
        ["Literature with moments"];
        literature_path_list=[cm_lit_path],
        show_moments=true,
        is_save=false,
        is_show=false
    )
    @test fig_moments !== nothing
    safe_rm(cm_lit_path)

    # Test show_moments=false (default)
    no_cm_path = tempname() * "_lit_no_cm.csv"
    open(no_cm_path, "w") do io_no_cm
        write(io_no_cm,
            "alpha,cl,cd\n" *
            "0.0,0.1,0.01\n5.0,0.5,0.02\n")
    end
    fig_no_moments = plot_polars(
        Any[],
        Any[],
        ["Literature no moments"];
        literature_path_list=[no_cm_path],
        show_moments=false,
        is_save=false,
        is_show=false
    )
    @test fig_no_moments !== nothing
    safe_rm(no_cm_path)

    # Tests for save_plot function
    @testset "_active_backend_prefers_vector_output" begin
        @test makie_ext !== nothing

        active_backend_prefers_vector_output =
            getfield(makie_ext, :_active_backend_prefers_vector_output)


        @test active_backend_prefers_vector_output(FakeMakieNoCurrentBackend) == false
        @test active_backend_prefers_vector_output(FakeMakieCurrentBackendThrows) == false

        @test active_backend_prefers_vector_output(FakeMakieReturnsCairoModule) == true
        @test active_backend_prefers_vector_output(FakeMakieReturnsOtherModule) == false
    end

    body_aero = create_body_aero()
    fig = plot_geometry(
        body_aero,
        "save_plot_test";
        is_save=false,
        is_show=false)
    @test fig isa Figure

    active_backend_prefers_vector_output =
        getfield(makie_ext, :_active_backend_prefers_vector_output)

    save_test_dir = tempdir()
    
    # Test 1: save_plot with explicit data_type (".png")
    VortexStepMethod.save_plot(fig, save_test_dir, "test_explicit_png", data_type=".png")
    @test isfile(joinpath(save_test_dir, "test_explicit_png.png"))
    safe_rm(joinpath(save_test_dir, "test_explicit_png.png"))

    # Test 2: save_plot with explicit data_type (".pdf")
    VortexStepMethod.save_plot(fig, save_test_dir, "test_explicit_pdf", data_type=".pdf")
    @test isfile(joinpath(save_test_dir, "test_explicit_pdf.pdf"))
    safe_rm(joinpath(save_test_dir, "test_explicit_pdf.pdf"))

    # Test 3: save_plot with data_type=nothing (backend-aware detection)
    backend_aware_dir = mktempdir()
    try
        VortexStepMethod.save_plot(fig, backend_aware_dir, "test_backend_aware", data_type=nothing)
        pdf_path = joinpath(backend_aware_dir, "test_backend_aware.pdf")
        png_path = joinpath(backend_aware_dir, "test_backend_aware.png")
        expected_ext = active_backend_prefers_vector_output(Makie) ? ".pdf" : ".png"

        @test xor(isfile(pdf_path), isfile(png_path))
        @test isfile(joinpath(backend_aware_dir, "test_backend_aware" * expected_ext))
    finally
        safe_rm(joinpath(backend_aware_dir, "test_backend_aware.pdf"))
        safe_rm(joinpath(backend_aware_dir, "test_backend_aware.png"))
        rm(backend_aware_dir; force=true, recursive=true)

        # Test 4: save_plot with title containing spaces (should be sanitized to underscores)
        VortexStepMethod.save_plot(fig, save_test_dir, "test with spaces", data_type=".png")
        @test isfile(joinpath(save_test_dir, "test_with_spaces.png"))
        safe_rm(joinpath(save_test_dir, "test_with_spaces.png"))

        # Test 5: save_plot with title containing percent signs (should be sanitized to "pct")
        VortexStepMethod.save_plot(fig, save_test_dir, "test%efficiency", data_type=".png")
        @test isfile(joinpath(save_test_dir, "testpctefficiency.png"))
        safe_rm(joinpath(save_test_dir, "testpctefficiency.png"))

        # Test 6: save_plot with title containing both spaces and percent signs
        VortexStepMethod.save_plot(fig, save_test_dir, "test %efficiency metric", data_type=".png")
        @test isfile(joinpath(save_test_dir, "test_pctefficiency_metric.png"))
        safe_rm(joinpath(save_test_dir, "test_pctefficiency_metric.png"))

        # Test 7: save_plot creates directory if it doesn't exist
        nested_dir = joinpath(save_test_dir, "nested_save_plot_dir")
        !isdir(nested_dir) && @test !isdir(nested_dir)
        VortexStepMethod.save_plot(fig, nested_dir, "test_nested_dir", data_type=".png")
        @test isdir(nested_dir)
        @test isfile(joinpath(nested_dir, "test_nested_dir.png"))
        safe_rm(joinpath(nested_dir, "test_nested_dir.png"))
        rm(nested_dir; force=true)

        # Test 8: save_plot raises error when save_path is nothing
        @test_throws ArgumentError VortexStepMethod.save_plot(fig, nothing, "test_title", data_type=".png")
    end
end

"""
Build a body whose refined sections carry a [`SectionAero`](@ref) surface table, so
the lofted airfoil-skin plotting path (`plot!(...; airfoils=true)`) has contours.
"""
function create_body_aero_with_skin(; n_panels=4)
    xs = [1.0, 0.75, 0.25, 0.0, 0.25, 0.75, 1.0]
    ys = [0.0, 0.05, 0.05, 0.0, -0.05, -0.05, 0.0]
    n_node = length(xs)
    alpha_range = deg2rad.([-2.0, 0.0, 2.0])
    x = reshape(xs, n_node, 1)
    y = reshape(ys, n_node, 1)
    cp = zeros(n_node, length(alpha_range), 1)
    cf = fill(0.003, n_node, length(alpha_range), 1)
    section_aero = SectionAero(alpha_range, [0.0], x, y, cp, cf)

    # POLAR_VECTORS, so the panels carry the interpolations a live polar rewrites.
    polar = (alpha_range, [0.0, 0.5, 1.0], fill(0.02, 3), fill(-0.05, 3))
    wing = Wing(n_panels, spanwise_distribution=LINEAR)
    add_section!(wing, [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], POLAR_VECTORS, polar, section_aero)
    add_section!(wing, [0.0, -2.0, 0.0], [1.0, -2.0, 0.0], POLAR_VECTORS, polar,
                 section_aero)
    refine!(wing)
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, [20.0, 0.0, 1.0])
    return body_aero, n_node
end

@testset "Airfoil skin (Makie)" begin
    airfoil_skin_geometry = getfield(makie_ext, :airfoil_skin_geometry)
    skin_observables = getfield(makie_ext, :AIRFOIL_SKIN_OBSERVABLES)
    skin_observables[] = nothing

    body_aero, n_node = create_body_aero_with_skin(; n_panels=4)
    n_sections = length(body_aero.wings[1].refined_sections)

    vertices, faces, ribs = airfoil_skin_geometry(body_aero)
    @test length(ribs) == n_sections
    @test all(rib -> length(rib) == n_node, ribs)
    @test !isempty(vertices)
    # Each of the (n_sections - 1) lofted strips triangulates (n_node - 1) quads.
    @test length(faces) == 2 * (n_sections - 1) * (n_node - 1)

    # Static airfoil-skin plot.
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plots = Makie.plot!(ax, body_aero; airfoils=true)
    @test !isempty(plots)

    # Observable airfoil-skin plot registers the body for pose updates.
    fig_obs = Figure()
    ax_obs = Axis3(fig_obs[1, 1])
    Makie.plot!(ax_obs, body_aero; airfoils=true, use_observables=true)
    @test !isnothing(skin_observables[])
    @test haskey(skin_observables[], objectid(body_aero))

    # In-place update refreshes the registered skin observables without error.
    @test_nowarn Makie.plot!(body_aero)

    # Live polars: the skin has to draw the panel's stored deformed shape, not the
    # tabulated contour, or a deformation bug is invisible in the picture.
    basis = VortexStepMethod.AirfoilAero.KulfanBasis()
    base = VortexStepMethod.KulfanParameters(fill(0.15, 8), fill(-0.05, 8), 0.0, 0.0)
    cambered = VortexStepMethod.AirfoilAero.deform_kulfan(basis, base,
        @. 0.08 * basis.x * (1 - basis.x))
    for panel in body_aero.panels
        VortexStepMethod.set_polar!(panel, deg2rad.([-4.0, 0.0, 4.0]),
            [0.3, 0.5, 0.7], [0.02, 0.02, 0.03], [-0.05, -0.05, -0.05];
            shape=cambered)
    end
    @test body_aero.panels[1].live_shape === cambered
    v_live, f_live, ribs_live = airfoil_skin_geometry(body_aero)
    live_nodes = length(VortexStepMethod.AirfoilAero.kulfan_to_coordinates(cambered)[1])
    @test live_nodes != n_node                       # the two routes are distinguishable
    @test all(rib -> length(rib) == live_nodes, ribs_live)
    @test length(f_live) == 2 * (n_sections - 1) * (live_nodes - 1)
    @test_nowarn Makie.plot!(Axis3(Figure()[1, 1]), body_aero; airfoils=true)
    for panel in body_aero.panels
        panel.live_shape = nothing
    end
    @test all(rib -> length(rib) == n_node, airfoil_skin_geometry(body_aero)[3])

    # A body without any surface aero yields no skin geometry, but still plots.
    plain_body = create_body_aero()
    v_plain, f_plain, ribs_plain = airfoil_skin_geometry(plain_body)
    @test isempty(v_plain)
    @test isempty(f_plain)
    @test isempty(ribs_plain)

    # Update path is a no-op (no error) when the body was never registered.
    skin_observables[] = nothing
    @test_nowarn Makie.plot!(plain_body)

    # border_linewidth flows through the standard (non-airfoil) panel plot.
    fig_lw = Figure()
    ax_lw = Axis3(fig_lw[1, 1])
    @test_nowarn Makie.plot!(ax_lw, plain_body; border_linewidth=3.0)
end
nothing

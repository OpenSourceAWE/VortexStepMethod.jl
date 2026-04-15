backend = if "plot-controlplots" in ARGS
    using ControlPlots
    import ControlPlots: plt
    "ControlPlots"
else
    using CairoMakie
    "Makie"
end

using VortexStepMethod
using Test

# Resolve repo data directory for ram air kite assets
_ram_data_dir = joinpath(dirname(dirname(@__DIR__)),
                         "data", "ram_air_kite")

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

let
    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
    body_src = joinpath(_ram_data_dir, "ram_air_kite_body.obj")
    foil_src = joinpath(_ram_data_dir, "ram_air_kite_foil.dat")
    cp(body_src, body_path; force=true)
    cp(foil_src, foil_path; force=true)
    global ram_wing = ObjWing(body_path, foil_path;
                              alpha_range=deg2rad.(-1:1),
                              delta_range=deg2rad.(-1:1))
end

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

@testset "Plotting ($backend)" begin
    save_dir = tempdir()
    body_aero = create_body_aero()

    fig = plot_geometry(
        body_aero,
        "Rectangular_wing_geometry";
        data_type=".png",
        save_path=save_dir,
        is_save=true,
        is_show=false)
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end

    if backend == "Makie"
        @test hasmethod(VortexStepMethod.show_plot, Tuple{Figure})
        @test_throws MethodError VortexStepMethod.show_plot(nothing)
        @test_nowarn VortexStepMethod.show_plot(fig)
    else
        FigType = plt.Figure
        @test hasmethod(VortexStepMethod.show_plot, Tuple{FigType})
        @test_throws MethodError VortexStepMethod.show_plot(nothing)
    end

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
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end

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
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end
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
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end

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
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end

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
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end

    # Test polar data plotting
    body_aero = BodyAerodynamics([ram_wing])
    fig = plot_polar_data(body_aero; is_show=false)
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end

    if backend == "ControlPlots"
        body_aero_empty = create_body_aero()
        empty!(body_aero_empty.panels)
        @test_throws ArgumentError plot_geometry(
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
        fig_dpi = plt.figure()
        default_dpi = fig_dpi.get_dpi()
        @test default_dpi != 173

        show_plot(fig_dpi; dpi=173)
        @test fig_dpi.get_dpi() == 173

        # Also verify the default keyword value is applied.
        show_plot(fig_dpi)
        @test fig_dpi.get_dpi() == 130
        plt.close(fig_dpi)

        ext = Base.get_extension(VortexStepMethod, :VortexStepMethodControlPlotsExt)
        @test ext !== nothing

        # Unit-test tuple parsing branch (e.g. readdlm(...; header=true) shape).
        tuple_table = [0.0 0.10 0.010; 5.0 0.20 0.020]
        tuple_header = [" AoA " "CL" "CD"]
        tuple_pd = ext._extract_literature_polar_data((tuple_table, tuple_header), "tuple.csv")
        @test tuple_pd[1] == tuple_table[:, 1]
        @test tuple_pd[2] == tuple_table[:, 2]
        @test tuple_pd[3] == tuple_table[:, 3]
        @test tuple_pd[4] == zeros(size(tuple_table, 1))

        # Unit-test matrix parsing branch (header in first row + explicit CS column).
        matrix_data = Any[
            "alpha" "cl" "cd" "cs";
            0.0 0.11 0.011 0.001;
            4.0 0.21 0.021 0.002
        ]
        matrix_pd = ext._extract_literature_polar_data(matrix_data, "matrix.csv")
        @test Float64.(matrix_pd[1]) == [0.0, 4.0]
        @test Float64.(matrix_pd[2]) == [0.11, 0.21]
        @test Float64.(matrix_pd[3]) == [0.011, 0.021]
        @test Float64.(matrix_pd[4]) == [0.001, 0.002]

        # Missing required columns should throw a clear ArgumentError.
        bad_data = Any[
            "aoa" "cl" "cs";
            0.0 0.1 0.0
        ]
        @test_throws ArgumentError ext._extract_literature_polar_data(bad_data, "bad.csv")

        # Integration: literature CSV with AoA alias and no CS should still plot.
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

        # Integration: missing CD column should fail.
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
    end
end
nothing

backend = if "plot-controlplots" in ARGS
    using ControlPlots
    "ControlPlots"
else
    using CairoMakie
    "Makie"
end

using VortexStepMethod
using Test

# Resolve repo data directory for ram air kite assets
const _ram_data_dir = joinpath(dirname(dirname(@__DIR__)),
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

if !@isdefined ram_wing
    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
    body_src = joinpath(_ram_data_dir, "ram_air_kite_body.obj")
    foil_src = joinpath(_ram_data_dir, "ram_air_kite_foil.dat")
    cp(body_src, body_path; force=true)
    cp(foil_src, foil_path; force=true)
    ram_wing = ObjWing(body_path, foil_path;
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

    # Test polar data plotting
    body_aero = BodyAerodynamics([ram_wing])
    fig = plot_polar_data(body_aero; is_show=false)
    if backend == "Makie"
        @test fig isa Figure
    else
        @test fig !== nothing
    end
end
nothing

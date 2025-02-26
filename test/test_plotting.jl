using VortexStepMethod
using ControlPlots
using Test

function create_wa()
    # Step 1: Define wing parameters
    n_panels = 20          # Number of panels
    span = 20.0            # Wing span [m]
    chord = 1.0            # Chord length [m]
    v_a = 20.0            # Magnitude of inflow velocity [m/s]
    density = 1.225        # Air density [kg/m³]
    alpha_deg = 30.0       # Angle of attack [degrees]
    alpha = deg2rad(alpha_deg)

    # Step 2: Create wing geometry with linear panel distribution
    wing = Wing(n_panels, spanwise_panel_distribution=:linear)

    # Add wing sections - defining only tip sections with inviscid airfoil model
    add_section!(wing, 
        [0.0, span/2, 0.0],   # Left tip LE 
        [chord, span/2, 0.0],  # Left tip TE
        :inviscid)
    add_section!(wing, 
        [0.0, -span/2, 0.0],  # Right tip LE
        [chord, -span/2, 0.0], # Right tip TE
        :inviscid)

    # Step 3: Initialize aerodynamics
    wa = BodyAerodynamics([wing])
    # Set inflow conditions
    vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
    set_va!(wa, vel_app)
    wa
end

plt.ioff()
@testset "Plotting" begin
    fig = plt.figure(figsize=(14, 14))
    res = plt.plot([1,2,3])
    @test fig isa plt.PyPlot.Figure
    @test res isa Vector{plt.PyObject}
    save_plot(fig, "/tmp", "plot")
    @test isfile("/tmp/plot.pdf")
    rm("/tmp/plot.pdf")
    show_plot(fig)
    wa = create_wa()
    if Sys.islinux()
        fig = plot_geometry(
            wa,
            "Rectangular_wing_geometry";
            data_type=".pdf",
            save_path="/tmp",
            is_save=true,
            is_show=false)
        @test fig isa plt.PyPlot.Figure
        @test isfile("/tmp/Rectangular_wing_geometry_angled_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_angled_view.pdf")
        @test isfile("/tmp/Rectangular_wing_geometry_front_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_front_view.pdf")
        @test isfile("/tmp/Rectangular_wing_geometry_side_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_side_view.pdf")
        @test isfile("/tmp/Rectangular_wing_geometry_top_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_top_view.pdf")

        # Step 5: Initialize the solvers
        vsm_solver = Solver(aerodynamic_model_type=VSM)
        llt_solver = Solver(aerodynamic_model_type=LLT)

        # Step 6: Solve the VSM and LLT
        results_vsm = solve(vsm_solver, wa)
        results_llt = solve(llt_solver, wa)

        # Step 7: Plot spanwise distributions
        y_coordinates = [panel.aero_center[2] for panel in wa.panels]

        fig = plot_distribution(
            [y_coordinates, y_coordinates],
            [results_vsm, results_llt],
            ["VSM", "LLT"],
            title="Spanwise Distributions"
        )
        @test fig isa plt.PyPlot.Figure

        # Step 8: Plot polar curves
        v_a = 20.0            # Magnitude of inflow velocity [m/s]
        angle_range = range(0, 20, 20)
        fig = plot_polars(
            [llt_solver, vsm_solver],
            [wa, wa],
            ["VSM", "LLT"],
            angle_range=angle_range,
            angle_type="angle_of_attack",
            v_a=v_a,
            title="Rectangular Wing Polars",
            data_type=".pdf",
            save_path="/tmp",
            is_save=true,
            is_show=false
        )
        @test fig isa plt.PyPlot.Figure
        @test isfile("/tmp/Rectangular_Wing_Polars.pdf")
        rm("/tmp/Rectangular_Wing_Polars.pdf")
    end
end
plt.ion()
nothing
using VortexStepMethod
using ControlPlots
using Test

function result_vsm()
    # Step 1: Define wing parameters
    n_panels = 20          # Number of panels
    span = 20.0            # Wing span [m]
    chord = 1.0            # Chord length [m]
    Umag = 20.0            # Magnitude of inflow velocity [m/s]
    density = 1.225        # Air density [kg/m³]
    alpha_deg = 30.0       # Angle of attack [degrees]
    alpha = deg2rad(alpha_deg)

    # Step 2: Create wing geometry with linear panel distribution
    wing = Wing(n_panels, spanwise_panel_distribution="linear")

    # Add wing sections - defining only tip sections with inviscid airfoil model
    add_section!(wing, 
        [0.0, span/2, 0.0],   # Left tip LE 
        [chord, span/2, 0.0],  # Left tip TE
        "inviscid")
    add_section!(wing, 
        [0.0, -span/2, 0.0],  # Right tip LE
        [chord, -span/2, 0.0], # Right tip TE
        "inviscid")

    # Step 3: Initialize aerodynamics
    wa = WingAerodynamics([wing])
    # Set inflow conditions
    vel_app = [cos(alpha), 0.0, sin(alpha)] .* Umag
    set_va!(wa, (vel_app, 0.0))  # Second parameter is yaw rate
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
    wa = result_vsm()
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
    end
end
plt.ion()
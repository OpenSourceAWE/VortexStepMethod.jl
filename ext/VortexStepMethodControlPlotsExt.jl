module VortexStepMethodControlPlotsExt
using ControlPlots, LaTeXStrings, VortexStepMethod, LinearAlgebra, Statistics, DelimitedFiles
import ControlPlots: plt
import VortexStepMethod: calculate_filaments_for_plotting

export plot_wing, plot_circulation_distribution, plot_geometry, plot_distribution,
       plot_polars, save_plot, show_plot, plot_polar_data, plot_combined_analysis

"""
    set_plot_style(titel_size=16; use_tex=false)

Set the default style for plots using LaTeX.
``
# Arguments:
- `titel_size`: size of the plot title in points (default: 16)
- `ùse_tex`: if the external `pdflatex` command shall be used
"""
function set_plot_style(titel_size=16; use_tex=false)
    rcParams = plt.PyDict(plt.matplotlib."rcParams")
    rcParams["text.usetex"] = use_tex
    rcParams["font.family"] = "serif"
    if use_tex
        rcParams["font.serif"] = ["Computer Modern Roman"]
    end
    rcParams["axes.titlesize"] = titel_size
    rcParams["axes.labelsize"] = 12
    rcParams["axes.linewidth"] = 1
    rcParams["lines.linewidth"] = 1
    rcParams["lines.markersize"] = 6
    rcParams["xtick.labelsize"] = 10
    rcParams["ytick.labelsize"] = 10
    rcParams["legend.fontsize"] = 10
    rcParams["figure.titlesize"] = 16
    if use_tex
        rcParams["pgf.texsystem"] = "pdflatex"  # Use pdflatex
    end
    rcParams["pgf.rcfonts"] = false
    rcParams["figure.figsize"] = (10, 6)  # Default figure size
end


"""
    save_plot(fig, save_path, title; data_type=".pdf")

Save a plot to a file.

# Arguments
- `fig`: Plots figure object
- `save_path`: Path to save the plot
- `title`: Title of the plot

# Keyword arguments
- `data_type`: File extension (default: ".pdf")
"""
function VortexStepMethod.save_plot(fig, save_path, title; data_type=".pdf")
    isnothing(save_path) && throw(ArgumentError("save_path should be provided"))

    !isdir(save_path) && mkpath(save_path)
    full_path = joinpath(save_path, title * data_type)

    @debug "Attempting to save figure to: $full_path"
    @debug "Current working directory: $(pwd())"

    try
        hasproperty(fig, :savefig) || throw(ArgumentError(
            "Figure object of type $(typeof(fig)) does not support savefig()."
        ))
        getproperty(fig, :savefig)(full_path)
        @debug "Figure saved as $data_type"

        if isfile(full_path)
            @debug "File successfully saved to $full_path"
            @debug "File size: $(filesize(full_path)) bytes"
        else
            @info "File does not exist after save attempt: $full_path"
        end
    catch e
        @error "Error saving figure: $e"
        @error "Error type: $(typeof(e))"
        rethrow(e)
    end
end

"""
    show_plot(fig::plt.Figure; dpi=130)

Display a plot at specified DPI.

# Arguments
- `fig`: Plots figure object

# Keyword arguments
- `dpi`: Dots per inch for the figure (default: 130)
"""
function VortexStepMethod.show_plot(fig::plt.Figure; dpi=130)
    fig.set_dpi(dpi)
    plt.display(fig)
end

"""
    plot_line_segment!(ax, segment, color, label; width=3)

Plot a line segment in 3D with arrow.

# Arguments
- `ax`: Plot axis
- `segment`: Array of two points defining the segment
- `color`: Color of the segment
- `label`: Label for the legend

# Keyword Arguments
- `width`: Line width (default: 3)
"""
function plot_line_segment!(ax, segment, color, label; width=3)
    ax.plot(
        [segment[1][1], segment[2][1]],
        [segment[1][2], segment[2][2]],
        [segment[1][3], segment[2][3]],
        color=color, label=label, linewidth=width
    )

    dir = segment[2] - segment[1]
    ax.quiver(
        [segment[1][1]], [segment[1][2]], [segment[1][3]],
        [dir[1]], [dir[2]], [dir[3]],
        color=color
    )
end

"""
    set_axes_equal!(ax; zoom=1.8)

Set 3D plot axes to equal scale.

# Arguments
- ax: 3D plot axis

# Keyword arguments
zoom: zoom factor (default: 1.8)
"""
function set_axes_equal!(ax; zoom=1.8)
    x_lims = ax.get_xlim3d() ./ zoom
    y_lims = ax.get_ylim3d() ./ zoom
    z_lims = ax.get_zlim3d() ./ zoom

    x_range = abs(x_lims[2] - x_lims[1])
    y_range = abs(y_lims[2] - y_lims[1])
    z_range = abs(z_lims[2] - z_lims[1])

    max_range = max(x_range, y_range, z_range)

    x_mid = mean(x_lims)
    y_mid = mean(y_lims)
    z_mid = mean(z_lims)

    ax.set_xlim3d((x_mid - max_range / 2, x_mid + max_range / 2))
    ax.set_ylim3d((y_mid - max_range / 2, y_mid + max_range / 2))
    ax.set_zlim3d((z_mid - max_range / 2, z_mid + max_range / 2))
end

"""
    create_geometry_plot(body_aero::BodyAerodynamics, title, view_elevation, view_azimuth;
                         zoom=1.8, use_tex=false)

Create a 3D plot of wing geometry including panels and filaments.

# Arguments
- body_aero: struct of type BodyAerodynamics
- title: plot title
- view_elevation: initial view elevation angle [°]
- view_azimuth: initial view azimuth angle [°]

# Keyword arguments
- zoom: zoom factor (default: 1.8)
"""
function create_geometry_plot(body_aero::BodyAerodynamics, title, view_elevation, view_azimuth;
                              zoom=1.8, use_tex=false)
    set_plot_style(28; use_tex)

    panels = body_aero.panels
    isempty(panels) && throw(ArgumentError("Cannot plot geometry: body_aero.panels is empty."))

    va = if body_aero.has_distributed_va
        body_aero._va
    else
        isa(body_aero.va, Tuple) ? body_aero.va[1] : body_aero.va
    end

    # Extract geometric data
    control_points = [panel.control_point for panel in panels]
    aero_centers = [panel.aero_center for panel in panels]

    # Create plot
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title(title)

    # Plot panels
    legend_used = Dict{String,Bool}()
    for (i, panel) in enumerate(panels)
        # Plot panel edges and surfaces
        corners = Matrix{Float64}(panel.corner_points)
        x_corners = corners[1, :]
        y_corners = corners[2, :]
        z_corners = corners[3, :]

        push!(x_corners, x_corners[1])
        push!(y_corners, y_corners[1])
        push!(z_corners, z_corners[1])

        ax.plot(x_corners,
            y_corners,
            z_corners,
            color=:grey,
            linewidth=1,
            label=i == 1 ? "Panel Edges" : "")

        # Plot control points and aerodynamic centers
        ax.scatter([control_points[i][1]], [control_points[i][2]], [control_points[i][3]],
            color=:green, label=i == 1 ? "Control Points" : "")
        ax.scatter([aero_centers[i][1]], [aero_centers[i][2]], [aero_centers[i][3]],
            color=:blue, label=i == 1 ? "Aerodynamic Centers" : "")

        # Plot filaments
        filaments = calculate_filaments_for_plotting(panel)
        legends = ["Bound Vortex", "side1", "side2", "wake_1", "wake_2"]

        for (filament, legend) in zip(filaments, legends)
            x1, x2, color = filament
            @debug "Legend: $legend"
            show_legend = !get(legend_used, legend, false)
            plot_line_segment!(ax, [x1, x2], color, show_legend ? legend : "")
            legend_used[legend] = true
        end
    end

    # Plot velocity vector
    max_chord = maximum(panel.chord for panel in panels)
    va_mag = norm(va)
    va_vector_begin = -2 * max_chord * va / va_mag
    va_vector_end = va_vector_begin + 1.5 * va / va_mag
    plot_line_segment!(ax, [va_vector_begin, va_vector_end], :lightblue, "va")

    # Add legends for the first occurrence of each label
    # by_label = Dict(zip(labels, handles))
    # ax.legend(values(by_label), keys(by_label), bbox_to_anchor=(0, 0, 1.1, 1))

    # Set labels and make axes equal
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    set_axes_equal!(ax; zoom)

    # Set the initial view
    ax.view_init(elev=view_elevation, azim=view_azimuth)

    # Ensure the figure is fully rendered
    # fig.canvas.draw()
    plt.tight_layout(rect=(0, 0, 1, 0.97))

    return fig
end

"""
    plot_geometry(body_aero::BodyAerodynamics, title;
                  data_type=".pdf", save_path=nothing,
                  is_save=false, is_show=false,
                  view_elevation=15, view_azimuth=-120, use_tex=false)

Plot wing geometry from different viewpoints and optionally save/show plots.

# Arguments:
- `body_aero`: the [BodyAerodynamics](@ref) to plot
- title: plot title

# Keyword arguments:
- `data_type``: string with the file type postfix (default: ".pdf")
- `save_path`: path for saving the graphic (default: `nothing`)
- `is_save`: boolean value, indicates if the graphic shall be saved (default: `false`)
- `is_show`: boolean value, indicates if the graphic shall be displayed (default: `false`)
- `view_elevation`: initial view elevation angle (default: 15) [°]
- `view_azimuth`: initial view azimuth angle (default: -120) [°]
- `use_tex`: if the external `pdflatex` command shall be used (default: false)

"""
function VortexStepMethod.plot_geometry(body_aero::BodyAerodynamics, title;
    data_type=".pdf",
    save_path=nothing,
    is_save=false,
    is_show=false,
    view_elevation=15,
    view_azimuth=-120,
    use_tex=false)

    if is_save
        plt.ioff()
        # Angled view
        fig = create_geometry_plot(body_aero, "$(title)_angled_view", 15, -120; use_tex)
        save_plot(fig, save_path, "$(title)_angled_view", data_type=data_type)

        # Top view
        fig = create_geometry_plot(body_aero, "$(title)_top_view", 90, 0; use_tex)
        save_plot(fig, save_path, "$(title)_top_view", data_type=data_type)

        # Front view
        fig = create_geometry_plot(body_aero, "$(title)_front_view", 0, 0; use_tex)
        save_plot(fig, save_path, "$(title)_front_view", data_type=data_type)

        # Side view
        fig = create_geometry_plot(body_aero, "$(title)_side_view", 0, -90; use_tex)
        save_plot(fig, save_path, "$(title)_side_view", data_type=data_type)
    end

    if is_show
        plt.ion()
        fig = create_geometry_plot(body_aero, title, view_elevation, view_azimuth; use_tex)
        plt.display(fig)
    else
        fig = create_geometry_plot(body_aero, title, view_elevation, view_azimuth; use_tex)
    end
    fig
end

"""
    plot_distribution(y_coordinates_list, results_list, label_list;
                      title="spanwise_distribution", data_type=".pdf",
                      save_path=nothing, is_save=false, is_show=true, use_tex=false)

Plot spanwise distributions of aerodynamic properties.

# Arguments
- `y_coordinates_list`: List of spanwise coordinates
- `results_list`: List of result dictionaries
- `label_list`: List of labels for different results

# Keyword arguments
- `title`: Plot title (default: "spanwise_distribution")
- `data_type`: File extension for saving (default: ".pdf")
- `save_path`: Path to save plots (default: nothing)
- `is_save`: Whether to save plots (default: false)
- `is_show`: Whether to display plots (default: true)
- `use_tex`: if the external `pdflatex` command shall be used
"""
function VortexStepMethod.plot_distribution(y_coordinates_list, results_list, label_list;
    title="spanwise_distribution",
    data_type=".pdf",
    save_path=nothing,
    is_save=false,
    is_show=true,
    use_tex=false)

    length(results_list) == length(label_list) || throw(ArgumentError(
        "Number of results ($(length(results_list))) must match number of labels ($(length(label_list)))"
    ))

    # Set the plot style
    set_plot_style(; use_tex)

    # Initializing plot
    fig, axs = plt.subplots(3, 3, figsize=(16, 10))
    fig.suptitle(title, fontsize=16)

    # CL plot
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        value = "$(round(result_i["cl"], digits=2))"
        if label_i == "LLT"
            label = label_i * L" $~C_\mathrm{L}$: " * value
        else
            label = label_i * L" $C_\mathrm{L}$: " * value
        end
        axs[1, 1].plot(
            y_coordinates_i,
            result_i["cl_distribution"],
            label=label
        )
    end
    axs[1, 1].set_title(L"$C_\mathrm{L}$ Distribution", size=16)
    axs[1, 1].set_xlabel(L"Spanwise Position $y/b$")
    axs[1, 1].set_ylabel(L"Lift Coefficient $C_\mathrm{L}$")
    axs[1, 1].legend()

    # CD plot
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        value = "$(round(result_i["cl"], digits=2))"
        if label_i == "LLT"
            label = label_i * L" $~C_\mathrm{D}$: " * value
        else
            label = label_i * L" $C_\mathrm{D}$: " * value
        end
        axs[1, 2].plot(
            y_coordinates_i,
            result_i["cd_distribution"],
            label=label
        )
    end
    axs[1, 2].set_title(L"$C_\mathrm{D}$ Distribution", size=16)
    axs[1, 2].set_xlabel(L"Spanwise Position $y/b$")
    axs[1, 2].set_ylabel(L"Drag Coefficient $C_\mathrm{D}$")
    axs[1, 2].legend()

    # Gamma Distribution
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[1, 3].plot(
            y_coordinates_i,
            result_i["gamma_distribution"],
            label=label_i
        )
    end
    axs[1, 3].set_title(L"\Gamma~Distribution", size=16)
    axs[1, 3].set_xlabel(L"Spanwise Position $y/b$")
    axs[1, 3].set_ylabel(L"Circulation~\Gamma")
    axs[1, 3].legend()

    # Geometric Alpha
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[2, 1].plot(
            y_coordinates_i,
            result_i["alpha_geometric"],
            label=label_i
        )
    end
    axs[2, 1].set_title(L"$\alpha$ Geometric", size=16)
    axs[2, 1].set_xlabel(L"Spanwise Position $y/b$")
    axs[2, 1].set_ylabel(L"Angle of Attack $\alpha$ (deg)")
    axs[2, 1].legend()

    # Calculated/ Corrected Alpha
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[2, 2].plot(
            y_coordinates_i,
            result_i["alpha_at_ac"],
            label=label_i
        )
    end
    axs[2, 2].set_title(L"$\alpha$ result (corrected to aerodynamic center)", size=16)
    axs[2, 2].set_xlabel(L"Spanwise Position $y/b$")
    axs[2, 2].set_ylabel(L"Angle of Attack $\alpha$ (deg)")
    axs[2, 2].legend()

    # Uncorrected Alpha plot
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[2, 3].plot(
            y_coordinates_i,
            result_i["alpha_uncorrected"],
            label=label_i
        )
    end
    axs[2, 3].set_title(L"$\alpha$ Uncorrected (if VSM, at the control point)", size=16)
    axs[2, 3].set_xlabel(L"Spanwise Position $y/b$")
    axs[2, 3].set_ylabel(L"Angle of Attack $\alpha$ (deg)")
    axs[2, 3].legend()

    # Force Components
    for (idx, component) in enumerate(["x", "y", "z"])
        axs[3, idx].set_title("Force in $component direction", size=16)
        axs[3, idx].set_xlabel(L"Spanwise Position $y/b$")
        axs[3, idx].set_ylabel(raw"$F_\mathrm" * "{$component}" * raw"$")
        for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
            # Extract force components for the current direction (idx)
            forces = results["F_distribution"][idx, :]
            # Verify dimensions match
            if length(y_coords) != length(forces)
                @warn "Dimension mismatch in force plotting" length(y_coords) length(forces) component
                continue  # Skip this component instead of throwing error
            end
            space = ""
            if label == "LLT"
                space = "~"
            end
            axs[3, idx].plot(
                y_coords,
                forces,
                label="$label" * space * raw"$~\Sigma~F_\mathrm" * "{$component}:~" *
                      raw"$" * "$(round(results["F$component"], digits=2)) N"
            )
            axs[3, idx].legend()
        end
    end

    fig.tight_layout()

    # Save and show plot
    if is_save
        save_plot(fig, save_path, title, data_type=data_type)
    end

    if is_show
        show_plot(fig)
    end

    return fig
end

"""
    generate_polar_data(solver, body_aero::BodyAerodynamics, angle_range;
                        angle_type="angle_of_attack", angle_of_attack=0.0,
                        side_slip=0.0, v_a=10.0, use_latex=false)

Generate polar data for aerodynamic analysis over a range of angles.

# Arguments
- `solver`: Aerodynamic solver object
- `body_aero`: Wing aerodynamics struct
- `angle_range`: Range of angles to analyze

# Keyword arguments
- `angle_type`: Type of angle variation ("angle_of_attack" or "side_slip")
- `angle_of_attack`: Initial angle of attack [rad]
- `side_slip`: Initial side slip angle in [rad]
- `v_a`: norm of apparent wind speed [m/s]

# Returns
- Tuple of polar data array and Reynolds number
"""

"""
    plot_polars(solver_list, body_aero_list, label_list;
                literature_path_list=String[],
                angle_range=range(0, 20, 2), angle_type="angle_of_attack",
                angle_of_attack=0.0, side_slip=0.0, v_a=10.0,
                title="polar", data_type=".pdf", save_path=nothing,
                is_save=true, is_show=true, use_tex=false)

Plot polar data comparing different solvers and configurations.

# Arguments
- `solver_list`: List of aerodynamic solvers
- `body_aero_list`: List of wing aerodynamics objects
- `label_list`: List of labels for each configuration

# Keyword arguments
- `literature_path_list`: Optional paths to literature data files
- `angle_range`: Range of angles to analyze [°]
- `angle_type`: "`angle_of_attack`" or "`side_slip`"; (default: `angle_of_attack`)
- `angle_of_attack:` AoA to be used for plotting the polars (default: 0.0) [rad]
- `side_slip`: side slip angle (default: 0.0) [rad]
- v_a: norm of apparent wind speed (default: 10.0) [m/s]
- title: plot title
- `data_type`: File extension for saving (default: ".pdf")
- `save_path`: Path to save plots (default: nothing)
- `is_save`: Whether to save plots (default: true)
- `is_show`: Whether to display plots (default: true)
- `use_tex`: if the external `pdflatex` command shall be used (default: false)
- `cl_over_cd`: Plot CL/CD vs angle instead of CL vs CD (default: true)
"""
function VortexStepMethod.plot_polars(
    solver_list,
    body_aero_list,
    label_list;
    literature_path_list=String[],
    angle_range=range(0, 20, 2),
    angle_type="angle_of_attack",
    angle_of_attack=0.0,
    side_slip=0.0,
    v_a=10.0,
    title="polar",
    data_type=".pdf",
    save_path=nothing,
    is_save=true,
    is_show=true,
    use_tex=false,
    cl_over_cd=true,
    show_moments=false,
)
    # Validate inputs
    total_cases = length(body_aero_list) + length(literature_path_list)
    if total_cases != length(label_list) || length(solver_list) != length(body_aero_list)
        throw(ArgumentError("Mismatch in number of solvers ($(length(solver_list))), " *
                            "cases ($total_cases), and labels ($(length(label_list)))"))
    end
    main_title = replace(title, " " => "_")
    set_plot_style(; use_tex)

    # Generate polar data
    polar_data_list = []
    cm_data_list = []
    for (i, (solver, body_aero)) in enumerate(zip(solver_list, body_aero_list))
        result = VortexStepMethod.generate_polar_data(
            solver, body_aero, angle_range;
            angle_type,
            angle_of_attack,
            side_slip,
            v_a
        )
        push!(polar_data_list, result.polar_data)
        push!(cm_data_list, (cmx=result.cmx, cmy=result.cmy,
                             cmz=result.cmz))
        # Update label with Reynolds number
        label_list[i] = "$(label_list[i]) Re = $(round(Int64, result.rey*1e-5))e5"
    end
    # Load literature data if provided
    if !isempty(literature_path_list)
        for path in literature_path_list
            raw_data = readdlm(path, ',')
            lit = VortexStepMethod.extract_literature_polar_data(
                raw_data, path; angle_type)
            push!(polar_data_list, lit.polar_data)
            push!(cm_data_list, (cmx=lit.cmx, cmy=lit.cmy,
                                 cmz=lit.cmz))
        end
    end

    # Number of computational results (excluding literature)
    n_solvers = length(solver_list)

    # Helper: format label and line style
    function format_label(label, i, n_solvers)
        if i < n_solvers
            linestyle, marker, markersize = "-", "*", 7
        else
            linestyle, marker, markersize = "-", ".", 5
        end
        if use_tex
            if contains(label, "LLT")
                label = replace(label, "e5" => raw"\cdot10^5")
                label = replace(label, " " => raw"~")
                label = replace(label,
                    "LLT" => raw"\mathrm{LLT}{~\,}")
                label = raw"$" * label * raw"$"
            else
                label = replace(label, "e5" => raw"\cdot10^5")
                label = replace(label, " " => "~")
                label = replace(label,
                    "VSM" => raw"\mathrm{VSM}")
                label = raw"$" * label * raw"$"
            end
        end
        return label, linestyle, marker, markersize
    end

    if show_moments
        # 2x3 layout: CL, CD, CS, CMx, CMy, CMz
        fig, axs = plt.subplots(2, 3, figsize=(21, 14))
        coeff_specs = [
            (1, raw"$C_\mathrm{L}$", 2, nothing),
            (2, raw"$C_\mathrm{D}$", 3, nothing),
            (3, raw"$C_\mathrm{S}$", 4, nothing),
            (4, raw"$C_\mathrm{Mx}$", nothing, :cmx),
            (5, raw"$C_\mathrm{My}$", nothing, :cmy),
            (6, raw"$C_\mathrm{Mz}$", nothing, :cmz),
        ]
        for (ax_idx, ylabel, pd_col, cm_field) in coeff_specs
            row = (ax_idx - 1) ÷ 3 + 1
            col = (ax_idx - 1) % 3 + 1
            ax = axs[row, col]
            for (i, (polar_data, cm, label)) in enumerate(
                    zip(polar_data_list, cm_data_list,
                        label_list))
                label, ls, mk, ms = format_label(
                    label, i, n_solvers)
                if pd_col !== nothing
                    y_data = polar_data[pd_col]
                else
                    y_data = Float64.(getfield(cm, cm_field))
                    all(isnan, y_data) && continue
                end
                ax.plot(polar_data[1], y_data;
                    label=label, linestyle=ls,
                    marker=mk, markersize=ms)
            end
            ax.set_title(
                ylabel * " vs $angle_type [°]")
            ax.set_xlabel("$angle_type [°]")
            ax.set_ylabel(ylabel)
            ax.legend()
        end
    else
        # 2x2 layout: CL, CD, CS, CL/CD or CL-vs-CD
        fig, axs = plt.subplots(2, 2, figsize=(14, 14))
        coeff_specs = [
            ((1, 1), raw"$C_\mathrm{L}$", 2),
            ((1, 2), raw"$C_\mathrm{D}$", 3),
            ((2, 1), raw"$C_\mathrm{S}$", 4),
        ]
        for (pos, ylabel, pd_col) in coeff_specs
            ax = axs[pos...]
            for (i, (polar_data, label)) in enumerate(
                    zip(polar_data_list, label_list))
                label, ls, mk, ms = format_label(
                    label, i, n_solvers)
                ax.plot(polar_data[1], polar_data[pd_col];
                    label=label, linestyle=ls,
                    marker=mk, markersize=ms)
                if maximum(polar_data[2]) > 10
                    ax.set_ylim([-0.5, 2])
                end
            end
            ax.set_title(
                ylabel * " vs $angle_type [°]")
            ax.set_xlabel("$angle_type [°]")
            ax.set_ylabel(ylabel)
            ax.legend()
        end
        # Fourth panel: CL/CD or CL-vs-CD
        ax4 = axs[2, 2]
        for (i, (polar_data, label)) in enumerate(
                zip(polar_data_list, label_list))
            label, ls, mk, ms = format_label(
                label, i, n_solvers)
            if cl_over_cd
                cl_cd = polar_data[2] ./ polar_data[3]
                ax4.plot(polar_data[1], cl_cd;
                    label=label, linestyle=ls,
                    marker=mk, markersize=ms)
            else
                ax4.plot(polar_data[3], polar_data[2];
                    label=label, linestyle=ls,
                    marker=mk, markersize=ms)
            end
        end
        if cl_over_cd
            ax4.set_title(
                raw"$C_\mathrm{L}/C_\mathrm{D}$" *
                " vs $angle_type [°]")
            ax4.set_xlabel("$angle_type [°]")
            ax4.set_ylabel(
                raw"$C_\mathrm{L}/C_\mathrm{D}$")
        else
            ax4.set_title(
                raw"$C_\mathrm{L}$" * " vs " *
                raw"$C_\mathrm{D}$")
            ax4.set_xlabel(raw"$C_\mathrm{D}$")
            ax4.set_ylabel(raw"$C_\mathrm{L}$")
        end
        ax4.legend()
    end

    fig.tight_layout(h_pad=3.5, rect=(0.01, 0.01, 0.99, 0.99))

    # Save and show plot
    if is_save && !isnothing(save_path)
        save_plot(fig, save_path, main_title; data_type)
    end

    if is_show
        show_plot(fig)
    end

    return fig
end

"""
    plot_polar_data(body_aero::BodyAerodynamics; alphas=collect(deg2rad.(-5:0.3:25)), delta_tes=collect(deg2rad.(-5:0.3:25)))

Plot polar data (Cl, Cd, Cm) as 3D surfaces against alpha and delta_te angles. delta_te is the trailing edge deflection angle
relative to the 2d airfoil or panel chord line.

# Arguments
- `body_aero`: Wing aerodynamics struct

# Keyword arguments
- `alphas`: Range of angle of attack values in radians (default: -5° to 25° in 0.3° steps)
- `delta_tes`: Range of trailing edge angles in radians (default: -5° to 25° in 0.3° steps)
- `is_show`: Whether to display plots (default: true)
- `use_tex`: if the external `pdflatex` command shall be used
"""
function VortexStepMethod.plot_polar_data(body_aero::BodyAerodynamics;
        alphas=collect(deg2rad.(-5:0.3:25)),
        delta_tes = collect(deg2rad.(-5:0.3:25)),
        is_show = true,
        use_tex = false
        )
    if body_aero.panels[1].aero_model == POLAR_MATRICES
    set_plot_style(; use_tex)

        # Create figure with subplots
        fig = plt.figure(figsize=(15, 6))

        # Get interpolation functions and labels
        interp_data = [
            (body_aero.panels[1].cl_interp, L"$C_l$"),
            (body_aero.panels[1].cd_interp, L"$C_d$"),
            (body_aero.panels[1].cm_interp, L"$C_m$")
        ]

        # Create each subplot
        for (idx, (interp, label)) in enumerate(interp_data)
            ax = fig.add_subplot(1, 3, idx, projection="3d")

            # Create interpolation matrix
            interp_matrix = zeros(length(alphas), length(delta_tes))
            interp_matrix .= [interp(alpha, delta_te) for alpha in alphas, delta_te in delta_tes]
            X = collect(delta_tes) .+ zeros(length(alphas))'
            Y = collect(alphas)' .+ zeros(length(delta_tes))

            # Plot surface
            ax.plot_wireframe(X, Y, interp_matrix,
                            edgecolor="blue",
                            lw=0.5,
                            rstride=5,
                            cstride=5,
                            alpha=0.6)

            # Set labels and title
            ax.set_xlabel(L"$\delta$ [rad]")
            ax.set_ylabel(L"$\alpha$ [rad]")
            ax.set_zlabel(label)
            ax.set_title(label * L" vs $\alpha$ and $\delta$")
            ax.grid(true)
        end

        # Adjust layout and display
        plt.tight_layout(rect=(0.01, 0.01, 0.99, 0.99))
        if is_show
            show_plot(fig)
        end
        return fig
    else
        throw(ArgumentError("Plotting polar data for $(body_aero.panels[1].aero_model) is not implemented."))
    end
end

"""
    plot_combined_analysis(solver, body_aero, results; kwargs...)

Create combined analysis by calling plot_geometry, plot_distribution,
and plot_polars sequentially. Each creates a separate matplotlib window.

# Arguments
- `solver`: Solver or array of solvers
- `body_aero`: BodyAerodynamics object or array
- `results`: Results dict or array of results dicts

See individual functions for detailed parameter descriptions.
"""
function VortexStepMethod.plot_combined_analysis(
    solver,
    body_aero,
    results;
    solver_label="VSM",
    labels=nothing,
    angle_range=range(0, 20, length=20),
    angle_type="angle_of_attack",
    angle_of_attack=0.0,
    side_slip=0.0,
    v_a=10.0,
    title="Combined Analysis",
    data_type=".pdf",
    save_path=nothing,
    is_save=false,
    is_show=true,
    view_elevation=15,
    view_azimuth=-120,
    use_tex=false,
    literature_path_list=String[],
    cl_over_cd=true,
    angle_of_attack_for_spanwise_distribution=5.0,
)
    # Normalize inputs to arrays for consistent handling
    solvers = solver isa Vector ? solver : [solver]
    body_aeros = body_aero isa Vector ? body_aero : [body_aero]
    results_list = results isa Vector ? results : [results]
    n_solvers = length(solvers)
    n_literature = length(literature_path_list)

    # Label normalization (matches Makie version)
    label_source = isnothing(labels) ? solver_label : labels
    labels_in = label_source isa AbstractVector ?
        string.(label_source) : [string(label_source)]
    solver_labels = length(labels_in) == 1 ?
        fill(labels_in[1], n_solvers) :
        labels_in[1:n_solvers]

    # Compute spanwise results at specified AoA
    results_spanwise = copy(results_list)
    if !isnothing(angle_of_attack_for_spanwise_distribution)
        α_span = deg2rad(
            angle_of_attack_for_spanwise_distribution)
        β_span = deg2rad(side_slip)
        for (i, (s, ba)) in enumerate(
                zip(solvers, body_aeros))
            va_old = copy(getfield(ba, :_va))
            omega_old = copy(ba.omega)
            set_va!(ba, [cos(α_span) * cos(β_span),
                sin(β_span), sin(α_span)] * v_a)
            results_spanwise[i] = solve(s, ba,
                s.sol.gamma_distribution)
            set_va!(ba, va_old, omega_old)
        end
    end

    # Extract y-coordinates for distribution plot
    y_coords_list = [
        [p.aero_center[2] for p in ba.panels]
        for ba in body_aeros]

    # Plot geometry (first body_aero only)
    plot_geometry(
        body_aeros[1],
        title;
        data_type, save_path, is_save, is_show,
        view_elevation, view_azimuth, use_tex
    )

    # Plot spanwise distributions
    plot_distribution(
        y_coords_list,
        results_spanwise,
        solver_labels;
        title=title * " - Distributions",
        data_type, save_path, is_save, is_show, use_tex
    )

    # Plot polars (include literature labels)
    polar_labels = if n_literature > 0 &&
            length(labels_in) == n_solvers + n_literature
        labels_in
    else
        solver_labels
    end
    plot_polars(
        solvers,
        body_aeros,
        polar_labels;
        literature_path_list, angle_range, angle_type,
        angle_of_attack, side_slip, v_a,
        title=title * " - Polars",
        data_type, save_path, is_save, is_show,
        use_tex, cl_over_cd
    )
end

end
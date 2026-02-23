module VortexStepMethodMakieExt
using Makie, VortexStepMethod, LinearAlgebra, Statistics, DelimitedFiles
import VortexStepMethod: calculate_filaments_for_plotting

export plot_geometry, plot_distribution, plot_polars, save_plot, show_plot,
    plot_polar_data, plot_combined_analysis, plot_airfoil_slices

# Global storage for panel mesh observables (for dynamic plotting)
const PANEL_MESH_OBSERVABLES = Ref{Union{Nothing,Dict}}(nothing)

"""
    plot!(ax, panel::VortexStepMethod.Panel; use_observables=false, kwargs...)

Plot a single `Panel` as a `mesh`.
The corner points are ordered as: LE1, TE1, TE2, LE2.
This creates two triangles: (LE1, TE1, TE2) and (LE1, TE2, LE2).

If `use_observables=true`, creates observables for dynamic updates.
"""
function Makie.plot!(ax, panel::VortexStepMethod.Panel; color=(:red, 0.2), R_b_w=nothing, T_b_w=nothing,
    use_observables=false, kwargs...)
    plots = []
    points = [Point3f(panel.corner_points[:, i]) for i in 1:4]
    if !isnothing(R_b_w) && !isnothing(T_b_w)
        points = [Point3f(R_b_w * p + T_b_w) for p in points]
    end

    if use_observables
        # Create observables for dynamic updates
        vertices_obs = Observable(points)
        faces_obs = Observable([Makie.GLTriangleFace(1, 2, 3), Makie.GLTriangleFace(1, 3, 4)])
        border_obs = Observable([points..., points[1]])

        p = mesh!(ax, vertices_obs, faces_obs; color, transparency=true, kwargs...)
        push!(plots, p)
        p = lines!(ax, border_obs; color=:black, transparency=true, kwargs...)
        push!(plots, p)

        # Note: Observables are stored at the body level, not individual panel level
        # Individual panels need their parent body for proper tracking
    else
        # Static plotting (original behavior)
        faces = [Makie.GLTriangleFace(1, 2, 3), Makie.GLTriangleFace(1, 3, 4)]
        p = mesh!(ax, points, faces; color, transparency=true, kwargs...)
        push!(plots, p)
        border_points = [points..., points[1]]
        p = lines!(ax, border_points; color=:black, transparency=true, kwargs...)
        push!(plots, p)
    end

    return plots
end

"""
    plot!(ax, body::VortexStepMethod.BodyAerodynamics; use_observables=false, kwargs...)

Plot a `BodyAerodynamics` object by plotting each of its panels.

If `use_observables=true`, creates observables for dynamic updates keyed by (body_id, panel_index).
Otherwise, creates static plots (original behavior).
"""
function Makie.plot!(ax, body::VortexStepMethod.BodyAerodynamics; color=(:red, 0.2), R_b_w=nothing, T_b_w=nothing,
    use_observables=false, kwargs...)
    plots = []

    if use_observables
        # Initialize global storage if needed
        if isnothing(PANEL_MESH_OBSERVABLES[])
            PANEL_MESH_OBSERVABLES[] = Dict()
        end

        body_id = objectid(body)

        # Create observables for each panel
        for (panel_idx, panel) in enumerate(body.panels)
            # Compute initial points
            points = [Point3f(panel.corner_points[:, i]) for i in 1:4]
            if !isnothing(R_b_w) && !isnothing(T_b_w)
                points = [Point3f(R_b_w * p + T_b_w) for p in points]
            end

            # Create observables
            vertices_obs = Observable(points)
            faces_obs = Observable([Makie.GLTriangleFace(1, 2, 3), Makie.GLTriangleFace(1, 3, 4)])
            border_obs = Observable([points..., points[1]])

            # Plot using observables
            p = mesh!(ax, vertices_obs, faces_obs; color, transparency=true, kwargs...)
            push!(plots, p)
            p = lines!(ax, border_obs; color=:black, transparency=true, kwargs...)
            push!(plots, p)

            # Store observables with stable key
            PANEL_MESH_OBSERVABLES[][(body_id, panel_idx)] = (
                vertices=vertices_obs,
                border=border_obs,
                faces=faces_obs
            )
        end
    else
        # Static plotting (original behavior)
        for panel in body.panels
            p = Makie.plot!(ax, panel; color, R_b_w, T_b_w, use_observables=false, kwargs...)
            push!(plots, p)
        end
    end

    return plots
end


"""
    plot!(body::VortexStepMethod.BodyAerodynamics; R_b_w=nothing, T_b_w=nothing)

Update existing body aerodynamics plot observables with current geometry.
This updates all panels in the body using their current corner_points.

Requires that `plot(body; use_observables=true)` or `plot!(ax, body; use_observables=true)`
was called first to create the observables.
"""
function Makie.plot!(body::VortexStepMethod.BodyAerodynamics; R_b_w=nothing, T_b_w=nothing, kwargs...)
    # Check if observables exist
    if isnothing(PANEL_MESH_OBSERVABLES[])
        error("No panel observables found. Call plot(body; use_observables=true) first.")
    end

    body_id = objectid(body)

    # Update each panel using stable (body_id, panel_idx) key
    for (panel_idx, panel) in enumerate(body.panels)
        key = (body_id, panel_idx)
        if !haskey(PANEL_MESH_OBSERVABLES[], key)
            error("No observables found for body $body_id panel $panel_idx. " *
                  "Call plot(body; use_observables=true) first.")
        end

        # Get observables for this panel
        obs = PANEL_MESH_OBSERVABLES[][key]

        # Recompute vertices from current panel.corner_points
        points = [Point3f(panel.corner_points[:, i]) for i in 1:4]
        if !isnothing(R_b_w) && !isnothing(T_b_w)
            points = [Point3f(R_b_w * p + T_b_w) for p in points]
        end

        # Update observables
        obs.vertices[] = points
        obs.border[] = [points..., points[1]]
    end

    return nothing
end

function Makie.plot(panel::VortexStepMethod.Panel; size=(1200, 800),
    R_b_w=nothing, T_b_w=nothing, color=(:red, 0.2), kwargs...)
    fig = Figure(; size)
    ax = Axis3(fig[1, 1]; aspect=:data,
        xlabel="X", ylabel="Y", zlabel="Z",
        azimuth=9 / 8 * π, zoommode=:cursor, viewmode=:fit,
    )

    # Create observables for panel geometry
    points = [Point3f(panel.corner_points[:, i]) for i in 1:4]
    if !isnothing(R_b_w) && !isnothing(T_b_w)
        points = [Point3f(R_b_w * p + T_b_w) for p in points]
    end

    vertices_obs = Observable(points)
    faces_obs = Observable([Makie.GLTriangleFace(1, 2, 3), Makie.GLTriangleFace(1, 3, 4)])

    # Plot mesh using observables
    mesh!(ax, vertices_obs, faces_obs; color, transparency=true, kwargs...)

    # Plot border
    border_obs = Observable([points..., points[1]])
    lines!(ax, border_obs; color=:black, transparency=true, kwargs...)

    # Store observables globally for updates
    panel_id = objectid(panel)
    if isnothing(PANEL_MESH_OBSERVABLES[])
        PANEL_MESH_OBSERVABLES[] = Dict()
    end
    PANEL_MESH_OBSERVABLES[][panel_id] = (
        vertices=vertices_obs,
        border=border_obs,
        faces=faces_obs
    )

    return fig
end

function Makie.plot(body_aero::VortexStepMethod.BodyAerodynamics; size=(1200, 800),
    limitmargin=0.1, R_b_w=nothing, T_b_w=nothing, color=(:red, 0.2),
    kwargs...)
    fig = Figure(; size)
    ax = Axis3(fig[1, 1]; aspect=:data,
        xlabel="X", ylabel="Y", zlabel="Z",
        azimuth=9 / 8 * π, zoommode=:cursor, viewmode=:fit,
        xautolimitmargin=(limitmargin, limitmargin),
        yautolimitmargin=(limitmargin, limitmargin),
        zautolimitmargin=(limitmargin, limitmargin),
    )

    # Initialize global storage if needed
    if isnothing(PANEL_MESH_OBSERVABLES[])
        PANEL_MESH_OBSERVABLES[] = Dict()
    end

    body_id = objectid(body_aero)

    # Create observables for each panel using stable (body_id, panel_idx) key
    for (panel_idx, panel) in enumerate(body_aero.panels)
        # Compute initial points
        points = [Point3f(panel.corner_points[:, i]) for i in 1:4]
        if !isnothing(R_b_w) && !isnothing(T_b_w)
            points = [Point3f(R_b_w * p + T_b_w) for p in points]
        end

        # Create observables
        vertices_obs = Observable(points)
        faces_obs = Observable([Makie.GLTriangleFace(1, 2, 3), Makie.GLTriangleFace(1, 3, 4)])
        border_obs = Observable([points..., points[1]])

        # Plot using observables
        mesh!(ax, vertices_obs, faces_obs; color, transparency=true, kwargs...)
        lines!(ax, border_obs; color=:black, transparency=true, kwargs...)

        # Store observables with stable key
        PANEL_MESH_OBSERVABLES[][(body_id, panel_idx)] = (
            vertices=vertices_obs,
            border=border_obs,
            faces=faces_obs
        )
    end

    return fig
end

"""
    save_plot(fig, save_path, title; data_type=".png")

Save a Makie figure to a file.

# Arguments
- `fig`: Makie Figure object
- `save_path`: Path to save the plot
- `title`: Title of the plot

# Keyword arguments
- `data_type`: File extension (default: ".png", also supports ".jpeg")
"""
function VortexStepMethod.save_plot(fig, save_path, title; data_type=".png")
    isnothing(save_path) && throw(ArgumentError("save_path should be provided"))

    !isdir(save_path) && mkpath(save_path)
    full_path = joinpath(save_path, title * data_type)

    @debug "Attempting to save figure to: $full_path"
    @debug "Current working directory: $(pwd())"

    try
        save(full_path, fig)
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
    show_plot(fig; dpi=130)

Display a Makie figure.

# Arguments
- `fig`: Makie Figure object

# Keyword arguments
- `dpi`: Dots per inch for the figure (default: 130) - currently unused in Makie
"""
function VortexStepMethod.show_plot(fig; dpi=130)
    display(fig)
end

"""
    plot_line_segment_makie!(ax, segment, color, label; width=3)

Plot a line segment in 3D with arrow using Makie.

# Arguments
- `ax`: Makie Axis3
- `segment`: Array of two points defining the segment
- `color`: Color of the segment
- `label`: Label for the legend

# Keyword Arguments
- `width`: Line width (default: 3)
"""
function plot_line_segment_makie!(ax, segment, color, label; width=3)
    # Plot line
    lines!(ax, [Point3f(segment[1]), Point3f(segment[2])];
        color=color, linewidth=width, label=label)

    # Plot arrow
    dir = segment[2] - segment[1]
    arrows3d!(ax, [Point3f(segment[1])], [Point3f(dir)];
        color=color, shaftradius=0.01, tipradius=0.03, tiplength=0.1)
end

"""
    set_axes_equal_makie!(ax, panels; zoom=1.8)

Set 3D Makie axis to equal scale based on panel data.

# Arguments
- `ax`: Makie Axis3
- `panels`: Array of panels
- `zoom`: zoom factor (default: 1.8)
"""
function set_axes_equal_makie!(ax, panels; zoom=1.8)
    # Calculate bounds from all panels
    all_x = Float64[]
    all_y = Float64[]
    all_z = Float64[]

    for panel in panels
        for i in 1:4
            push!(all_x, panel.corner_points[1, i])
            push!(all_y, panel.corner_points[2, i])
            push!(all_z, panel.corner_points[3, i])
        end
    end

    x_range = (maximum(all_x) - minimum(all_x)) / zoom
    y_range = (maximum(all_y) - minimum(all_y)) / zoom
    z_range = (maximum(all_z) - minimum(all_z)) / zoom

    max_range = max(x_range, y_range, z_range)

    x_mid = mean([maximum(all_x), minimum(all_x)])
    y_mid = mean([maximum(all_y), minimum(all_y)])
    z_mid = mean([maximum(all_z), minimum(all_z)])

    limits!(ax,
        x_mid - max_range / 2, x_mid + max_range / 2,
        y_mid - max_range / 2, y_mid + max_range / 2,
        z_mid - max_range / 2, z_mid + max_range / 2)
end

"""
    create_geometry_plot_makie(body_aero::BodyAerodynamics, title,
                               view_elevation, view_azimuth; zoom=1.8)

Create a 3D Makie plot of wing geometry including panels and filaments.

# Arguments
- `body_aero`: struct of type BodyAerodynamics
- `title`: plot title
- `view_elevation`: initial view elevation angle [°]
- `view_azimuth`: initial view azimuth angle [°]

# Keyword arguments
- `zoom`: zoom factor (default: 1.8)
"""
function create_geometry_plot_makie(body_aero::BodyAerodynamics, title,
    view_elevation, view_azimuth; zoom=0.5)
    panels = body_aero.panels
    va = isa(body_aero.va, Tuple) ? body_aero.va[1] : body_aero.va

    # Create figure
    fig = Figure(size=(1400, 1400))
    ax = Axis3(fig[1, 1];
        title=title,
        xlabel="x", ylabel="y", zlabel="z",
        aspect=:data,
        azimuth=deg2rad(view_azimuth),
        elevation=deg2rad(view_elevation))

    # Plot panels
    legend_used = Dict{String,Bool}()
    for (i, panel) in enumerate(panels)
        # Panel edges
        corners = [Point3f(panel.corner_points[:, j]) for j in 1:4]
        push!(corners, corners[1])
        lines!(ax, corners; color=:grey, linewidth=1,
            label=i == 1 ? "Panel Edges" : nothing)

        # Control points
        scatter!(ax, [Point3f(panel.control_point)];
            color=:green, markersize=10,
            label=i == 1 ? "Control Points" : nothing)

        # Aerodynamic centers
        scatter!(ax, [Point3f(panel.aero_center)];
            color=:blue, markersize=10,
            label=i == 1 ? "Aerodynamic Centers" : nothing)

        # Plot filaments
        filaments = calculate_filaments_for_plotting(panel)
        legends = ["Bound Vortex", "side1", "side2", "wake_1", "wake_2"]

        for (filament, legend) in zip(filaments, legends)
            x1, x2, color = filament
            show_legend = !get(legend_used, legend, false)
            plot_line_segment_makie!(ax, [x1, x2], color,
                show_legend ? legend : nothing)
            legend_used[legend] = true
        end
    end

    # Plot velocity vector
    max_chord = maximum(panel.chord for panel in panels)
    va_mag = norm(va)
    va_vector_begin = -2 * max_chord * va / va_mag
    va_vector_end = va_vector_begin + 1.5 * va / va_mag
    plot_line_segment_makie!(ax, [va_vector_begin, va_vector_end], :lightblue, "va")

    # Set equal axes
    set_axes_equal_makie!(ax, panels; zoom)

    # Add legend
    axislegend(ax; position=:lt)

    return fig
end

"""
    plot_geometry(body_aero::BodyAerodynamics, title;
                  data_type=".png", save_path=nothing,
                  is_save=false, is_show=false,
                  view_elevation=15, view_azimuth=-120, use_tex=false)

Plot wing geometry from different viewpoints using Makie.

# Arguments:
- `body_aero`: the BodyAerodynamics to plot
- `title`: plot title

# Keyword arguments:
- `data_type`: File extension (default: ".png", also supports ".jpeg")
- `save_path`: Path for saving (default: nothing)
- `is_save`: Whether to save (default: false)
- `is_show`: Whether to display (default: false)
- `view_elevation`: View elevation angle [°] (default: 15)
- `view_azimuth`: View azimuth angle [°] (default: -120)
- `use_tex`: Ignored for Makie (default: false)
"""
function VortexStepMethod.plot_geometry(body_aero::BodyAerodynamics, title;
    data_type=".png",
    save_path=nothing,
    is_save=false,
    is_show=false,
    view_elevation=15,
    view_azimuth=-120,
    use_tex=false)

    if is_save
        # Angled view
        fig = create_geometry_plot_makie(body_aero, "$(title)_angled_view", 15, -120)
        save_plot(fig, save_path, "$(title)_angled_view", data_type=data_type)

        # Top view
        fig = create_geometry_plot_makie(body_aero, "$(title)_top_view", 90, 0)
        save_plot(fig, save_path, "$(title)_top_view", data_type=data_type)

        # Front view
        fig = create_geometry_plot_makie(body_aero, "$(title)_front_view", 0, 0)
        save_plot(fig, save_path, "$(title)_front_view", data_type=data_type)

        # Side view
        fig = create_geometry_plot_makie(body_aero, "$(title)_side_view", 0, -90)
        save_plot(fig, save_path, "$(title)_side_view", data_type=data_type)
    end

    fig = create_geometry_plot_makie(body_aero, title, view_elevation, view_azimuth)

    if is_show
        display(fig)
    end

    return fig
end

"""
    plot_distribution(y_coordinates_list, results_list, label_list;
                      title="spanwise_distribution", data_type=".png",
                      save_path=nothing, is_save=false, is_show=true, use_tex=false)

Plot spanwise distributions of aerodynamic properties using Makie.

# Arguments
- `y_coordinates_list`: List of spanwise coordinates
- `results_list`: List of result dictionaries
- `label_list`: List of labels for different results

# Keyword arguments
- `title`: Plot title (default: "spanwise_distribution")
- `data_type`: File extension (default: ".png", also supports ".jpeg")
- `save_path`: Path to save plots (default: nothing)
- `is_save`: Whether to save (default: false)
- `is_show`: Whether to display (default: true)
- `use_tex`: Ignored for Makie (default: false)
"""
function VortexStepMethod.plot_distribution(y_coordinates_list, results_list, label_list;
    title="spanwise_distribution",
    data_type=".png",
    save_path=nothing,
    is_save=false,
    is_show=true,
    use_tex=false)

    length(results_list) == length(label_list) || throw(ArgumentError(
        "Number of results ($(length(results_list))) must match labels ($(length(label_list)))"
    ))

    # Create figure with 3x3 grid
    fig = Figure(size=(1600, 1000))
    Label(fig[0, :], title, fontsize=20)

    # Row 1: CL, CD, Gamma
    ax_cl = Axis(fig[1, 1], title="CL Distribution",
        xlabel="Spanwise Position y/b", ylabel="Lift Coefficient CL")
    ax_cd = Axis(fig[1, 2], title="CD Distribution",
        xlabel="Spanwise Position y/b", ylabel="Drag Coefficient CD")
    ax_gamma = Axis(fig[1, 3], title="Γ Distribution",
        xlabel="Spanwise Position y/b", ylabel="Circulation Γ")

    # Row 2: Alpha geometric, alpha at ac, alpha uncorrected
    ax_alpha_geo = Axis(fig[2, 1], title="α Geometric",
        xlabel="Spanwise Position y/b", ylabel="Angle of Attack α (deg)")
    ax_alpha_ac = Axis(fig[2, 2], title="α result (corrected to aerodynamic center)",
        xlabel="Spanwise Position y/b", ylabel="Angle of Attack α (deg)")
    ax_alpha_unc = Axis(fig[2, 3], title="α Uncorrected (if VSM, at control point)",
        xlabel="Spanwise Position y/b", ylabel="Angle of Attack α (deg)")

    # Row 3: Force components
    ax_fx = Axis(fig[3, 1], title="Force in x direction",
        xlabel="Spanwise Position y/b", ylabel="Fx")
    ax_fy = Axis(fig[3, 2], title="Force in y direction",
        xlabel="Spanwise Position y/b", ylabel="Fy")
    ax_fz = Axis(fig[3, 3], title="Force in z direction",
        xlabel="Spanwise Position y/b", ylabel="Fz")

    # Plot CL
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        value = round(results["cl"], digits=2)
        lines!(ax_cl, Vector(y_coords), Vector(results["cl_distribution"]),
            label="$label CL: $value")
    end
    axislegend(ax_cl, position=:lt)

    # Plot CD
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        value = round(results["cd"], digits=2)
        lines!(ax_cd, Vector(y_coords), Vector(results["cd_distribution"]),
            label="$label CD: $value")
    end
    axislegend(ax_cd, position=:lt)

    # Plot Gamma
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_gamma, Vector(y_coords), Vector(results["gamma_distribution"]),
            label=label)
    end
    axislegend(ax_gamma, position=:lt)

    # Plot alpha geometric
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_alpha_geo, Vector(y_coords), Vector(results["alpha_geometric"]),
            label=label)
    end
    axislegend(ax_alpha_geo, position=:lt)

    # Plot alpha at ac
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_alpha_ac, Vector(y_coords), Vector(results["alpha_at_ac"]),
            label=label)
    end
    axislegend(ax_alpha_ac, position=:lt)

    # Plot alpha uncorrected
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_alpha_unc, Vector(y_coords), Vector(results["alpha_uncorrected"]),
            label=label)
    end
    axislegend(ax_alpha_unc, position=:lt)

    # Plot force components
    force_axes = [ax_fx, ax_fy, ax_fz]
    components = ["x", "y", "z"]
    for (idx, (ax, comp)) in enumerate(zip(force_axes, components))
        for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
            forces = results["F_distribution"][idx, :]
            if length(y_coords) != length(forces)
                @warn "Dimension mismatch" length(y_coords) length(forces) comp
                continue
            end
            total_force = round(results["F$comp"], digits=2)
            lines!(ax, Vector(y_coords), Vector(forces),
                label="$label ΣF$comp: $total_force N")
        end
        axislegend(ax, position=:lt)
    end

    # Save and show
    if is_save
        save_plot(fig, save_path, title, data_type=data_type)
    end

    if is_show
        display(fig)
    end

    return fig
end

"""
    generate_polar_data(solver, body_aero::BodyAerodynamics, angle_range;
                        angle_type="angle_of_attack", angle_of_attack=0.0,
                        side_slip=0.0, v_a=10.0)

Generate polar data for aerodynamic analysis over a range of angles.

# Arguments
- `solver`: Aerodynamic solver object
- `body_aero`: Wing aerodynamics struct
- `angle_range`: Range of angles to analyze

# Keyword arguments
- `angle_type`: Type of angle variation ("angle_of_attack" or "side_slip")
- `angle_of_attack`: Initial angle of attack [rad]
- `side_slip`: Initial side slip angle [rad]
- `v_a`: norm of apparent wind speed [m/s]

# Returns
- Tuple of polar data array and Reynolds number
"""
function generate_polar_data_makie(
    solver,
    body_aero::BodyAerodynamics,
    angle_range;
    angle_type="angle_of_attack",
    angle_of_attack=0.0,
    side_slip=0.0,
    v_a=10.0
)
    n_panels = length(body_aero.panels)
    n_angles = length(angle_range)

    # Initialize arrays
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    gamma_distribution = zeros(n_angles, n_panels)
    reynolds_number = zeros(n_angles)

    for (i, angle_i) in enumerate(angle_range)
        # Set angle based on type
        if angle_type == "angle_of_attack"
            α = deg2rad(angle_i)
            β = side_slip
        elseif angle_type == "side_slip"
            α = angle_of_attack
            β = deg2rad(angle_i)
        else
            throw(ArgumentError("angle_type must be 'angle_of_attack' or 'side_slip'"))
        end

        # Update inflow conditions
        set_va!(
            body_aero,
            [
                cos(α) * cos(β),
                sin(β),
                sin(α)
            ] * v_a
        )

        # Solve and store results
        results = solve(solver, body_aero, gamma_distribution[i, :])

        cl[i] = results["cl"]
        cd[i] = results["cd"]
        cs[i] = results["cs"]
        gamma_distribution[i, :] = results["gamma_distribution"]
        reynolds_number[i] = results["Rey"]
    end

    polar_data = [angle_range, cl, cd, cs]
    return polar_data, reynolds_number[1]
end

"""
    compute_polar_with_cmy(solver, body_aero, angle_range; angle_type=\"angle_of_attack\",
                           angle_of_attack=0.0, side_slip=0.0, v_a=10.0)

Compute CL/CD/CS polars and optional CMy for a sweep, reusing the previous gamma
as initial guess for faster convergence. Returns a named tuple with angle,
cl, cd, cs, cmy (NaN when unavailable) and per-sample Reynolds numbers.
"""
function compute_polar_with_cmy(
    solver,
    body_aero,
    angle_range;
    angle_type::String="angle_of_attack",
    angle_of_attack::Float64=0.0,
    side_slip::Float64=0.0,
    v_a::Float64=10.0
)
    n_angles = length(angle_range)
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    cmy = fill(NaN, n_angles)  # moment coefficient about body y (pitch)
    reynolds_number = zeros(n_angles)

    gamma_prev = nothing
    for (i, angle_i) in enumerate(angle_range)
        if angle_type == "angle_of_attack"
            α = deg2rad(angle_i)
            β = deg2rad(side_slip)
        elseif angle_type == "side_slip"
            α = deg2rad(angle_of_attack)
            β = deg2rad(angle_i)
        else
            throw(ArgumentError("angle_type must be 'angle_of_attack' or 'side_slip'"))
        end

        set_va!(body_aero, [cos(α) * cos(β), sin(β), sin(α)] * v_a)
        results = solve(solver, body_aero, gamma_prev)

        cl[i] = results["cl"]
        cd[i] = results["cd"]
        cs[i] = results["cs"]
        cmy[i] = get(results, "cmy", NaN)
        reynolds_number[i] = results["Rey"]
        gamma_prev = results["gamma_distribution"]
    end

    return (
        angle=collect(angle_range),
        cl=cl,
        cd=cd,
        cs=cs,
        cmy=cmy,
        rey=reynolds_number,
    )
end

"""
    plot_polars(solver_list, body_aero_list, label_list;
                literature_path_list=String[],
                angle_range=range(0, 20, 2), angle_type="angle_of_attack",
                angle_of_attack=0.0, side_slip=0.0, v_a=10.0,
                title="polar", data_type=".png", save_path=nothing,
                is_save=true, is_show=true, use_tex=false)

Plot polar data comparing different solvers using Makie.

# Arguments
- `solver_list`: List of aerodynamic solvers
- `body_aero_list`: List of wing aerodynamics objects
- `label_list`: List of labels for each configuration

# Keyword arguments
- `literature_path_list`: Optional paths to literature data files
- `angle_range`: Range of angles [°]
- `angle_type`: "angle_of_attack" or "side_slip" (default: angle_of_attack)
- `angle_of_attack`: AoA [rad] (default: 0.0)
- `side_slip`: Side slip angle [rad] (default: 0.0)
- `v_a`: Wind speed [m/s] (default: 10.0)
- `title`: Plot title
- `data_type`: File extension (default: ".png", also supports ".jpeg")
- `save_path`: Path to save (default: nothing)
- `is_save`: Whether to save (default: true)
- `is_show`: Whether to display (default: true)
- `use_tex`: Ignored for Makie (default: false)
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
    data_type=".png",
    save_path=nothing,
    is_save=true,
    is_show=true,
    use_tex=false
)
    # Validate inputs
    total_cases = length(body_aero_list) + length(literature_path_list)
    if total_cases != length(label_list) || length(solver_list) != length(body_aero_list)
        throw(ArgumentError("Mismatch in solvers/cases/labels"))
    end
    main_title = replace(title, " " => "_")

    # Generate polar data
    polar_data_list = []
    labels_with_re = copy(label_list)
    for (i, (solver, body_aero)) in enumerate(zip(solver_list, body_aero_list))
        polar_data, rey = generate_polar_data_makie(
            solver, body_aero, angle_range;
            angle_type, angle_of_attack, side_slip, v_a
        )
        push!(polar_data_list, polar_data)
        labels_with_re[i] = "$(label_list[i]) Re = $(round(Int64, rey*1e-5))e5"
    end

    # Load literature data if provided
    if !isempty(literature_path_list)
        for path in literature_path_list
            data = readdlm(path, ',')
            header = lowercase.(string.(data[1, :]))
            alpha_idx = findfirst(x -> occursin("alpha", x), header)
            cl_idx = findfirst(x -> occursin("cl", x), header)
            cd_idx = findfirst(x -> occursin("cd", x), header)
            cs_idx = findfirst(x -> occursin("cs", x), header)
            cs_col = cs_idx === nothing ? zeros(size(data, 1) - 1) : data[2:end, cs_idx]
            push!(polar_data_list, [
                data[2:end, alpha_idx],
                data[2:end, cl_idx],
                data[2:end, cd_idx],
                cs_col
            ])
        end
    end

    # Create figure with 2x2 grid
    fig = Figure(size=(1400, 1400))

    ax_cl = Axis(fig[1, 1], title="CL vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CL")
    ax_cd = Axis(fig[1, 2], title="CD vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CD")
    ax_cs = Axis(fig[2, 1], title="CS vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CS")
    ax_polar = Axis(fig[2, 2], title="CL vs CD",
        xlabel="CD", ylabel="CL")

    # Number of computational results
    n_solvers = length(solver_list)

    # Plot CL vs angle
    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, labels_with_re))
        marker = i <= n_solvers ? :star5 : :circle
        markersize = i <= n_solvers ? 12 : 8
        scatterlines!(ax_cl, polar_data[1], polar_data[2];
            label=label, marker=marker, markersize=markersize)
        if maximum(polar_data[2]) > 10
            ylims!(ax_cl, -0.5, 2)
        end
    end
    axislegend(ax_cl, position=:lt)

    # Plot CD vs angle
    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, labels_with_re))
        marker = i <= n_solvers ? :star5 : :circle
        markersize = i <= n_solvers ? 12 : 8
        scatterlines!(ax_cd, polar_data[1], polar_data[3];
            label=label, marker=marker, markersize=markersize)
        if maximum(polar_data[2]) > 10
            ylims!(ax_cd, -0.5, 2)
        end
    end
    axislegend(ax_cd, position=:lt)

    # Plot CS vs angle
    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, labels_with_re))
        marker = i <= n_solvers ? :star5 : :circle
        markersize = i <= n_solvers ? 12 : 8
        scatterlines!(ax_cs, polar_data[1], polar_data[4];
            label=label, marker=marker, markersize=markersize)
        if maximum(polar_data[2]) > 10
            ylims!(ax_cs, -0.5, 2)
        end
    end
    axislegend(ax_cs, position=:lt)

    # Plot CL vs CD
    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, labels_with_re))
        marker = i <= n_solvers ? :star5 : :circle
        markersize = i <= n_solvers ? 12 : 8
        scatterlines!(ax_polar, polar_data[3], polar_data[2];
            label=label, marker=marker, markersize=markersize)
        if maximum(polar_data[2]) > 10 || maximum(polar_data[3]) > 10
            ylims!(ax_polar, -0.5, 2)
            xlims!(ax_polar, -0.5, 2)
        end
    end
    axislegend(ax_polar, position=:lt)

    # Save and show
    if is_save && !isnothing(save_path)
        save_plot(fig, save_path, main_title; data_type)
    end

    if is_show
        display(fig)
    end

    return fig
end

"""
    plot_polar_data(body_aero::BodyAerodynamics;
                    alphas=collect(deg2rad.(-5:0.3:25)),
                    delta_tes=collect(deg2rad.(-5:0.3:25)),
                    is_show=true, use_tex=false)

Plot polar data (Cl, Cd, Cm) as 3D surfaces using Makie.

# Arguments
- `body_aero`: Wing aerodynamics struct

# Keyword arguments
- `alphas`: Range of AoA values [rad] (default: -5° to 25° in 0.3° steps)
- `delta_tes`: Range of trailing edge angles [rad] (default: -5° to 25° in 0.3° steps)
- `is_show`: Whether to display (default: true)
- `use_tex`: Ignored for Makie (default: false)
"""
function VortexStepMethod.plot_polar_data(body_aero::BodyAerodynamics;
    alphas=collect(deg2rad.(-5:0.3:25)),
    delta_tes=collect(deg2rad.(-5:0.3:25)),
    is_show=true,
    use_tex=false)

    if body_aero.panels[1].aero_model == POLAR_MATRICES
        # Create figure with 3 subplots
        fig = Figure(size=(1500, 600))

        # Get interpolation functions
        interp_data = [
            (body_aero.panels[1].cl_interp, "Cl"),
            (body_aero.panels[1].cd_interp, "Cd"),
            (body_aero.panels[1].cm_interp, "Cm")
        ]

        # Create each subplot
        for (idx, (interp, label)) in enumerate(interp_data)
            ax = Axis3(fig[1, idx];
                title="$label vs α and δ",
                xlabel="δ [rad]",
                ylabel="α [rad]",
                zlabel=label,
                azimuth=1.275 * π)

            # Create interpolation matrix
            interp_matrix = [interp(alpha, delta_te)
                             for alpha in alphas, delta_te in delta_tes]

            # Create wireframe
            wireframe!(ax, delta_tes, alphas, interp_matrix;
                color=:blue, linewidth=0.5, transparency=true)
        end

        if is_show
            display(fig)
        end
        return fig
    else
        throw(ArgumentError(
            "Plotting polar data for $(body_aero.panels[1].aero_model) not implemented."
        ))
    end
end

"""
    plot_combined_analysis(solver, body_aero, results;
                          solver_label="VSM",
                          angle_range=range(0,20,length=20),
                          angle_type="angle_of_attack",
                          angle_of_attack=0.0, side_slip=0.0, v_a=10.0,
                          title="Combined Analysis",
                          view_elevation=15, view_azimuth=-120,
                          is_show=true, use_tex=false,
                          literature_path_list=String[],
                          data_type=".png", save_path=nothing, is_save=false)

Create combined multi-panel figure with geometry, polar data, distributions, and polars.

# Arguments
- `solver`: Aerodynamic solver
- `body_aero`: BodyAerodynamics object
- `results`: Solution dictionary from solve()

# Keyword arguments
- `labels`: Optional label string or label vector. If a vector with length
  `length(solver)+length(literature_path_list)` is provided, the tail is used
  for literature labels; otherwise literature labels default to file basenames.
- `solver_label`: Backward-compatible alias for `labels`.
- `angle_range`: Range of angles for polars (default: range(0,20,length=20))
- `angle_type`: "angle_of_attack" or "side_slip" (default: "angle_of_attack")
- `angle_of_attack`: AoA in degrees (default: 0.0)
- `side_slip`: Side slip in degrees (default: 0.0)
- `v_a`: Wind speed in m/s (default: 10.0)
- `title`: Overall figure title (default: "Combined Analysis")
- `view_elevation`: Geometry view elevation [°] (default: 15)
- `view_azimuth`: Geometry view azimuth [°] (default: -120)
- `is_show`: Display figure (default: true)
- `use_tex`: Ignored for Makie (default: false)
- `literature_path_list`: Paths to literature CSV files (default: String[])
- `data_type`: File extension (default: ".png", also supports ".jpeg")
- `save_path`: Directory path to save files (default: nothing)
- `is_save`: Save plots to files (default: false)
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
    view_elevation=15,
    view_azimuth=-120,
    is_show=true,
    use_tex=false,
    literature_path_list=String[],
    data_type=".png",
    save_path=nothing,
    is_save=false,
    angle_of_attack_for_spanwise_distribution=5.0,
)
    # Normalize inputs to arrays for consistent handling
    solvers = solver isa Vector ? solver : [solver]
    body_aeros = body_aero isa Vector ? body_aero : [body_aero]
    results_list = results isa Vector ? results : [results]
    n_solvers = length(solvers)
    n_literature = length(literature_path_list)

    length(body_aeros) == n_solvers || throw(ArgumentError(
        "Length mismatch: got $(length(body_aeros)) body_aero entries for $n_solvers solvers"
    ))
    length(results_list) == n_solvers || throw(ArgumentError(
        "Length mismatch: got $(length(results_list)) results entries for $n_solvers solvers"
    ))

    label_source = isnothing(labels) ? solver_label : labels
    labels_in = label_source isa AbstractVector ? string.(label_source) : [string(label_source)]
    valid_label_lengths = (1, n_solvers, n_solvers + n_literature)
    length(labels_in) in valid_label_lengths || throw(ArgumentError(
        "labels (or solver_label) must be a string or a vector of length 1, $n_solvers, or $(n_solvers + n_literature); got $(length(labels_in))"
    ))

    solver_labels = length(labels_in) == 1 ? fill(labels_in[1], n_solvers) : labels_in[1:n_solvers]
    literature_labels = if n_literature == 0
        String[]
    elseif length(labels_in) == n_solvers + n_literature
        labels_in[(n_solvers+1):end]
    else
        [splitext(basename(path))[1] for path in literature_path_list]
    end

    # Auto-detect screen size and use 80% of it
    fig = try
        screen_size = Makie.primary_resolution()
        fig_width = round(Int, screen_size[1] * 0.8)
        fig_height = round(Int, screen_size[2] * 0.8)
        Figure(size=(fig_width, fig_height))
    catch
        # Fallback if screen detection fails
        Figure(size=(1800, 1200))
    end
    Label(fig[0, :], title, fontsize=20, font=:bold)

    # Use first body_aero for geometry and polar data display
    first_body = body_aeros[1]
    panels = first_body.panels
    va = isa(first_body.va, Tuple) ? first_body.va[1] : first_body.va

    # Compute spanwise results for each solver
    results_spanwise_list = copy(results_list)
    if !isnothing(angle_of_attack_for_spanwise_distribution)
        α_span = deg2rad(angle_of_attack_for_spanwise_distribution)
        β_span = deg2rad(side_slip)
        for (i, (s, ba)) in enumerate(zip(solvers, body_aeros))
            va_old = copy(ba.va)
            set_va!(ba, [cos(α_span) * cos(β_span), sin(β_span),
                sin(α_span)] * v_a)
            results_spanwise_list[i] = solve(s, ba,
                s.sol.gamma_distribution)
            set_va!(ba, va_old)
        end
    end

    # [1,1] Wing Geometry
    ax_geo = Axis3(fig[1, 1];
        title="Wing Geometry",
        xlabel="x", ylabel="y", zlabel="z",
        aspect=:data,
        azimuth=deg2rad(view_azimuth),
        elevation=deg2rad(view_elevation))

    legend_used = Dict{String,Bool}()
    for (i, panel) in enumerate(panels)
        corners = [Point3f(panel.corner_points[:, j]) for j in 1:4]
        push!(corners, corners[1])
        lines!(ax_geo, corners; color=:grey, linewidth=1,
            label=i == 1 ? "Panel Edges" : nothing)

        scatter!(ax_geo, [Point3f(panel.control_point)];
            color=:green, markersize=10,
            label=i == 1 ? "Control Points" : nothing)

        scatter!(ax_geo, [Point3f(panel.aero_center)];
            color=:blue, markersize=10,
            label=i == 1 ? "Aerodynamic Centers" : nothing)

        filaments = calculate_filaments_for_plotting(panel)
        legends = ["Bound Vortex", "side1", "side2", "wake_1", "wake_2"]

        for (filament, legend) in zip(filaments, legends)
            x1, x2, color = filament
            show_legend = !get(legend_used, legend, false)
            plot_line_segment_makie!(ax_geo, [x1, x2], color,
                show_legend ? legend : nothing)
            legend_used[legend] = true
        end
    end

    max_chord = maximum(panel.chord for panel in panels)
    va_mag = norm(va)
    va_vector_begin = -2 * max_chord * va / va_mag
    va_vector_end = va_vector_begin + 1.5 * va / va_mag
    plot_line_segment_makie!(ax_geo, [va_vector_begin, va_vector_end],
        :lightblue, "va")

    set_axes_equal_makie!(ax_geo, panels; zoom=0.5)
    axislegend(ax_geo; position=:lt)

    # [1,2] Polar Data Surfaces or Curves
    if first_body.panels[1].aero_model == POLAR_MATRICES
        alphas = collect(deg2rad.(-5:0.3:20))
        delta_tes = collect(deg2rad.(-5:0.3:20))

        interp_data = [
            (first_body.panels[1].cl_interp, "Cl"),
            (first_body.panels[1].cd_interp, "Cd"),
            (first_body.panels[1].cm_interp, "Cm")
        ]

        for (idx, (interp, label)) in enumerate(interp_data)
            ax = Axis3(fig[1, 2][1, idx];
                title="$label vs α and δ",
                xlabel="δ [rad]",
                ylabel="α [rad]",
                zlabel=label,
                azimuth=1.275 * π)

            interp_matrix = [interp(alpha, delta_te)
                             for alpha in alphas, delta_te in delta_tes]

            wireframe!(ax, delta_tes, alphas, interp_matrix;
                color=:blue, linewidth=0.5, transparency=true)
        end
    elseif first_body.panels[1].aero_model == POLAR_VECTORS
        alphas_deg = collect(-5:0.5:20)
        alphas = deg2rad.(alphas_deg)

        ax_cl_curve = Axis(fig[1, 2][1, 1];
            title="2D Cl vs α",
            xlabel="α [°]",
            ylabel="Cl")
        ax_cd_curve = Axis(fig[1, 2][1, 2];
            title="2D Cd vs α",
            xlabel="α [°]",
            ylabel="Cd")
        ax_cm_curve = Axis(fig[1, 2][1, 3];
            title="2D Cm vs α",
            xlabel="α [°]",
            ylabel="Cm")

        cl_vals = [first_body.panels[1].cl_interp(a) for a in alphas]
        cd_vals = [first_body.panels[1].cd_interp(a) for a in alphas]
        cm_vals = [first_body.panels[1].cm_interp(a) for a in alphas]

        lines!(ax_cl_curve, alphas_deg, cl_vals; color=:blue, linewidth=2)
        lines!(ax_cd_curve, alphas_deg, cd_vals; color=:red, linewidth=2)
        lines!(ax_cm_curve, alphas_deg, cm_vals; color=:green, linewidth=2)
    end

    # [2,1] Spanwise Distributions (3×3 grid)
    y_coords = [panel.aero_center[2] for panel in first_body.panels]

    aoa_span = round(angle_of_attack_for_spanwise_distribution, digits=1)

    ax_cl = Axis(fig[2, 1][1, 1], title="CL Distribution (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="CL")
    ax_cd = Axis(fig[2, 1][1, 2], title="CD Distribution (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="CD")
    ax_gamma = Axis(fig[2, 1][1, 3], title="Γ Distribution (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="Γ")
    ax_alpha_geo = Axis(fig[2, 1][2, 1], title="α Geometric (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="α (deg)")
    ax_alpha_ac = Axis(fig[2, 1][2, 2], title="α at 1/4c aero center (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="α (deg)")
    ax_alpha_unc = Axis(fig[2, 1][2, 3], title="α at 3/4c uncorrected (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="α (deg)")

    ax_fx = Axis(fig[2, 1][3, 1], title="Force x (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="Fx")
    ax_fy = Axis(fig[2, 1][3, 2], title="Force y(α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="Fy")
    ax_fz = Axis(fig[2, 1][3, 3], title="Force z (α=$aoa_span)",
        xlabel="Spanwise position y[m]", ylabel="Fz")

    colors = Makie.wong_colors()
    for (si, (lbl, rs)) in enumerate(zip(solver_labels, results_spanwise_list))
        color = colors[mod1(si, length(colors))]
        y_si = [panel.aero_center[2] for panel in body_aeros[si].panels]

        cl_val = round(rs["cl"], digits=2)
        lines!(ax_cl, Vector(y_si), Vector(rs["cl_distribution"]);
            label="$lbl CL: $cl_val", color)

        cd_val = round(rs["cd"], digits=2)
        lines!(ax_cd, Vector(y_si), Vector(rs["cd_distribution"]);
            label="$lbl CD: $cd_val", color)

        lines!(ax_gamma, Vector(y_si),
            Vector(rs["gamma_distribution"]); label=lbl, color)

        lines!(ax_alpha_geo, Vector(y_si),
            rad2deg.(Vector(rs["alpha_geometric"])); label=lbl, color)

        lines!(ax_alpha_ac, Vector(y_si),
            rad2deg.(Vector(rs["alpha_at_ac"])); label=lbl, color)

        lines!(ax_alpha_unc, Vector(y_si),
            rad2deg.(Vector(rs["alpha_uncorrected"]));
            label=lbl, color)

        force_axes = [ax_fx, ax_fy, ax_fz]
        components = ["x", "y", "z"]
        for (idx, (ax, comp)) in enumerate(zip(force_axes, components))
            forces = rs["F_distribution"][idx, :]
            total_force = round(rs["F$comp"], digits=2)
            lines!(ax, Vector(y_si), Vector(forces);
                label="$lbl ΣF$comp: $total_force N", color)
        end
    end

    for ax in [ax_cl, ax_cd, ax_gamma, ax_alpha_geo, ax_alpha_ac,
        ax_alpha_unc, ax_fx, ax_fy, ax_fz]
        axislegend(ax, position=:lt)
    end

    # force y-limits
    ylims!(ax_alpha_geo, -0.25 * aoa_span, 1.1 * aoa_span)

    # [2,2] Polars (2×2 grid)
    polar_series = Tuple[]
    for (si, (s, ba, lbl)) in enumerate(zip(solvers, body_aeros, solver_labels))
        polar_solver = compute_polar_with_cmy(
            s, ba, angle_range;
            angle_type=angle_type,
            angle_of_attack=angle_of_attack,
            side_slip=side_slip,
            v_a=v_a
        )
        label_re = "$lbl Re = $(round(Int64,
                     first(polar_solver.rey) * 1e-5))e5"
        push!(polar_series, (polar_solver, label_re))
    end

    # Load literature data (if any)
    if !isempty(literature_path_list)
        for (path, lit_label) in zip(literature_path_list, literature_labels)
            data = readdlm(path, ',')
            header_raw = string.(data[1, :])
            header = lowercase.(strip.(header_raw))
            alpha_idx = findfirst(x -> occursin("alpha", x) || occursin("aoa", x), header)
            cl_idx = findfirst(x -> occursin("cl", x), header)
            cd_idx = findfirst(x -> occursin("cd", x), header)
            cs_idx = findfirst(x -> occursin("cs", x), header)
            cmy_idx = findfirst(x -> occursin("cmy", x), header)

            parse_col(col) = begin
                vals = Float64[]
                for v in col
                    if v isa Real
                        push!(vals, Float64(v))
                    else
                        s = strip(String(v))
                        y = tryparse(Float64, s)
                        push!(vals, isnothing(y) ? NaN : y)
                    end
                end
                vals
            end

            angles = parse_col(data[2:end, alpha_idx])
            cl_vals = parse_col(data[2:end, cl_idx])
            cd_vals = parse_col(data[2:end, cd_idx])
            cs_vals = cs_idx === nothing ? zeros(size(data, 1) - 1) : parse_col(data[2:end, cs_idx])
            cmy_vals = cmy_idx === nothing ? fill(NaN, length(angles)) : parse_col(data[2:end, cmy_idx])

            push!(polar_series, ((angle=angles, cl=cl_vals, cd=cd_vals, cs=cs_vals, cmy=cmy_vals, rey=fill(NaN, length(angles))), lit_label))
        end
    end

    ax_cl_polar = Axis(fig[2, 2][1, 1], title="CL vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CL")
    ax_cd_polar = Axis(fig[2, 2][1, 2], title="CD vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CD")
    ax_cs_polar = Axis(fig[2, 2][2, 1], title="CS vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CS")
    ax_polar = Axis(fig[2, 2][2, 2], title="CL vs CD",
        xlabel="CD", ylabel="CL")

    for (idx, (pd, lbl)) in enumerate(polar_series)
        color = colors[mod1(idx, length(colors))]
        marker = idx == 1 ? :star5 : :circle
        markersize = idx == 1 ? 12 : 8

        scatterlines!(ax_cl_polar, pd.angle, pd.cl; label=lbl, marker, markersize, color)
        scatterlines!(ax_cd_polar, pd.angle, pd.cd; label=lbl, marker, markersize, color)
        scatterlines!(ax_cs_polar, pd.angle, pd.cs; label=lbl, marker, markersize, color)
        scatterlines!(ax_polar, pd.cd, pd.cl; label=lbl, marker, markersize, color)
    end
    # axislegend(ax_cl_polar, position=:lt)
    # axislegend(ax_cd_polar, position=:lt)
    axislegend(ax_cs_polar, position=:lb)
    # axislegend(ax_polar, position=:lt)

    # Set column widths: left column wider for 3x3 grid
    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    if is_show
        display(fig)
    end

    return fig
end

# ============================================================================
# Airfoil Slice and Kulfan Fit Visualization
# ============================================================================

"""
    plot_airfoil_slices(obj_path::String; n_slices=5, is_show=true,
                        save_path=nothing, data_type=".png")

Plot airfoil cross-sections extracted from OBJ mesh with Kulfan CST fits overlaid.

# Arguments
- `obj_path`: Path to OBJ file

# Keyword Arguments
- `n_slices`: Number of spanwise slices (default 5)
- `is_show`: Display figure (default true)
- `save_path`: Directory to save figure (default nothing)
- `data_type`: Image format (default ".png")

# Returns
- Makie Figure object
"""
function VortexStepMethod.plot_airfoil_slices(obj_path::String; n_slices::Int=5,
                                              is_show::Bool=true,
                                              save_path=nothing,
                                              data_type::String=".png")
    # Slice the mesh
    airfoils, y_positions = VortexStepMethod.slice_obj_wing(obj_path, n_slices)

    # Determine grid layout
    n_cols = min(3, n_slices)
    n_rows = ceil(Int, n_slices / n_cols)

    fig = Figure(size=(400 * n_cols, 350 * n_rows))
    Label(fig[0, :], "Airfoil Slices from OBJ with Kulfan CST Fits", fontsize=16)

    for (i, ((x_slice, y_slice), y_pos)) in enumerate(zip(airfoils, y_positions))
        row = div(i - 1, n_cols) + 1
        col = mod1(i, n_cols)

        ax = Axis(fig[row, col];
                  title="Slice $i (y=$(round(y_pos, digits=3)))",
                  xlabel="x/c",
                  ylabel="y/c",
                  aspect=DataAspect())

        if isempty(x_slice) || length(x_slice) < 5
            text!(ax, 0.5, 0.5, text="No data", align=(:center, :center))
            continue
        end

        # Get normalized coordinates for visualization
        x_norm, y_norm, _ = VortexStepMethod.normalize_airfoil(x_slice, y_slice)

        # Plot slice points
        scatter!(ax, x_norm, y_norm;
                color=:blue, markersize=4,
                label="$(length(x_norm)) pts")

        # Try Kulfan fit
        try
            params = VortexStepMethod.fit_kulfan_parameters(x_slice, y_slice)

            # Check if fit is reasonable (TE thickness should be small)
            if abs(params.TE_thickness) < 1.0
                # Generate fitted curve
                x_fit, y_fit = VortexStepMethod.kulfan_to_coordinates(params; n_points=100)

                # Plot fitted curve
                lines!(ax, x_fit, y_fit;
                      color=:red, linewidth=2, label="Kulfan fit")

                # Add TE thickness annotation
                text!(ax, 0.02, 0.98;
                     text="TE=$(round(params.TE_thickness, digits=4))",
                     align=(:left, :top),
                     fontsize=10,
                     space=:relative)
            end
        catch e
            # Silently skip fit - just show points
        end

        # Add legend only to first plot
        if i == 1
            axislegend(ax; position=:rt)
        end
    end

    if !isnothing(save_path)
        VortexStepMethod.save_plot(fig, save_path, "airfoil_slices"; data_type)
    end

    if is_show
        display(fig)
    end

    return fig
end

"""
    plot_airfoil_fit(x::Vector, y::Vector; title="Airfoil Fit", is_show=true)

Plot a single airfoil with its Kulfan CST fit.

# Arguments
- `x, y`: Airfoil coordinates

# Keyword Arguments
- `title`: Plot title
- `is_show`: Display figure

# Returns
- Makie Figure object
- KulfanParameters
"""
function plot_airfoil_fit(x::Vector, y::Vector; title::String="Airfoil Fit",
                          is_show::Bool=true)
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1, 1];
              title=title,
              xlabel="x/c",
              ylabel="y/c",
              aspect=DataAspect())

    # Normalize coordinates
    x_norm, y_norm, _ = VortexStepMethod.normalize_airfoil(x, y)

    # Plot original points
    scatter!(ax, x_norm, y_norm;
            color=:blue, markersize=6, label="Input ($(length(x)) pts)")

    # Fit Kulfan parameters
    params = VortexStepMethod.fit_kulfan_parameters(x, y)

    # Generate fitted curve
    x_fit, y_fit = VortexStepMethod.kulfan_to_coordinates(params; n_points=150)

    # Plot fitted curve
    lines!(ax, x_fit, y_fit;
          color=:red, linewidth=2, linestyle=:dash, label="Kulfan fit")

    # Add parameter info
    info_text = """
    Upper: $(round.(params.upper_weights[1:4], digits=3))...
    Lower: $(round.(params.lower_weights[1:4], digits=3))...
    LE: $(round(params.leading_edge_weight, digits=3))
    TE: $(round(params.TE_thickness, digits=5))
    """

    ax_info = Axis(fig[1, 2]; title="Kulfan Parameters")
    hidedecorations!(ax_info)
    hidespines!(ax_info)
    text!(ax_info, 0.1, 0.9;
         text=info_text,
         align=(:left, :top),
         fontsize=12)

    axislegend(ax; position=:rt)

    if is_show
        display(fig)
    end

    return fig, params
end

end

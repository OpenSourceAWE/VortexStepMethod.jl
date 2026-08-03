module VortexStepMethodMakieExt
using MakieControlPlots.Makie, VortexStepMethod, LinearAlgebra, Statistics, DelimitedFiles
import MakieControlPlots
import VortexStepMethod: calculate_filaments_for_plotting
import VortexStepMethod: ObjAdapter, AirfoilAero

export plot_geometry, plot_distribution, plot_polars, save_plot, show_plot,
    plot_polar_data, plot_combined_analysis,
    plot_section_polars

# Global storage for panel mesh observables (for dynamic plotting)
const PANEL_MESH_OBSERVABLES = Ref{Union{Nothing,Dict}}(nothing)
# Global storage for airfoil-skin observables, keyed by body objectid.
const AIRFOIL_SKIN_OBSERVABLES = Ref{Union{Nothing,Dict}}(nothing)

const PLATE_FACES = [Makie.GLTriangleFace(1, 2, 5), Makie.GLTriangleFace(1, 5, 6),
    Makie.GLTriangleFace(2, 3, 4), Makie.GLTriangleFace(2, 4, 5)]
const PLATE_BORDER_IDX = [1, 2, 3, 4, 5, 6, 1]

"""
    panel_normal(panel) -> Point3f

Body-frame airfoil-upper-surface unit normal of a panel, from its `corner_points`
and sign-aligned to `z_airf` (both body-frame, so the sign relation survives the
kite's attitude changes). Falls back to `+z` for a degenerate quad.
"""
function panel_normal(panel)
    c = panel.corner_points
    up = cross(Point3f(c[:, 2]) - Point3f(c[:, 1]), Point3f(c[:, 1]) - Point3f(c[:, 4]))
    n = norm(up)
    up = n < 1e-9 ? Point3f(0, 0, 1) : Point3f(up / n)
    dot(up, Point3f(panel.z_airf)) < 0 && (up = -up)
    return up
end

"""
    plate_hinge_local(crease_frac, delta) -> (lx, ly)

Chordwise/thickness fractions (chord = 1, LE at 0) of the flap hinge for a plate
whose trailing edge is deflected by `delta` [rad] about `crease_frac` and then
re-pinned so the TE stays at the panel TE corner. A downward deflection therefore
lifts the hinge (`ly > 0`, the "bulge up"). `crease_frac` outside `(0, 1)` disables
the kink (hinge stays on the chord line).
"""
function plate_hinge_local(crease_frac, delta)
    (crease_frac <= 0 || crease_frac >= 1) && return (Float64(crease_frac), 0.0)
    tex = crease_frac + (1 - crease_frac) * cos(delta)
    tey = -(1 - crease_frac) * sin(delta)
    den = tex^2 + tey^2
    den < 1e-12 && return (Float64(crease_frac), 0.0)
    return (crease_frac * tex / den, -crease_frac * tey / den)
end

"""
    panel_plate_geometry(panel; R_b_w=nothing, T_b_w=nothing) -> Vector{Point3f}

The 6 vertices `[LE_1, hinge_1, TE_1, TE_2, hinge_2, LE_2]` of a panel's flat-plate
skin, kinked at [`plate_hinge_local`](@ref) by the panel's `delta`/`crease_frac`
(TE pinned to the corners, hinge bulging up). Triangulated by [`PLATE_FACES`](@ref);
with `delta == 0` the two quads are coplanar (the original flat quad). Transformed to
world by `R_b_w`/`T_b_w`.
"""
function panel_plate_geometry(panel; R_b_w=nothing, T_b_w=nothing)
    to_world(p) = (isnothing(R_b_w) || isnothing(T_b_w)) ? Point3f(p) :
        Point3f(R_b_w * p + T_b_w)
    c = panel.corner_points
    c1 = Point3f(c[:, 1]); c2 = Point3f(c[:, 2]); c3 = Point3f(c[:, 3]); c4 = Point3f(c[:, 4])
    up = panel_normal(panel)
    lx, ly = plate_hinge_local(panel.crease_frac, panel.delta)
    chord1 = c2 - c1; chord2 = c3 - c4
    hinge1 = c1 + lx * chord1 + (ly * norm(chord1)) * up
    hinge2 = c4 + lx * chord2 + (ly * norm(chord2)) * up
    return [to_world(c1), to_world(hinge1), to_world(c2),
            to_world(c3), to_world(hinge2), to_world(c4)]
end

"""
    plot!(ax, panel::VortexStepMethod.Panel; use_observables=false, kwargs...)

Plot a single `Panel` as a flat-plate `mesh`, kinked at the flap hinge by the
panel's `delta`/`crease_frac` (see [`panel_plate_geometry`](@ref)); with `delta == 0`
this is the flat quad LE1-TE1-TE2-LE2.

If `use_observables=true`, creates observables for dynamic updates.
"""
function Makie.plot!(ax, panel::VortexStepMethod.Panel; color=(:red, 0.2), R_b_w=nothing, T_b_w=nothing,
    use_observables=false, border_linewidth=1.5, transparency=true, kwargs...)
    plots = []
    points = panel_plate_geometry(panel; R_b_w, T_b_w)

    if use_observables
        # Create observables for dynamic updates
        vertices_obs = Observable(points)
        faces_obs = Observable(copy(PLATE_FACES))
        border_obs = Observable(points[PLATE_BORDER_IDX])

        p = mesh!(ax, vertices_obs, faces_obs; color, transparency, kwargs...)
        push!(plots, p)
        p = lines!(ax, border_obs; color=:black, linewidth=border_linewidth,
                   transparency, kwargs...)
        push!(plots, p)

        # Note: Observables are stored at the body level, not individual panel level
        # Individual panels need their parent body for proper tracking
    else
        # Static plotting (original behavior)
        p = mesh!(ax, points, PLATE_FACES; color, transparency, kwargs...)
        push!(plots, p)
        p = lines!(ax, points[PLATE_BORDER_IDX]; color=:black, linewidth=border_linewidth,
                   transparency, kwargs...)
        push!(plots, p)
    end

    return plots
end

"""
    airfoil_skin_geometry(body; R_b_w=nothing, T_b_w=nothing) -> (vertices, faces, ribs)

Lofted airfoil skin of a `BodyAerodynamics`: each section's deflected contour
(`section_surface` at the panel's `delta`) is fitted between the panel's `corner_points`
by a 2D similarity, TE pinned so a deflection bulges the fore body up. The skin reflects
`delta` only when the geometry carries per-`delta` slices (`obj_to_yaml` with a
`delta_range`); with δ=0-only data it renders undeflected, a deliberate cue that the
deflected slices are missing. Transformed to world by `R_b_w`/`T_b_w`.
`vertices`/`faces` triangulate the skin between consecutive equal-node sections; `ribs` is
one closed contour polyline per section. Sections without contour data are skipped.
"""
function airfoil_skin_geometry(body; R_b_w=nothing, T_b_w=nothing)
    to_world(p) = (isnothing(R_b_w) || isnothing(T_b_w)) ? Point3f(p) :
        Point3f(R_b_w * p + T_b_w)
    ribs = Vector{Point3f}[]
    panel_offset = 0
    for wing in body.wings
        sections = wing.refined_sections
        n = length(sections)
        n_panels = n - 1
        n_panels < 1 && continue
        spanwise = Point3f(wing.spanwise_direction)
        increasing = dot(Point3f(sections[n].LE_point) -
                         Point3f(sections[1].LE_point), spanwise) > 0
        for (i, section) in enumerate(sections)
            isnothing(section.section_aero) && continue
            panel_idx = panel_offset + min(i, n_panels)
            panel_idx <= length(body.panels) || continue
            panel = body.panels[panel_idx]
            corners = panel.corner_points
            corner1 = Point3f(corners[:, 1]); corner3 = Point3f(corners[:, 3])
            corner2 = Point3f(corners[:, 2]); corner4 = Point3f(corners[:, 4])
            plus_edge = i <= n_panels ? !increasing : increasing
            leading = plus_edge ? corner1 : corner4
            trailing = plus_edge ? corner2 : corner3
            chord = trailing - leading
            chord_len = norm(chord)
            chord_len < 1e-9 && continue
            up = panel_normal(panel)
            xs, ys, _, _ = VortexStepMethod.section_surface(section.section_aero,
                                                            0.0, panel.delta)
            le_i = argmin(xs)
            le_x = xs[le_i]; le_y = ys[le_i]
            te_x = 0.5 * (xs[1] + xs[end]); te_y = 0.5 * (ys[1] + ys[end])
            dx = te_x - le_x; dy = te_y - le_y
            den = dx^2 + dy^2
            den < 1e-18 && continue
            rib = Vector{Point3f}(undef, length(xs))
            for k in eachindex(xs)
                px = xs[k] - le_x; py = ys[k] - le_y
                lx = (px * dx + py * dy) / den
                ly = (py * dx - px * dy) / den
                rib[k] = to_world(leading + lx * chord + (ly * chord_len) * up)
            end
            push!(ribs, rib)
        end
        panel_offset += n_panels
    end
    vertices = Point3f[]
    faces = Makie.GLTriangleFace[]
    for i in 1:length(ribs) - 1
        length(ribs[i]) == length(ribs[i + 1]) || continue
        base = length(vertices)
        append!(vertices, ribs[i])
        append!(vertices, ribs[i + 1])
        node_count = length(ribs[i])
        for k in 1:node_count - 1
            a1, a2 = base + k, base + k + 1
            b1, b2 = base + node_count + k, base + node_count + k + 1
            push!(faces, Makie.GLTriangleFace(a1, a2, b2))
            push!(faces, Makie.GLTriangleFace(a1, b2, b1))
        end
    end
    return vertices, faces, ribs
end

"""
    plot!(ax, body::VortexStepMethod.BodyAerodynamics; use_observables=false,
          airfoils=false, kwargs...)

Plot a `BodyAerodynamics` object. By default draws each panel as a flat quad; with
`airfoils=true` instead draws the lofted airfoil skin (one see-through wing-shaped
mesh with a contour rib line per section, see [`airfoil_skin_geometry`](@ref)).

If `use_observables=true`, creates observables for dynamic updates keyed by (body_id, panel_index).
Otherwise, creates static plots (original behavior).
"""
function Makie.plot!(ax, body::VortexStepMethod.BodyAerodynamics; color=(:red, 0.2), R_b_w=nothing, T_b_w=nothing,
    use_observables=false, airfoils=false,
    airfoil_color=:deepskyblue, airfoil_opacity=0.2, rib_color=:black,
    border_linewidth=1.5, transparency=true, kwargs...)
    plots = []

    if airfoils
        vertices, faces, ribs = airfoil_skin_geometry(body; R_b_w, T_b_w)
        skin_color = (airfoil_color, airfoil_opacity)
        if use_observables
            isnothing(AIRFOIL_SKIN_OBSERVABLES[]) && (AIRFOIL_SKIN_OBSERVABLES[] = Dict())
            vertices_obs = Observable(vertices)
            rib_obs = [Observable(rib) for rib in ribs]
            isempty(faces) || push!(plots, mesh!(ax, vertices_obs, faces;
                color=skin_color, transparency))
            for rib in rib_obs
                push!(plots, lines!(ax, rib; color=rib_color,
                    linewidth=border_linewidth, transparency))
            end
            AIRFOIL_SKIN_OBSERVABLES[][objectid(body)] =
                (vertices=vertices_obs, ribs=rib_obs)
        else
            isempty(faces) || push!(plots, mesh!(ax, vertices, faces;
                color=skin_color, transparency))
            for rib in ribs
                push!(plots, lines!(ax, rib; color=rib_color,
                    linewidth=border_linewidth, transparency))
            end
        end
        return plots
    end

    if use_observables
        # Initialize global storage if needed
        if isnothing(PANEL_MESH_OBSERVABLES[])
            PANEL_MESH_OBSERVABLES[] = Dict()
        end

        body_id = objectid(body)

        # Create observables for each panel
        for (panel_idx, panel) in enumerate(body.panels)
            # Compute initial points
            points = panel_plate_geometry(panel; R_b_w, T_b_w)

            # Create observables
            vertices_obs = Observable(points)
            faces_obs = Observable(copy(PLATE_FACES))
            border_obs = Observable(points[PLATE_BORDER_IDX])

            # Plot using observables
            p = mesh!(ax, vertices_obs, faces_obs; color, transparency, kwargs...)
            push!(plots, p)
            p = lines!(ax, border_obs; color=:black, linewidth=border_linewidth,
                       transparency, kwargs...)
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
            p = Makie.plot!(ax, panel; color, R_b_w, T_b_w, use_observables=false,
                            border_linewidth, transparency, kwargs...)
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
    body_id = objectid(body)

    # Panel meshes (if plotted with use_observables): refresh from corner_points.
    if !isnothing(PANEL_MESH_OBSERVABLES[])
        for (panel_idx, panel) in enumerate(body.panels)
            key = (body_id, panel_idx)
            haskey(PANEL_MESH_OBSERVABLES[], key) || continue
            obs = PANEL_MESH_OBSERVABLES[][key]
            points = panel_plate_geometry(panel; R_b_w, T_b_w)
            obs.vertices[] = points
            obs.border[] = points[PLATE_BORDER_IDX]
        end
    end

    # Airfoil skin (if plotted with use_observables): refresh from the current pose.
    if !isnothing(AIRFOIL_SKIN_OBSERVABLES[]) &&
            haskey(AIRFOIL_SKIN_OBSERVABLES[], body_id)
        skin = AIRFOIL_SKIN_OBSERVABLES[][body_id]
        vertices, _, ribs = airfoil_skin_geometry(body; R_b_w, T_b_w)
        skin.vertices[] = vertices
        for (i, rib) in enumerate(ribs)
            i <= length(skin.ribs) && (skin.ribs[i][] = rib)
        end
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
    points = panel_plate_geometry(panel; R_b_w, T_b_w)

    vertices_obs = Observable(points)
    faces_obs = Observable(copy(PLATE_FACES))

    # Plot mesh using observables
    mesh!(ax, vertices_obs, faces_obs; color, transparency=true, kwargs...)

    # Plot border
    border_obs = Observable(points[PLATE_BORDER_IDX])
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
        points = panel_plate_geometry(panel; R_b_w, T_b_w)

        # Create observables
        vertices_obs = Observable(points)
        faces_obs = Observable(copy(PLATE_FACES))
        border_obs = Observable(points[PLATE_BORDER_IDX])

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

function _active_backend_prefers_vector_output(makie=Makie)
    isdefined(makie, :current_backend) || return false

    backend = try
        makie.current_backend()
    catch
        return false
    end

    return nameof(backend) == :CairoMakie
end

"""
    save_plot(fig, save_path, title; data_type=nothing)

Save a Makie figure to a file.

# Arguments
- `fig`: Makie Figure object
- `save_path`: Path to save the plot
- `title`: Title of the plot

# Keyword arguments
- `data_type`: File extension. If `nothing`, defaults to `".pdf"` when the
    active Makie backend is CairoMakie and `".png"` otherwise.
"""
function VortexStepMethod.save_plot(fig::Makie.Figure, save_path, title; data_type=nothing)
    if isnothing(data_type)
        data_type = _active_backend_prefers_vector_output() ? ".pdf" : ".png"
    end
    isnothing(save_path) && throw(ArgumentError("save_path should be provided"))

    !isdir(save_path) && mkpath(save_path)
    sanitized_title = replace(replace(String(title), ' ' => '_'), '%' => "pct")
    full_path = joinpath(save_path, sanitized_title * data_type)
    fallback_path = joinpath(save_path, sanitized_title * ".png")

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
        # GLMakie cannot export vector formats such as PDF/SVG directly.
        # If that happens, save as PNG so batch example runs keep working.
        if e isa MethodError && lowercase(data_type) in (".pdf", ".svg")
            @warn "Vector export format $data_type is not supported by the active Makie backend; falling back to PNG" requested_path=full_path fallback_path=fallback_path
            save(fallback_path, fig)
            @debug "Figure saved as PNG fallback"
            return nothing
        end
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
function VortexStepMethod.show_plot(fig::Makie.Figure; dpi=130)
    isinteractive() && display(fig)
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
    va = getfield(body_aero, :_va)

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
                  data_type=nothing, save_path=nothing,
                  is_save=false, is_show=false,
                  view_elevation=15, view_azimuth=-120, use_tex=false)

Makie implementation of [`plot_geometry`](@ref).

# Arguments:
- `body_aero`: the BodyAerodynamics to plot
- `title`: plot title

# Keyword arguments:
- `data_type`: File extension (default: `nothing`; delegated to `save_plot` backend-aware default)
- `save_path`: Path for saving (default: nothing)
- `is_save`: Whether to save (default: false)
- `is_show`: Whether to display (default: false)
- `view_elevation`: View elevation angle in degrees (default: 15)
- `view_azimuth`: View azimuth angle in degrees (default: -120)
- `use_tex`: Ignored for Makie (default: false)
"""
function VortexStepMethod.plot_geometry(body_aero::BodyAerodynamics, title;
    data_type=nothing,
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

    if is_show && isinteractive()
        display(fig)
    end

    return fig
end

"""
    plot_distribution(y_coordinates_list, results_list, label_list;
                      title="spanwise_distribution", data_type=nothing,
                      save_path=nothing, is_save=false, is_show=true, use_tex=false)

Makie implementation of [`plot_distribution`](@ref).

# Arguments
- `y_coordinates_list`: List of spanwise coordinates
- `results_list`: List of result dictionaries
- `label_list`: List of labels for different results

# Keyword arguments
- `title`: Plot title (default: "spanwise_distribution")
- `data_type`: File extension (default: `nothing`; delegated to `save_plot` backend-aware default)
- `save_path`: Path to save plots (default: nothing)
- `is_save`: Whether to save (default: false)
- `is_show`: Whether to display (default: true)
- `use_tex`: Ignored for Makie (default: false)
"""
function VortexStepMethod.plot_distribution(y_coordinates_list, results_list, label_list;
    title="spanwise_distribution",
    data_type=nothing,
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

    # Plot CD
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        value = round(results["cd"], digits=2)
        lines!(ax_cd, Vector(y_coords), Vector(results["cd_distribution"]),
            label="$label CD: $value")
    end

    # Plot Gamma
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_gamma, Vector(y_coords), Vector(results["gamma_distribution"]),
            label=label)
    end

    # Plot alpha geometric
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_alpha_geo, Vector(y_coords), rad2deg.(Vector(results["alpha_geometric"])),
            label=label)
    end

    # Plot alpha at ac
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_alpha_ac, Vector(y_coords), rad2deg.(Vector(results["alpha_at_ac"])),
            label=label)
    end

    # Plot alpha uncorrected
    for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
        lines!(ax_alpha_unc, Vector(y_coords), rad2deg.(Vector(results["alpha_uncorrected"])),
            label=label)
    end

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
    end

    # Shared legend at bottom of grid
    Legend(fig[4, :], ax_gamma;
        orientation=:horizontal, tellwidth=false, tellheight=true)

    # Save and show
    if is_save
        save_plot(fig, save_path, title, data_type=data_type)
    end

    if is_show && isinteractive()
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
- `angle_of_attack`: Initial angle of attack [°]
- `side_slip`: Initial side slip angle [°]
- `v_a`: Norm of apparent wind speed [m/s]

# Returns
- Tuple of polar data array and Reynolds number
"""

"""
    plot_polars(solver_list, body_aero_list, label_list;
                literature_path_list=String[],
                angle_range=range(0, 20, 2), angle_type="angle_of_attack",
                angle_of_attack=0.0, side_slip=0.0, v_a=10.0,
                title="polar", data_type=nothing, save_path=nothing,
                is_save=true, is_show=true, use_tex=false)

Makie implementation of [`plot_polars`](@ref).

# Arguments
- `solver_list`: List of aerodynamic solvers
- `body_aero_list`: List of wing aerodynamics objects
- `label_list`: List of labels for each configuration

# Keyword arguments
- `literature_path_list`: Optional paths to literature data files
- `angle_range`: Range of angles in degrees
- `angle_type`: "angle_of_attack" or "side_slip" (default: angle_of_attack)
- `angle_of_attack`: AoA [°] (default: 0.0)
- `side_slip`: Side slip angle [°] (default: 0.0)
- `v_a`: Wind speed [m/s] (default: 10.0)
- `title`: Plot title
- `data_type`: File extension (default: `nothing`; delegated to `save_plot` backend-aware default)
- `save_path`: Path to save (default: nothing)
- `is_save`: Whether to save (default: true)
- `is_show`: Whether to display (default: true)
- `use_tex`: Ignored for Makie (default: false)
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
    data_type=nothing,
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
        throw(ArgumentError("Mismatch in solvers/cases/labels"))
    end
    main_title = replace(title, " " => "_")

    # Generate polar data
    polar_data_list = []
    cm_data_list = []
    labels_with_re = copy(label_list)
    for (i, (solver, body_aero)) in enumerate(zip(solver_list, body_aero_list))
        result = VortexStepMethod.generate_polar_data(
            solver, body_aero, angle_range;
            angle_type, angle_of_attack, side_slip, v_a
        )
        push!(polar_data_list, result.polar_data)
        push!(cm_data_list, (cmx=result.cmx, cmy=result.cmy,
                             cmz=result.cmz))
        labels_with_re[i] = "$(label_list[i]) Re = $(round(Int64, result.rey*1e-5))e5"
    end

    # Load literature data if provided
    if !isempty(literature_path_list)
        for path in literature_path_list
            lit = VortexStepMethod.extract_literature_polar_data(
                readdlm(path, ','), path; angle_type)
            push!(polar_data_list, lit.polar_data)
            push!(cm_data_list, (cmx=lit.cmx, cmy=lit.cmy,
                                 cmz=lit.cmz))
        end
    end

    # Number of computational results
    n_solvers = length(solver_list)

    if show_moments
        # 2x3 layout: CL, CD, CS, CMx, CMy, CMz
        fig = Figure(size=(1800, 1000))
        ax_cl = Axis(fig[1, 1], title="CL vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CL")
        ax_cd = Axis(fig[1, 2], title="CD vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CD")
        ax_cs = Axis(fig[1, 3], title="CS vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CS")
        ax_cmx = Axis(fig[2, 1], title="CMx vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CMx")
        ax_cmy = Axis(fig[2, 2], title="CMy vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CMy")
        ax_cmz = Axis(fig[2, 3], title="CMz vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CMz")

        for ax in (ax_cl, ax_cd, ax_cs, ax_cmx, ax_cmy, ax_cmz)
            ax.yticklabelspace = 36.0
            ax.xticklabelspace = 24.0
        end

        for (i, (polar_data, cm, label)) in enumerate(
                zip(polar_data_list, cm_data_list, labels_with_re))
            marker = i <= n_solvers ? :star5 : :circle
            markersize = i <= n_solvers ? 12 : 8
            linestyle = i <= n_solvers ? :solid : :dash
            angles = polar_data[1]
            scatterlines!(ax_cl, angles, polar_data[2];
                label=label, marker=marker, markersize=markersize, linestyle)
            scatterlines!(ax_cd, angles, polar_data[3];
                label=label, marker=marker, markersize=markersize, linestyle)
            scatterlines!(ax_cs, angles, polar_data[4];
                label=label, marker=marker, markersize=markersize, linestyle)
            if !all(isnan, cm.cmx)
                scatterlines!(ax_cmx, angles,
                    Float64.(cm.cmx);
                    label=label, marker=marker,
                    markersize=markersize, linestyle)
            end
            if !all(isnan, cm.cmy)
                scatterlines!(ax_cmy, angles,
                    Float64.(cm.cmy);
                    label=label, marker=marker,
                    markersize=markersize, linestyle)
            end
            if !all(isnan, cm.cmz)
                scatterlines!(ax_cmz, angles,
                    Float64.(cm.cmz);
                    label=label, marker=marker,
                    markersize=markersize, linestyle)
            end
        end
        Legend(fig[3, 1:3], ax_cl;
            orientation=:horizontal, tellwidth=false,
            tellheight=true)
    else
        # 2x2 layout: CL, CD, CS, CL/CD or CL-vs-CD
        fig = Figure(size=(1400, 1400))
        ax_cl = Axis(fig[1, 1], title="CL vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CL")
        ax_cd = Axis(fig[1, 2], title="CD vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CD")
        ax_cs = Axis(fig[2, 1], title="CS vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CS")
        ax_fourth = if cl_over_cd
            Axis(fig[2, 2], title="CL/CD vs $angle_type [°]",
                xlabel="$angle_type [°]", ylabel="CL/CD")
        else
            Axis(fig[2, 2], title="CL vs CD",
                xlabel="CD", ylabel="CL")
        end

        for (i, (polar_data, label)) in enumerate(
                zip(polar_data_list, labels_with_re))
            marker = i <= n_solvers ? :star5 : :circle
            markersize = i <= n_solvers ? 12 : 8
            linestyle = i <= n_solvers ? :solid : :dash
            scatterlines!(ax_cl, polar_data[1], polar_data[2];
                label=label, marker=marker,
                markersize=markersize, linestyle)
            scatterlines!(ax_cd, polar_data[1], polar_data[3];
                label=label, marker=marker,
                markersize=markersize, linestyle)
            scatterlines!(ax_cs, polar_data[1], polar_data[4];
                label=label, marker=marker,
                markersize=markersize, linestyle)
            if cl_over_cd
                cl_cd = polar_data[2] ./ polar_data[3]
                scatterlines!(ax_fourth, polar_data[1], cl_cd;
                    label=label, marker=marker,
                    markersize=markersize, linestyle)
            else
                scatterlines!(ax_fourth, polar_data[3],
                    polar_data[2];
                    label=label, marker=marker,
                    markersize=markersize, linestyle)
            end
        end
        Legend(fig[3, :], ax_cl;
            orientation=:horizontal, tellwidth=false,
            tellheight=true)
    end

    # Save and show
    if is_save && !isnothing(save_path)
        save_plot(fig, save_path, main_title; data_type)
    end

    if is_show && isinteractive()
        display(fig)
    end

    return fig
end

"""
    plot_polar_data(body_aero::BodyAerodynamics;
                    alphas=collect(deg2rad.(-5:0.3:25)),
                    delta_tes=collect(deg2rad.(-5:0.3:25)),
                    is_show=true, use_tex=false)

Makie implementation of [`plot_polar_data`](@ref).

# Arguments
- `body_aero`: Wing aerodynamics struct

# Keyword arguments
- `alphas`: Range of AoA values in radians (default: `deg2rad.(-5:0.3:25)`)
- `delta_tes`: Range of trailing edge angles in radians (default: `deg2rad.(-5:0.3:25)`)
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
                             for delta_te in delta_tes, alpha in alphas]

            # Create wireframe
            wireframe!(ax, delta_tes, alphas, interp_matrix;
                color=:blue, linewidth=0.5, transparency=true)
        end

        if is_show && isinteractive()
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

Makie implementation of [`plot_combined_analysis`](@ref).

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
- `view_elevation`: Geometry view elevation in degrees (default: 15)
- `view_azimuth`: Geometry view azimuth in degrees (default: -120)
- `is_show`: Display figure (default: true)
- `use_tex`: Ignored for Makie (default: false)
- `literature_path_list`: Paths to literature CSV files (default: String[])
- `data_type`: File extension (default: ".png", also supports ".jpeg")
- `save_path`: Directory path to save files (default: nothing)
- `is_save`: Save plots to files (default: false)
- `cl_over_cd`: Plot CL/CD vs angle instead of CL vs CD (default: true)
- `angle_of_attack_for_spanwise_distribution`: AoA for spanwise plots (default: 5.0)
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
    cl_over_cd=true,
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
    va = getfield(first_body, :_va)

    # Compute spanwise results for each solver
    results_spanwise_list = copy(results_list)
    if !isnothing(angle_of_attack_for_spanwise_distribution)
        α_span = deg2rad(angle_of_attack_for_spanwise_distribution)
        β_span = deg2rad(side_slip)
        for (i, (s, ba)) in enumerate(zip(solvers, body_aeros))
            va_old = copy(getfield(ba, :_va))
            omega_old = copy(ba.omega)
            set_va!(ba, [cos(α_span) * cos(β_span), sin(β_span),
                sin(α_span)] * v_a)
            results_spanwise_list[i] = solve(s, ba,
                s.sol.gamma_distribution)
            set_va!(ba, va_old, omega_old)
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

    Legend(fig[2, 1][4, :], ax_gamma;
        orientation=:horizontal, tellwidth=false, tellheight=true)

    # force y-limits
    ylims!(ax_alpha_geo, -0.25 * aoa_span, 1.1 * aoa_span)

    # [2,2] Polars (2×2 grid)
    polar_series = Tuple[]
    for (si, (s, ba, lbl)) in enumerate(
            zip(solvers, body_aeros, solver_labels))
        result = VortexStepMethod.generate_polar_data(
            s, ba, angle_range;
            angle_type, angle_of_attack, side_slip, v_a)
        pd = result.polar_data
        label_re = "$lbl Re = $(round(Int64,
                     result.rey * 1e-5))e5"
        push!(polar_series, (
            (angle=pd[1], cl=pd[2], cd=pd[3], cs=pd[4],
             cmy=result.cmy,
             rey=pd[9]),
            label_re))
    end

    # Load literature data (if any)
    if !isempty(literature_path_list)
        for (path, lit_label) in zip(
                literature_path_list, literature_labels)
            lit = VortexStepMethod.extract_literature_polar_data(
                readdlm(path, ','), path; angle_type)
            pd = lit.polar_data
            push!(polar_series, (
                (angle=pd[1], cl=pd[2], cd=pd[3], cs=pd[4],
                 cmy=lit.cmy,
                 rey=fill(NaN, length(pd[1]))),
                lit_label))
        end
    end

    ax_cl_polar = Axis(fig[2, 2][1, 1], title="CL vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CL")
    ax_cd_polar = Axis(fig[2, 2][1, 2], title="CD vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CD")
    ax_cs_polar = Axis(fig[2, 2][2, 1], title="CS vs $angle_type [°]",
        xlabel="$angle_type [°]", ylabel="CS")
    ax_fourth_polar = if cl_over_cd
        Axis(fig[2, 2][2, 2], title="CL/CD vs $angle_type [°]",
            xlabel="$angle_type [°]", ylabel="CL/CD")
    else
        Axis(fig[2, 2][2, 2], title="CL vs CD",
            xlabel="CD", ylabel="CL")
    end

    for (idx, (pd, lbl)) in enumerate(polar_series)
        color = colors[mod1(idx, length(colors))]
        marker = idx == 1 ? :star5 : :circle
        markersize = idx == 1 ? 12 : 8

        scatterlines!(ax_cl_polar, pd.angle, pd.cl; label=lbl, marker, markersize, color)
        scatterlines!(ax_cd_polar, pd.angle, pd.cd; label=lbl, marker, markersize, color)
        scatterlines!(ax_cs_polar, pd.angle, pd.cs; label=lbl, marker, markersize, color)
        if cl_over_cd
            cl_cd = pd.cl ./ pd.cd
            scatterlines!(ax_fourth_polar, pd.angle, cl_cd;
                label=lbl, marker, markersize, color)
        else
            scatterlines!(ax_fourth_polar, pd.cd, pd.cl;
                label=lbl, marker, markersize, color)
        end
    end
    Legend(fig[2, 2][3, :], ax_cl_polar;
        orientation=:horizontal, nbanks=2,
        tellwidth=false, tellheight=true)

    # Set column widths: left column wider for 3x3 grid
    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    if is_show && isinteractive()
        display(fig)
    end

    return fig
end

"""
    plot_section_polars(body_aero, coefficient=:cl; is_show=true,
                        is_save=false, save_path=nothing, data_type=".png")

Implementation of [`plot_section_polars`](@ref); rendered through `MakieControlPlots`.
"""
function VortexStepMethod.plot_section_polars(body_aero::BodyAerodynamics,
    coefficient::Symbol=:cl; is_show::Bool=true, is_save::Bool=false,
    save_path=nothing, data_type::String=".png")

    coefficient in (:cl, :cd, :cm) ||
        throw(ArgumentError("coefficient must be :cl, :cd, or :cm, got :$coefficient"))
    idx = coefficient === :cl ? 2 : coefficient === :cd ? 3 : 4
    label = uppercasefirst(string(coefficient))

    alphas_deg = nothing
    series = Vector{Float64}[]
    labels = String[]
    for wing in body_aero.wings
        for (s, section) in enumerate(wing.unrefined_sections)
            section.aero_model == POLAR_VECTORS || continue
            aero = section.aero_data
            aero === nothing && continue
            section_alphas = rad2deg.(aero[1])
            if isnothing(alphas_deg)
                alphas_deg = collect(section_alphas)
            elseif length(section_alphas) != length(alphas_deg)
                @warn "section $s has a different α grid; plotting against the first section's α"
            end
            push!(series, Float64.(aero[idx]))
            push!(labels, "section $s")
        end
    end
    isempty(series) && error("No POLAR_VECTORS sections found in body")

    plt = MakieControlPlots.plot(alphas_deg, series;
        xlabel="α [deg]", ylabel=label, title="$label per section",
        labels=labels, disp=(is_show || is_save))

    if is_save && !isnothing(save_path)
        isdir(save_path) || mkpath(save_path)
        MakieControlPlots.savefig(
            joinpath(save_path, "section_polars_$(coefficient)$(data_type)"))
    end
    return plt
end

# --- OBJ mesh / airfoil plotting (formerly ObjAdapterMakieExt) ---

"""
    fitted_airfoil_3d(section, wrap_method; delta=0.0, crease_frac=0.75) -> 3×N or nothing

Shrink-wrap a section's sliced contour with `wrap_method` and map it back into 3D
through the section's local airfoil frame, for overlaying on the 3D slice diagnostic.
A nonzero `delta` (degrees) deflects the trailing edge and re-wraps
([`deform_section`](@ref)), showing the geometry the solvers consume.
Returns `nothing` for a degenerate slice.
"""
function fitted_airfoil_3d(s, wrap_method; delta=0.0, crease_frac=0.75)
    frame = ObjAdapter.airfoil_frame(s.LE_point, s.TE_point, s.span_dir)
    frame === nothing && return nothing
    x_af, _, z_af = frame
    px = [dot(p .- s.LE_point, x_af) for p in s.contour3d]
    pz = [dot(p .- s.LE_point, z_af) for p in s.contour3d]
    chord = maximum(px) - minimum(px)
    chord < 1e-9 && return nothing
    x0 = minimum(px)
    xf, yf = AirfoilAero.shrink_wrap(collect((px .- x0) ./ chord),
                                     collect(pz ./ chord), wrap_method)
    maximum(abs, yf) > 1.0 && return nothing
    if !iszero(delta)
        def = AirfoilAero.deform_section(xf, yf, deg2rad(delta); crease_frac)
        xf, yf = def.x, def.y
    end
    return reduce(hcat, [s.LE_point .+ x_af .* (x0 + xf[i] * chord) .+ z_af .* (yf[i] * chord)
                         for i in eachindex(xf)])
end

"""
    map_airfoil_3d(le, te, tangent, x, y) -> 3×N or nothing

Map normalized airfoil coordinates into 3D through the airfoil frame of the station
given by its LE/TE points and leading-edge tangent (chord scale = `|TE - LE|`).
"""
function map_airfoil_3d(le, te, tangent, x, y)
    frame = ObjAdapter.airfoil_frame(le, te, tangent)
    (frame === nothing || isempty(x)) && return nothing
    x_af, _, z_af = frame
    chord = norm(te .- le)
    return reduce(hcat, [le .+ x_af .* (x[i] * chord) .+ z_af .* (y[i] * chord)
                         for i in eachindex(x)])
end

"""
    generated_slices(out_dir, delta, fit_pts) -> (slices, le, te)

Read the stations of a generated [`obj_to_yaml`](@ref) output directory and their
written `.dat` airfoils — raw slice, wrap, and the `delta`-degree deformed wrap when
it was generated — assembled for [`plot_slices_3d`](@ref). Nothing is re-sliced or
re-wrapped; only the Kulfan fits of the stored coordinates are recomputed (via
`fit_pts`), exactly as the polar pipeline fits them.
"""
function generated_slices(out_dir, delta, fit_pts)
    geom = VortexStepMethod.YAML.load_file(joinpath(out_dir, "geometry.yaml"))
    info = Dict(r[1] => r[3] for r in geom["wing_airfoils"]["data"])
    rows = geom["wing_sections"]["data"]
    les = [Float64.(r[2:4]) for r in rows]
    tes = [Float64.(r[5:7]) for r in rows]
    n = length(rows)
    deg = round(float(delta); digits=1)
    tag = "_d" * (deg == round(deg) ? string(Int(deg)) : string(deg)) * ".dat"
    missing_deltas = String[]
    slices = map(1:n) do i
        id = rows[i][1]
        tangent = normalize(les[min(i + 1, n)] .- les[max(i - 1, 1)])
        xw, yw = AirfoilAero.read_dat_coordinates(joinpath(out_dir,
                                                           info[id]["dat_file"]))
        xr, yr = AirfoilAero.read_dat_coordinates(joinpath(out_dir,
                                                           info[id]["raw_dat_file"]))
        d2 = (; raw=Point2f.(xr, yr), fit=Point2f.(xw, yw),
              fit_kulfan=fit_pts(xw, yw), def=Point2f[], def_kulfan=Point2f[])
        def3d = nothing
        if !iszero(delta)
            dpath = joinpath(out_dir, "airfoils", "$(id)$(tag)")
            if isfile(dpath)
                xd, yd = AirfoilAero.read_dat_coordinates(dpath)
                def3d = map_airfoil_3d(les[i], tes[i], tangent, xd, yd)
                d2 = (; d2..., def=Point2f.(xd, yd), def_kulfan=fit_pts(xd, yd))
            else
                push!(missing_deltas, basename(dpath))
            end
        end
        (; centroid=Point3f((les[i] .+ tes[i]) ./ 2), label_y=les[i][2],
         cloud3d=map_airfoil_3d(les[i], tes[i], tangent, xr, yr),
         wrap3d=map_airfoil_3d(les[i], tes[i], tangent, xw, yw), def3d, d2)
    end
    isempty(missing_deltas) ||
        @warn "No generated .dat for delta=$(delta)° ($(join(missing_deltas, ", ")));" *
              " generated deflections are named airfoils/<i>_d<degrees>.dat."
    return slices, reduce(hcat, les), reduce(hcat, tes)
end

"""
    plot_slices_3d(path; n_slices=10, rotation=I, n_bins=60, wingtip_distance=0.0,
                   wrap_method=ShrinkWrap(), delta=0.0, crease_frac=0.75,
                   obj_path=nothing, is_show=true)

3D slice diagnostic with a hover 2D airfoil panel: raw slice points (green), the
shrink-wrapped airfoil (crimson) and its Kulfan fit (black dashed), and, for nonzero
`delta` (degrees), the deflected re-wrapped airfoil (purple — the
[`deform_section`](@ref) geometry the solvers consume) and its Kulfan fit (orange
dashed). **Hovering a slice** updates the 2D panel.

`path` selects the source:
- a mesh `.obj` file: live preview — slices and wraps here with `wrap_method`
  (`n_slices`, `n_bins`, `wingtip_distance`, `crease_frac`).
- a generated [`obj_to_yaml`](@ref) output directory: audit mode — stations and
  airfoils are read from `geometry.yaml` and the written `.dat` files, so the plot
  shows exactly what the polar pipeline analysed. `delta` must then match a
  generated deflection value; pass `obj_path` to also draw the mesh (with the same
  `rotation` used at generation).

Pass a `3×3` rotation matrix to reorient the mesh first.
"""
function ObjAdapter.plot_slices_3d(path::String; n_slices::Int=10, rotation=I,
        n_bins::Int=60, wingtip_distance=0.0,
        wrap_method=AirfoilAero.ShrinkWrap(), delta=0.0, crease_frac=0.75,
        obj_path=nothing, is_show::Bool=true)
    kulfan_pts(k) = Point2f.(AirfoilAero.kulfan_to_coordinates(k; n_points=150)...)
    fit_pts(x, y) = kulfan_pts(AirfoilAero.fit_kulfan_parameters(
        x, y, AirfoilAero.LeastSquaresFit()))
    mesh_path = isdir(path) ? obj_path : path
    vertices = faces = nothing
    if mesh_path !== nothing
        vertices, faces = ObjAdapter.read_faces(mesh_path)
        rotation === I || (vertices = [rotation * v for v in vertices])
    end

    if isdir(path)
        slices, le, te = generated_slices(path, delta, fit_pts)
    else
        span = maximum(v -> v[2], vertices) - minimum(v -> v[2], vertices)
        m = ObjAdapter.march_edges(vertices, faces; step=span / n_bins)
        le = reduce(hcat, m.le)
        te = reduce(hcat, m.te)
        idx = ObjAdapter.station_indices(m.arclen, n_slices; wingtip_distance)
        secs = filter(!isnothing,
                      [ObjAdapter.build_section(vertices, faces, m.le[i], m.te[i],
                                                m.point[i], m.tangent[i]) for i in idx])
        slices = map(secs) do s
            xf, yf = AirfoilAero.shrink_wrap(collect(Float64, s.x_airfoil),
                                             collect(Float64, s.y_airfoil), wrap_method)
            def = iszero(delta) ? nothing :
                  AirfoilAero.deform_section(xf, yf, deg2rad(delta); crease_frac)
            d2 = (; raw=Point2f.(s.x_airfoil, s.y_airfoil), fit=Point2f.(xf, yf),
                  fit_kulfan=fit_pts(xf, yf),
                  def=def === nothing ? Point2f[] : Point2f.(def.x, def.y),
                  def_kulfan=def === nothing ? Point2f[] : kulfan_pts(def.kulfan))
            (; centroid=Point3f((s.LE_point .+ s.TE_point) ./ 2),
             label_y=s.LE_point[2], cloud3d=reduce(hcat, s.contour3d),
             wrap3d=fitted_airfoil_3d(s, wrap_method),
             def3d=def === nothing ? nothing :
                   fitted_airfoil_3d(s, wrap_method; delta, crease_frac), d2)
        end
    end
    show_delta = any(!isempty(s.d2.def) for s in slices)

    fig = Figure(size=(1500, 850))
    ax = Axis3(fig[1, 1]; aspect=:data, title="slices (hover a slice →)")
    if vertices !== nothing
        coords = permutedims(reduce(hcat, vertices))
        tri = permutedims(reduce(hcat, [Int.(f) for f in faces]))
        mesh!(ax, coords, tri; color=(:gray, 0.15), transparency=true)
    end
    lines!(ax, le[1, :], le[2, :], le[3, :]; color=:dodgerblue, linewidth=3, label="LE")
    lines!(ax, te[1, :], te[2, :], te[3, :]; color=:orange, linewidth=3, label="TE")
    for s in slices
        s.cloud3d === nothing ||
            scatter!(ax, s.cloud3d[1, :], s.cloud3d[2, :], s.cloud3d[3, :];
                     color=:seagreen, markersize=3)
        s.wrap3d === nothing ||
            lines!(ax, s.wrap3d[1, :], s.wrap3d[2, :], s.wrap3d[3, :];
                   color=:crimson, linewidth=1.5)
        s.def3d === nothing ||
            lines!(ax, s.def3d[1, :], s.def3d[2, :], s.def3d[3, :];
                   color=:purple, linewidth=1.5)
    end
    axislegend(ax; position=:rt)

    sel = Observable(max(1, cld(length(slices), 2)))
    title2 = @lift("slice $($sel)  (y = $(round(slices[$sel].label_y, digits=2)))")
    gl2 = GridLayout(fig[1, 2])
    ax2 = Axis(gl2[2, 1]; title=title2, xlabel="x/c", ylabel="y/c", aspect=DataAspect())
    colsize!(fig.layout, 1, Relative(0.6))
    raw2 = scatter!(ax2, @lift(slices[$sel].d2.raw); color=:seagreen, markersize=6)
    fit2 = lines!(ax2, @lift(slices[$sel].d2.fit); color=:crimson, linewidth=2)
    fitk2 = lines!(ax2, @lift(slices[$sel].d2.fit_kulfan); color=:black,
                   linewidth=1.5, linestyle=:dash)
    handles, labels = [raw2, fit2, fitk2], ["slice", "wrap", "kulfan"]
    if show_delta
        def2 = lines!(ax2, @lift(slices[$sel].d2.def); color=:purple, linewidth=2)
        defk2 = lines!(ax2, @lift(slices[$sel].d2.def_kulfan); color=:darkorange,
                       linewidth=1.5, linestyle=:dash)
        push!(handles, def2, defk2)
        push!(labels, "δ=$(delta)° wrap", "δ kulfan")
    end
    ylims!(ax2, -1, 1)   # fixed y/c so airfoil thickness is comparable across slices
    Legend(gl2[1, 1], handles, labels; orientation=:horizontal, framevisible=false)

    on(events(ax.scene).mouseposition) do mp
        is_mouseinside(ax.scene) || return
        best, best_d = 1, Inf
        for (i, s) in enumerate(slices)
            pr = Makie.project(ax.scene, s.centroid)
            d = (pr[1] - mp[1])^2 + (pr[2] - mp[2])^2
            d < best_d && (best_d = d; best = i)
        end
        sel[] = best
    end

    is_show && display(fig)
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
- Makie Figure object and the fitted `KulfanParameters`
"""
function ObjAdapter.plot_airfoil_fit(x::Vector, y::Vector; title::String="Airfoil Fit",
                                     is_show::Bool=true)
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1, 1];
              title=title,
              xlabel="x/c",
              ylabel="y/c",
              aspect=DataAspect())

    x_norm, y_norm, _ = AirfoilAero.normalize_airfoil(x, y)
    scatter!(ax, x_norm, y_norm;
            color=:blue, markersize=6, label="Input ($(length(x)) pts)")

    params = AirfoilAero.fit_kulfan_parameters(x, y)
    x_fit, y_fit = AirfoilAero.kulfan_to_coordinates(params; n_points=150)
    lines!(ax, x_fit, y_fit;
          color=:red, linewidth=2, linestyle=:dash, label="Kulfan fit")

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

    is_show && display(fig)
    return fig, params
end

"""
    plot_airfoils(geometry_file; overlay=nothing, symmetric=false,
                  idxs=nothing, n_cols=3, is_show=true, is_save=false,
                  save_path=nothing, data_type=".png")

Makie implementation of [`plot_airfoils`](@ref).

With `overlay=false` each airfoil gets its own subplot; with `overlay=true` all
airfoils are drawn on one axis coloured by section. By default (`overlay=nothing`)
the overlay is used when there are more than 12 sections, where a subplot grid
becomes unreadable. Both modes show the raw `_raw.dat` slice points as dots with
the fitted airfoil as a line. Pass `idxs` (e.g. `idxs=[1]`) to plot only those
airfoils by position; otherwise `symmetric=true` shows just the first half of the
sections (the wing is mirror-symmetric, so the other half is redundant).
"""
function ObjAdapter.plot_airfoils(geometry_file::String;
    overlay=nothing, symmetric::Bool=false, idxs=nothing, n_cols::Int=3,
    is_show::Bool=true, is_save::Bool=false, save_path=nothing,
    data_type::String=".png")

    airfoils = ObjAdapter.airfoils_from_yaml(geometry_file)
    isempty(airfoils) && error("No airfoils with a dat_file found in $geometry_file")
    if idxs !== nothing
        airfoils = airfoils[idxs]
    elseif symmetric
        airfoils = airfoils[1:cld(length(airfoils), 2)]
    end

    n = length(airfoils)
    use_overlay = overlay === nothing ? n > 12 : overlay
    title = "Airfoils: $(basename(geometry_file))"

    if use_overlay
        ids = [af.id for af in airfoils]
        crange = (minimum(ids), maximum(ids))
        fig = Figure(size=(900, 600))
        ax = Axis(fig[1, 1]; title, xlabel="x/c", ylabel="y/c", aspect=DataAspect())
        for af in airfoils
            isempty(af.x_raw) || scatter!(ax, af.x_raw, af.y_raw;
                color=af.id, colorrange=crange, colormap=:viridis, markersize=3)
            lines!(ax, af.x, af.y; color=af.id, colorrange=crange, colormap=:viridis)
        end
        Colorbar(fig[1, 2]; colormap=:viridis, limits=crange, label="section")
    else
        ncol = min(n_cols, n)
        nrow = ceil(Int, n / ncol)
        fig = Figure(size=(380 * ncol, 320 * nrow))
        Label(fig[0, :], title, fontsize=16)
        for (i, af) in enumerate(airfoils)
            ax = Axis(fig[div(i - 1, ncol) + 1, mod1(i, ncol)];
                title="Airfoil $(af.id)", xlabel="x/c", ylabel="y/c",
                aspect=DataAspect())
            isempty(af.x_raw) ||
                scatter!(ax, af.x_raw, af.y_raw; color=(:gray, 0.5), markersize=4)
            lines!(ax, af.x, af.y; color=:black)
        end
    end

    if is_save && !isnothing(save_path)
        VortexStepMethod.save_plot(fig, save_path, "airfoils"; data_type)
    end
    is_show && display(fig)
    return fig
end

end

module ObjAdapterMakieExt

using Makie
using LinearAlgebra
import ObjAdapter
import AirfoilAero
import VortexStepMethod

"""
    fitted_airfoil_3d(section, fit_method) -> 3×N matrix or nothing

Fit a Kulfan airfoil to a section's sliced contour with `fit_method` and map it back
into 3D through the section's local airfoil frame, for overlaying on the 3D slice
diagnostic. Returns `nothing` for a degenerate slice or a blown-up fit.
"""
function fitted_airfoil_3d(s, fit_method)
    frame = ObjAdapter.airfoil_frame(s.LE_point, s.TE_point, s.span_dir)
    frame === nothing && return nothing
    x_af, _, z_af = frame
    px = [dot(p .- s.LE_point, x_af) for p in s.contour3d]
    pz = [dot(p .- s.LE_point, z_af) for p in s.contour3d]
    chord = maximum(px) - minimum(px)
    chord < 1e-9 && return nothing
    x0 = minimum(px)
    params = AirfoilAero.fit_kulfan_parameters(collect((px .- x0) ./ chord),
                                               collect(pz ./ chord), fit_method)
    maximum(abs, [params.upper_weights; params.lower_weights]) > 50 && return nothing
    xf, yf = AirfoilAero.kulfan_to_coordinates(params)
    return reduce(hcat, [s.LE_point .+ x_af .* (x0 + xf[i] * chord) .+ z_af .* (yf[i] * chord)
                         for i in eachindex(xf)])
end

"""
    plot_slices_3d(obj_path; n_slices=10, rotation=I, n_bins=60,
                   fit_method=EnvelopeFit(min_distance=0.01), is_show=true)

3D diagnostic: the OBJ surface mesh with the LE/TE edge curves, each station's
sliced contour (green), and the `fit_method` Kulfan fit (crimson). **Hovering a
slice** updates the 2D airfoil panel on the right (raw slice points + the fit).
Pass `rotation=:auto` (or a 3×3 matrix) to reorient the mesh first.
"""
function ObjAdapter.plot_slices_3d(obj_path::String; n_slices::Int=10, rotation=I,
        n_bins::Int=60, fit_method=AirfoilAero.EnvelopeFit(min_distance=0.01),
        is_show::Bool=true)
    vertices, faces = ObjAdapter.read_faces(obj_path)
    rotation === :auto && (rotation = ObjAdapter.auto_rotation(vertices))
    rotation === I || (vertices = [rotation * v for v in vertices])

    fig = Figure(size=(1500, 850))
    ax = Axis3(fig[1, 1]; aspect=:data, title="OBJ slices (hover a slice →)")

    coords = permutedims(reduce(hcat, vertices))
    tri = permutedims(reduce(hcat, [Int.(f) for f in faces]))
    mesh!(ax, coords, tri; color=(:gray, 0.15), transparency=true)

    span = maximum(v -> v[2], vertices) - minimum(v -> v[2], vertices)
    m = ObjAdapter.march_edges(vertices, faces; step=span / n_bins)
    le = reduce(hcat, m.le)
    te = reduce(hcat, m.te)
    lines!(ax, le[1, :], le[2, :], le[3, :]; color=:dodgerblue, linewidth=3, label="LE")
    lines!(ax, te[1, :], te[2, :], te[3, :]; color=:orange, linewidth=3, label="TE")

    idx = ObjAdapter.station_indices(m.arclen, n_slices)
    secs = filter(!isnothing, [ObjAdapter.build_section(vertices, faces, m.le[i], m.te[i],
                                                        m.point[i], m.tangent[i]) for i in idx])
    centroids = [Point3f((s.LE_point .+ s.TE_point) ./ 2) for s in secs]
    fit2d(s) = begin
        p = AirfoilAero.fit_kulfan_parameters(collect(Float64, s.x_airfoil),
                                              collect(Float64, s.y_airfoil), fit_method)
        maximum(abs, [p.upper_weights; p.lower_weights]) > 50 && return Point2f[]
        Point2f.(AirfoilAero.kulfan_to_coordinates(p)...)
    end
    data2d = [(; raw=Point2f.(s.x_airfoil, s.y_airfoil), fit=fit2d(s)) for s in secs]

    for s in secs
        c = reduce(hcat, s.contour3d)
        scatter!(ax, c[1, :], c[2, :], c[3, :]; color=:seagreen, markersize=3)
        f = fitted_airfoil_3d(s, fit_method)
        f === nothing || lines!(ax, f[1, :], f[2, :], f[3, :]; color=:crimson, linewidth=1.5)
    end
    axislegend(ax; position=:rt)

    sel = Observable(1)
    title2 = @lift("slice $($sel)  (y = $(round(secs[$sel].LE_point[2], digits=2)))")
    ax2 = Axis(fig[1, 2]; title=title2, xlabel="x/c", ylabel="y/c", aspect=DataAspect())
    colsize!(fig.layout, 1, Relative(0.6))
    scatter!(ax2, @lift(data2d[$sel].raw); color=:seagreen, markersize=6, label="slice")
    lines!(ax2, @lift(data2d[$sel].fit); color=:crimson, linewidth=2, label="fit")
    axislegend(ax2; position=:rt)

    on(events(ax.scene).mouseposition) do mp
        best, best_d = 1, Inf
        for (i, ctr) in enumerate(centroids)
            pr = Makie.project(ax.scene, ctr)
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
    plot_airfoils(geometry_file, ::MakieBackend; overlay=nothing, symmetric=false,
                  idxs=nothing, n_cols=3, is_show=true, is_save=false,
                  save_path=nothing, data_type=".png")

Makie backend implementation of [`plot_airfoils`](@ref).

With `overlay=false` each airfoil gets its own subplot; with `overlay=true` all
airfoils are drawn on one axis coloured by section. By default (`overlay=nothing`)
the overlay is used when there are more than 12 sections, where a subplot grid
becomes unreadable. Both modes show the raw `_raw.dat` slice points as dots with
the fitted airfoil as a line. Pass `idxs` (e.g. `idxs=[1]`) to plot only those
airfoils by position; otherwise `symmetric=true` shows just the first half of the
sections (the wing is mirror-symmetric, so the other half is redundant).
"""
function ObjAdapter.plot_airfoils(geometry_file::String,
    ::VortexStepMethod.MakieBackend;
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

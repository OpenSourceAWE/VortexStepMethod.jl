module VortexStepMethodMakieExt
using Makie, VortexStepMethod

"""
    plot!(ax::Makie.Axis3, panel::VortexStepMethod.Panel; kwargs...)

Plot a single `Panel` as a `mesh`.
The corner points are ordered as: LE1, TE1, TE2, LE2.
This creates two triangles: (LE1, TE1, TE2) and (LE1, TE2, LE2).
"""
function Makie.plot!(ax::Makie.Axis3, panel::VortexStepMethod.Panel; kwargs...)
    points = [Point3f(panel.corner_points[:, i]) for i in 1:4]
    faces = [Makie.GLTriangleFace(1, 2, 3), Makie.GLTriangleFace(1, 3, 4)]
    mesh!(ax, points, faces; kwargs...)
end

"""
    plot!(ax::Makie.Axis3, body::VortexStepMethod.BodyAerodynamics; kwargs...)

Plot a `BodyAerodynamics` object by plotting each of its panels.
"""
function Makie.plot!(ax::Makie.Axis3, body::VortexStepMethod.BodyAerodynamics; kwargs...)
    for panel in body.panels
        Makie.plot!(ax, panel; kwargs...)
    end
end

function Makie.plot(panel::VortexStepMethod.Panel; size = (1200, 800), kwargs...)
    fig = Figure(; size)
    ax = Axis3(fig[1, 1]; aspect = :data,
               xlabel = "X", ylabel = "Y", zlabel = "Z",
               azimuth = 9/8*π, zoommode = :cursor, viewmode = :fit,
               xautolimitmargin=(1.0, 1.0),
               yautolimitmargin=(10.0, 10.0),
           )

    plot!(ax, panel; kwargs...)
    return fig
end

function Makie.plot(body_aero::VortexStepMethod.BodyAerodynamics; size = (1200, 800), kwargs...)
    fig = Figure(; size)
    ax = Axis3(fig[1, 1]; aspect = :data,
               xlabel = "X", ylabel = "Y", zlabel = "Z",
               azimuth = 9/8*π, zoommode = :cursor, viewmode = :fit,
               # xautolimitmargin=(1.0, 1.0),
               # yautolimitmargin=(10.0, 10.0),
           )

    plot!(ax, body_aero; kwargs...)
    return fig
end

end

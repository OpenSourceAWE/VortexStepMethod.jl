"""
    XFoilSolver

XFoil backend. Loads the deformed coordinates (optionally repaneling them) and reads
the pressure distribution from `Xfoil.cpdump`. Defaults match NeuralFoil's training
runs (`ncrit=9`, `max_iter=100`, incompressible) so the two backends are comparable.

# Fields
- `npan`: XFoil panel count when repaneling (default 160, XFoil's `PANE` default).
- `max_iter`: viscous iterations per angle (default 100).
- `xtrip`: forced transition `(upper, lower)` as `x/c` (default `(0.05, 0.05)`). Tripping
  near the leading edge is a valid NeuralFoil transition input and, unlike free
  transition, converges on the shrink-wrapped hinge of a deflected section.
- `ncrit`: e^N transition criticality (default 9.0, the standard clean-tunnel value).
- `mach`: Mach number (default 0.0, incompressible).
- `repanel`: repanel the input coordinates with `Xfoil.pane` before solving (default
  `false`). The [`shrink_wrap`](@ref) already emits smoothed cosine panels, and XFoil's
  own curvature-attracted repaneling re-clusters a deflected section's hinge crease —
  which NeuralFoil (analysing the smooth Kulfan fit) never sees — so `false` tracks
  NeuralFoil more closely. Set `true` to let XFoil repanel anyway.
"""
@with_kw struct XFoilSolver <: AbstractAirfoilSolver
    npan::Int = 160
    max_iter::Int = 100
    xtrip::Tuple{Float64,Float64} = (0.05, 0.05)
    ncrit::Float64 = 9.0
    mach::Float64 = 0.0
    repanel::Bool = false
end

"""
    analyze_sweep(solver::XFoilSolver, def, alpha_range, Re) -> Vector{SectionSolution}

Set the deformed coordinates once (repaneling if `solver.repanel`), then solve
`alpha_range` (radians) sweeping negative and positive angles outward from zero with a
reinit at each side for convergence. Non-converged angles yield empty Cp and `NaN`
confidence.
"""
function analyze_sweep(solver::XFoilSolver, def::DeformedSection, alpha_range, Re)
    Xfoil.set_coordinates(def.x, def.y)
    solver.repanel && Xfoil.pane(npan=solver.npan)
    sols = Vector{SectionSolution}(undef, length(alpha_range))
    neg = sort(findall(<(0.0), alpha_range); rev=true)
    pos = sort(findall(>=(0.0), alpha_range))
    for order in (neg, pos)
        reinit = true
        for ia in order
            cl, cd, _, cm, converged = Xfoil.solve_alpha(rad2deg(alpha_range[ia]), Re;
                iter=solver.max_iter, reinit=reinit, mach=solver.mach,
                xtrip=solver.xtrip, ncrit=solver.ncrit)
            reinit = false
            if converged
                xc, cp = Xfoil.cpdump()
                xu, cu, xl, cl2 = split_surfaces(xc, cp)
                sols[ia] = SectionSolution(alpha_range[ia], cl, cd, cm, 1.0, xu, cu, xl, cl2)
            else
                sols[ia] = SectionSolution(alpha_range[ia], NaN, NaN, NaN, NaN,
                                           Float64[], Float64[], Float64[], Float64[])
            end
        end
    end
    return sols
end

"""
    analyze_section(solver::XFoilSolver, def, alpha, Re) -> SectionSolution

Single-angle convenience; delegates to [`analyze_sweep`](@ref).
"""
function analyze_section(solver::XFoilSolver, def::DeformedSection, alpha, Re)
    return analyze_sweep(solver, def, [alpha], Re)[1]
end

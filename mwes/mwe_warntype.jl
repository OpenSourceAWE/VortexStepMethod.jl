# Run with: julia --project mwes/mwe_warntype.jl 2>&1 | less
# Look for red/yellow "Any" or "Union" in output

using VortexStepMethod
using VortexStepMethod: calculate_AIC_matrices!, gamma_loop!, calculate_results,
                       calculate_velocity_induced_single_ring_semiinfinite!,
                       calculate_velocity_induced_bound_2D!,
                       velocity_3D_bound_vortex!,
                       velocity_3D_trailing_vortex!,
                       velocity_3D_trailing_vortex_semiinfinite!,
                       cross3!, calc_norm_array!,
                       Panel, reinit!, solve_base!
using LinearAlgebra
using StaticArrays
using InteractiveUtils

# Setup
n_panels = 20
span = 20.0
chord = 1.0
alpha = deg2rad(30.0)

wing = Wing(n_panels, spanwise_distribution=LINEAR)
add_section!(wing, [0.0, span/2, 0.0], [chord, span/2, 0.0], INVISCID)
add_section!(wing, [0.0, -span/2, 0.0], [chord, -span/2, 0.0], INVISCID)
refine!(wing)
body_aero = BodyAerodynamics([wing])
vel_app = [cos(alpha), 0.0, sin(alpha)] .* 20.0
set_va!(body_aero, vel_app)
solver = Solver(body_aero)

va_norm_array = ones(n_panels)
va_unit_array = ones(n_panels, 3)

# Prepare args for individual functions
panel = body_aero.panels[1]
filaments = panel.filaments
ep = panel.control_point
velind = @MVector zeros(3)
tempvel = @MVector zeros(3)
va_unit = @MVector zeros(3)
va_unit .= 1.0 / sqrt(3.0)
work = body_aero.work_vectors

sep = "=" ^ 80
printstyled("\n$sep\n cross3!\n$sep\n"; color=:cyan)
a = MVector{3,Float64}(1.0, 2.0, 3.0)
b = MVector{3,Float64}(4.0, 5.0, 6.0)
r = MVector{3,Float64}(0.0, 0.0, 0.0)
@code_warntype cross3!(r, a, b)

printstyled("\n$sep\n velocity_3D_bound_vortex!\n$sep\n"; color=:cyan)
@code_warntype velocity_3D_bound_vortex!(
    velind, filaments[1], ep, 1.0, 0.001, work)

printstyled("\n$sep\n velocity_3D_trailing_vortex!\n$sep\n"; color=:cyan)
@code_warntype velocity_3D_trailing_vortex!(
    tempvel, filaments[2], ep, 1.0, 20.0, work)

printstyled("\n$sep\n velocity_3D_trailing_vortex_semiinfinite!\n$sep\n";
    color=:cyan)
@code_warntype velocity_3D_trailing_vortex_semiinfinite!(
    tempvel, filaments[4], va_unit, ep, 1.0, 20.0, work)

printstyled("\n$sep\n calculate_velocity_induced_single_ring_semiinfinite!\n$sep\n";
    color=:cyan)
@code_warntype calculate_velocity_induced_single_ring_semiinfinite!(
    velind, tempvel, filaments, ep, false, 20.0, va_unit, 1.0, 0.001, work)

printstyled("\n$sep\n calculate_velocity_induced_bound_2D!\n$sep\n"; color=:cyan)
@code_warntype calculate_velocity_induced_bound_2D!(
    velind, panel, ep, work)

printstyled("\n$sep\n calculate_AIC_matrices!\n$sep\n"; color=:cyan)
@code_warntype calculate_AIC_matrices!(
    body_aero, VSM, 0.001, va_norm_array, va_unit_array)

printstyled("\n$sep\n gamma_loop!\n$sep\n"; color=:cyan)
@code_warntype gamma_loop!(
    solver, body_aero, body_aero.panels, 0.5; log=false)

printstyled("\n$sep\n solve_base!\n$sep\n"; color=:cyan)
@code_warntype solve_base!(solver, body_aero, nothing)

printstyled("\n$sep\n calc_norm_array!\n$sep\n"; color=:cyan)
@code_warntype calc_norm_array!(solver.br.va_norm_dist, solver.sol._va_dist)

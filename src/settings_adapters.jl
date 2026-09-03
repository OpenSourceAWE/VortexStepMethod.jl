# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Turning settings into the arguments the slicer, the polar generator and the live
# polars take, so a caller passes the settings and not a dozen loose numbers — and
# so the tables and the live polars cannot be built off different ones.

using .AirfoilAero: ShrinkWrap, NeuralFoilSolver, XFoilSolver, LivePolarSettings,
                    AbstractAirfoilSolver

"""
    ShrinkWrap(mesh::MeshSettings)

The shrink wrap a mesh's slices are cleaned with.
"""
AirfoilAero.ShrinkWrap(mesh::MeshSettings) =
    ShrinkWrap(clearance = mesh.clearance,
               min_concave_radius = mesh.min_concave_radius)

"""
    NeuralFoilSolver(airfoil::AirfoilSettings)

NeuralFoil set up the way the settings ask for it.
"""
AirfoilAero.NeuralFoilSolver(airfoil::AirfoilSettings) =
    NeuralFoilSolver(model_size = airfoil.model_size, n_crit = airfoil.n_crit,
                     xtr_upper = airfoil.xtr_upper, xtr_lower = airfoil.xtr_lower)

"""
    XFoilSolver(airfoil::AirfoilSettings)

XFoil set up the way the settings ask for it: the same forced transition and the
same `n_crit` as [`NeuralFoilSolver`](@ref), so a comparison between the two
measures the network rather than a difference in how the boundary layer was posed.
"""
AirfoilAero.XFoilSolver(airfoil::AirfoilSettings) =
    XFoilSolver(ncrit = airfoil.n_crit,
                xtrip = (airfoil.xtr_upper, airfoil.xtr_lower))

"""
    LivePolarSettings(airfoil::AirfoilSettings)

The live-polar sampling the settings ask for. Live polars are a batched network
pass per solve, so they are NeuralFoil whatever `solver` names — but they read the
same network size and transition criticality the tables were generated at, which
is the point of their sharing a block.
"""
AirfoilAero.LivePolarSettings(airfoil::AirfoilSettings) =
    LivePolarSettings(offsets = deg2rad.(airfoil.live_offsets),
                      model_size = airfoil.model_size, n_crit = airfoil.n_crit)

"""
    airfoil_solver(airfoil::AirfoilSettings) -> AbstractAirfoilSolver

The section backend `airfoil.solver` names.
"""
airfoil_solver(airfoil::AirfoilSettings) =
    airfoil.solver == "xfoil" ? XFoilSolver(airfoil) : NeuralFoilSolver(airfoil)

"""
    settings_range(triple) -> StepRangeLen

A settings file's `[first, step, last]` as a range.
"""
settings_range(triple) = triple[1]:triple[2]:triple[3]

"""
    alpha_range(airfoil::AirfoilSettings) -> StepRangeLen

The angle-of-attack sweep [deg] the polars are tabulated over.
"""
alpha_range(airfoil::AirfoilSettings) = settings_range(airfoil.alpha_range)

"""
    delta_range(airfoil::AirfoilSettings) -> Union{Nothing, StepRangeLen}

The flap-deflection sweep [deg], or `nothing` for a dataset without one.
"""
delta_range(airfoil::AirfoilSettings) =
    isnothing(airfoil.delta_range) ? nothing : settings_range(airfoil.delta_range)

"""
    reynolds(set::VSMSettings, wing::WingSettings) -> Float64

`density * v_app * chord_ref / mu`, the definition the solver uses. The air comes
from `solver_settings`, which is the only place it is stated, so polars generated
at this number cannot drift from the solve that reads them.
"""
reynolds(set::VSMSettings, wing::WingSettings) =
    set.solver_settings.density * wing.airfoil.v_app * wing.airfoil.chord_ref /
    set.solver_settings.mu

"""
    rotation_matrix(mesh::MeshSettings) -> Matrix{Float64}

The mesh-to-slicer rotation, as the 3×3 the slicer takes.
"""
rotation_matrix(mesh::MeshSettings) = permutedims(reduce(hcat, mesh.rotation))

"""
    slice_args(wing::WingSettings) -> NamedTuple

The keyword arguments `obj_to_yaml` and `plot_slices_3d` share: how the mesh is
oriented, how far the outer sections stop short of the tips, and how a section is
wrapped. Splat it into either so the picture and the dataset cannot disagree.
"""
slice_args(wing::WingSettings) = (rotation = rotation_matrix(wing.mesh),
    wrap_method = ShrinkWrap(wing.mesh),
    wingtip_distance = wing.mesh.wingtip_distance)

"""
    preview_args(wing::WingSettings) -> NamedTuple

[`slice_args`](@ref) plus the two `plot_slices_3d` also takes: the section count
and the leading-edge marching resolution. `obj_to_yaml` accepts neither by that
name, which is why they are not in `slice_args`.
"""
preview_args(wing::WingSettings) = (n_slices = wing.mesh.n_sections,
    n_bins = wing.mesh.n_bins, slice_args(wing)...)

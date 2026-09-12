# The settings file

A `vsm_settings.yaml` configures a whole run in one place: the flight condition, each
wing and how it is discretised, and the numerical solver. [`VSMSettings`](@ref) reads
it, and [`Wing`](@ref), [`Solver`](@ref) and [`set_va!`](@ref) are each built from the
object it returns.

```julia
settings = VSMSettings("ram_air_kite/vsm_settings.yaml")           # under data/
settings = VSMSettings("my/vsm_settings.yaml"; data_prefix=false)  # as written

wing = Wing(settings)
body_aero = BodyAerodynamics([wing])
solver = Solver(body_aero, settings)
set_va!(body_aero, settings)
```

The file has three top-level blocks — `condition:`, `wings:` and `solver_settings:` —
and each may be left out, in which case its defaults apply. `wings:` is a list, so a
multi-wing configuration repeats the entry.

## An annotated file

Every key below is optional except a wing's `name`, `n_panels`,
`spanwise_panel_distribution`, `spanwise_direction` and `remove_nan`, and — whenever
`solver_settings:` is present — its `aerodynamic_model_type` and
`type_initial_gamma_distribution`. An omitted key keeps its default: the `condition:`
values shown are those defaults, and the docstrings linked below carry the rest.

```yaml
condition:
  wind_speed: 10.0                # free-stream velocity magnitude [m/s]
  alpha: 5.0                      # angle of attack [°]
  beta: 0.0                       # sideslip angle [°]
  yaw_rate: 0.0                   # yaw rate [°/s]

wings:
  - name: main_wing               # label the wing carries into plots and output
    # sections and polars, resolved against the working directory
    geometry_file: data/ram_air_kite/geometry.yaml
    n_panels: 50                  # panels over the span; two sections make a panel
    # LINEAR, COSINE, SPLIT_PROVIDED, UNCHANGED or BILLOWING
    spanwise_panel_distribution: LINEAR
    spanwise_direction: [0, 1, 0] # unit vector along the span, in the CAD frame
    remove_nan: true              # interpolate over NaN entries in the polar tables
    use_prior_polar: false        # reuse polars on disk instead of regenerating them
    billowing_percentage: 0.0     # trailing-edge billow, as % of arc length
    crease_frac: 0.75             # chordwise position of the deflection hinge [-]

    mesh:                         # how the sections were sliced from a CAD mesh
      obj_file: data/ram_air_kite/ram_air_kite.obj  # the mesh sections come from
      n_sections: 45              # sections sliced from the mesh
      n_bins: 60                  # leading-edge stations marched across the span
      # rows of the mesh-to-slicer rotation, whose x = chord, y = span, z = up
      rotation: [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
      wingtip_distance: 0.0       # arc length the outermost sections stop short [m]
      clearance: 0.006            # shrink-wrap offset outside the cloud [chord fraction]
      min_concave_radius: 0.02    # shrink-wrap rolling-ball radius [chord fraction]

    airfoil:                      # the 2D backend and the polars it tabulates
      solver: xfoil               # section backend: neuralfoil or xfoil
      model_size: large           # NeuralFoil network size
      n_crit: 9.0                 # e^N transition criticality; lower transitions earlier
      xtr_upper: 0.05             # forced upper-surface transition [chord fraction]
      xtr_lower: 0.05             # forced lower-surface transition [chord fraction]
      alpha_range: [-180, 1, 180] # angle-of-attack sweep [°] as [first, step, last]
      delta_range: [-40, 10, 40]  # flap-deflection sweep [°]; null for no flap sweep
      # angles off the reference angle a live polar is re-solved at [°]
      live_offsets: [-12, -9, -6, -3, 0, 3, 6, 9, 12]
      v_app: 25.0                 # apparent wind the Reynolds number is taken at [m/s]
      chord_ref: 1.0              # reference (maximum panel) chord [m]
      table_format: arrow         # per-node table format: csv or arrow

solver_settings:
  aerodynamic_model_type: VSM     # VSM or LLT
  type_initial_gamma_distribution: ELLIPTIC  # ELLIPTIC or ZEROS
  solver_type: LOOP               # LOOP or NONLIN
  density: 1.225                  # air density [kg/m³]
  mu: 1.81e-5                     # dynamic viscosity [N·s/m²]
  rtol: 1e-6                      # relative tolerance on the circulation residual [-]
  relaxation_factor: 0.01         # under-relaxation of the circulation update [-]
```

[`SolverSettings`](@ref) lists the rest of `solver_settings:`; a key left out keeps its
default. `n_panels` given there is ignored — the total is summed from the wings.

## `mesh:`

[`MeshSettings`](@ref) — which `.obj` mesh a wing's sections were cut from, and how.
The block's `obj_file` is the mesh the sections are *generated from*, an input to the
`geometry_file` the wing then flies; a wing's own top-level `obj_file` is a different
route — straight from a mesh and a `.dat`, with no polar generation — and cannot be
given alongside `geometry_file`.

[`rotation_matrix`](@ref), [`slice_args`](@ref) and [`preview_args`](@ref) turn the
block into the arguments
[`obj_to_yaml`](@ref VortexStepMethod.ObjAdapter.obj_to_yaml) and
[`plot_slices_3d`](@ref) take, and [`ShrinkWrap`](@ref) is built from `clearance` and
`min_concave_radius`. A wing naming no `mesh:` block slices exactly as an
unconfigured `obj_to_yaml` call does. [From CAD mesh to aerodynamic model](@ref)
walks through what those arguments do.

## `airfoil:`

[`AirfoilSettings`](@ref) — the 2D section backend and the polars it tabulates.
`solver:` chooses the viscous panel code (`xfoil`) or the neural surrogate
(`neuralfoil`) for a whole dataset from the file, and [`airfoil_solver`](@ref)
returns the [`XFoilSolver`](@ref) or [`NeuralFoilSolver`](@ref) it names.

One block answers for both the tables a mesh is sliced into and the live polars a
deformed section is re-solved on, so the two cannot be generated at different
transition settings or off different networks. [`alpha_range`](@ref),
[`delta_range`](@ref) and [`reynolds`](@ref) turn the sweeps and the
`density * v_app * chord_ref / mu` reference into what the polar generator takes.

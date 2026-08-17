# Changelog

## Unreleased

### Added
- `table_format` keyword on `write_section_aero`, `generate_airfoils`, `obj_to_yaml`
  and `surfplan_to_aero_yaml`: `:csv` (default, readable) or `:arrow` (binary, ~40×
  faster to load and 2.6× smaller). `read_section_aero` detects the format from the
  file suffix, so a geometry YAML can reference either.
- `convert_node_table` and `write_node_rows` rewrite a per-node table in the format
  the destination suffix names. `obj_to_yaml` migrates an existing dataset with
  them when `table_format` differs from what the directory holds, so a dataset
  changes format without re-running the airfoil solver that produced it.
- `flow_curvature` solver setting (default `false`). When enabled, each section
  gets the thin-airfoil pitch-rate moment increment `Δcm = -(π/4) q̂` with
  `q̂ = q c / (2 v_rel)`. A section rotating about its own spanwise axis sees an
  incidence that varies linearly along the chord, which is equivalent to
  parabolic camber and produces a quarter-chord moment that a single control
  point cannot represent. The lift response to `q` was already exact because the
  inflow is sampled at the three-quarter-chord point, so only the moment was
  missing.
- `pitch_rate_dist` field on `BodyAerodynamics`, holding the `q` above per panel.
  `set_va!(body_aero, va, omega)` fills it by projecting `omega` onto every
  panel's `y_airf`, so panels at different dihedral see different rates from one
  body rate. The distributed `set_va!(body_aero, va_distribution;
  pitch_rate_dist)` takes it directly, so twist and flapping rates of a deforming
  wing — which no single body rate can express — reach the moment. Omitting the
  keyword zeroes it rather than reusing a stale `omega`, which the distributed
  form never sets.
- `section_pitch_rate(velocity_leading, velocity_trailing, z_airf, chord)` builds
  one entry of that distribution from a section's edge velocities, and reduces to
  `ω ⋅ y_airf` for rigid motion.
- `bin/update_default_manifests` regenerates both checked-in default manifests in
  one command, and `bin/install` takes `+X.Y` / `--version X.Y` to pick a Julia
  channel without prompting. Selecting a version no longer moves the juliaup
  default, so regenerating the 1.11 manifest leaves the shell on whatever channel
  it was using.

### Fixed
- The `Solver` docstring quoted the `SolverSettings` defaults for
  `type_initial_gamma_distribution` and `core_radius_fraction` (`ELLIPTIC`, `1e-20`)
  instead of its own (`ZEROS`, `0.05`), and never said what `core_radius_fraction`
  measures. It now documents the `Solver` defaults and cites Damiani et al. (2019) for
  the 0.05 cut-off.

### Changed
- `SolverSettings` now defaults to the same values as `Solver`: `core_radius_fraction`
  `1e-20` → `0.05` and `type_initial_gamma_distribution` `ELLIPTIC` → `ZEROS`. The 0.05
  bound vortex core cut-off follows Damiani et al. (2019), "A Vortex Step Method for
  Nonlinear Airfoil Polar Data as Implemented in KiteAeroDyn", and matches the upstream
  `awegroup/Vortex-Step-Method` default; at `1e-20` the Biot-Savart singularity guard
  never engaged. Coefficients are unchanged for well-separated geometry, since the guard
  only triggers where a control point falls within 5% of a filament length of a bound
  vortex, and every settings file shipped in `data/` sets both keys explicitly.
- `BodyAerodynamics.AIC` is stored as `(n_panels, n_panels, 3)` instead of
  `(3, n_panels, n_panels)`, so each component slice `AIC[:, :, k]` is contiguous and
  the induced-velocity products reach BLAS `gemv` instead of the generic fallback.
  `solve!` is 3.1–5.3× faster (n=120, VSM: 26.1 ms → 4.9 ms inviscid, 21.1 ms →
  6.0 ms with polars). Code reading `AIC[k, i, j]` must become `AIC[i, j, k]`.
- The filament induced-velocity kernels reuse the `r0`/`length` each `BoundFilament`
  already stores, take the core-radius cutoff from `|r1.r0|/|r0|` without forming the
  perpendicular vector, and defer the cross products that only one branch reads. Output
  is bit-identical; `solve!` is a further 1.47-1.55x faster.
- `read_node_table` parses into a preallocated matrix instead of `reduce(vcat, …)`
  over a generator, which was quadratic in the row count: ~21× faster on a 16 MB
  surface table (2.49 s → 0.12 s), benefiting every existing dataset.

## VortexStepMethod v4.0.0 2026-08-03

### Added
- `NeuralFoil`-based airfoil polar generation via the new `AirfoilAero` submodule
  (`NeuralFoilSolver`, `XFoilSolver`, `analyze_section`, `analyze_sweep`,
  `fit_kulfan_parameters`, `shrink_wrap`, `ShrinkWrap`)
- `ObjAdapter` submodule: converts a 3D wing `.obj` mesh to the native YAML/CSV
  geometry format (`obj_to_yaml`, `perpendicular_sections`,
  `write_yaml`, `plot_slices_3d`, `plot_airfoils`)
- `SurfplanAdapter` submodule: `surfplan_to_aero_yaml` turns a SurfplanAdapter aero
  export into the native pressure-ready YAML/CSV geometry (shared `generate_airfoils`
  core with `ObjAdapter`)
- `ObjWing(obj_path[, dat_path]; Re, n_panels, aero_solver, remake, ...)` convenience
  constructor restored for backward compatibility — internally calls
  `ObjAdapter.obj_to_yaml` then `Wing`; `aero_solver` selects the polar backend
  (default `NeuralFoilSolver()`; pass `XFoilSolver()` for old behavior);
  `remake=false` (default) reuses an existing `geometry.yaml` in `output_dir` to skip
  expensive polar generation when only `n_panels` changes; set `remake=true` to force
  regeneration
- Per-section surface aero tables: `SectionAero` (full contour + surface pressure
  `cp` + skin friction `cf` per node over an `(α, δ)` grid), `section_surface`,
  `read_section_aero` / `write_section_aero`; surface aero propagates through
  `Section` and `Wing` and is spanwise-interpolated during `refine!`
- `POLY` aero model for polynomial cl/cd/cm (exported)
- `lei_poly_coeffs(tube_diameter, camber)` (exported from `AirfoilAero`) returning
  the Breukels α-polynomial cl/cd/cm coefficients
- Li/Gaunaa spanwise artificial viscosity (Li, Gaunaa, Pirrung & Lønbæk,
  TORQUE 2026) for the LOOP solver, stabilizing post-stall circulation
  distributions that otherwise develop non-physical sawtooth oscillations;
  opt-in via `is_with_artificial_viscosity` (default `false`) and
  `artificial_viscosity_factor` (default `0.035`) on the solver settings
- `crease_frac` wing setting (chordwise flap-hinge fraction, default `0.75`)
  for drawing the δ-deflected plate/skin, readable from YAML
- `examples/V3_neuralfoil.jl` and `examples/obj_to_yaml_kite.jl`

### Changed
- BREAKING - `ObjWing` polar generation now uses NeuralFoil (via
  `ObjAdapter.obj_to_yaml`) instead of XFoil + a user-supplied `.dat` file. The
  constructor signature is backward-compatible (`dat_path` is accepted and
  silently ignored), but the resulting aerodynamic polars will differ
  numerically from v3.x outputs. To reproduce old XFoil-based polars pass
  `aero_solver=XFoilSolver()` to `ObjAdapter.obj_to_yaml` and call `Wing`
  directly.
- OBJ-based wings are now built via `ObjAdapter.obj_to_yaml` + `Wing(yaml_path)`
  instead of the old `ObjWing` pipeline (XFoil + single `.dat` file);
  `ObjWing` is kept as a shim that accepts but ignores `dat_path`
- BREAKING - `LEI_AIRFOIL_BREUKELS` is now a deprecated alias of `POLY`, and the
  built-in Breukels regression no longer runs at solve time. Sections must carry
  the α-polynomial coefficients directly instead of `(tube_diameter, camber)`;
  compute them up front with `lei_poly_coeffs(tube_diameter, camber)`
- `Section` and `add_section!` accept an optional `section_aero` argument
- BREAKING - plotting is now Makie-only; the `VortexStepMethodMakieExt` extension
  loads once a Makie backend and
  [`MakieControlPlots`](https://github.com/OpenSourceAWE/MakieControlPlots.jl) are
  available, so plotting code must now load `MakieControlPlots` (and drop
  `set_plot_backend!`); `plot_section_polars` is rendered through `MakieControlPlots`

### Removed
- `ObjWing` as a standalone pipeline (replaced by `ObjAdapter`); the name is
  re-exported as a compatibility wrapper
- `auto_rotation` helper (internal, removed from public API)
- `PanelGroupingMethod` enum (already removed in v3.0.0 — stale docs entry cleaned up)
- ControlPlots test run removed from CI (`plot-controlplots` arg) due to
  a `libraqm`/HarfBuzz symbol conflict in the GitHub Actions environment
- BREAKING - the `ControlPlots` plotting extension, `examples_cp/`, and the
  `PythonCall`/Matplotlib setup (`bin/install_controlplots`, CondaPkg
  `LocalPreferences` defaults)
- BREAKING - the `PlotBackend`/`MakieBackend`/`ControlPlotsBackend` types and
  `set_plot_backend!`; plotting works as soon as a Makie backend and
  `MakieControlPlots` are loaded
- the never-implemented `plot_circulation_distribution`

## VortexStepMethod v3.3.6 2026-06-13

### Added
- `calc_forces!` and `solve_base!` (both exported): `solve!` is now
  `solve_base!` followed by `calc_forces!`, so a frozen circulation can be
  mapped to forces without re-running the nonlinear gamma solve (#245)
- `calculate_cd` and `calculate_cm`, splitting the combined `calculate_cd_cm`
  into separate drag- and moment-coefficient functions; `calculate_cd_cm` is
  kept as a thin wrapper (#246)

### Changed
- `calc_forces!` is now allocation-free in the per-step hot path
  (preallocated `panel_area_dist` and `unrefined_count_dist` buffers) (#245)

### Fixed
- 3D polar plotting (#245)
- flaky Aqua `persistent_tasks` test now actually disabled via
  `persistent_tasks=false` (`()` did not disable it) (#246)

## VortexStepMethod v3.3.5 2026-06-05

### Added
- `moment_coeff_unrefined_dist` field in `VSMSolution`: the summed
  `moment_frac`-referenced pitching-moment coefficient per unrefined section [-]

## VortexStepMethod v3.3.4 2026-05-31

### Added
- `PlotBackend`, `MakieBackend`, `ControlPlotsBackend`, and
  `set_plot_backend!` so applications can explicitly choose which plotting
  extension the backend-agnostic plotting API should use
- `PythonCall` added as a weak dependency to support the `ControlPlots`
  backend with PythonPlot

### Changed
- backend-agnostic plotting wrappers now route through the active plotting
  backend, and each plotting extension initializes itself as the default only
  when no backend has been selected yet
- relaxed `ControlPlots` compatibility to include both `0.2.5` and `0.3`
- improved `bin/install` and `bin/install_controlplots` scripts

### Fixed
- corrected projection onto core radius in `velocity_3D_bound_vortex!` and
  semi-infinite trailing vortex projection so that the radial direction is
  always measured from the filament axis, not from the origin (#241)
- fixed 0-based subplot indexing in `ControlPlotsExt` for PythonPlot
  compatibility (`plot_distribution` no longer errors with PythonPlot backend)
- fixed missing initialization of `damp` in `solver.jl`

## VortexStepMethod v3.3.3 2026-05-21

### Fixed
- `MakieExt` and `ControlPlotsExt` no longer both define
  `VortexStepMethod.plot_geometry` for the same type, resolving a method
  ambiguity when both extensions were loaded (#236)
- `menu()` and `menu_cp()` now always set the active backend before dispatching
  to a plot function, preventing stale-backend errors

## VortexStepMethod v3.3.2 2026-05-18

### Changed
- use 2-arg version of atan to avoid possible NaN

## VortexStepMethod v3.3.1 2026-05-13

### Changed
- `unrefined_deform!` linearly interpolates twist and TE deflection between
  unrefined sections and rotates each refined section about the average of its
  adjacent local airfoil normals (#234)

### Fixed
- `smooth_sqrt` in the solver hot loop keeps gradients defined at zero
  velocity magnitude

## VortexStepMethod v3.3.0 2026-05-05

### Added
- `ForwardDiff` compatibility, used by default in `linearize` (#232)
- `backend` keyword argument for `linearize`
- example `linearize_check.jl` comparing FiniteDiff and ForwardDiff tangents

### Changed
- core structs are parameterized on the scalar type `T` so dual numbers can
  propagate through them; public constructors are unchanged

## VortexStepMethod v3.2.0 2026-05-02

### Added
- support for both CairoMakie and GLMakie in the example `menu()` via new
  `CairoMakie.activate()` / `GLMakie.activate()` entries (#222)
- plots are saved as PDF when CairoMakie is active
- `bin/install` now creates the `output` folder used by examples to store
  generated plots

### Changed
- examples write plots into the shared `output` folder instead of the working
  directory
- example plot file names are sanitized: spaces replaced with `_` and `%`
  replaced with `pct`

### Fixed
- NONLIN solver no longer returns stale `gamma` on repeated `solve!` calls
  (#228)

## VortexStepMethod v3.1.3 2026-04-23

### Fixed
- bug in `linearize` where `body_aero.va` was used instead of `body_aero._va`,
  causing incorrect initial velocity storage (#227)
- `set_va!` no longer overwrites `_va` with a computed reference velocity when
  omega is nonzero; this simplifies the function and fixes linearization with
  turn rates
- linearize now correctly preserves and restores `omega` across perturbations

## VortexStepMethod v3.1.2 2025-04-20
### Added
- add back compat entry v2 for SciMLBase

## VortexStepMethod v3.1.1 2025-04-20
### Added
- add back compat entry v3 for RecursiveArrayTools

## VortexStepMethod v3.1.0 2025-04-19

### Breaking
- `billowing_angle` replaced by `billowing_percentage` on `Wing` and
  `WingSettings` (percentage of arc length, not radians)
- `billowing_angle_from_percentage()` removed
- `BILLOWING` distribution now uses sinusoidal rotation instead of circular arc

### Added
- `billowing.jl` example comparing flat vs billowed V3 kite
- Coarse V3 kite geometry, settings, and combined CFD polar data
- `cl_over_cd` keyword for `plot_polars` and `plot_combined_analysis`
- the function `menu_cp()` can now be used to run the examples with the ControlPlots backend
- the script `bin/install` can and should be used to instantiate the project and all sub-projects after git checkout
- the files `Manifest-v1.11.toml.default` and `Manifest-v12.toml.default` for enhanced reproducibility
- the scripts `bin/jetls` and `bin/jetls_examples` for running a static check of the source code
- the script `test_bench.jl` for benchmarking the refinement method and measuring allocations

### Changed
- Plot legends moved to shared horizontal legend at bottom of grid layouts
- the script bin/run_julia now can be called with a script name as parameter
- fixed all JETLS warnings in the source code for improved performance and stability

### Fixed errors
- Domain error in elliptical gamma distribution when control points lie
  outside the nominal span envelope

## VortexStepMethod v3.0.1 2025-04-04

### Changed
- the file `CITATION.cff`
- compat entry for RecursiveArrayTools
### Added
- the file `.zenodo.json`

## VortexStepMethod v3.0.0

### Breaking Changes

#### Unified Wing type

- `Wing` and `RamAirWing` merged into a single `Wing` struct
- `RamAirWing` constructor replaced by `ObjWing` (returns a `Wing` with
  OBJ-specific fields populated)
- Renamed `ram_geometry.jl` → `obj_geometry.jl`

#### Unrefined sections replace panel grouping

- `n_groups` replaced by `n_unrefined_sections` in `Wing`, `WingSettings`, and
  `SolverSettings`
- `PanelGroupingMethod` enum removed
- Each refined panel maps back to its unrefined section via
  `refined_panel_mapping`
- `n_groups` and `grouping_method` in YAML settings files emit a deprecation
  warning and are ignored

#### Renamed VSMSolution fields (`_array` → `_dist`)

- `panel_width_array` → `width_dist`
- `alpha_array` → `alpha_dist`
- `cl_array` → `cl_dist`, `cd_array` → `cd_dist`, `cm_array` → `cm_dist`
- `panel_lift` → `lift_dist`, `panel_drag` → `drag_dist`
- `panel_moment` → `panel_moment_dist`
- `group_moment_dist` → `moment_unrefined_dist`
- Type parameter `Solver{P,G}` / `VSMSolution{P,G}` changed to `Solver{P,U}` /
  `VSMSolution{P,U}` (U = unrefined sections)

#### Mesh refinement separated from reinit

- `reinit!` no longer re-refines the mesh; it only updates panel properties from
  existing `refined_sections`
- New `refine!` function handles mesh refinement explicitly

### Changed

- Consistent naming convention: `_dist` = per-panel, `_unrefined_dist` =
  per-unrefined-section
- `correct_aoa` solver setting added (default `false`)
- `alpha_geometric` is now the angle between apparent wind and chord
- Panel `width` calculation uses sum instead of average
- `src/plotting.jl` removed; all plotting moved to Makie and ControlPlots
  extensions
- `TestEnv` no longer required as a global dependency

### Added

- `n_unrefined_sections` field in `Wing` for tracking pre-refinement sections
- `refined_panel_mapping` for automatic panel-to-section association
- Unrefined distribution fields in `VSMSolution`: `cl_unrefined_dist`,
  `cd_unrefined_dist`, `cm_unrefined_dist`, `alpha_unrefined_dist`,
  `moment_unrefined_dist`, `width_unrefined_dist`
- `alpha_geometric_dist` field in `VSMSolution`
- `unrefined_deform!` for applying twist and TE deflection to unrefined
  sections; `non_deformed_sections` preserved for reset/re-deformation
- `refine!` function for explicit mesh refinement
- `plot_combined_analysis` in the Makie extension
- Literature comparison plotting in the Makie extension
- `obj_file` and `dat_file` fields in `WingSettings`
- `sort_sections` kwarg for section ordering control
- Separate `test/Project.toml` and `examples/Project.toml`
- New tests: `test_unrefined_dist.jl`, `test_refinement_validation.jl`,
  `test_yaml_wing_deformation.jl`

### Removed

- `RamAirWing` type (use `ObjWing` or `Wing` directly)
- `PanelGroupingMethod` enum
- `n_groups` from `WingSettings`, `SolverSettings`, and YAML files
- `src/plotting.jl` (920 lines moved to extensions)

## VortexStepMethod v2.3.0 2025-10-16

### Added

- A Makie plotting extension.

## VortexStepMethod v2.2.1 2025-10-07

### Fixed

- Reference in vsm_settings.yaml now points to aero_geometry.yaml.

### Added

- None.

### Removed

- Bridle configuration data (nodes, lines, connections).

### Changed

- Renamed wing_geometry_polars_CFD.yaml → aero_geometry.yaml.

## VortexStepMethod v2.2.0 2025-08-28

### Added

- The kwarg `aero_coeffs` to the function `linearize`: if true the linearization
  will output normalized coefficients instead of moments and forces.

## VortexStepMethod v2.1.0 2025-08-11

### Changed

#### 1. Core New Functionality: YAML Geometry Support

- **New file:** `yaml_geometry.jl` (290+ lines) — Complete YAML-based wing
  geometry loading.
- **New structs:** `WingSectionData`, `WingAirfoilData`, `WingAirfoilInfo` (with
  `@with_kw` macros).
- **New function:** `load_polar_data()` — Robust CSV polar data loading with
  error handling.
- **New constructors:**
  - `Wing(geometry_file::String)` — Create wings from YAML files.
  - `Wing(settings::VSMSettings)` — Create wings from settings.

#### 2. Enhanced Settings System

- Renamed `vs()` to `VSMSettings()` constructor.
- Added convenience `Solver(body_aero, settings)` constructor.
- Improved settings structure and validation.

#### 3. Comprehensive Test Infrastructure

- **Split tests:** Reorganized from monolithic files to modular structure
  (`test/module_name/test_*.jl`).
- **New test module:** `yaml_geometry` with `test_load_polar_data.jl` and
  `test_wing_constructor.jl`.
- **Test utilities:** `test_data_utils.jl` with shared helper functions.
- **Test data:** Extensive YAML and CSV test files in `data`.

#### 4. Data and Examples

- **New data sets:** Complete `TUDELFT_V3_KITE` with CFD polars and literature
  results.
- **Enhanced examples:** Updated examples to use YAML geometry (e.g.,
  `pyramid_model.jl`, `V3_kite.jl`).
- **Real-world configs:** Production-ready YAML geometry files for various kite
  configurations.

#### 5. Core Module Improvements

- **Path handling:** Robust file path resolution for relative/absolute paths.
- **Error handling:** Comprehensive validation and graceful fallback to INVISCID
  mode.
- **Memory management:** Improved file I/O and cleanup in tests.
- **Documentation:** Added comprehensive docstrings and examples.

## VortexStepMethod v2.0.0 2025-07-11

### Changed

- bump Interpolations to 0.16
- breaking: rename `init!` to `reinit!`

## VortexStepMethod v1.2.6 2025-04-30

### Changed

- update NonlinearSolve.jl

## VortexStepMethod v1.2.5 2025-04-18

### Changed

- suppress `@info` messages when creating a RamAirWing
- improve examples `ram_air_kite.jl` and `bench.jl`

## VortexStepMethod v1.2.4 2025-04-17

### Changed

- implement export of `solve` and `solve!` correctly
- do not export `menu()` because KiteUtils exports it
- the polars of the ram air kite are now saved in `.csv` files
- update CI.yml

## VortexStepMethod v1.2.3 2025-04-13

### Changed

- expose the angle of attack `alpha_array` in the `VSMSolution` #167

## VortexStepMethod v1.2.2 2025-04-13

### Added

- the parameter `prn`to the constructor `RamAirWing` which allows it suppress
  info messages

## VortexStepMethod v1.2.1 2025-04-08

### Added

- Add back `bench2.jl` and rename it to `bench_solve.jl` #150

### Removed

- Remove problematic parallel Xfoil computing and use single thread instead #161

## VortexStepMethod v1.2.0 2025-03-27

### Changed

- The fields that had as type a `Matrix of size Px1` have now the type `Vector`
- Many new fields of the type `VSMSolution` documented
- `reinit!(body_aero)` is now a public function

### Added

- Added the option to use nonlinear solve to calculate the gamma distribution
  #140
- New page `Tips and tricks` added to the documentation
- Fast and modular linearization added around an operating point #140

## VortexStepMethod v1.1.2 2025-03-23

### Added

- The function `install_examples()` which allows to easily install the examples
  without using `git`
- The function `solve!` returns a struct now. The function `solve`that returns a
  dict is still available.
- The moment coefficients distribution in `solve!`
- The script `install` to the `bin` folder for users who checked out this git
  repository
- The script `bench2.jl` was added for allocation testing of the `solve!`
  function

### Changed

- Read the y-coordinates in the correct direction from the
  `ram_air_kite_body.obj` file
- In the `menu.jl`, changed `help` to `help_me`. It works better now, no more
  warnings on Linux, it should also work on MacOS now
- The coordinate frames of the panels now use the same convention as the kite
  body frame
- The page "Glossary" of the documentation is quite complete now
- The center of mass field of the `RamAirWing` is removed, and the geometry is
  created such that `[0, 0, 0]` is the center of mass
- The enumeration `WingType` was added and replaces the symbols, used before
- The allocations of the function `solve!` where reduced by a factor of 11 to 85
  allocations
- `align_to_principal` option added to `RamAirWing`
- `deform!` by a distribution instead of just a left and right angle

### Fixed

- The function `calculate_circulation_distribution_elliptical_wing()` was never
  called
- Fix the calculation of force coefficients in `solve!`
- The continues integration scripts (CI.yml) use now separate runs for the test
  coverage and for the allocation tests.

## VortexStepMethod v1.1.0 2025-03-04

### Added

- Dynamically deform the RamAirWing by twisting the left side and right side,
  and deforming the trailing edges using deform! #19
- Set turn rate `omega = [omega_x, omega_y, omega_z]` in kite body frame using
  set_va! #49
- Add moment coefficient calculations around specified point to `solve` #17
- Add moment distribution of the moment around the local panel y-axes around
  user-defined points on the panels to `solve!` #90
- Add function `solve!()` which returns a `VSMSolution` struct #87
- Add the option to remove the NaNs in `aero_data` vectors or matrices using the
  `remove_nan` keyword in the `Wing` and `RamAirWing` constructors #98

### Changed

- Add origin argument to `BodyAerodynamics` constructor #66
- Improve documentation

## Initial Release

This project is based on version 1.0 of the Python project
[Vortex-Step-Method](https://github.com/ocayon/Vortex-Step-Method)

## Noteworthy Differences of v1.0.1 to the Python version

- implemented in Julia, therefore about 50 times faster
- an importer for `.obj` wing geometry files was added (see:
  [#10](https://github.com/OpenSourceAWE/VortexStepMethod.jl/issues/10))
- automatic creation of polars using Xfoil.jl was added (see:
  [#43](https://github.com/OpenSourceAWE/VortexStepMethod.jl/pull/43))
- a ram-air kite example was added
- `Umag` was replaced with `v_a` as variable name for the norm of the apparent
  wind speed
- memory allocations were significantly reduced
- a menu (examples/menu.jl) for running the examples was added
- plotting was moved to an extension #55
- added improved online documentation

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

# Architecture: `lib/` sub-packages + standard-format core

## Principle

The VSM **core** ingests exactly one format — a standard geometry YAML plus standard
polar CSVs — and depends on nothing heavy. Everything that *produces* that format is
a separate `lib/` sub-package. Two generators sit above core in the DAG:

```
ObjAdapter        obj (+ future filetypes) → standard geometry YAML,
   │  depends on   calls AirfoilAero for the polar CSVs
   ▼
AirfoilAero       2D airfoil solvers → cl/cd/cm/cp tables + poly coeffs
   │  writes         (xfoil, neuralfoil, poly/breukels; future cfd, masure)
   ▼  csv / yaml
              ─────►  VSM core   (reads interp-or-poly from csv; solves)
```

Core depends on **neither** at runtime. `Xfoil`, the NeuralFoil weights, and `Kulfan`
leave core's `Project.toml` entirely. No `Wing(obj_file)` / on-the-fly generation
convenience in core — those live in the generators.

## Repository layout (SciML `lib/` monorepo pattern)

Real sub-packages (own `Project.toml`, UUID, `src/`, `test/`), unregistered, wired
with `[sources]` (Julia ≥ 1.11). Extraction later = move the folder and register.

```
VortexStepMethod.jl/
  Project.toml            # core — sheds Xfoil, NeuralFoil, Kulfan
  src/                    # solver, wing_geometry, yaml/csv readers, interp+poly eval
  lib/
    AirfoilAero/
      Project.toml        # owns Xfoil, NeuralFoil deps
      src/AirfoilAero.jl
      test/
    ObjAdapter/
      Project.toml        # [sources] AirfoilAero = {path=../AirfoilAero}
      src/ObjAdapter.jl
      test/
```

Top `Project.toml` gains only for tests/examples that exercise generators:

```toml
[sources]
AirfoilAero = {path = "lib/AirfoilAero"}
ObjAdapter  = {path = "lib/ObjAdapter"}
```

Core `src` itself has no dependency on either.

## What goes where

### VSM core (`src/`)
- `solver.jl`, `body_aerodynamics.jl`, `wake.jl`, `filament.jl`, `panel.jl`,
  `wing_geometry.jl`, `settings.jl`.
- Standard-format **readers**: `load_polar_data`, `read_aero_matrix`, `read_cp_data`.
- Cp runtime: `cp_types.jl` (`CpData`, `CpPolar`), `CpPolar(::CpData)`,
  `cp_distribution`, `delta_cp`.
- Shared util: `interpolate_matrix_nans!` (moved out of `polars.jl` into `utils.jl`).
- Aero-model evaluation for the two core models: **INTERP** (tables) and **POLY**
  (evaluate polynomial from stored coeffs). Generic; core does not know who made the
  coeffs or tables.

### AirfoilAero (`lib/AirfoilAero`)
General airfoil-aero I/O: **`.dat` in → up to four CSVs out** (`cl`/`cd`/`cm`/`cp`).
No `cp.csv` for the poly/breukels backend (a polynomial model has no meaningful
chord-slice Cp). Backends are interchangeable: `xfoil`, `neuralfoil`, `poly`
(future `cfd`, `masure`).
- `kulfan.jl` (CST parameterization), `neuralfoil.jl` (+ weights),
  `airfoil_solvers/` (`AbstractAirfoilSolver`, `XFoilSolver`, `NeuralFoilSolver`),
  the deform helpers (`turn_trailing_edge!`, `get_lower_upper`, `deform_section`).
- Polar **generation**: `create_polars` (XFoil sweep), `generate_cp_polar`; writes
  via core's `write_aero_matrix` / `write_cp_data`.
- **poly** (formerly Breukels): `compute_lei_coeffs` becomes a *generator* that emits
  the α-polynomial coefficients for a section's tube-diameter/camber to CSV.
- Owns `Xfoil`, NeuralFoil deps. Depends on core only for the CSV/`CpData` IO types.

### ObjAdapter (`lib/ObjAdapter`)
- `.obj` (and future filetypes) mesh reading (`read_faces`,
  `create_interpolations`), `obj_slice.jl` (`perpendicular_sections`),
  `obj_to_yaml.jl`.
- **Generates** `.dat` airfoil files from the mesh slices (it does *not* take `.dat`
  as input) plus the standard geometry YAML; the `.dat` files then feed AirfoilAero
  to produce the polar CSVs.
- Orchestrates: mesh → perpendicular slices → `.dat` → AirfoilAero polars/cp/poly →
  standard YAML + CSVs. Output is files, not `Wing` objects.
- Depends on AirfoilAero. Does **not** import core `Wing`/`Section` internals; the
  live `ObjWing` runtime type is dropped (convert-then-load via `Wing(yaml)`).

## Generated artifacts out of git

Pipeline-generated polars/dat/cp CSVs (e.g. `data/*/polars_neuralfoil/`,
`data/*/polars_CFD_NF_combined/`, generated `.dat`) are **regeneratable** and should
not be tracked: untrack + `.gitignore`, and scrub from history (destructive — needs
explicit confirmation and a force-push). Reference/experimental data
(`literature_results/`, measured `polars_CFD/`) stays tracked.

## Breukels → poly

Core keeps two aero models: **INTERP** and **POLY**. `LEI_AIRFOIL_BREUKELS` becomes
`POLY`. The *derivation* `(tube_diameter, camber) → α-polynomial coeffs`
(`compute_lei_coeffs`, today in `panel.jl`) moves to AirfoilAero. A section then
stores the **coefficients** (a short vector per cl/cd/cm), not `(d_tube, camber)`.
`init_aero!`/`calculate_cl/cd/cm` in core evaluate `poly(α)` from those coeffs — a
tiny generic evaluator, unchanged in spirit but fed from data instead of computed
in-solver.

`AeroData` union in core collapses toward: `Nothing` (inviscid), poly-coeff tuple
(POLY), vector tables (INTERP/`POLAR_VECTORS`), matrix tables (`POLAR_MATRICES`).
Naming: rename the enum values to `INTERP_*`/`POLY` in a follow-up; keep the existing
`POLAR_VECTORS`/`POLAR_MATRICES` names working during migration to limit churn.

## Cp fit into the split

Already implemented and unchanged by the split:
- **Core reader/eval**: `cp_types.jl`, `read_cp_data`, `Section.cp_data`,
  `Panel.cp_polar`, `init_aero!` build, `VSMSolution.cp_*`, `prepare_cp_output!`,
  `delta_cp(sol)`, all-or-none + uniform-`n_chord` validation.
- **AirfoilAero generator**: `generate_cp_polar`, `write_cp_data`, the airfoil solvers
  and `deform_section` that feed them.

Only relocation is needed: move `generate_cp_polar`/`write_cp_data` and the airfoil
solvers into AirfoilAero; keep `read_cp_data`/`CpData`/`CpPolar`/`cp_distribution` in
core.

## `polars.jl` disposition

- `interpolate_matrix_nans!` → core `utils.jl` (core `wing_geometry` depends on it).
- `read_aero_matrix` → core (a reader).
- `create_polars`, `write_aero_matrix`, `turn_trailing_edge!`, `get_lower_upper`,
  XFoil helpers → AirfoilAero (generation).

## Migration order (staged, verify load after each)

1. **Core util split.** Move `interpolate_matrix_nans!` into `src/utils.jl`; leave
   `read_aero_matrix` in core. Core still loads.
2. **AirfoilAero package.** Create `lib/AirfoilAero` (Project.toml + module); move
   `kulfan.jl`, `neuralfoil.jl`, `airfoil_solvers/`, generation half of `polars.jl`,
   `generate_cp_polar`/`write_cp_data`. Drop them from core includes/exports/deps.
   Verify AirfoilAero loads standalone and core loads without Xfoil/NeuralFoil/Kulfan.
3. **Breukels → poly.** Move `compute_lei_coeffs` to AirfoilAero; add core POLY eval
   from stored coeffs; migrate section data + loaders + `LEI_AIRFOIL_BREUKELS` name.
4. **ObjAdapter package.** Create `lib/ObjAdapter`; move mesh/slice/dat/obj_to_yaml;
   complete obj→standard-yaml round-trip through the core loader; drop `ObjWing`.
5. **Tests/examples.** Repoint `ObjWing`/generation callers to the sub-packages;
   move generator tests under each `lib/*/test`.

## Open risks

- **obj→yaml completeness.** `obj_to_yaml` must emit a YAML the core loader fully
  reconstructs (unrefined sections at slice resolution, polar/cp/poly CSV refs). Any
  gap surfaces in stage 4.
- **Enum rename churn.** Deferred; keep old names as aliases through migration.
- **`[sources]` + CI.** Ensure test/docs projects declare the `[sources]` paths.

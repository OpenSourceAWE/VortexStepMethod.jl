# Plan: align the standard geometry format to awesIO

## Source

`github.com/awegroup/awesIO` (IEA Wind Task 48, Airborne Wind Europe glossary) is a
JSON-Schema YAML validation standard for AWE systems. Formal schemas exist for
system / power_curves / wind_resource / operational_constraints. Geometry is not yet
a formal schema, but `examples/original_config_files/TUDELFT_V3_KITE/aero_geometry.yaml`
is the reference format — and it already matches our `wing_sections`/`wing_airfoils`
`headers`/`data` layout. We adopt it instead of a parallel format.

## Key conceptual point: type = solver, not airfoil

An airfoil is a **shape** (`.dat`), optionally plus regression parameters. The
awesIO `wing_airfoils` `type` field (`neuralfoil | breukels_regression |
masure_regression | polars | inviscid`) is really **which solver/method produces the
polar** — an `AbstractAirfoilSolver` concern in AirfoilAero, not an intrinsic airfoil
property. The **core never sees the solver type**: a preprocessing/resolve step turns
the rich types into one of three core data representations:

| awesIO `type` | resolves to (core) | resolver |
|---|---|---|
| `polars` (csv: α[rad],cl,cd,cm) | `POLAR_VECTORS` | (already resolved) |
| `inviscid` | `INVISCID` | (already resolved) |
| `breukels_regression` (`t`, `kappa`) | `POLY` coeffs | `AirfoilAero.lei_poly_coeffs` |
| `neuralfoil` (`dat_file_path`, `model_size`, `n_crit`, `xtr_*`) | `POLAR_VECTORS` csv | `AirfoilAero` NeuralFoil |
| `masure_regression` (`t,eta,kappa,delta,lambda,phi`) | `POLY`/csv (future) | `AirfoilAero` (future) |

So: **awesIO rich yaml = converter input; core reads the resolved polars/poly/inviscid
subset.** This is exactly the converter/core split already built.

## Format deltas to adopt

### `wing_sections`
awesIO headers: `[airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z, VUP_x, VUP_y, VUP_z]`.
- **Add the per-section up-vector `VUP`** (section orientation / through-thickness
  direction). Ours currently stops at TE_z.
- Core loader: parse `VUP` when the header includes it; store on `Section` (new
  optional field, e.g. `up_vector`). Backward compatible — absent ⇒ derive as today.
- Relevance: the up-vector is the natural anchor for the panel frame and for the
  `thickness_frac` through-thickness pivot in `deform_section`.

### `wing_airfoils`
awesIO headers: `[airfoil_id, type, info_dict]` + top-level `alpha_range`, `reynolds`.
- info_dict fields by type (see table): `dat_file_path`, `csv_file_path`, `chord`,
  `is_strut`, and the regression params.
- Our current loader reads `type` + `info_dict.csv_file_path` → `POLAR_VECTORS`. Extend
  to the full taxonomy, but keep the resolve in the converter layer.

## Work items

1. **Core loader (`yaml_geometry.jl`)**: accept the awesIO `wing_sections` header with
   `VUP_*`; parse `alpha_range`/`reynolds` if present; read `type ∈ {polars, inviscid}`
   and `poly` (resolved). Add optional `Section.up_vector`.
2. **AirfoilAero/ObjAdapter resolve step**: `resolve_aero_geometry(yaml) -> yaml'` that
   turns `neuralfoil`/`breukels_regression`/`masure_regression` entries into `polars`
   (csv) or `poly` (coeffs), writing the CSVs/`.dat` and a core-ready yaml. Reuse
   `generate_polar_from_dat`, `lei_poly_coeffs`.
3. **ObjAdapter `obj_to_yaml`**: emit the awesIO layout (VUP per section from the slice
   normal; `type: neuralfoil` or resolved `polars`).
4. **Cp**: awesIO has no surface-pressure schema. Our chord-slice `cp.csv` stays our
   own extension for now; propose upstreaming (see `docs/awesio_cp_issue_proposal.md`).
5. **Validation**: round-trip the TUDELFT_V3_KITE `aero_geometry.yaml` from awesIO
   through the resolve step + core `Wing(yaml)` and solve.

## Open questions

- Should core store `up_vector` per section and use it in the panel frame now, or just
  carry it for round-trip fidelity until the frame refactor?
- `chord`/`is_strut` in info_dict — carry through, ignore in core, or use `is_strut`
  to skip lift on struts?
- Where does the resolve step live — AirfoilAero (airfoil-level) or ObjAdapter
  (whole-geometry)? Likely AirfoilAero exposes per-airfoil resolve; ObjAdapter/CLI
  drives it over a geometry file.

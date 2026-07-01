# Draft issue for awegroup/awesIO — surface-pressure (Cp) distribution support

> Not submitted. Draft for review before opening on `github.com/awegroup/awesIO`.

**Title:** Add a surface-pressure (Cp) distribution format for aero-structural coupling

## Summary

awesIO's `aero_geometry.yaml` covers integrated section coefficients (cl/cd/cm) via
`polars` CSVs and the regression/`neuralfoil` types. For **fluid–structure (aero-
elastic) coupling** we also need the **chordwise surface-pressure distribution** per
section, which the current schemas don't cover. This proposes a small, optional
`Cp(alpha, delta)` format that fits the existing `wing_airfoils` structure.

## Motivation

Membrane/soft kites deform under load; coupling an aerodynamic solver to a structural
solver needs a spanwise-and-chordwise pressure field, not just section totals.
Integrated cl/cd/cm cannot reconstruct the chordwise load `ΔCp(x/c)` that drives local
deflection. A standard, tool-agnostic Cp table would let structural codes consume
aero output the same way they consume `polars`.

## Proposed format

A per-airfoil surface-pressure table, referenced from `wing_airfoils` info_dict, e.g.
`cp_file_path: polars/1_cp.csv`, evaluated at **N equal-x/c vertical chord slices**
(so upper/lower pressures share an x station and `ΔCp = Cp_lower - Cp_upper` is direct):

```
alpha, delta, up@0.025, up@0.075, …, lo@0.025, lo@0.075, …
```

- `alpha` (deg), `delta` (trailing-edge deflection, deg) — one row per operating point.
- `up@<x/c>` / `lo@<x/c>` — upper/lower `Cp` at chord fraction `x/c`; `N` = number of
  slices, read from the header (no fixed spacing assumed).
- Self-describing: `N`, the `x/c` stations, and the `(alpha, delta)` grid all come from
  the file, so no companion metadata is needed.

Rationale for chord slices (vs raw panel nodes): a fixed x/c grid is what a structural
mesh wants for interpolation, and it is directly differenceable into `ΔCp`.

## Alternatives considered

- Panel-node `(x, Cp)` pairs (XFoil-style): variable count per angle, not aligned
  across α/δ, awkward for structural interpolation.
- Storing full 2D pressure fields: overkill for section-based lifting-line/VSM coupling.

## Scope

- Optional; purely additive to `wing_airfoils` (a new info_dict key + CSV).
- A JSON schema entry mirroring the `polars` one, plus an example under a kite config.
- Producible by NeuralFoil (edge-velocity → Cp) and XFoil (cpdump); already prototyped
  in VortexStepMethod.jl (`read_cp_data`/`write_cp_data`, `CpData`).

## Reference implementation

VortexStepMethod.jl generates and consumes this format today
(`generate_cp_polar` → `write_cp_data` → `read_cp_data` → per-panel `Cp(alpha, delta)`
in the solver output). Happy to contribute the schema + example.

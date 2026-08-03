# From CAD mesh to aerodynamic model

A kite or wing usually starts life as a 3D `.obj` mesh from CAD, not as a set of
clean 2D airfoils with polars. [`obj_to_yaml`](@ref VortexStepMethod.ObjAdapter.obj_to_yaml)
bridges that gap: it turns the mesh into the package's native YAML geometry route
— per-section airfoil `.dat` files, polar CSVs, and a `geometry.yaml` that ties them
together — so the rest of the solver never has to know the input came from CAD.

The conversion runs four stages per spanwise station:

```
 .obj mesh
    │  1. slice          perpendicular_sections
    ▼
 raw 2D point cloud  (ribs, spars, membrane — noisy)
    │  2. shrink-wrap    shrink_wrap / ShrinkWrap
    ▼
 clean closed airfoil contour  (cosine-clustered panels)
    │  3. 2D aero solver
    ├── NeuralFoil:  fit Kulfan CST  →  neural network      (default, fast)
    └── XFoil:       coordinates  →  panel code  (± repanel)
    │  4. write
    ▼
 airfoils/*.dat  +  polars/*.csv  +  geometry.yaml   →   Wing(geometry.yaml)
```

## 1. Slice

[`perpendicular_sections`](@ref VortexStepMethod.ObjAdapter.perpendicular_sections)
places `n_sections` stations at equal leading-edge arc-length intervals and cuts the
mesh *perpendicular to the local span*, rather than along a fixed global plane. On a
curved kite tip a fixed-plane cut would smear the profile out and exaggerate the
chord; a perpendicular cut keeps each airfoil undistorted. Each slice comes back as a
raw 2D point cloud together with the station's leading- and trailing-edge positions in
3D. The cloud is messy: for a ram-air kite it contains the outer membrane plus every
interior rib, spar and seam the cutting plane happened to cross.

## 2. Shrink-wrap

[`shrink_wrap`](@ref), configured by [`ShrinkWrap`](@ref), turns that noisy cloud into
a single clean closed airfoil. It builds a distance field on a grid, thresholds it at a
rolling-ball radius (bridging gaps between points and ignoring interior structure),
flood-fills the outside, erodes the boundary back to a small `clearance`, and traces the
resulting level set with marching squares. The contour is parameterised by arclength, so
the leading edge comes out genuinely round and the blunt trailing edge is capped by an
arc — and the output points are cosine-clustered toward both edges. The same wrapped
airfoil is what *both* backends analyse, so the geometry the polar is generated for
always matches the `.dat` written to disk.

This is the robustness win: an arbitrary, self-intersecting, structurally-detailed slice
becomes one well-formed airfoil with no manual cleanup.

## 3. The 2D aerodynamic solver

`obj_to_yaml`'s `aero_solver` argument selects the backend; both consume the *same*
wrapped contour but along different paths.

**NeuralFoil** ([`NeuralFoilSolver`](@ref), the default). The wrapped airfoil is fitted
to Kulfan CST parameters with [`fit_kulfan_parameters`](@ref) (a
[`LeastSquaresFit`](@ref)), and those parameters are fed through a small pre-trained
neural network. NeuralFoil was trained on XFoil results over Kulfan-parametrised
airfoils, so the CST fit is not just a convenience — it puts the airfoil into exactly the
representation the network expects. It is fast, differentiable, and returns an
`analysis_confidence` alongside the coefficients. This is the path used everywhere in the
test suite and the recommended default for generating full polars.

**XFoil** ([`XFoilSolver`](@ref)). The wrapped coordinates go straight into the XFoil
panel code — no Kulfan fit. Because the shrink-wrap already emits smooth cosine panels,
XFoil needs no internal repaneling, and `repanel=false` is the default: on a clean
airfoil XFoil's own curvature-attracted repaneling gives an identical result, and on a
trailing-edge-deflected shape it would re-cluster the hinge crease into panels NeuralFoil
never sees and drift away from it. Set `repanel=true` only if you deliberately want
XFoil's paneling. XFoil is the viscous reference; it is slower and can fail to converge
at some angles (those points come back as `NaN`).

The two backends agree closely with matched settings (transition, `n_crit`, Reynolds),
which is what the `test_backend_comparison.jl` suite checks — the residual difference on a
clean airfoil is NeuralFoil's own model accuracy, not a geometry artefact.

### Trailing-edge deflection

Passing a `delta_range` sweeps trailing-edge deflection as a second polar axis. For each
deflection angle [`deform_section`](@ref) pivots the trailing edge about a crease,
*re-shrink-wraps* the deflected shape, and re-fits Kulfan parameters. Re-normalising means
the chord is always measured leading-edge to trailing-edge, so a deflected section becomes
a cambered, asymmetric airfoil rather than simply a rotated one. The result is written as
a `POLAR_MATRICES` table (coefficients over `alpha × delta`); with no `delta_range` a plain
`POLAR_VECTORS` table (coefficients over `alpha`) is written instead.

## 4. Write to YAML

For each unique airfoil id `j`, `obj_to_yaml` writes into `output_dir`:

- `airfoils/{j}.dat` — the shrink-wrapped airfoil (matches the polar)
- `airfoils/{j}_raw.dat` — the raw sliced points the wrap enclosed (for inspection)
- `airfoils/{j}_d{tag}.dat` — each deflected shape, when a `delta_range` is given. The
  `{tag}` encodes the deflection in degrees (`m` for a minus sign, `p` for the decimal
  point), e.g. `_d5.dat` for 5°, `_dm3.dat` for −3°, `_d2p5.dat` for 2.5°
- `polars/{j}.csv` — the generated polar (`POLAR_VECTORS` or `POLAR_MATRICES`)
- `geometry.yaml` — `wing_sections` (leading/trailing-edge points) plus `wing_airfoils`
  (each section's `type` and the `.dat`/`.csv` paths above)

A near-vanishing wingtip slice can shrink-wrap to an implausibly thick blob; such a
degenerate section reuses its nearest valid neighbour's airfoil and polar while keeping its
own edge positions, and a warning lists the reuse. All floats are rounded to millimetre
precision by the single [`write_yaml`](@ref VortexStepMethod.ObjAdapter.write_yaml) writer,
so generated geometry files stay diff-friendly and consistent.

## Why this is useful

- **CAD in, solver-ready model out.** No hand-drawing airfoils or manually pairing them
  with polars — the mesh alone is enough.
- **Robust to messy geometry.** Ribs, spars and self-intersections in a slice are absorbed
  by the shrink-wrap instead of breaking the airfoil.
- **One geometry, either backend.** Fast NeuralFoil for full-envelope polars, or XFoil as a
  viscous cross-check, both from the identical wrapped contour.
- **Generate once, reuse cheaply.** The `.dat`/`.csv`/`geometry.yaml` bundle is a plain-text
  artifact you check in and load with `Wing(geometry.yaml)`; the expensive slicing,
  wrapping and polar generation happen only when the geometry changes.

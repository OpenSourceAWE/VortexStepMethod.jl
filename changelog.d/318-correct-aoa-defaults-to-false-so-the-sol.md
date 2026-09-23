### Changed

- `correct_aoa` defaults to `true`: the polars are still read at the 3/4-chord control point, and each section force is turned by the flow at its 1/4-chord aerodynamic centre (Gaunaa, Li & Pirrung, TORQUE 2026). Every VSM result moves, the induced drag most. That flow now leaves out only the section's own bound vortex, not every panel's, which also moves `LLT` results on swept or curved wings.

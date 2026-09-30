### Fixed

- The `LOOP` solver caps each panel's relaxation at `1 / (1 + π c |z_airf ⋅ AIC_ii|)`, the diagonal step of its self-induced stiffness, so narrow panels (cosine-spaced tips) no longer make it diverge to NaN; the curved reference wing and finely refined wings now converge at the default `relaxation_factor`, with unchanged iteration counts on wings that already converged. A `LOOP` solve that still does not converge keeps halving `relaxation_factor`, warning each time, while it is above 1e-3, and each attempt stops at the first non-finite residual.

### Changed

- `cl_unrefined_dist`, `cd_unrefined_dist`, `cm_unrefined_dist` and `alpha_unrefined_dist` are area-weighted means over the section's panels instead of plain means, and `chord_unrefined_dist` is the section area over its width, so a coefficient times `chord_unrefined_dist * width_unrefined_dist` is the section's load. `moment_unrefined_dist` is the section's moment, summed over its panels, instead of their mean; `linearize` with `aero_coeffs=false` returns it.

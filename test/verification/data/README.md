# Verification data

Copied from `tests/verification_cases` of
[ocayon/Vortex-Step-Method](https://github.com/ocayon/Vortex-Step-Method) at
commit `c82d521` (MIT licence). Angles in degrees.

| File | Contents | Python source |
| --- | --- | --- |
| `clarky_polar.csv` | Clark Y 2D section polar, CFD to 18°, polynomial fit above | `curved_wing/polars/clarky_maneia.csv`, first of its three identical 20–32° blocks |
| `curved_wing_rans.csv` | curved Clark Y wing, 3D RANS (Maneia) | `curved_wing/polars/curved_wing_polars_maneia.csv` |
| `naca4415_cfd_polar.csv` | NACA 4415 2D section lift, CFD at Re 3e6 | `swept_wing/polars/NACA4415_CFD_Re3e6.csv` |
| `rectangular_wing_ar12_cfd.csv` | unswept AR 12 NACA 4415 wing, 3D CFD | `swept_wing/polars/0sweepAR12_CFD.csv` |
| `v3_kite_rans_cl.csv`, `v3_kite_rans_cd.csv` | TU Delft V3 kite with struts, RANS (Lebesque) | `TUDELFT_V3_KITE/CFD_data/RANS_C{L,D}_alpha_struts.csv` |

## Initial Release
This project is based on version 1.0 of the Python project [Vortex-Step-Method](https://github.com/ocayon/Vortex-Step-Method)

## Noteworthy Differences
- implemented in Julia, therefore about 50 times faster
- an importer for `.obj` wing geometry files was added (see: [#10](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/issues/10))
- automatic creation of polars using Xfoil.jl was added (see: [#43](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/pull/43))
- a ram-air kite example was added
- `Umag` was replaced with `v_a` as variable name for the norm of the apparent wind speed

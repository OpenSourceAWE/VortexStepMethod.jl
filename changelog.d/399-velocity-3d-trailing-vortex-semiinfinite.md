### Changed

- Inside the core of a semi-infinite trailing vortex the induced velocity scales linearly with the distance to the axis, as it does for bound and finite trailing filaments, rather than holding the core-boundary value.

### Fixed

- A semi-infinite trailing vortex evaluated at its start point returns zero rather than NaN, and ForwardDiff sees the velocity's slope on its axis rather than zero.

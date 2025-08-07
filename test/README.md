# Test Data Organization

This directory contains organized test data for the VortexStepMethod.jl test suite.

## Running Tests

You can run tests in different ways:
```bash
# Run all tests
julia --project=. test/runtests.jl
```
```bash
# Run specific test files
julia --project=. test/yaml_geometry/test_yaml_geometry.jl
julia --project=. test/body_aerodynamics/test_body_aerodynamics.jl
julia --project=. test/solver/test_solver.jl
# etc.
```
## Adding New Test Data

When adding test data for a new source module:

1. Create test files in the appropriate module directory under `test/`
2. Add any required data files (YAML, CSV) in the same directory
3. Use `test_data_path()` helper function to access files with proper paths
4. Follow the existing naming conventions


## Structure

The test data is organized by source module being tested, following the structure of the `src/` directory:

```
test/
├── body_aerodynamics/          # Tests for src/body_aerodynamics.jl
│   ├── test_body_aerodynamics.jl
│   ├── test_results.jl
│   └── test_wing.yaml         # Wing geometry for body aero tests
├── yaml_geometry/              # Tests for src/yaml_geometry.jl
│   ├── test_yaml_geometry.jl
│   ├── simple_wing.yaml       # Basic 2-section wing
│   ├── complex_wing.yaml      # Multi-section wing
│   ├── standard_airfoil.csv   # Standard polar data
│   └── alternate_airfoil.csv  # Alternative airfoil data
├── solver/                     # Tests for src/solver.jl
│   ├── test_solver.jl
│   ├── solver_settings.yaml   # Settings for solver tests
│   ├── solver_test_polar.csv  # Polar data for solver tests
│   └── solver_test_wing.yaml  # Wing for solver tests
├── settings/                   # Tests for src/settings.jl
│   └── test_settings.jl
├── filament/                   # Tests for src/filament.jl
│   ├── test_bound_filament.jl
│   └── test_semi_infinite_filament.jl
├── panel/                      # Tests for src/panel.jl
│   └── test_panel.jl
├── plotting/                   # Tests for src/plotting.jl
│   └── test_plotting.jl
├── polars/                     # Tests for src/polars.jl
│   └── test_polars.jl
├── ram_geometry/               # Tests for src/ram_geometry.jl
│   └── test_kite_geometry.jl
├── wake/                       # Tests for src/wake.jl
│   └── test_wake.jl
├── wing_geometry/              # Tests for src/wing_geometry.jl
│   └── test_wing_geometry.jl
├── VortexStepMethod/           # Tests for main module
│   └── test_VortexStepMethod.jl
├── test_data_utils.jl          # Shared utilities for test data access
├── runtests.jl                 # Main test runner
└── README.md                   # This file
```

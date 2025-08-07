# Test Data Organization

This directory contains organized test data for the VortexStepMethod.jl test suite.

## Structure

The test data is organized by source module being tested, following the structure of the `src/` directory:

```
test/data/
├── body_aerodynamics/          # Tests for src/body_aerodynamics.jl
│   ├── wings/
│   │   └── test_wing.yaml
│   ├── polars/
│   │   └── simple_polar.csv
│   └── complete_settings.yaml
├── yaml_geometry/              # Tests for src/yaml_geometry.jl
│   ├── wings/
│   │   ├── simple_wing.yaml
│   │   └── complex_wing.yaml
│   └── polars/
│       ├── standard_airfoil.csv
│       └── alternate_airfoil.csv
├── solver/                     # Tests for src/solver.jl
│   ├── wings/
│   │   └── solver_test_wing.yaml
│   ├── polars/
│   │   └── solver_test_polar.csv
│   └── solver_settings.yaml
├── settings/                   # Common settings files
│   ├── basic_vsm.yaml
│   └── basic_llt.yaml
└── wing_geometry/              # Tests for src/wing_geometry.jl (future)
```

## Usage

Use the helper functions in `test_data_utils.jl` to access test data:

### Basic Path Access
```julia
include("test_data_utils.jl")

# Get paths to test data files
wing_file = test_data_path("yaml_geometry", "wings", "simple_wing.yaml")
polar_file = test_data_path("yaml_geometry", "polars", "standard_airfoil.csv")
settings_file = test_data_path("settings", "basic_vsm.yaml")
```

### Module-Specific Convenience Functions
```julia
# Get standard wing file for a module
wing_file = get_standard_wing_file("body_aerodynamics")

# Get complete settings file for a module
settings_file = get_complete_settings_file("solver")

# Create temporary settings with custom parameters
temp_settings = create_temp_wing_settings("yaml_geometry", "simple_wing.yaml"; 
                                         alpha=15.0, wind_speed=25.0)
```

## Design Principles

1. **Module Isolation**: Each source module has its own test data directory
2. **Consistent Structure**: All modules follow the same subdirectory pattern (wings/, polars/, etc.)
3. **Reusable Components**: Common files (like basic settings) are shared in the `settings/` directory
4. **Clear Naming**: File names clearly indicate their purpose and content
5. **Self-Contained**: Each module's data is independent and can be used in isolation

## Adding New Test Data

When adding test data for a new source module:

1. Create a new directory under `test/data/` matching the source file name
2. Add subdirectories as needed (`wings/`, `polars/`, etc.)
3. Update `test_data_utils.jl` to include helper functions for the new module
4. Follow the existing naming conventions

## File Descriptions

### Wing Geometry Files (.yaml)
- `simple_wing.yaml`: Basic 2-section wing for common tests
- `complex_wing.yaml`: Multi-section wing with multiple airfoil types
- `test_wing.yaml`: Simple wing for body aerodynamics tests
- `solver_test_wing.yaml`: Wing optimized for solver testing

### Polar Data Files (.csv)
- `standard_airfoil.csv`: Standard 6-point polar data
- `alternate_airfoil.csv`: Alternative airfoil with different coefficients
- `simple_polar.csv`: Minimal 3-point polar for basic tests
- `solver_test_polar.csv`: Polar data for solver tests

### Settings Files (.yaml)
- `basic_vsm.yaml`: Basic VSM solver settings
- `basic_llt.yaml`: Basic LLT solver settings
- `complete_settings.yaml`: Full wing + solver configuration
- `solver_settings.yaml`: Solver-specific configuration

This organization makes it easy to:
- Find test data relevant to specific source modules
- Avoid conflicts between different test requirements
- Maintain and update test data independently
- Add new test cases without affecting existing ones

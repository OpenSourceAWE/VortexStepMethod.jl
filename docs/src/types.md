```@meta
CurrentModule = VortexStepMethod
```
## Enumerations
```@docs
Model
WingType
AeroModel
PanelDistribution
PanelGroupingMethod
InitialGammaDistribution
SolverType
SolverStatus
```

## Basic Vectors
```@docs
MVec3
PosVector
VelVector
```

## Aerodynamic data
```@docs
AeroData
```

## Settings
```@docs
VSMSettings
WingSettings
SolverSettings
```

## Wing Geometry, Panel and Aerodynamics
A body is constructed of one or more abstract wings. All wings are of type Wing.
A Wing has one or more sections and can be created from YAML files or OBJ geometry.
```@docs
Section
Section(LE_point, TE_point, aero_model)
Wing
Wing(n_panels::Int; spanwise_distribution::PanelDistribution=LINEAR,
     spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]))
ObjWing
ObjWing(obj_path, dat_path; alpha=0.0, crease_frac=0.75, wind_vel=10., mass=1.0,
         n_panels=54, n_sections=n_panels+1, spanwise_distribution=UNCHANGED,
         spanwise_direction=[0.0, 1.0, 0.0])
BodyAerodynamics
```

## The Solver and its results
```@docs
Solver
VSMSolution
```

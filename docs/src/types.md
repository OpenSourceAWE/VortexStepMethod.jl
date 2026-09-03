```@meta
CurrentModule = VortexStepMethod
```
## Enumerations
```@docs
Model
WingType
AeroModel
PanelDistribution
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
MeshSettings
AirfoilSettings
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
BodyAerodynamics
```

## The Solver and its results
```@docs
Solver
VSMSolution
```

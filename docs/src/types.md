```@meta
CurrentModule = VortexStepMethod
```
## Enumerations
```@docs
Model
AeroModel
PanelDistribution
InitialGammaDistribution
```

## Basic Vectors
```@docs
   MVec3
   PosVector
   VelVector
```

## Wing Geometry, Panel and Aerodynamics
A body is constructed of one or more abstract wings. An abstract wing can be a Wing or a KiteWing. 
A Wing/ KiteWing has one or more sections.
```@docs
    Section
    Wing
    KiteWing
    BodyAerodynamics
```

## The Solver
```@docs
Solver
```

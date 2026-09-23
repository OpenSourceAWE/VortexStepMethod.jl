## Reference Frames

### Introduction
Reference frames are needed for following purposes:
- for creating a CAD model of the wing (or the wings)
- for defining the apparent wind speed vector `va_vec`
- for calculating the lift and drag and side force coefficients
- for calculating the resulting forces and moments

### CAD reference frame (CAD)
A geometric model is always created using the CAD reference frame.
It can have any origin (with respect to the kite), but usually either the center of gravity of the body or the bridle point/ Kite Control Unit is used. 

- Y is defined spanwise, towards the right tip (seen from the front, the left side).
- X is defined chord wise, from LE to TE, positive.
- Z is defined as positive upwards.

### Kite aero body frame (KA)
The body-fixed frame of the kite is the **KA** frame of
[KiteUtils.jl](https://github.com/OpenSourceAWE/KiteUtils.jl), aft-right-up:
- X is defined chord wise, from LE to TE, positive.
- Y is defined spanwise, towards the right tip (seen from the front, the left side).
- Z is defined as the cross product of X and Y, so positive upwards.

The body-axis force components are Fx, Fy and Fz, with the matching coefficients `cfx`, `cfy` and `cfz`.

The origin of the KA frame can be defined by the user by passing the keyword argument `kite_body_origin = ...` to the `BodyAerodynamics` constructor.

## The turn rates
The turn rates $\mathrm{omega} = [\mathrm{omega_x}, \mathrm{omega_y} ,\mathrm{omega_z}]$ are defined in the **KA** frame. The unit of the components is $\mathrm{rad}~\mathrm{s^{-1}}$.

## Input and output
- when running a simulation, the turnrate of the kite must be provided on each time step
- the apparent wind speed vector `va_vec` is defined in the **KA** frame
- the resulting forces are defined in the **KA** frame
- the moments and moment coefficients are defined in the **KA** frame
- **CD** is along the apparent wind, **CL** along `va × y` and **CS** completes the wind axes, so they match Fx, Fz and Fy only at α = β = 0
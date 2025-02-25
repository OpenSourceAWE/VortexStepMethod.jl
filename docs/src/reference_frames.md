## Reference Frames

### Introduction
Reference frames are needed for following purposes:
- for creating a CAD model of the wing (or the wings)
- for defining the apparent wind speed vector $v_a$
- for calculating the lift and drag and side force coefficients
- for calculating the resulting forces and moments

### CAD reference frame (CAD)
A geometric model is always created using the CAD reference frame.
It can have any origin (with respect to the kite), but usually either the center of gravity of the body or the bridle point/ Kite Control Unit is used. 

- Y defined spanwise, looking at the kite from the front (so seeing the LE first) the front left is positive.
- X is defined chord wise, from LE to TE, positive.
- Z is defined as positive upwards.

### Kite reference frame (KP)
This is a body-fixed reference frame.
- Y defined spanwise, looking at the kite from the front (so seeing the LE first) the front left is positive.
- X is defined chord wise, from LE to TE, positive.
- Z is defined as the cross product of Y and X

The origin of the kite reference frame can be defined by the user by calling the function `init_kp(origin::MVec3)` where the origin must be defined in the `CAD` reference frame.

This reference frame is different from the kite reference frame **K** used in `KiteModels.jl` and `KiteUtils.jl`.

## The turn rates
The turn rates $\mathrm{omega} = [\mathrm{omega_x}, \mathrm{omega_y} ,\mathrm{omega_z}]$ are defined in the **KP** reference frame. The unit of the components is $\mathrm{rad}~\mathrm{s^{-1}}$.

## Input and output
- when running a simulation, the turnrate of the kite must be provided on each time step
- the apparent wind speed vector `v_a` is defined in the **KP** reference frame
- the resulting forces are defined in the **KP** reference frame
- the **CL**, **CD**, **CS** and the resulting moments and moment coefficients are defined in the **KP** reference frame
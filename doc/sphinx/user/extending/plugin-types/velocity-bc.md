(sec:extending:plugin-types:velocity-bc)=
# Velocity and traction boundary conditions

## Prescribed velocity boundary conditions

Most of the time, one chooses relatively simple boundary values for the
velocity: either a zero boundary velocity, a tangential flow model in which
the tangential velocity is unspecified but the normal velocity is zero at the
boundary, or one in which all components of the velocity are unspecified
(i.e., for example, an outflow or inflow condition where the total stress in
the fluid is assumed to be zero). However, sometimes we want to choose a
velocity model in which the velocity on the boundary equals some prescribed
value. A typical example is one in which plate velocities are known, for
example their current values or historical reconstructions. In that case, one
needs a model in which one needs to be able to evaluate the velocity at
individual points at the boundary. This can be implemented via plugins.

To implement a new boundary velocity model, you need to overload the
[aspect::VelocityBoundaryConditions::Interface](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1BoundaryVelocity_1_1Interface.html)
class and use the
`ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL` macro to register your new
class. The implementation of the new class should be in namespace
`aspect::BoundaryVelocity`.

In essence, the main function you have to implement for this plugin system is one
that, given a point returns the prescribed velocity at that point.
There are also member functions in the base class you can overload that are called
at the beginning of each time step; this is useful if one needs to perform some
operation once for each time step; examples are models that need to call an
external program to obtain plate velocities for the current time, or from
historical records, in which case it is far cheaper to do so only once at the
beginning of the time step than for every boundary point separately. See, for
example, the `aspect::VelocityBoundaryConditions::GPlates` class.

The remaining functions are discussed in the
documentation of this interface class at
[aspect::VelocityBoundaryConditions::Interface](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1BoundaryVelocity_1_1Interface.html).


## Prescribed traction boundary conditions

Alternatively, at a boundary one can prescribe the traction (i.e., a force density)
that *drives* the velocity, rather than the velocity itself. As for prescribed
velocities, a plugin system allows to describe this information, with the key
function being one that for a given point returns the prescribed traction.

To implement a new boundary traction model, you need to overload the
[aspect::BoundaryTraction::Interface](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1BoundaryTraction_1_1Interface.html)
class and use the
`ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL` macro to register your new
class. The implementation of the new class should be in namespace
`aspect::BoundaryTraction`.
The member functions that can be overloaded are discussed in the
documentation of this interface class at
[aspect::BoundaryTraction::Interface](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1BoundaryTraction_1_1Interface.html).

# Prescribed velocity boundary conditions

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
`ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS` macro to register your new
class. The implementation of the new class should be in namespace
`aspect::VelocityBoundaryConditions`.

Specifically, your new class needs to implement the following basic interface:

```{code-block} c++
template <int dim>
    class aspect::VelocityBoundaryConditions::Interface
    {
      public:
        virtual
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const = 0;

        virtual
        void
        initialize (const GeometryModel::Interface<dim> &geometry_model);

        virtual
        void
        update ();

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
```

The first of these functions needs to provide the velocity at the given point.
The next two are other member functions that can (but need not) be overloaded
if a model wants to do initialization steps at the beginning of the program or
at the beginning of each time step. Examples are models that need to call an
external program to obtain plate velocities for the current time, or from
historical records, in which case it is far cheaper to do so only once at the
beginning of the time step than for every boundary point separately. See, for
example, the `aspect::VelocityBoundaryConditions::GPlates` class.

The remaining functions are obvious, and are also discussed in the
documentation of this interface class at
[aspect::VelocityBoundaryConditions::Interface](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1BoundaryVelocity_1_1Interface.html).
The purpose of the last two
functions has been discussed in the general overview of plugins above.

# Initial conditions

The initial conditions model is responsible for describing the initial
temperature distribution throughout the domain. It essentially has to provide
a function that for each point can return the initial temperature. Note that
the model {math:numref}`eq:stokes-1`2]&ndash;{math:numref}`eq:temperature`3] does not require
initial values for the pressure or velocity. However, if coefficients are
nonlinear, one can significantly reduce the number of initial nonlinear
iterations if a good guess for them is available; consequently, <span
ASPECT initializes the pressure with the
adiabatically computed hydrostatic pressure, and a zero velocity. Neither of
these two has to be provided by the objects considered in this section.

To implement a new initial conditions model, you need to overload the
[aspect::InitialConditions::Interface][] class and use the
`ASPECT_REGISTER_INITIAL_CONDITIONS` macro to register your new class. The
implementation of the new class should be in namespace
`aspect::InitialConditions`.

Specifically, your new class needs to implement the following basic interface:

``` c++
template <int dim>
    class aspect::InitialConditions::Interface
    {
      public:
        void
        initialize (const GeometryModel::Interface<dim>       &geometry_model,
                    const BoundaryTemperature::Interface<dim> &boundary_temperature,
                    const AdiabaticConditions<dim>            &adiabatic_conditions);

        virtual
        double
        initial_temperature (const Point<dim> &position) const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
```

The meaning of the first class should be clear. The purpose of the last two
functions has been discussed in the general overview of plugins above.

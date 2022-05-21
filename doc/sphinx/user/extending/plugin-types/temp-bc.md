# Temperature boundary conditions

The boundary conditions are responsible for describing the temperature values
at those parts of the boundary at which the temperature is fixed (see
{ref}`sec:1.4.3][] for how it is determined which parts of the boundary
this applies to).

To implement a new boundary conditions model, you need to overload the
[aspect::BoundaryTemperature::Interface][] class and use the
`ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::BoundaryTemperature`.

Specifically, your new class needs to implement the following basic interface:

```{code-block} c++
template <int dim>
    class aspect::BoundaryTemperature::Interface
    {
      public:
        virtual
        double
        temperature (const GeometryModel::Interface<dim> &geometry_model,
                     const unsigned int                   boundary_indicator,
                     const Point<dim>                    &location) const = 0;

        virtual
        double minimal_temperature () const = 0;

        virtual
        double maximal_temperature () const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
```

The first of these functions needs to provide the fixed temperature at the
given point. The geometry model and the boundary indicator of the particular
piece of boundary on which the point is located is also given as a hint in
determining where this point may be located; this may, for example, be used to
determine if a point is on the inner or outer boundary of a spherical shell.
The remaining functions are obvious, and are also discussed in the
documentation of this interface class at
[aspect::BoundaryTemperature::Interface][]. The purpose of the last two
functions has been discussed in the general overview of plugins above.

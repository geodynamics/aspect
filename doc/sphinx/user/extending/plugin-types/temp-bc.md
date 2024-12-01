
# Temperature and compositional boundary conditions

Boundary conditions for the temperature (and compositional fields, if
present) are responsible for describing the fields' values at those
parts of the boundary at which the temperature or composition is
fixed, as well as at inflow boundaries.

To implement a new boundary conditions model for the temperature, you need to overload the
`aspect::BoundaryTemperature::Interface` class and use the
`ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::BoundaryTemperature`.

Correspondingly, for a boundary conditions model for compositional fields, you need to overload the
`aspect::BoundaryComposition::Interface` class and use the
`ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::BoundaryComposition`.

In either case,
the principle function of the respective `Interface` class that you need to overload
provides the prescribed temperature or composition at a
given point. The boundary indicator of the particular
piece of boundary on which the point is located is also given as a hint in
determining where this point may be located; this may, for example, be used to
determine if a point is on the inner or outer boundary of a spherical shell.
Models may also need to query the geometry model (via the `SimulatorAccess` class)
to understand what a specific boundary indicator is supposed to mean &ndash; e.g.,
if the given boundary indicator corresponds to an inner or outer surface of a shell
geometry.

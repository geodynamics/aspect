
# Temperature boundary conditions

The boundary conditions are responsible for describing the temperature values
at those parts of the boundary at which the temperature is fixed.

To implement a new boundary conditions model, you need to overload the
`aspect::BoundaryTemperature::Interface` class and use the
`ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::BoundaryTemperature`.

The principle function you need to overload needs to provide the fixed temperature at a
given point. The boundary indicator of the particular
piece of boundary on which the point is located is also given as a hint in
determining where this point may be located; this may, for example, be used to
determine if a point is on the inner or outer boundary of a spherical shell.
Models may also need to query the geometry model (via the `SimulatorAccess` class)
to understand what a specific boundary indicator is supposed to mean &ndash; e.g.,
if the given boundary indicator corresponds to an inner or outer surface of a shell
geometry.

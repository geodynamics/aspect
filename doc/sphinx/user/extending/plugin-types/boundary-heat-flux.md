
# Boundary heat flux

In contrast to prescribing the actual temperature at a boundary, it is
also possible to prescribe the *heat flux* across a boundary. A number
of heat flux models are already implemented, but if you want to
implement a new boundary heat flux model, you need to overload the
`aspect::BoundaryHeatFlux::Interface` class and use the
`ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::BoundaryTemperature`.

The `Interface` class that you need to overload provides the
prescribed heat flux at a given set of points via the `heat_flux()`
function. The function receives information about material parameters
at these points as well. For historical reasons, the function is asked
to provide the heat flux as a vector, even though the place where the
heat flux is used only uses the component of this vector that is
*normal* to the boundary (which it computes by taking the dot product
between the returned vector and the normal vector). Because there are
situations where all you can do is compute the normal heat flux as a
scalar, the `heat_flux()` function also receives the normal vector as
an input argument. As a consequence, one way for the function to
compute the required heat flux vector is to compute the scalar heat
flux and multiply it by the normal vector.

Finally, the function also receives the boundary indicator of the particular
piece of boundary on which the point is located, as a hint in
determining where this point may be located; this may, for example, be used to
determine if a point is on the inner or outer boundary of a spherical shell.
Models may also need to query the geometry model (via the `SimulatorAccess` class)
to understand what a specific boundary indicator is supposed to mean &ndash; e.g.,
if the given boundary indicator corresponds to an inner or outer surface of a shell
geometry.

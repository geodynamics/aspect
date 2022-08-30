# Switching off pressure normalization

In most practically relevant cases, the Stokes equations
{math:numref}`eq:stokes-1`-{math:numref}`eq:stokes-2` only determine the pressure up
to a constant because only the pressure gradient appears in the equations, not
the actual value of it. However, unlike this "mathematical"
pressure, we have a very specific notion of the "physical"
pressure: namely a well-defined quantity that at the surface of Earth equals
the air pressure, which compared to the hydrostatic pressure inside Earth is
essentially zero.

As a consequence, the default in ASPECT is to
normalize the computed "mathematical" pressure in such a way that
either the mean pressure at the surface is zero (where the geometry model
describes where the "surface" is, see
{ref}`sec:extending:plugin-types:geometry-models`), or that the mean pressure in the
domain is zero. This normalization is important if your model describes
densities, viscosities and other quantities in dependence of the pressure
because you almost certainly had the "physical" pressure
in mind, not some unspecified "mathematical" one. On the other
hand, if you have a material model in which the pressure does not enter, then
you don't need to normalize the pressure at all - simply go with
whatever the solver provides. In that case, you can switch off pressure
normalization by looking at the `Pressure normalization` parameter at the top
level of the input file, see {ref}`parameters:global`.

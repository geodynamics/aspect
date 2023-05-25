(sec:methods:pressure-norm)=
# Pressure normalization

The equations described above, {math:numref}`eq:stokes-1`-{math:numref}`eq:temperature`, only determine the pressure $p$ up to an additive constant.
On the other hand, since the pressure appears in the definition of many of the coefficients, we need a pressure that has some sort of *absolute* definition.
A physically useful definition would be to normalize the pressure in such a way that the average pressure along the "surface" has a prescribed value where the geometry description (see {ref}`sec:extending:plugin-types:geometry-models`) has to determine which part of the boundary of the domain is the "surface" (we call a part of the boundary the "surface" if its depth is "close to zero").

Typically, one will choose this average pressure to be zero, but there is a parameter "`Surface pressure`" in the input file (see {ref}`parameters:global`) to set it to a different value.
One may want to do that, for example, if one wants to simulate the Earth's mantle without the overlying lithosphere.
In that case, the "surface" would be the interface between mantle and lithosphere, and the average pressure at the surface to which the solution of the equations will be normalized should in this case be the hydrostatic pressure at the bottom of the lithosphere.

An alternative is to normalize the pressure in such a way that the *average* pressure throughout the domain is zero or some constant value.
This is not a useful approach for most geodynamics applications but is common in benchmarks for which analytic solutions are available.
Which kind of normalization is chosen is determined by the "`Pressure normalization`" flag in the input file, see {ref}`parameters:global`.

(sec:cookbooks:platelike-boundary)=
# Convection in a box with prescribed, variable velocity boundary conditions

A similarly simple setup to the ones considered in the previous subsections is
to equip the model we had with a different set of boundary conditions. There,
we used slip boundary conditions, i.e., the fluid can flow tangentially along
the four sides of our box but this tangential velocity is unspecified. On the
other hand, in many situations, one would like to actually prescribe the
tangential flow velocity as well. A typical application would be to use
boundary conditions at the top that describe experimentally determined
velocities of plates. This cookbook shows a simple version of something like
this. To make it slightly more interesting, we choose a $2\times 1$ domain in
2d.

Like for many other things, has a set of plugins for prescribed velocity
boundary values (see
Sections&nbsp;{ref}`parameters:Boundary_20velocity_20model` and
{ref}`sec:extending:plugin-types:velocity-bc`). These plugins allow one
to write sophisticated models for the boundary velocity on parts or all of the
boundary, but there is also one simple implementation that just takes a
formula for the components of the velocity.

To illustrate this, let us consider the
[cookbooks/platelike-boundary.prm](https://github.com/geodynamics/aspect/blob/main/cookbooks/platelike-boundary/platelike-boundary.prm) input file. It
essentially extends the input file considered in the previous example. The
part of this file that we are particularly interested in in the current
context is the selection of the kind of velocity boundary conditions on the
four sides of the box geometry, which we do using a section like this:

``` {literalinclude} boundary.part.prm
```

We use tangential flow at boundaries named left, right and bottom.
Additionally, we specify a comma separated list (here with only a single
element) of pairs consisting of the name of a boundary and the name of a
prescribed velocity boundary model. Here, we use the `function` model on the
`top` boundary, which allows us to provide a function-like notation for the
components of the velocity vector at the boundary.

The second part we need is that we actually describe the function that sets
the velocity. We do this in the subsection `Function`. The first of these
parameters gives names to the components of the position vector (here, we are
in 2d and we use $x$ and $z$ as spatial variable names) and the time. We could
have left this entry at its default, `x,y,t`, but since we often think in
terms of &ldquo;depth&rdquo; as the vertical direction, let us use `z` for the
second coordinate. In the second parameter we define symbolic constants that
can be used in the formula for the velocity that is specified in the last
parameter. This formula needs to have as many components as there are space
dimensions, separated by semicolons. As stated, this means that we prescribe
the (horizontal) $x$-velocity and set the vertical velocity to zero. The
horizontal component is here either $1$ or $-1$, depending on whether we are
to the right or the left of the point $1+\sin(\pi t/2)$ that is moving back
and forth with time once every four time units. The `if` statement understood
by the parser we use for these formulas has the syntax
`if(condition, value-if-true, value-if-false)`.

The remainder of the setup is described in the following, complete input file:

``` {literalinclude} platelike.prm
```

This model description yields a setup with a Rayleigh number of 200 (taking
into account that the domain has size 2). It would, thus, be dominated by heat
conduction rather than convection if the prescribed velocity boundary
conditions did not provide a stirring action. Visualizing the results of this
simulation[^footnote1] yields images like the ones shown in {numref}`fig:platelike`.

```{figure-md} fig:platelike
<img src="platelike.*" style="width:30.0%" />

Variable velocity boundary conditions: Temperature and velocity fields at the initial time (top left) and at various other points in time during the simulation.
```

[^footnote1]: In fact, the pictures are generated using a twice more refined mesh to
provide adequate resolution. We keep the default setting of five global
refinements in the parameter file as documented above to keep compute time
reasonable when using the default settings.

(sec:cookbooks:overview)=
# How to set up computations

<span class="smallcaps">ASPECT</span>'s computations are controlled by
input parameter files such as those we will discuss in the following
sections.[^footnote1] Basically, these are just regular text files you can edit with
programs like `gedit`, `kwrite` or `kate` when working on Linux, or something
as simple as `NotePad` on Windows. When setting up these input files for a
model you have in mind, you have to describe everything that characterizes the
situation you are considering. In particular, this includes the following:

-   What internal forces act on the medium (the left-hand side of the equation)?

-   What external forces act on the medium (the right-hand of the equation side)?

-   What is the domain (geometry)?

-   What happens at the boundary for each variable involved (boundary
    conditions)?

-   How did it look at the beginning (initial conditions)?

For each of these questions, there are one or more input parameters (sometimes
grouped into sections) that allow you to specify what you want. For example,
to choose a geometry, you will typically have a block like this in your input
file:

```{literalinclude} ../../../manual/cookbooks/overview/doc/geometry.part.prm
```

This indicates that you want to do a computation in 2d, using a rectangular
geometry (a &ldquo;box&rdquo;) with the edge length ("extent") equal to one in both the $x$-
and $y$-directions. Of course, there are other geometries you can choose from
for the `Model name` parameter, and consequently other subsections that
specify the details of these geometries.

Similarly, you describe boundary conditions using parameters such as this:

```{literalinclude} ../../../manual/cookbooks/overview/doc/boundary-conditions.part.prm
```

This snippet describes which of the four boundaries of the two-dimensional box
we have selected above should have a prescribed temperature or an insulating
boundary, and at which parts of the boundary we want zero, tangential or
prescribed velocities.[^footnote2]

If you go down the list of questions about the setup above, you have already
done the majority of the work describing your computation. The remaining
parameters you will typically want to specify have to do with the computation
itself. For example, what variables do you want to output and how often? What
statistics do you want to compute. The following sections will give ample
examples for all of this, but using the questions above as a guideline is
already a good first step.

:::{note}
It is of course possible to set up input files for computations completely from scratch.
However, in practice, it is often simpler to go through the list of cookbooks already provided and
find one that comes close to what you want to do. You would then modify this cookbook until it
does what you want to do. The advantage is that you can start with something you already know
works, and you can inspect how each change you make – changing the details of the geometry,
changing the material model, or changing what is being computed at the end of each time step –
affects what you get.
:::

[^footnote1]: You can also extend ASPECT using plugins - i.e., pieces of code you compile separately and either link into the ASPECT executable itself, or reference from the input file. This is discussed in {ref}`cha:extending`.

[^footnote2]: Internally, the geometry models ASPECT uses label every part of the boundary with what is called a *boundary indicator*
– a number that identifies pieces of the boundary. If you know which number each piece has, you can list these numbers on
the right hand sides of the assignments of boundary types above. For example, the left boundary of the box has boundary
indicator zero (see {ref}`parameters:Geometry_20model`), and using this number instead of the `left` would have been equally valid. However, numbers
are far more difficult to remember than names, and consequently every geometry model provides string aliases such as “`left`”
for each boundary indicator describing parts of the boundary. These symbolic aliases are specific to the geometry – for the box,
they are “`left`”, “`right`”, “`bottom`”, etc., whereas for a spherical shell they are “`inner`” and “`outer`” – but are described in the
documentation of every geometry model, see {ref}`parameters:Geometry_20model`.

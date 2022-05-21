### How to set up computations

<span class="smallcaps">ASPECT</span>'s computations are controlled by
input parameter files such as those we will discuss in the following
sections.[1] Basically, these are just regular text files you can edit with
programs like `gedit`, `kwrite` or `kate` when working on Linux, or something
as simple as `NotePad` on Windows. When setting up these input files for a
model you have in mind, you have to describe everything that characterizes the
situation you are considering. In particular, this includes the following:

-   What internal forces act on the medium (the equation)?

-   What external forces do we have (the right hand side)

-   What is the domain (geometry)?

-   What happens at the boundary for each variable involved (boundary
    conditions)?

-   How did it look at the beginning (initial conditions)?

For each of these questions, there are one or more input parameters (sometimes
grouped into sections) that allow you to specify what you want. For example,
to choose a geometry, you will typically have a block like this in your input
file:

``` prmfile
```

This indicates that you want to do a computation in 2d, using a rectangular
geometry (a &ldquo;box&rdquo;) with edge length equal to one in both the $x$-
and $y$-directions. Of course, there are other geometries you can choose from
for the `Model name` parameter, and consequently other subsections that
specify the details of these geometries.

Similarly, you describe boundary conditions using parameters such as this:

``` prmfile
```

This snippet describes which of the four boundaries of the two-dimensional box
we have selected above should have a prescribed temperature or an insulating
boundary, and at which parts of the boundary we want zero, tangential or
prescribed velocities.[2]

If you go down the list of questions about the setup above, you have already
done the majority of the work describing your computation. The remaining
parameters you will typically want to specify have to do with the computation
itself. For example, what variables do you want to output and how often? What
statistics do you want to compute. The following sections will give ample
examples for all of this, but using the questions above as a guideline is
already a good first step.

<div class="center">

</div>

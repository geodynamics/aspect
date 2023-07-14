(sec:methods:particles)=
# Particles

ASPECT can, optionally, also deal with
particles (sometimes called "tracers"). Particles can be thought
of as point-like objects that are simply advected along with the flow. In
other words, if $\mathbf u(\mathbf x,t)$ is the flow field that results from
solving equations {math:numref}`eq:stokes-1`-{math:numref}`eq:stokes-2`, then the
$k$th particle's position satisfies the equations
```{math}
\begin{aligned}
  \frac{d}{dt} \mathbf x_k(t)
  = \mathbf u(\mathbf x_k(t),t).
\end{aligned}
```
The initial positions of all
particles also need to be given and are usually either chosen randomly, based
on a fixed pattern, or are read from a file.

Particles are typically used to track visually where material that starts
somewhere ends up after some time of a simulation. It can also be used to
track the *history* of the volume of the fluid that surrounds a particle, for
example by tracking how much strain has accumulated, or what the minimal or
maximal temperature may have been in the medium along the trajectory of a
particle. To this end, particles can carry *properties*. These are scalar- or
vector-valued quantities that are attached to each particle, that are
initialized at the beginning of a simulation, and that are then updated at
each time step. In other words, if we denote by $\mathbf p_{k,m}(t)$ the value
of the $m$th property attached to the $k$th particle, then
$\mathbf p_{k,m}(t)$ will satisfy a differential equation of the form
```{math}
\begin{aligned}
\frac{\partial}{\partial t} \mathbf p_{k,m}(t)
= \mathbf g_m\left(\mathbf p_{k,m},
p(\mathbf x_k(t),t)), T(\mathbf x_k(t),t)),
\varepsilon(\mathbf u(\mathbf x_k(t),t)),
\mathfrak c(\mathbf x_k(t),t)\right).
\end{aligned}
```
The exact form of
$\mathbf g_m$ of course depends on what exactly a particular property
represents. Like with compositional fields (see {ref}`sec:methods:compositional-fields`), it is
possible to describe the right hand side $\mathbf g_m$ in ways that also
allows for impulse (delta) functions in time.

How particles are used in practice is probably best explained using examples.
To this end, see in particular {ref}`sec:cookbooks:using-particles`. All particle-related
input parameters are listed in
{ref}`parameters:Postprocess/Particles`. The implementation of
particles is discussed in great detail in {cite:t}`gassmoller:etal:2018`.

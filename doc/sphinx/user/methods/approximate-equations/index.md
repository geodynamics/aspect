(sec:methods:approximate-equations)=
# Approximate equations
There are a number of common variations to equations
{math:numref}`eq:stokes-1`-{math:numref}`eq:temperature` that are used in the geosciences.
For example, one frequently finds references to the anelastic liquid approximation (ALA), truncated anelastic liquid approximation (TALA), and the Boussinesq approximation (BA).
These can all be derived from the basic equations{math:numref}`eq:stokes-1`-{math:numref}`eq:temperature` via various approximations, and we will discuss them in the following sections.
Since they are typically only provided considering velocity, pressure and temperature, we will in the following omit the dependence on the compositional fields used in previous sections, though this dependence can easily be added back into the equations stated below.
A detailed discussion of the approximations introduced below can also be found in {cite:t}`schubert:etal:2001` and {cite:t}`king:etal:2010`; a theoretical and practical comparison of many of these formulations using ASPECT can be found in {cite:t}`gassmoller:etal:2020`.

:::{note}
Historically, the mantle convection community has typically used one or another of these simplified formulations for computer simulations - oftentimes the simplest of them, the Boussinesq Approximation (BA) discussed in {ref}`sec:methods:approximate-equations:ba`.
These kinds of approximations are appropriate in many contexts; for example, for crustal dynamics simulations, the hydrostatic pressures are never high enough to lead to noticeable compression effects and as a consequence the density really is more or less independent of the pressure - as assumed in several of the approximations below.
Yet it is worth pointing out that many older publications showing mantle convection simulations did not rely on these approximations because they describe the physical situation better than equations {math:numref}`eq:stokes-1`-{math:numref}`eq:temperature`, but *because simulation technology did not allow for anything else at the time*.
This has changed today, and ASPECT implements more realistic formulations as discussed in this section and in {ref}`sec:methods:choosing-a-formulation`.
As a consequence, you should evaluate which formulation is appropriate for what you want to do.
The fact that someone else in the past used a simplified formulation does not mean that you should do the same for a similar situation: it could just indicate that they did not have the technology to use a more complete formulation at the time.
:::

The three approximations mentioned all start by writing the pressure and temperature as the sum of a (possibly depth dependent) reference state plus a perturbation, i.e., we will write
```{math}
\begin{aligned}
  p(\mathbf x,t) &= \bar p(z) + p'(\mathbf x,t), \\
  T(\mathbf x,t) &= \bar T(z) + T'(\mathbf x,t).
\end{aligned}
```
Here, barred quantities are reference states and may depend on the depth $z$ (not necessarily the third component of $\mathbf x$) whereas primed quantities are the spatially and temporally variable deviations of the temperature and pressure fields from this reference state.
In particular, the reference pressure is given by solving the hydrostatic equation,
```{math}
:label: eq:hydrostatic-pressure
\begin{aligned}
  \nabla \bar p = \bar\rho \mathbf g,
\end{aligned}
```
where $\bar\rho=\rho(\bar p,\bar T)$ is a *reference density* that depends on depth and represents a typical change of material parameters and solution variables with depth.
$\bar T(z)$ is chosen as an adiabatic profile accounting for the fact that the temperature increases as the pressure increases.
With these definitions, equations {math:numref}`eq:stokes-1`-{math:numref}`eq:stokes-2` can equivalently be written as follows:
```{math}
\begin{aligned}
  -\nabla \cdot \left[2\eta \left(\varepsilon(\mathbf u)
                                  - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)
                \right] + \nabla p' &=  (\rho-\bar\rho) \mathbf g  & \qquad  & \textrm{in $\Omega$},  \\
  \nabla \cdot (\rho \mathbf u) &= 0
  & \qquad
  & \textrm{in $\Omega$}.
\end{aligned}
```
The temperature equation, when omitting entropic effects, still reads as
```{math}
:label: eq:temperature-decomposed
\begin{gathered}
  \rho C_p \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla T\right)
  - \nabla\cdot k\nabla T   \\
  =  \rho H + 2\eta  \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right):
  \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right) +\alpha T \left( \mathbf u \cdot \nabla p \right)  \quad  \textrm{in $\Omega$},
\end{gathered}
```
where the right-hand side includes radiogenic heat production, shear heating and adiabatic heating (in that order).

Starting from these equations, the approximations discussed in the next few subsections make use of the fact that for the flows for which these approximations are valid, the perturbations $p'$, $T'$ are much smaller than typical values of the reference quantities $\bar p$, $\bar T$.
The terms influenced by these approximations are $\nabla \cdot (\rho u) =0$ in the continuity equation, and all occurrences of $\rho(p,T)$ in the temperature equation, and we will discuss them separately below.
The equations for these approximations are almost always given in terms of non-dimensionalized quantities.
We will for now stick with the dimensional form because it expresses in a clearer way the approximations that are made.
The non-dimensionalization can then be done on each of the forms below separately.

:::{toctree}
ala.md
tala.md
ba.md
ica.md
:::

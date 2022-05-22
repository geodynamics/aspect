(sec:methods:approximate-equations:ba)=
# The Boussinesq approximation (BA)

If we further assume that the reference temperature and the reference density are constant, $\bar T(z)=T_0$, $\bar\rho(\bar p,\bar T)=\rho_0$, - in other words, density variations are so small that they are negligible everywhere except for in the right-hand side of the velocity equation (the buoyancy term), which describes the driving force of the flow, then we can further simplify the mass conservation equations of the TALA to $\nabla \cdot \mathbf u=0$.
This means that the density in all other parts of the equations is not only independent of the pressure variations $p'$ as assumed in the TALA, but also does not depend on the much larger hydrostatic pressure $\bar p$ nor on the reference temperature $\bar T$.
We then obtain the following set of equations that also uses the incompressibility in the definition of the strain rate:
```{math}
\begin{aligned}
  -\nabla \cdot \left[2\eta \varepsilon(\mathbf u) \right] + \nabla p' &=
  -\bar \alpha \bar\rho T' \mathbf g  & \qquad  & \textrm{in $\Omega$}, \\
  \nabla \cdot \mathbf u &= 0  & \qquad  & \textrm{in $\Omega$}.
\end{aligned}
```
In addition, as the reference temperature is constant, one needs to neglect the adiabatic and shear heating in the energy equation
```{math}
:label: eq:temperature-BA
  \bar\rho C_p \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla
  T\right) - \nabla\cdot k\nabla T  =   \bar\rho H   \quad \textrm{in $\Omega$}.
```

## On incompressibility

The Boussinesq approximation assumes that the density can be considered constant in all occurrences in the equations with the exception of the buoyancy term on the right hand side of {math:numref}`eq:stokes-1`.
The primary result of this assumption is that the continuity equation {math:numref}`eq:stokes-2` will now read
```{math}
\begin{gathered}
  \nabla \cdot \mathbf u = 0.
\end{gathered}
```
This makes the equations *much* simpler to solve: First, because the divergence operation in this equation is the transpose of the gradient of the pressure in the momentum equation {math:numref}`eq:stokes-1`, making the system of these two equations symmetric.
And secondly, because the two equations are now linear in pressure and velocity (assuming that the viscosity $\eta$ and the density $\rho$ are considered fixed).
In addition, one can drop all terms involving $\nabla \cdot \mathbf u$ from the left hand side of the momentum equation {math:numref}`eq:stokes-1`; while dropping these terms does not affect the solution of the equations, it makes assembly of linear systems faster.

From a physical perspective, the assumption that the density is constant in the continuity equation but variable in the momentum equation is of course inconsistent.
However, it is justified if the variation is small since the momentum equation can be rewritten to read
```{math}
\begin{gathered}
  -\nabla \cdot 2\eta \varepsilon(\mathbf u) + \nabla p' =
  (\rho-\rho_0) \mathbf g,
\end{gathered}
```
where $p'$ is the *dynamic* pressure and $\rho_0$ is the constant reference density.
This makes it clear that the true driver of motion is in fact the *deviation* of the density from its background value, however small this value is: the resulting velocities are simply proportional to the density variation, not to the absolute magnitude of the density.

As such, the Boussinesq approximation can be justified.
On the other hand, given the real pressures and temperatures at the bottom of the Earth's mantle, it is arguable whether the density can be considered to be almost constant.
Most realistic models predict that the density of mantle rocks increases from somewhere around 3300 at the surface to over 5000 kilogram per cubic meters at the core-mantle boundary, due to the increasing lithostatic pressure.
While this appears to be a large variability, if the density changes slowly with depth, this is not in itself an indication that the Boussinesq approximation will be wrong.
To this end, consider that the continuity equation can be rewritten as $\frac 1\rho \nabla \cdot (\rho \mathbf u)=0$, which we can multiply out to obtain
```{math}
\begin{gathered}
  \nabla \cdot \mathbf u  +  \frac 1\rho \mathbf u \cdot \nabla \rho = 0.
\end{gathered}
```
The question whether the Boussinesq approximation is valid is then whether the second term (the one omitted in the Boussinesq model) is small compared to the first.
To this end, consider that the velocity can change completely over length scales of maybe 10 km, so that $\nabla \cdot\mathbf u \approx \|u\| /10\text{ km}$.
On the other hand, given a smooth dependence of density on pressure, the length scale for variation of the density is the entire earth mantle, i.e., $\frac 1\rho \mathbf u \cdot \nabla\rho \approx \|u\| 0.5 / 3000 \text{ km}$ (given a variation between minimal and maximal density of 0.5 times the
density itself).
In other words, for a smooth variation, the contribution of the compressibility to the continuity equation is very small.
This may be different, however, for models in which the density changes rather abruptly, for example due to phase changes at mantle discontinuities.

## On almost linear models.

A further simplification can be obtained if one assumes that all coefficients with the exception of the density do not depend on the solution variables but are, in fact, constant.
In such models, one typically assumes that the density satisfies a relationship of the form $\rho=\rho(T)=\rho_0(1-\alpha(T-T_0))$ with a small thermal expansion coefficient $\alpha$ and a reference density $\rho_0$ that is attained at temperature $T_0$.
Since the thermal expansion is considered small, this naturally leads to the following variant of the Boussinesq model discussed above, with the replacement of $\bar \rho(z) \bar \alpha(z)$ with a constant ($\rho_0 \alpha$):
```{math}
:label: eq:stokes-1-Boussinesq-linear
  -\nabla \cdot \left[2\eta \varepsilon(\mathbf u)
                \right] + \nabla p' =
  -\rho_0 \alpha T \mathbf g   \qquad   \textrm{in $\Omega$},
```
```{math}
:label: eq:stokes-2-Boussinesq-linear
  \nabla \cdot \mathbf u = 0   \qquad   \textrm{in $\Omega$},
```
```{math}
:label: eq:temperature-Boussinesq-linear
\rho_0 C_p \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla T\right)
  - \nabla\cdot k\nabla T
  =  \rho H   \quad   \textrm{in $\Omega$}.
```
Note that the right hand side forcing term in {math:numref}`eq:stokes-1-Boussinesq-linear` is now only the deviation of the gravitational force from the force that would act if the material were at temperature $T_0$.

Under the assumption that all other coefficients are constant, one then arrives at equations in which the only nonlinear term is the advection term, $\mathbf u \cdot \nabla T$ in the temperature equation {math:numref}`eq:temperature-Boussinesq-linear`.
This facilitates the use of a particular class of time stepping schemes in which one does not solve the whole set of equations at once, iterating out nonlinearities as necessary, but instead in each time step solves first the Stokes system with the previous time step's temperature, and then uses the so-computed velocity to solve the temperature equation. These kind of time stepping schemes are often referred to as *operator splitting* methods.

:::{note}
ASPECT does not solve the equations in the way described in this paragraph, however, a particular operator splitting method was used in earlier ASPECT versions.
It first solves the Stokes equations and then uses a semi-explicit time stepping method for the temperature equation where diffusion is handled implicitly and advection explicitly. This algorithm is often called *IMPES* (it originated in the porous media flow community, where the acronym stands for *Im*plicit *P*ressure, *E*xplicit *S*aturation) and is explained in more detail in {cite:t}`kronbichler:etal:2012`. Since then the algorithm in ASPECT has been rewritten to use an implicit time stepping algorithm also for the temperature equation because this allows to use larger time steps.
:::

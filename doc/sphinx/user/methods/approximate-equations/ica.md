(sec:methods:approximate-equations:ica)=
# The isothermal/isentropic compression approximation (ICA)

In the compressible case and without the assumption of a reference state, the conservation of mass equation in equation {math:numref}`eq:stokes-2` is $\nabla \cdot \left( \rho \textbf{u} \right)= 0$, which is nonlinear and not symmetric to the $\nabla p$ term in the force balance equation {math:numref}`eq:stokes-1`, making solving and preconditioning the resulting linear and nonlinear systems difficult.
To make this work in ASPECT, we consequently reformulate this equation.
Dividing by $\rho$ and applying the product rule of differentiation gives
```{math}
\frac{1}{\rho} \nabla \cdot \left( \rho \textbf{u} \right) = \nabla \cdot \textbf{u} + \frac{1}{\rho} \nabla \rho \cdot  \textbf{u}.
```
We will now make two basic assumptions: First, the variation of the density $\rho(p,T,\mathbf x, \mathfrak c)$ is dominated by the dependence on the (total) pressure; in other words, $\nabla \rho \approx \frac{\partial \rho}{\partial   p}\nabla p$.
This assumption is primarily justified by the fact that, in the Earth's mantle, the density increases by at least 50% between Earth's crust and the core-mantle boundary due to larger pressure there.
Secondly, we assume that the pressure is dominated by the static pressure, which implies that $\nabla p \approx \nabla p_s \approx \rho \textbf{g}$.
This is justified, because the viscosity in the Earth is large and velocities are small, hence $\nabla p' \ll \nabla p_s$. This finally allows us to write
```{math}
\frac{1}{\rho} \nabla \rho \cdot \textbf{u} \approx \frac{1}{\rho} \frac{\partial \rho}{\partial p} \nabla p \cdot \textbf{u} \approx \frac{1}{\rho} \frac{\partial \rho}{\partial p} \nabla p_s \cdot \textbf{u} \approx \frac{1}{\rho} \frac{\partial \rho}{\partial p} \rho \textbf{g} \cdot \textbf{u}
```
so we get
```{math}
:label: eq:stokes-2-compressible
\nabla \cdot \textbf{u} = - \frac{1}{\rho} \frac{\partial \rho}{\partial p} \rho \textbf{g} \cdot \textbf{u}
```
where $\frac{1}{\rho} \frac{\partial \rho}{\partial p}$ is often referred to as the compressibility.
Note that we have not yet made any assumptions about the change in temperature with pressure; we need to do this in order to calculate the compressibility.
There are two simple choices we could make; either to ignore adiabatic heating and use the isothermal compressibility:
```{math}
\nabla \cdot \textbf{u} = -\frac{1}{\rho} \frac{\partial \rho}{\partial p}_T \rho \textbf{g} \cdot \textbf{u} = -\beta_T \rho \textbf{g} \cdot \textbf{u}
```
or to assume that heating is everywhere adiabatic and use the isentropic compressibility:
```{math}
\nabla \cdot \textbf{u} = -\frac{1}{\rho} \frac{\partial \rho}{\partial p}_S \rho \textbf{g} \cdot \textbf{u} = -\beta_S \rho \textbf{g} \cdot \textbf{u}
```
Both choices are possible in ASPECT, the user simply needs to specify their preferred compressibility in the material model.
The isentropic compressibility is likely to be the more accurate approximation in models of mantle convection.

For this approximation, Equation {math:numref}`eq:stokes-2-compressible` replaces Equation {math:numref}`eq:stokes-2`. It has the advantage that it retains the symmetry of the Stokes equations if we can treat the right hand side of {math:numref}`eq:stokes-2-compressible` as known.
We do so by evaluating $\rho$ and $\mathbf u$ using the solution from the last time step (or values extrapolated from previous time steps), or using a nonlinear solver scheme.

:::{note}
This is the default approximation ASPECT uses to model compressible convection, see {ref}`sec:methods:combined-formulations`. The approximation is named "isothermal compression" for historical reasons, but the compressibility can be either isentropic or isothermal.
:::

(sec:methods:approximate-equations:ala)=
# The anelastic liquid approximation (ALA)

The *anelastic liquid approximation (ALA)* is based on two assumptions.
First, that the density variations relative to the adiabatic reference state at any given depth $\rho(p,T)-\bar\rho (z)$ are small and in particular can be accurately described by a Taylor expansion in pressure and temperature {cite}`schubert:etal:2001`:
```{math}
\begin{aligned}
  \rho(p,T) &\approx  \bar\rho
  + \left( \frac{\partial \rho(\bar p,\bar T)}{\partial T} \right)_{p} T'
  + \left( \frac{\partial \rho(\bar p,\bar T)}{\partial P} \right)_{T} p' \\
  \left( \frac{\partial \rho(\bar p,\bar T)}{\partial T} \right)_{p} &= -\bar \alpha \bar \rho(\bar p,\bar T) \\
  \left( \frac{\partial \rho(\bar p,\bar T)}{\partial P} \right)_{T} &= \bar \beta_T \bar \rho(\bar p,\bar T)
\end{aligned}
```
where $\bar \alpha$ is the thermal expansion coefficient ($\alpha = -\frac{1}{\rho}\left(\frac{\partial \rho}{\partial T}\right)_p$) and $\bar \beta_T$ is the isothermal compressibility ($\beta_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial p}\right)_T$), both on the adiabatic reference curve.
The subscripts ($p$ or $T$) indicate the variable that is held fixed. The second assumption is that the variation of the density from the reference density can be neglected in the mass balance and temperature equations.
This yields the following system of equations for the velocity and pressure equations:
```{math}
:label: eq:stokes-ALA-1
  -\nabla \cdot \left[2\eta \left(\varepsilon(\mathbf u)
                                  - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)
                \right] + \nabla p' =
  \bar \rho \left(\bar \beta_T p' - \bar \alpha T' \right) \mathbf g
   \qquad   \textrm{in $\Omega$},
```
```{math}
:label: eq:stokes-ALA-2
  \nabla \cdot (\bar\rho \mathbf u) = 0
   \qquad  \textrm{in $\Omega$}.
```

For the temperature equation, using the definition of the hydrostatic pressure gradient {math:numref}`eq:hydrostatic-pressure`, we arrive at the following:

```{math}
:label: eq:temperature-ala
\begin{gathered}
  \bar\rho C_p \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla
  T\right) - \nabla\cdot k\nabla T \\
  = \bar\rho H + 2\eta
  \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right):
  \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)
  +\alpha \bar\rho T (\mathbf u \cdot \mathbf g)
  \quad   \textrm{in $\Omega$}.
\end{gathered}
```

:::{note}
Our energy equation is formed in terms of $T$, while in the literature, the equation has sometimes been formulated in terms of $T'$, which yields additional terms containing $\bar{T}$ on the right-hand side.
Both ways of writing the equation are equivalent.
:::

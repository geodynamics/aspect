(sec:methods:initial-conditions)=
# Initial conditions and the adiabatic pressure/temperature

Equations {math:numref}`eq:stokes-1`-{math:numref}`eq:temperature` require us to pose initial conditions for the temperature, and this is done by selecting one of the existing models for initial conditions in the input parameter file, see {ref}`parameters:Initial_20temperature_20model`.
The equations themselves do not require that initial conditions are specified for the velocity and pressure variables (since there are no time derivatives on these variables in the model).

Nevertheless, a nonlinear solver will have difficulty converging to the correct solution if we start with a completely unphysical pressure for models in which coefficients such as density $\rho$ and viscosity $\eta$ depend on the pressure and temperature.
To this end, ASPECT uses pressure and temperature fields $p_{\textrm{ad}}(z), T_{\textrm{ad}}(z)$ computed in the adiabatic conditions model (see {ref}`parameters:Adiabatic_20conditions_20model`).
By default, these fields satisfy adiabatic conditions:
```{math}
\begin{aligned}
  \rho C_p \frac{\textrm{d}}{\textrm{d}z} T_{\textrm{ad}}(z)
  &=
  \frac{\partial\rho}{\partial T} T_{\textrm{ad}}(z) g_z,
\\
  \frac{\textrm{d}}{\textrm{d}z} p_{\textrm{ad}}(z)
  &=
  \rho g_z,
\end{aligned}
```
where strictly speaking $g_z$ is the magnitude of the vertical component of the gravity vector field, but in practice we take the magnitude of the entire gravity vector.

These equations can be integrated numerically starting at $z=0$, using the depth dependent gravity field and values of the coefficients $\rho=\rho(p,T,z), C_p=C_p(p,T,z)$.
As starting conditions at $z=0$, we choose a pressure $p_{\textrm{ad}}(0)$ equal to the average surface pressure (often chosen to be zero, see {ref}`sec:methods:pressure-norm`), and an adiabatic surface temperature $T_{\textrm{ad}}(0)$ that is also selected in the input parameter file.

However, users can also supply their own adiabatic conditions models or define an arbitrary profile using the "function" plugin.

:::{note}
The adiabatic surface temperature is often chosen significantly higher than the actual surface temperature.
For example, on Earth, the actual surface temperature is on the order of 290 K, whereas a reasonable adiabatic surface temperature is maybe 1600 K.
The reason is that the bulk of the mantle is more or less in thermal equilibrium with a thermal profile that corresponds to the latter temperature, whereas the very low actual surface temperature and the very high bottom temperature at the core-mantle boundary simply induce a thermal boundary layer.
Since the temperature and pressure profile we compute using the equations above are simply meant to be good starting points for nonlinear solvers, it is important to choose this profile in such a way that it covers most of the mantle well; choosing an adiabatic surface temperature of 290 K would yield a temperature and pressure profile that is wrong almost throughout the entire mantle.
:::

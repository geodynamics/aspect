
# Reference state: The adiabatic profile

The reference temperature profile $\bar{T}$, reference density profile $\bar{\rho}$ and the reference pressure $\bar{p}$ are computed in the adiabatic conditions model (provided by the class `AdiabaticConditions`, see {ref}`sec:methods:initial-conditions`).
By default, these fields satisfy adiabatic conditions (if adiabatic heating is included in the model, see {ref}`parameters:Heating_20model/Adiabatic_20heating`):

```{math}
\begin{aligned}
  \frac{\textrm{d} \bar{T}(z)}{\textrm{d}z}  &=  \frac{\alpha \bar{T}(z) g_z}{C_p}, \\
  \frac{\textrm{d} \bar{p}(z)}{\textrm{d}z} &= \bar\rho g_z,\\
  \bar{\rho} &= \bar\rho (\bar{p}, \bar{T}, z) \qquad \text{(as defined by the material model)},
\end{aligned}
```
where strictly speaking $g_z$ is the magnitude of the vertical component of the gravity vector field, but in practice we take the magnitude of the entire gravity vector.
If there is no adiabatic heating in the model, $\bar{T}$ is constant by default and set to the adiabatic surface temperature.
The density gradient is always computed by a simple finite difference approximation of the depth derivative of $\bar{\rho}$.

However, users can also supply their own adiabatic conditions models or define an arbitrary profile using the "function" plugin, which allows the user to define arbitrary functions for $\bar{T}(z)$, $\bar{p}(z)$ and $\bar{\rho}(z)$, see {ref}`parameters:Adiabatic_20conditions_20model`.

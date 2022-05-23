# A comment on adiabatic heating

Other codes and texts sometimes make a simplification to the adiabatic heating term in the previous equation.
If you assume the vertical component of the gradient of the *dynamic* pressure to be small compared to the gradient of the *total* pressure (in other words, the gradient is dominated by the gradient of the hydrostatic pressure), then $-\rho \mathbf g \approx \nabla \mathbf{p}$, and we have the following relation (the negative sign is due to $\mathbf g$ pointing downwards)
```{math}
 \alpha T \left( \mathbf u \cdot \nabla \mathbf p \right) \approx -\alpha \rho T \mathbf u \cdot \mathbf g.
```
While this simplification is possible, it is not necessary if you have access to the total pressure.
ASPECT therefore by default implements the original term without this simplification, but allows to simplify this term by setting the "`Use simplified adiabatic heating`" parameter in {ref}`parameters:Heating_20model/Adiabatic_20heating`.

(sec:advection-stabilization)=
# Advection Stabilization

ASPECT implements several advection schemes for
the temperature and compositional field equations. Specifically, the parameter
{ref}`parameters:Discretization/Stabilization_20parameters/Stabilization_20method`
allows using one of the following methods:

-   Entropy Viscosity Stabilization

-   SUPG Stabilization

Both add additional terms to the temperature (or compositional field)
equation. We will discuss the case for the temperature equation here. The
compositional fields only differ in having a zero conductivity, fewer
right-hand side terms, and $\rho C_p=1$. The strong form of the temperature
equation reads
```{math}
\rho C_p \frac{\partial T}{\partial t} + \rho C_p \mathbf{u} \cdot \nabla T - \nabla \cdot k\nabla T = F,
```
where $F$ is the combination of source and reaction terms, while the weak
form - with test function $\varphi$ and L2 inner product $(\cdot,\cdot)$ - is
```{math}
:label: eqn:weak-form-for-advection
a(T,\varphi) =
 \left(\rho C_p \frac{\partial T}{\partial t}, \varphi \right)
 + \left(\rho C_p \mathbf{u} \cdot \nabla T, \varphi \right)
 + \left( k \nabla T, \nabla \varphi \right) = (F,\varphi) = f(\varphi).
```

:::{toctree}
supg.md
entropy-viscosity.md
:::


# Constitutive laws

Equation {math:numref}`eq:stokes-1` describes buoyancy-driven flow in an isotropic fluid where strain rate is related to stress by a scalar (possibly spatially variable) multiplier, $\eta$.
For some material models it is useful to generalize this relationship to anisotropic materials, or other exotic constitutive laws.
For these cases ASPECT can optionally include a generalized, fourth-order tensor field as a material model state variable which changes equation {math:numref}`eq:stokes-1` to
```{math}
:label: eq:stokes-1-anisotropic
  -\nabla \cdot \left[2\eta \left(C \varepsilon(\mathbf u) - \frac{1}{3}(tr(C \varepsilon(\mathbf u)))\mathbf 1\right) \right] + \nabla p = \rho \mathbf g   \qquad  \textrm{in $\Omega$}
```
and the shear heating term in equation {math:numref}`eq:temperature` to
```{math}
:label: eq:temperature-anisotropic
\begin{aligned}
  \dots  \notag
  \\ + 2 \eta \left(C \varepsilon(\mathbf u) - \frac{1}{3}(tr(C \varepsilon(\mathbf u)))\mathbf 1\right)  : \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)
  \\ \dots \notag
\end{aligned}
```
where $C = C_{ijkl}$ is defined by the material model.
For physical reasons, $C$ needs to be a symmetric rank-4 tensor: i.e., when multiplied by a symmetric (strain rate) tensor of rank 2 it needs to return another symmetric tensor of rank 2.
In mathematical terms, this means that $C_{ijkl}=C_{jikl}=C_{ijlk}=C_{jilk}$.
Energy considerations also require that $C$ is positive definite: i.e., for any $\varepsilon \neq 0$, the scalar $\varepsilon : (C \varepsilon)$ must be positive.

This functionality can be optionally invoked by any material model that chooses to define a $C$ field, and falls back to the default case ($C=\mathbb I$) if no such field is defined.
It should be noted that $\eta$ still appears in equations {math:numref}`eq:stokes-1-anisotropic` and {math:numref}`eq:temperature-anisotropic`.
$C$ is therefore intended to be thought of as a "director" tensor rather than a replacement for the viscosity field, although in practice either interpretation is okay.

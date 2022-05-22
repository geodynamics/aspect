
# SUPG Stabilization

For streamline upwind/Petrov-Galerkin (SUPG) (see for example
{cite:t}`john:knobloch:2006,clevenger:heister:2019`), we add to the weak form
$a(\cdot,\cdot)$ the cell-wise defined weak form
```{math}
a_{\text{SUPG}} (T, \varphi) =
 \sum_{K \in \mathcal{T}_h}
  \delta_K \left( \rho C_p \frac{\partial T}{\partial t} - k \triangle T + \mathbf{\beta} \cdot \nabla T - F, \mathbf{\beta} \cdot \nabla \varphi \right)_K,
```
where $K \in \mathcal{T}_h$ are the cells in the computation,
$\delta_K \geq 0$ is a stabilization coefficient defined on each cell,
$\mathbf{\beta} = \rho C_p \mathbf{u}$ is the effective advection velocity.
The standard literature about SUPG does not contain $\rho C_p$, so it makes
sense to include this in the velocity. The first argument in the inner product
is the strong form of the residual of PDE, which is tested with the expression
$\mathbf{\beta} \cdot \nabla \varphi$ representing the solution in streamline
direction. We have to assume $k$ to be constant per cell, as we can not
compute the spatial derivatives easily.

For the implementation, $\frac{\partial T}{\partial t}$ is replaced by the
BDF2 approximation, and its terms from older timesteps and $-F$, are moved to
the right-hand side of the PDE.

We use the parameter design presented in {cite:t}`john:knobloch:2006` for
$\delta_K$:
```{math}
\delta_K = \frac{h}{2d\|\mathbf{\beta}\|_{\infty,K}} \left( \coth(Pe)-\frac{1}{Pe} \right)
```
where the Peclet number is given by
```{math}
Pe = \frac{ h \| \mathbf{\beta} \|_{\infty,K}}{2 d k_\text{max}},
```
$d$ is the polynomial degree of the temperature or composition element (typically 2),
$\coth(x) = (1+\exp(-2x)) / (1-\exp(-2x)),$ and
$k_\text{max}=\| k \|_{\infty, K}$ is the maximum conductivity in the cell $K$.

If $Pe<1$, the equation is diffusion-dominated and no stabilization is needed,
so we set $\delta_K=0$. Care needs to be taken in the definition if
$\| \beta \|$ or $k$ become zero:

1.  If $k$ is zero, then $Pe=\infty$ and the right part of the product in the
    definition of $\delta_K$ is equal to one.

2.  If $\| \beta \|$ is zero, $Pe < 1$, so we set $\delta_K=0$.

3.  If both are zero, no stabilization is needed (the field remains constant).

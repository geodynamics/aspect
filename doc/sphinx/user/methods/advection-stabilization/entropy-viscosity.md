
# Entropy viscosity

The entropy viscosity method
({cite}`guermond:etal:2011,kronbichler:etal:2012`) adds an artificial diffusion $\nu_h$
to the weak form {math:numref}`eqn:weak-form-for-advection`, where the diffusion
term $\left (k\nabla T, \nabla \varphi \right)$ is replaced by
```{math}
\left(\max (k, \nu_h) \nabla T, \nabla \varphi \right).
```
The parameter $\nu_h$ is chosen as a constant per cell as
```{math}
v_h \vert_K = \min \left( v_h^\text{max} \vert_K, v_h^E \vert_K \right),
```
where $v_h^\text{max}$ is the maximum dissipation defined as
```{math}
v_h^\text{max} \vert_K = \alpha_\text{max} h \| \mathbf u \|_{\infty,K}
```
on each cell $K$ with parameter $\alpha_\text{max}$ (known as "beta"
in the parameter files, see
{ref}`parameters:Discretization/Stabilization_20parameters/beta`). By
itself, this is commonly known as a first-order viscosity stabilization
scheme, which is effective at stabilization, but too diffusive to be used by
itself. In fact, one can show that this reduces the convergence order of
smooth solutions to be only first order. This is avoided by taking the minimum with the entropy viscosity $v_h^E|_K$ above. It is defined as
```{math}
v_h^E \vert_K = \alpha_E \frac{h^2 \| r_E \|_{\infty, K}}{\| E - E_\text{avg} \|_{\infty, \Omega}}.
```
The constant $\alpha_E$ is given by "cR" in the parameter files,
see {ref}`parameters:Discretization/Stabilization_20parameters/cR`. In the
denominator, the entropy viscosity above is scaled by the maximum deviation of
the temperature entropy $E=\frac{1}{2}(T-T_m)^2$ with
$T_m = \frac{1}{2}(T_\text{min}+T_\text{max})$ from the spatial average
$E_\text{avg} = \frac{1}{| \Omega |}\int E \;\text{d}x$. The residual $r_E$ of
the entropy equation for $E$ is defined as
```{math}
r_E = \frac{\partial E}{\partial t} + (T-T_m)(\mathbf{u}\cdot \nabla T - k\triangle T - F).
```
This residual is defined in such a way, that it is zero for the exact
solution, large where the numerical approximation is poor (for example in
areas with strong gradients), and small in areas where the numerical
approximation is good.

The above definition assumes the entropy residual exponent
("alpha" in the parameter files, see
{ref}`parameters:Discretization/Stabilization_20parameters/alpha`38]) is set
to 2 (the default and recommended). For the choice of 1 for "alpha," the
entropy viscosity is defined as
```{math}
v_h^E \vert_K = \alpha_E \frac{h |\Omega| \cdot \| \mathbf u \|_{\infty,K} \cdot \| r_E \|_{\infty, K}}
 {\| \mathbf u \|_{\infty,\Omega} \cdot (T_\text{max} - T_\text{min})}.
 ```
instead.

An additional parameter is the strain rate scaling factor "gamma"
(see {ref}`parameters:Discretization/Stabilization_20parameters/gamma`),
which changes the definition of the maximum dissipation $\nu_h^\text{max}$ to
```{math}
v_h^\text{max} \vert_K = \alpha_\text{max} h \|\lvert\mathbf u\rvert + \gamma h_K \lvert\varepsilon (\mathbf u)\rvert\|_{\infty,K},
```
where $\gamma\geq 0$ is the aforementioned parameter in front of the strain
rate.

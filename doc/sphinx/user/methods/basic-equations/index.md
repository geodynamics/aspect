(sec:methods:basic-equations)=
# Basic equations

ASPECT solves a system of equations in a $d=2$- or $d=3$-dimensional domain $\Omega$ that describes the motion of a highly viscous fluid driven by differences in the gravitational force due to a density that depends on the temperature.
In the following, we largely follow the exposition of this material in {cite:t}`schubert:etal:2001`.

Specifically, we consider the following set of equations for velocity $\mathbf{u}$, pressure $p$ and temperature $T$, as well as a set of advected quantities $c_i$ that we call *compositional fields*:
```{math}
:label: eq:stokes-1
  -\nabla \cdot \left[2\eta \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right) \right] + \nabla p = \rho \mathbf g \textrm{  in $\Omega$},
```
```{math}
:label: eq:stokes-2
  \nabla \cdot (\rho \mathbf u) = 0 \textrm{  in $\Omega$},
```
```{math}
:label: eq:temperature
\begin{aligned}
  \rho C_p \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla T\right) - \nabla\cdot k\nabla T &= \rho H \notag \\
  &\quad + 2\eta \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right): \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right) \\
  & \quad +\alpha T \left( \mathbf u \cdot \nabla p \right) \notag \\
  &\quad + \rho T \Delta S \left(\frac{\partial X}{\partial t} + \mathbf u\cdot\nabla X\right) \textrm{in $\Omega$}, \notag
\end{aligned}
```

```{math}
:label: eq:compositional
  \frac{\partial c_i}{\partial t} + \mathbf u\cdot\nabla c_i = q_i \textrm{  in $\Omega$}, i=1\ldots C
```
where $\varepsilon(\mathbf u) = \frac{1}{2}(\nabla \mathbf u + \nabla\mathbf u^T)$ is the symmetric gradient of the velocity (often called the *strain rate*)[^footnote1].

In this set of equations, {math:numref}`eq:stokes-1` and {math:numref}`eq:stokes-2` represent the compressible Stokes equations in which $\mathbf u=\mathbf u(\mathbf x,t)$ is the velocity field and $p=p(\mathbf x,t)$ the pressure field.
Both fields depend on space $\mathbf x$ and time $t$.
Fluid flow is driven by the gravity force that acts on the fluid and that is proportional to both the density of the fluid and the strength of the gravitational pull.

Coupled to this Stokes system is equation {math:numref}`eq:temperature` for the temperature field $T=T(\mathbf x,t)$ that contains heat conduction terms as well as advection with the flow velocity $\mathbf u$.
The right hand side terms of this equation correspond to

* internal heat production for example due to radioactive decay;
* friction heating;
* adiabatic compression of material;
* phase change.

The last term of the temperature equation corresponds to the latent heat generated or consumed in the process of phase change of material.
The latent heat release is proportional to changes in the fraction of material $X$ that has already undergone the phase transition (also called phase function) and the change of entropy $\Delta S$.
This process applies both to solid-state phase transitions and to melting/solidification.
Here, $\Delta S$ is positive for exothermic phase transitions.
As the phase of the material, for a given composition, depends on the temperature and pressure, the latent heat term can be reformulated:
```{math}
\begin{aligned}
\frac{\partial X}{\partial t} + \mathbf u\cdot\nabla X = \frac{DX}{Dt} = \frac{\partial X}{\partial T} \frac{DT}{Dt} + \frac{\partial X}{\partial p} \frac{Dp}{Dt} = \frac{\partial X}{\partial T} \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla T \right) + \frac{\partial X}{\partial p} \mathbf u\cdot\nabla p.
\end{aligned}
```
The last transformation results from the assumption that the flow field is always in equilibrium and consequently $\partial p/\partial t=0$ (this is the same assumption that underlies the fact that equation {math:numref}`eq:stokes-1` does not have a term $\partial \mathbf u / \partial t$).
With this reformulation, we can rewrite {math:numref}`eq:temperature` in the following way in which it is in fact implemented:
```{math}
:label: eq:temperature-reformulated
\begin{aligned}
  \left(\rho C_p - \rho T \Delta S \frac{\partial X}{\partial T}\right) \left(\frac{\partial T}{\partial t} + \mathbf u\cdot\nabla T\right) - \nabla\cdot k\nabla T &= \rho H \notag \\
  &\quad +  2\eta \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right): \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right) \\
  &\quad +\alpha T \left( \mathbf u \cdot \nabla p \right) \notag  \\
  &\quad  + \rho T \Delta S \frac{\partial X}{\partial p} \mathbf u\cdot\nabla p \quad \textrm{in $\Omega$}.  \notag
\end{aligned}
```

The last of the equations above, equation {math:numref}`eq:compositional`, describes the evolution of additional fields that are transported along with the velocity field $\mathbf u$ and may react with each other and react to other features of the solution, but that do not diffuse.
We call these fields $c_i$ *compositional fields*, although they can also be used for other purposes than just tracking chemical compositions.
We will discuss this equation in more detail in {ref}`sec:methods:compositional-fields`.

:::{toctree}
adiabatic-heating.md
bc.md
2d-models.md
:::

[^footnote1]: There is no consensus in the sciences on the notation used for strain and strain rate.
The symbols $\varepsilon$, $\dot\varepsilon$,  $\varepsilon(\mathbf u)$, and $\dot\varepsilon(\mathbf u)$, can all be found.
In this manual, and in the code, we will consistently use $\varepsilon$ as an *operator*, i.e., the symbol is not used on its own but only as applied to a field.
In other words, if $\mathbf u$ is the velocity field, then $\varepsilon(\mathbf u) = \frac{1}{2}(\nabla \mathbf u + \nabla\mathbf u^T)$ will denote the strain rate.
On the other hand, if $\mathbf d$ is the displacement field, then $\varepsilon(\mathbf d) = \frac{1}{2}(\nabla \mathbf d + \nabla\mathbf d^T)$ will denote the strain.

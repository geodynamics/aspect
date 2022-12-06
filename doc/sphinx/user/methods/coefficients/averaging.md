(sec:methods:coefficients:averaging)=
# Coefficient averaging

In multiphase rocks, or multirock areas in convection simulations, properties must be averaged because the length scales at which the rock types vary is far smaller than the resolution of the mesh.
As a consequence, we need to use &ldquo;effective coefficients,&rdquo; i.e., coefficients that do not correspond to any particular rock, but that lead to a macroscopic response that is a good match to the response of the correct, but unresolvable mixture of rocks.
For viscosity and conductivity, there is no single expression that describes how averaging should be performed; indeed, these properties are dependent on rock texture and mineral alignment, both of which may change through time as strain accumulates, and chemical diffusion and reactions take place.
Some of the existing multicomponent material models in <span class="smallcaps">ASPECT</span> allow the user to choose from a range of averaging schemes for viscosity.

In the case of density, thermal expansivity, heat capacity and bulk compressibility, there is one correct way of averaging.
Here we must consider conservation of mass and composition in a multicomponent rock $r$.
If component $i$ has masses $M_i$ and densities $\rho_i$, we can consider the summation of volume fractions:
```{math}
\begin{aligned}
V_r &=& \frac{M_r}{\rho_r} = \sum_i \frac{M_i}{\rho_i} \\
\frac{1}{\rho_r} &=& \sum_i \frac{x_i}{\rho_i}
\end{aligned}
```
where $x_i$ are mass fractions of the components in the rock.

Similarly, we can obtain averaging formulae for the other thermodynamic properties:
```{math}
\begin{aligned}
  \frac{\alpha}{\rho} &=& \sum_i x_i \frac{\alpha_i}{\rho_{i}} \\
  \frac{\beta_T}{\rho} &=& \sum_i x_i \frac{\beta_{Ti}}{\rho_{i}} \\
  C_p &=& \sum_i x_i C_{pi}
\end{aligned}
```

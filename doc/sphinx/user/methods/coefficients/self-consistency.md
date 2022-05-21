(sec:coefficient_self_consistency)=
# Coefficient self-consistency

*This section was contributed by Bob Myhill.*

The coefficients in the previous section may at first appear independent.
However, there are thermodynamic relations between these properties which must be satisfied in any self-consistent material model.
The following section describes the relations required for thermodynamic consistency, and presents some suggested ways by which consistency can be assured.

In order to derive the relationships between different material properties, we must introduce a thermodynamic potential known as the [*specific Gibbs free energy*](https://en.wikipedia.org/wiki/Gibbs_free_energy) $\mathcal{G}(p, T)$ with units $\text{ J}/\text{ kg}$.
The word "specific" indicates that the energy is given per unit mass, rather than volume or number of atoms or molecules.
This potential is equal to the maximum amount of non-expansion work that can be extracted from a thermodynamically closed system.
At equilibrium conditions and fixed temperature and pressure, the Gibbs free energy is minimized.
The following equations provide the definitions and relationships between thermodynamic properties in terms of the specific Gibbs free energy:
```{math}
  S = - \left( \frac{\partial \mathcal{G}}{\partial T} \right)_{p},
```
```{math}
:label: eq:mm_density
  \frac{1}{\rho} = \left( \frac{\partial \mathcal{G}}{\partial p} \right)_{T},
```
```{math}
:label: eq:mm_alpha_g
  \frac{\alpha}{\rho} = \frac{\partial^2 \mathcal{G}}{\partial {p} \, \partial {T}},
```
```{math}
:label: eq:mm_betaT_g
  \beta_T = -\rho \left( \frac{\partial^2 \mathcal{G}}{\partial {p}^2}  \right)_{T},
```
```{math}
:label: eq:mm_isobaric_heat_capacity
  C_p = -T \left( \frac{\partial^2 \mathcal{G}}{\partial {T}^2}  \right)_{p},
```
```{math}
:label: eq:mm_isentropic_compressibility
  \beta_S = \beta_T - \frac{\alpha^2 T}{\rho C_p},
```
```{math}
  \frac{C_V}{C_p} = \frac{\beta_S}{\beta_T},
```
```{math}
  \gamma = \frac{\alpha }{\beta_T \rho C_V}.
```
where $S$ is the specific entropy, $C_p$ and $C_V$ are the specific isobaric and isochoric heat capacities, $\beta_T$ and $\beta_S$ are the isothermal and isotropic compressibilities, and $\gamma$ is the thermodynamic Gr&uuml;neisen parameter.
The subscript indicates the thermodynamic variable ($p$ or $T$) that is held constant.

Thermodynamically self-consistent material models must obey the explicit and
implicit relations between the different properties *at all pressures and
temperatures*. Explicit relations are here defined as those between properties
and their derivatives, such as that between density and thermal expansivity.
Implicit relations involve mixed pressure and temperature derivatives, and
derive from the symmetry of second derivatives. The following paragraphs list
the relations most relevant for the construction of
thermodynamically-consistent material models in <span
class="smallcaps">ASPECT</span>.

## Consistency in $\boldsymbol{\rho}$-$\boldsymbol{\alpha}$ and $\boldsymbol{\rho}$-$\boldsymbol{\beta_T}$

Using the chain rule to combine {math:numref}`eq:mm_density`, {math:numref}`eq:mm_alpha_g` and {math:numref}`eq:mm_betaT_g` yields the more familiar definitions of $\alpha$ and $\beta_T$:
```{math}
:label: eq:mm_thermal_expansivity
\begin{aligned}
  \alpha &=& -\frac{1}{\rho} \left( \frac{\partial \rho}{\partial T} \right)_{p},
\end{aligned}
```
```{math}
:label: eq:mm_isothermal_compressibility
\begin{aligned}
  \beta_T &=& \frac{1}{\rho} \left( \frac{\partial \rho}{\partial p} \right)_{T}.
\end{aligned}
```
## Isobaric heat capacity

We start by taking the partial derivative of the isobaric heat capacity {math:numref}`eq:mm_isobaric_heat_capacity` with respect to pressure at constant temperature:
```{math}
:label: eq:heat_capacity_p_dependence
\begin{aligned}
  \left( \frac{\partial C_{p}}{\partial p} \right)_{T} &=& -T \frac{\partial^3 \mathcal{G}}{\partial {T}^2 \, \partial {p}} \\
  &=& -T \left( \frac{\partial \left(\alpha / \rho \right)}{\partial T} \right)_{p}.
\end{aligned}
```
From this expression it becomes clear that if $\alpha / \rho$ has any temperature dependence, the heat capacity $C_p$ *cannot* be globally constant.
One way to solve this issue is to define heat capacity at constant pressure, and then integrate {math:numref}`eq:heat_capacity_p_dependence` with respect
to pressure:
```{math}
C_p(p, T) = C_p(p_{\textrm{ref}}, T) -T \int_{p_{\textrm{ref}}}^p \left(\frac{\partial \left(\alpha / \rho \right)}{\partial T} \right)_{p} \text{d}p.
```
There is no guarantee that this expression will have a form for which the integral can be found analytically.

## Isentropic gradient

The material properties also define the slope of the adiabat (the change in temperature with pressure at constant entropy) at all pressures and temperatures.
Using the cyclic relation, we can define this slope in terms of partial differentials of the entropy with respect to pressure and temperature:
```{math}
\begin{aligned}
  \left( \frac{\partial T}{\partial p} \right)_{S} &=& - \left( \frac{\partial T}{\partial S} \right)_{p} \left( \frac{\partial S}{\partial p} \right)_{T} \\
  &=& - \left( \frac{T}{C_p} \right) \left( - \frac{\alpha}{\rho} \right) \\
  &=& \frac{\alpha T}{\rho C_p} \label{eq:mm_isentropic_gradient}
\end{aligned}
```
This expression does not pose a constraint on the material properties, but in order to be self-consistent, the adiabat must be computed following this relation.

For complex material models, obtaining analytical functions which obey all these relations may be a non-trivial exercise.
Furthermore, it is often not immediately clear when a given formulation is thermodynamically inconsistent.
Indeed, both the thermodynamic and the geodynamic literature contain many equations of state and material parameterizations which do not obey these relations!
This may not invalidate the results obtained with these models, but it is a point worth keeping in mind as the geodynamics community moves to more complicated and more realistic parameterizations.

:::{warning}
Some compressible formulations in ASPECT ({ref}`sec:methods:choosing-a-formulation:mass-conservation-approx`) use the isothermal compressibility, while others use the isentropic compressibility.
Fully self-consistent material models must either specify what approximation of the compressible equations they are consistent with (see {ref}`sec:methods:approximate-equations`), or have a switch so that they use the correct compressibility for each of the different approximations.
The conversion between isothermal and isentropic compressibilities is given in {math:numref}`eq:mm_isentropic_compressibility`.
:::

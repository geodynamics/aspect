# Boundary conditions

Having discussed {math:numref}`eq:temperature`, let us come to the last one of the original set of equations, {math:numref}`eq:compositional`.
It describes the motion of a set of advected quantities $c_i(\mathbf x,t),i=1\ldots C$.
We call these *compositional fields* because we think of them as spatially and temporally varying concentrations of different elements, minerals, or other constituents of the composition of the material that convects.
As such, these fields participate actively in determining the values of the various coefficients of these equations.
On the other hand, ASPECT also allows the definition of material models that are independent of these compositional fields, making them passively advected quantities.
Several of the cookbooks in {ref}`cha:cookbooks` consider compositional fields in this way, i.e., essentially as tracer quantities that only keep track of where material came from.

These equations are augmented by boundary conditions that can either be of Dirichlet, Neumann, or tangential type on subsets of the boundary $\Gamma=\partial\Omega$:
```{math}
\begin{aligned}
  \mathbf u &= 0 & \qquad &\textrm{on $\Gamma_{0,\mathbf u}$}, \\
  \mathbf u &= \mathbf u_{\text{prescribed}} & \qquad &\textrm{on $\Gamma_{\text{prescribed},\mathbf u}$}, \\
  \mathbf n \cdot \mathbf u &= 0 & \qquad &\textrm{on $\Gamma_{\parallel,\mathbf u}$}, \\
  (2\eta \varepsilon(\mathbf u) -p I)\mathbf n  &= \mathbf t & \qquad &\textrm{on $\Gamma_{\text{traction},\mathbf u}$}, \\
  T &= T_{\text{prescribed}} & \qquad &\textrm{on $\Gamma_{D,T}$}, \\
  \mathbf n \cdot k\nabla T &= 0  & \qquad &\textrm{on $\Gamma_{N,T}$}. \\
\end{aligned}
```
```{math}
:label: eq:gamma-in-composition
  c_i = c_{i,\text{prescribed}}  \qquad \textrm{on $\Gamma_{\text{in}}=\{\mathbf x: \mathbf u\cdot\mathbf n<0\}$}.
```
Here, the boundary conditions for velocity and temperature are subdivided into disjoint parts:

-   $\Gamma_{0,\mathbf u}$ corresponds to parts of the boundary on which the velocity is fixed to be zero.

-   $\Gamma_{\text{prescribed},\mathbf u}$ corresponds to parts of the boundary on which the velocity is prescribed to some value (which could also be zero).
    It is possible to restrict prescribing the velocity to only  certain components of the velocity vector.

-   $\Gamma_{\parallel,\mathbf u}$ corresponds to parts of the boundary on which the velocity may be nonzero but must be parallel to the boundary, with the tangential component undetermined.

-   $\Gamma_{\text{traction},\mathbf u}$ corresponds to parts of the boundary on which the traction is prescribed to some surface force density (a common application being $\mathbf t=-p\mathbf n$ if one just wants to prescribe a pressure component).
    It is possible to restrict prescribing the traction to only certain vector components.

-   $\Gamma_{D,T}$ corresponds to places where the temperature is prescribed (for example at the inner and outer boundaries of the Earth's mantle).

-   $\Gamma_{N,T}$ corresponds to places where the temperature is unknown but the heat flux across the boundary is zero (for example on symmetry surfaces if only a part of the shell that constitutes the domain the Earth's mantle occupies is simulated).

We require that one of these boundary conditions hold at each point for both velocity and temperature, i.e., $\Gamma_{0,\mathbf u}\cup\Gamma_{{\text{prescribed}}\mathbf u}\cup\Gamma_{\parallel,\mathbf u}\cup\Gamma_{{\text{traction}}\mathbf u}=\Gamma$ and $\Gamma_{D,T}\cup\Gamma_{N,T}=\Gamma$.

Boundary conditions have to be imposed for the compositional fields only at those parts of the boundary where flow points inward, see equation {math:numref}`eq:gamma-in-composition`, but not where it is either tangential to the boundary or points outward.
The difference in treatment between temperature and compositional boundary conditions is due to the fact that the temperature equation contains a (possibly small) diffusion component, whereas the compositional equations do not.

There are other equations that ASPECT can optionally solve.
For example, it can deal with free surfaces (see {ref}`sec:methods:freesurface`), melt generation and transport (see {ref}`sec:methods:melt-transport`), and it can advect along particles (see {ref}`sec:methods:particles`).
These optional models are discussed in more detail in the indicated sections.


# Arbitrary Lagrangian-Eulerian implementation

The question of how to handle the motion of the mesh with a free surface is
challenging. Eulerian meshes are well behaved, but they do not move with the
fluid motions, which makes them difficult for use with free surfaces.
Lagrangian meshes do move with the fluid, but they quickly become so distorted
that remeshing is required. ASPECT implements
an Arbitrary Lagrangian-Eulerian (ALE) framework for handling motion of the
mesh. The ALE approach tries to retain the benefits of both the Lagrangian
and the Eulerian approaches by allowing the mesh motion $\textbf{u}_m$ to be
largely independent of the fluid. The mass conservation condition requires
that $\textbf{u}_m \cdot \textbf{n} = \textbf{u} \cdot \textbf{n}$ on the free
surface, but otherwise the mesh motion is unconstrained, and should be chosen
to keep the mesh as well behaved as possible.

ASPECT uses a Laplacian scheme for calculating
the mesh velocity. The mesh velocity is calculated by solving
```{math}
\begin{aligned}
-\Delta \textbf{u}_m &= 0 & \qquad & \textrm{in } \Omega, \\
\textbf{u}_m &= \left( \textbf{u} \cdot \textbf{n} \right) \textbf{n} & \qquad & \textrm{on } \partial \Omega_{\textrm{free surface}}, \\
\textbf{u}_m \cdot \textbf{n} &= 0 & \qquad & \textrm{on } \partial \Omega_{\textrm{free slip}}, \\
\textbf{u}_m &= 0 & \qquad & \textrm{on } \partial \Omega_{\textrm{Dirichlet}}.
\end{aligned}
```
After this mesh velocity is calculated, the mesh vertices are time-stepped
explicitly. This scheme has the effect of choosing a minimally distorting
perturbation to the mesh. Because the mesh velocity is no longer zero in the
ALE approach, we must then correct the Eulerian advection terms in the
advection system with the mesh velocity (see, e.g., {cite}`donea:etal:2004`).
For instance, the temperature equation
{math:numref}`eq:temperature-Boussinesq-linear` becomes

```{math}
\rho C_p \left(\frac{\partial T}{\partial t} + \left(\mathbf u - \mathbf u_m \right) \cdot\nabla T\right)
  - \nabla\cdot k\nabla T
  =
  \rho H
   \quad
   \textrm{in $\Omega$}.
```

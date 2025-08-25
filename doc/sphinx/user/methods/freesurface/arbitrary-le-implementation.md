
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
\textbf{u}_m &= \sum_i \textbf{u}_{m,i} & \qquad & \textrm{on } \partial \Omega_{\textrm{deformation}}, \\
\textbf{u}_m \cdot \textbf{n} &= 0 & \qquad & \textrm{on } \partial \Omega_{\textrm{tangential}}, \\
\textbf{u}_m &= 0 & \qquad & \textrm{on } \partial \Omega_{\textrm{remaining}}.
\end{aligned}
```

Here, $\partial \Omega_{\textrm{deformation}}$ are all boundaries with active
mesh deformation processes (specified in {ref}`parameters:Mesh_20deformation/Mesh_20deformation_20boundary_20indicators`),
$\partial \Omega_{\textrm{tangential}}$ are all boundaries
along which tangential mesh deformation is allowed (all boundaries in {ref}`parameters:Boundary_20velocity_20model/Tangential_20velocity_20boundary_20indicators`
and {ref}`parameters:Mesh_20deformation/Additional_20tangential_20mesh_20velocity_20boundary_20indicators`),
and $\partial \Omega_{\textrm{remaining}}$ are
all remaining boundaries where no mesh deformation is allowed.

$\textbf{u}_{m,i}$ are the individual mesh deformation processes that deform the
boundary of the mesh on $\partial \Omega_{\textrm{deformation}}$, whose
individual contributions are added to compute the final mesh deformation along these
boundaries. Examples for such deformation processes are a free surface
($\textbf{u}_{m,i} = \left( \textbf{u} \cdot \textbf{n} \right) \textbf{n}$), an analytically
prescribed mesh velocity ($\textbf{u}_{m,i} = \textbf{f}(\textbf{x},t)$), or more complicated
expressions like the solution of a diffusion equation for topography or a
landscape evolution model.

After this mesh velocity is calculated, the mesh vertices are time-stepped
explicitly. This scheme has the effect of choosing a minimally distorting
perturbation to the mesh. Because the mesh velocity is no longer zero in the
ALE approach, we must then correct all velocities that act on quantities
defined relative to the mesh position with the mesh velocity, e.g. the
Eulerian advection terms in the
advection system (see, e.g., {cite}`donea:etal:2004`).
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

This correction ensures that material only enters or leaves the domain through a deforming boundary
if the flow velocity $\mathbf{u}$ is not equal to the mesh velocity $\mathbf{u}_m$ at the boundary.
Please note that this also means prescribed velocity boundary conditions are defined in a stationary
coordinate system, not a coordinate system that moves with the mesh. For example an inflow velocity
of 5 cm/yr at a boundary that moves with 7 cm/yr into the model domain actually leads to an
outflow of material with a relative velocity of 2 cm/yr.

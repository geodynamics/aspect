(sec:methods:freesurface)=
# Free surface calculations

In reality the boundary conditions of a convecting Earth are not no-slip or
free slip (i.e., no normal velocity). Instead, we expect that a free surface
is a more realistic approximation, since air and water should not prevent the
flow of rock upward or downward. This means that we require zero stress on the
boundary, or $\sigma \cdot \textbf{n} = 0$, where
$\sigma = 2 \eta \varepsilon (\textbf{u})$. In general there will be flow
across the boundary with this boundary condition. To conserve mass we must
then advect the boundary of the domain in the direction of fluid flow. Thus,
using a free surface necessitates that the mesh be dynamically deformable.

:::{toctree}
arbitrary-le-implementation.md
stabilization.md
:::

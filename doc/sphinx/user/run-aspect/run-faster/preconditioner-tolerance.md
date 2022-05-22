# Adjusting solver preconditioner tolerances

To solve the Stokes equations it is necessary to lower the condition number of
the Stokes matrix by preconditioning it. In
ASPECT a right preconditioner
```{math}
Y^{-1} =
\begin{pmatrix}
\widetilde{A^{-1}} & -\widetilde{A^{-1}}B^{T}\widetilde{S^{-1}} \\
0 & \widetilde{S^{-1}}
\end{pmatrix}
```
 is used to precondition the system, where $\widetilde{A^{-1}}$
is the approximate inverse of the A block and $\widetilde{S^{-1}}$ is the
approximate inverse of the Schur complement matrix. Matrix
$\widetilde{A^{-1}}$ and $\widetilde{S^{-1}}$ are calculated through a CG
solve, which requires a tolerance to be set. In comparison with the solver
tolerances of the previous section, these parameters are relatively safe to
use, since they only change the preconditioner, but can speed up or slow down
solving the Stokes system considerably.

In practice $\widetilde{A^{-1}}$ takes by far the most time to compute, but is
also very important in conditioning the system. The accuracy of the
computation of $\widetilde{A^{-1}}$ is controlled by the parameter
`Linear solver A block tolerance` which has a default value of $1e-2$. Setting
this tolerance to a less strict value will result in more outer iterations,
since the preconditioner is not as good, but the amount of time to compute
$\widetilde{A^{-1}}$ can drop significantly resulting in a reduced total solve
time. The cookbook {ref}`sec:cookbooks:crustal-deformation` for example can be
computed much faster by setting the `Linear solver A block tolerance` to
$5e-1$. The calculation of $\widetilde{S^{-1}}$ is usually much faster and the
conditioning of the system is less sensitive to the parameter
`Linear solver S block tolerance`, but for some problems it might be worth it
to investigate.

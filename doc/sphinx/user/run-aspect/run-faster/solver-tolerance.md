# Adjusting solver tolerances

At the heart of every time step lies the solution of linear systems for the
Stokes equations, the temperature field, and possibly for compositional
fields. In essence, each of these steps requires us to solve a linear system
of the form $Ax=b$ which we do through iterative solvers, i.e., we try to find
a sequence of approximations $x^{(k)}$ where $x^{(k)}\rightarrow x=A^{-1}b$.
This iteration is terminated at iteration $k$ if the approximation is
"close enough" to the exact solution. The solvers we use determine
this by testing after every iteration whether the *residual*,
$r^{(k)}=A(x-x^{(k)})=b-Ax^{(k)}$, satisfies
$\|r^{(k)}\|\le\varepsilon\|r^{(0)}\|$ where $\varepsilon$ is called the
(relative) *tolerance*.

Obviously, the smaller we choose $\varepsilon$, the more accurate the
approximation $x^{(k)}$ will be. On the other hand, it will also take more
iterations and, consequently, more CPU time to reach the stopping criterion
with a smaller tolerance. The default value of these tolerances are chosen so
that the approximation is typically sufficient. You can make
ASPECT run faster if you choose these tolerances
larger. The parameters you can adjust are all listed in
{ref}`parameters:Solver_20parameters` and are located in the
`Solver parameters` subsection of the input file. In particular, the
parameters you want to look at are `Linear solver tolerance`,
`Temperature solver tolerance` and `Composition solver tolerance`.

All this said, it is important to understand the consequences of choosing
tolerances larger. In particular, if you choose tolerances too large, then the
difference between the exact solution of a linear system $x$ and the
approximation $x^{(k)}$ may become so large that you do not get an accurate
output of your model any more. A rule of thumb in choosing tolerances is to
start with a small value and then increase the tolerance until you come to a
point where the output quantities start to change significantly. This is the
point where you will want to stop.

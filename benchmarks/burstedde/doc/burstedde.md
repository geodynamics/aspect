# The Burstedde variable viscosity benchmark

*This section was contributed by Iris van Zelst.*

This benchmark is intended to test solvers for variable viscosity Stokes
problems. It begins with postulating a smooth exact polynomial solution to the
Stokes equation for a unit cube, first proposed by {cite:t}`DohrmannBochev2004`
and also described by {cite:t}`burstedde:etal:2013`:
```{math}
:label: eq:burstedde-velocity
  {\mathbf u} &= \left( \begin{array}{c}
      x+x^2+xy+x^3y \\
      y + xy + y^2 + x^2 y^2\\
      -2z - 3xz - 3yz - 5x^2 yz
      \end{array}
  \right)
```
```{math}
:label: eq:burstedde-pressure
  p = xyz + x^3 y^3z - \frac{5}{32}.
```
It is then trivial to verify that the velocity field is divergence-free. The
constant $-\frac{5}{32}$ has been added to the expression of $p$ to ensure
that the volume pressure normalization of can be used in this benchmark (in
other words, to ensure that the exact pressure has mean value zero and,
consequently, can easily be compared with the numerically computed pressure).
Following {cite:t}`burstedde:etal:2013`, the viscosity $\mu$ is given by the
smoothly varying function
```{math}
:label: eq:burstedde-mu
\mu = \exp\left\{1 - \beta\left[x (1-x) + y(1-y) + z(1-z)\right]\right\}.
```
The maximum of this function is $\mu = e$, for
example at $(x,y,z)=(0,0,0)$, and the minimum of this function is
$\mu = \exp \Big( 1-\frac{3\beta}{4}\Big)$ at $(x,y,z) = (0.5,0.5,0.5)$. The
viscosity ratio $\mu^\ast$ is then given by
```{math}
\mu^\ast = \frac{\exp\Big(1-\frac{3\beta}{4}\Big)}{\exp(1)} = \exp\Big(\frac{-3\beta}{4}\Big).
```
Hence, by varying $\beta$ between 1 and 20, a difference of up to 7 orders of
magnitude viscosity is obtained. $\beta$ will be one of the parameters that
can be selected in the input file that accompanies this benchmark.

The corresponding body force of the Stokes equation can then be computed by
inserting this solution into the momentum equation,
```{math}
:label: eq:burstedde-momentum
{\nabla} p - \nabla \cdot (2  \mu {\epsilon(\mathbf u)}) = \rho \mathbf g.
```
Using equations {math:numref}`eq:burstedde-velocity`, {math:numref}`eq:burstedde-pressure` and
{math:numref}`eq:burstedde-mu` in the momentum equation
{math:numref}`eq:burstedde-momentum`, the following expression for the body force
$\rho\mathbf g$ can be found:
```{math}
\begin{gathered}
  {\rho\mathbf g}
  =
  \left(
    \begin{array}{c}
      yz+3x^2y^3z\\
      xz +3x^3y^2z \\
      xy+x^3y^3
    \end{array}
  \right)
  -\mu
  \left(
    \begin{array}{c}
      2+6xy  \\
      2 + 2x^2 +  2y^2 \\
      -10yz
    \end{array}
  \right)  \\
  +
  (1-2x)\beta \mu
  \left(
    \begin{array}{c}
      2+4x+2y+6x^2y \\
      x+y+2xy^2+x^3 \\
      -3z -10xyz
    \end{array}
  \right)
  +(1-2y)\beta \mu
  \left(
    \begin{array}{c}
      x+y+2xy^2+x^3 \\
      2+2x+4y+4x^2y \\
      -3z-5x^2z \\
    \end{array}
  \right)
  \\
  +(1-2z)\beta \mu
  \left(
    \begin{array}{c}
      -3z -10xyz \\
      -3z-5x^2z \\
      -4-6x-6y-10x^2y
    \end{array}
  \right)\end{gathered}
  ```
 Assuming $\rho = 1$, the above expression translates
into an expression for the gravity vector $\mathbf g$. This expression for the
gravity (even though it is completely unphysical), has consequently been
incorporated into the `BursteddeGravity` gravity model that is described in
the `benchmarks/burstedde/burstedde.cc` file that accompanies this benchmark.

We will use the input file `benchmarks/burstedde/burstedde.prm` as input,
which is very similar to the input file `benchmarks/inclusion/adaptive.prm`
discussed above in {ref}`sec:benchmarks:inclusion`. The major
changes for the 3D polynomial Stokes benchmark are listed below:

```{literalinclude} burstedde.prm
```

The boundary conditions that are used are simply the velocities from equation
{math:numref}`eq:burstedde-velocity` prescribed on each boundary. The viscosity
parameter in the input file is $\beta$. Furthermore, in order to compute the
velocity and pressure $L_1$ and $L_2$ norm, the postprocessor
`BursteddePostprocessor` is used. Please note that the linear solver tolerance
is set to a very small value (deviating from the default value), in order to
ensure that the solver can solve the system accurately enough to make sure
that the iteration error is smaller than the discretization error.

Expected analytical solutions at two locations are summarized in
{numref}`tab:burstedde-table` and can be deduced from equations
{math:numref}`eq:burstedde-velocity` and {math:numref}`eq:burstedde-pressure`.
{numref}`fig:burstedde-benchmark` shows that the analytical
solution is indeed retrieved by the model.

```{figure-md} fig:burstedde-benchmark
<img src="burstedde.*" alt = "Burstedde results" width="100%"/>

Burstedde benchmark: Results for the 3D polynomial Stokes benchmark, obtained with a resolution of $16\times 16$ elements, with $\beta = 10$.
```

```{table} Analytical solutions
:name: tab:burstedde-table

| Quantity       | $\mathbf{r} = (0,0,0)$ | $\mathbf{r} = (1,1,1)$ |
|:---------------|:----------------------:|:----------------------:|
| $p$            |       $-0.15625$       |       $1.84375$        |
| $\mathbf{u}$   |       $(0,0,0)$        |      $(4,4,-13)$       |
| $\|\mathbf{u}\|$ |          $0$           |        $14.177$        |

```

The convergence of the numerical error of this benchmark has been analyzed by
playing with the mesh refinement level in the input file, and results can be
found in {numref}`fig:burstedde:errors`. The velocity shows cubic error convergence, while
the pressure shows quadratic convergence in the $L_1$ and $L_2$ norms, as one
would hope for using $Q_2$ elements for the velocity and $Q_1$ elements for
the pressure.

```{figure-md} fig:burstedde:errors
<img src="errors.*" alt = "Error convergence" width="100%"/>

Burstedde benchmark: Error convergence for the 3D polynomial Stokes benchmark.
```

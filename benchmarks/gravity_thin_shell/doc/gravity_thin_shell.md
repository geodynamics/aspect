(sec:benchmarks:thin-shell-gravity)=
# Thin shell gravity benchmark

*This section was contributed by Cedric Thieulot, Bart Root and Paul Bremner.*

The idea behind this benchmark is to test the accuracy of the gravity
postprocessor: the domain is a thin shell of constant density somewhere in the
Earth mantle and we wish to compute the resulting gravity field and potential
at satellite height.

The domain is a spherical shell of 10&nbsp;km radius centered at a depth $D$,
i.e. the inner radius is $R_{inner}=6371-D-5~\text{ km}$ and the
outer radius $R_\text{outer}=6371-D+5~ \text{ km}$. It is filled with a
fluid of constant density $\rho_0=3300~\text{ kg}/\text{m}^3$ with total
mass
$M=\frac43 \pi (R_\text{outer}^3-R_\text{inner}^3)\rho_0$.

Both the Stokes and energy equations solve are bypassed:

```{literalinclude} thinshell_b.prm
```

We make use of the Custom mesh subdivision option to generate a mesh where a
single cell is used in the radial direction (parameterized with the number of
'slices') while 6 blocks (default value) of $2^5\times 2^5$ cells
are used in the lateral direction. This gives a total of 6,144 mesh cells. For
$D=100~\text{ km}$ the parameterization of the mesh is then as follows:

```{literalinclude} thinshell_a.prm
```

We wish to compute the gravity potential $U$ and the gravity vector
${\mathbf g}$ at a radius of $R_{s}=6371+250=6621~\text{ km}$ (i.e.
the gravity experienced by a satellite flying at a height of $250~\text{ km}$
above the surface of the Earth) on a regular longitude-latitude grid spanning
the whole Earth, with a $2^\circ\times 2^\circ$ resolution:

```{literalinclude} thinshell_c.prm
```

The Gauss-Legendre Quadrature (GLQ) algorithm is central to the code as it is
used to compute the elemental integrals stemming from the discretization of
the weak form of the PDEs. Logically, the gravity postprocessor is also based
on the GLQ since the potential $U$ at any location in space is given by the
integral
```{math}
U({\mathbf r}) = \iiint_\Omega \frac{G \rho({\mathbf r}')}{|{\mathbf r}-{\mathbf r}'|} d{\mathbf r}'
= \sum_K \iiint_{\Omega_K} \frac{G \rho({\mathbf r}')}{|{\mathbf r}-{\mathbf r}'|} d{\mathbf r}',
```
where the sum runs over all cells of the mesh and
${G}=6.67430\times 10^{-11}~\text{ m}^3/\text{kg}/\text{s}^2$ is the
gravitational constant. The gravity acceleration vector is obtained via
${\mathbf g}=-{\mathbf \nabla} U$, or
```{math}
{\mathbf g}({\mathbf r}) =
\iiint_\Omega  \frac{G \rho({\mathbf r}')}{|{\mathbf r}-{\mathbf r}'|^2} d{\mathbf r}'
=
\sum_K \iiint_{\Omega_K} \frac{G \rho({\mathbf r}')}{|{\mathbf r}-{\mathbf r}'|^2} d{\mathbf r}'.
```

The default number of quadrature points for each cell in the postprocessor is
$2^\text{ndim}$, where ndim is the number of dimensions. The
$n$-point GLQ allows exact integration of $(2n-1)$-order polynomials. However,
the integrand of the Newton integral is not a polynomial (it contains a term
$\sim r^{-m}$), so there is not an optimal number of quadrature points to use.
Therefore, the postprocessor allows the user to choose how many additional
points per dimension are used with an expectation that an increase in the
number of quadrature points inside the cells leads to a more accurate
calculation. This increase number $I$ is set to 1 in the example above,
although it can also be chosen to be 0, 1, 2 or even -1. The case $I=-1$
approximates the cell as a point mass since there is a single quadrature point
in the middle of the cell.

Given the symmetry of the problem, the values of the potential and
acceleration depend solely on the distance $r$ from the origin, see {cite:t}`turcotte:schubert:2014`. Their analytical values are given
by $U(r) = - GM/r$ and $|{\mathbf g}(r)|= g_r(r)=GM/r^2$. Minimum, maximum,
and average values of both the potential and the acceleration are printed in
the statistics file while measurements at all the desired points are written
in the `output_gravity` folder inside the output folder. We ran the input file
for $I\in \{-1,0,1,2\}$ for $D=0,100,500,1500,3000~\text{ km}$ and the results
are presented in Table&nbsp;[1] alongside the analytical values.

```{table} Thin shell gravity benchmark: $1\text{ mGal}=10^{-5}\text{ m}/\text{s}^2$. $n_q$ is the number of GLQ points per cell. 'a.v.' stands for analytical value.
:name: tab:thin_shell_gravity_benchmark

|      |      |       |           |              |           |              |                                    |              |
|:----:|:----:|:-----:|:---------:|:------------:|:---------:|:------------:|:----------------------------------:|:------------:|
|  D   | $I$  | $n_q$ |           | $g_r$ (mGal) |           |              | $U$ (J&nbsp;kg<sup>&minus;1</sup>) |              |
| (km) |      |       |   avrg.   |     min      |    max    |    avrg.     |                min                 |     max      |
|  0   |  -1  | $1^3$ | 2563.6541 |  2530.6859   | 2607.5210 | -169764.4978 |            -169832.3647            | -169744.3149 |
|      |  0   | $2^3$ | 2562.8680 |  2553.5683   | 2571.9157 | -169676.4060 |            -169681.1065            | -169671.6199 |
|      |  1   | $3^3$ | 2562.6878 |  2561.9847   | 2563.7170 | -169676.3142 |            -169676.8060            | -169675.9148 |
|      |  2   | $4^3$ | 2562.6993 |  2562.6215   | 2562.7395 | -169676.3211 |            -169676.3361            | -169676.2903 |
|      | a.v. |       | 2562.6993 |              |           | -169676.3210 |                                    |              |
| 100  |  -1  | $1^3$ | 2484.0821 |  2479.3619   | 2499.2301 | -164477.2574 |            -164520.0786            | -164472.4899 |
|      |  0   | $2^3$ | 2482.9027 |  2481.7033   | 2484.0834 | -164391.6151 |            -164392.2246            | -164390.9963 |
|      |  1   | $3^3$ | 2482.8799 |  2482.7731   | 2482.9959 | -164391.6031 |            -164391.6623            | -164391.5473 |
|      |  2   | $4^3$ | 2482.8819 |  2482.8759   | 2482.8865 | -164391.6041 |            -164391.6067            | -164391.6012 |
|      | a.v. |       | 2482.8818 |              |           | -164391.6041 |                                    |              |
| 500  |  -1  | $1^3$ | 2177.3562 |  2177.2415   | 2179.8737 | -144163.9859 |            -144178.2736            | -144162.2361 |
|      |  0   | $2^3$ | 2176.2392 |  2176.2253   | 2176.2404 | -144088.7939 |            -144088.7971            | -144088.7538 |
|      |  1   | $3^3$ | 2176.2391 |  2176.2391   | 2176.2392 | -144088.7937 |            -144088.7938            | -144088.7936 |
|      |  2   | $4^3$ | 2176.2391 |  2176.2391   | 2176.2391 | -144088.7937 |            -144088.7937            | -144088.7937 |
|      | a.v. |       | 2176.2391 |              |           | -144088.7937 |                                    |              |
| 1500 |  -1  | $1^3$ | 1498.7997 |  1498.7569   | 1499.0522 | -99236.0403  |            -99238.4122             | -99235.3167  |
|      |  0   | $2^3$ | 1498.0239 |  1498.0236   | 1498.0240 | -99184.1673  |            -99184.1678             | -99184.1647  |
|      |  1   | $3^3$ | 1498.0240 |  1498.0240   | 1498.0240 | -99184.1672  |            -99184.1672             | -99184.1672  |
|      |  2   | $4^3$ | 1498.0240 |  1498.0240   | 1498.0240 | -99184.1672  |            -99184.1672             | -99184.1672  |
|      | a.v. |       | 1498.0240 |              |           | -99184.1672  |                                    |              |
| 3000 |  -1  | $1^3$ | 717.8387  |   717.8315   | 717.8533  | -47528.1777  |            -47528.3488             | -47528.0725  |
|      |  0   | $2^3$ | 717.4641  |   717.4640   | 717.4641  | -47503.2981  |            -47503.2982             | -47503.2980  |
|      |  1   | $3^3$ | 717.4641  |   717.4641   | 717.4641  | -47503.2981  |            -47503.2981             | -47503.2981  |
|      |  2   | $4^3$ | 717.4641  |   717.4641   | 717.4641  | -47503.2981  |            -47503.2981             | -47503.2981  |
|      | a.v. |       | 717.4641  |              |           | -47503.2981  |                                    |              |

```



The accuracy of the calculations increases with $I$ but so does the time spent
in the postprocessor: for $I\in\{-1,0,1,2\}$ this time was about 18, 132, 440
and 1040&nbsp;s respectively. This is easily explained when one realizes that
increasing $I$ from 0 to 1 means that the number of GLQ points per cell goes
from $2^3=8$ to $3^3=27$, i.e. a $3.375$ increase in operations for the same
number of cells and measurement points. The time spent in the postprocessor
increases by a similar factor ($440/132\sim 3.4$). The accuracy obtained with
lower $I$ values increases with increasing anomaly depth (or increasing
distance from the observation point). Note that this benchmark has uniform
density so that the projection of the density from the nodes onto the
quadrature points is exact and it does not introduce any smoothing of data
which might occur in the presence of density discontinuities (e.g. air-water,
water-crust, moho, etc.) inside a cell. $I=1$ seems to be the best compromise
between accuracy (gravity acceleration errors are less than 0.01&nbsp;mGal for
all shells) and compute time for this experiment.

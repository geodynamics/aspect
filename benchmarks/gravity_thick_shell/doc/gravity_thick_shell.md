(sec:benchmarks:thick-shell-gravity)=
# Thick shell gravity benchmark

*This section was contributed by Cedric Thieulot.*

This benchmark tests the accuracy of the gravity field and gravitational
potential computed by the gravity postprocessor inside and outside an
Earth-sized planet (without its core) of constant density. The domain is a
spherical shell with inner radius $R_\text{inner}=3840~\text{ km}$ and
outer radius $R_\text{outer}=6371~\text{ km}$. The density is constant
in the domain and set to $\rho_0=3300~\text{ kg}/\text{m}^3$.

First, let us calculate the exact profile which we expect the benchmark to
reproduce. The gravitational potential $U$ of a spherically symmetric object
satisfies the Poisson equation $\Delta U = 4\pi G \rho(\mathbf r)$. For a
constant density shell, this equation can be solved analytically for the
gravitational acceleration and potential inside and outside the planet. Inside
($r<R_\text{inner}$) and outside ($r>R_\text{outer}$) the
spherical shell (i.e. where $\rho=0$) the Poisson equation simplifies to the
Laplace equation $\Delta U=0$:
```{math}
\frac{1}{r^2} \frac{\partial }{\partial r} \left(r^2 \frac{\partial U}{\partial r} \right) = 0.
```
The solution to this expression is:
```{math}
:label: eq:app1
g=\frac{\partial U}{\partial r} = \frac{C}{r^2},
```
where $C$
is a constant of integration. In order to avoid an infinite gravity field at
$r=0$ (where the density is also zero in this particular setup of a shell), we
need to impose $C=0$, i.e. the gravity is zero for
$r\leq R_\text{inner}$. Another way of arriving at the same conclusion
is to realize that $g$ is zero at the center of the body because the material
around it exerts an equal force in every direction. Inside the shell,
$\rho=\rho_0$, yielding
```{math}
g=\frac{\partial U}{\partial r} = \frac{4 \pi}{3} G \rho_0 r + \frac{A}{r^2},
```
where $A$ is another integration constant. At the inner boundary,
$r=R_\text{inner}$ and $g=0$, allowing $A$ to be computed.
Substituting in the value of $A$,
```{math}
:label: eqgin
g=\frac{\partial U}{\partial r} = \frac{4 \pi}{3} G \rho_0
\left(r - \frac{R_\text{inner}^3}{r^2} \right).
```
When
$r\geq R_\text{outer}$, the gravitational potential is given by {math:numref}`eq:app1`.
Requiring the gravity field to be continuous at
$r=R_\text{outer}$:
```{math}
:label: eq:gout
g(r) = \frac{G M}{r^2},
```
where $M=\frac{4 \pi}{3} \rho_0(R_\text{outer}^3-R_\text{inner}^3)$
is the mass contained in the shell. For $r\ge R_\text{outer}$, the
potential is obtained by integrating {math:numref}`eq:gout`):
```{math}
U(r)=-\frac{GM}{r} +D,
```
where $D$ is an integration constant which has to
be zero since we require the potential to vanish for $r\rightarrow \infty$.
The potential within the shell,
$R_\text{inner}\leq r \leq R_\text{outer}$, is found by
integrating {math:numref}`eqgin`:
```{math}
U(r)= \frac{4 \pi}{3} G \rho_0 \left(\frac{r^2}{2} + \frac{R_\text{inner}^3}{r} \right)  + E,
```
where $E$ is a constant. Continuity of the potential at
$r=R_\text{outer}$ requires that
$E=-2\pi\rho_0 G R_\text{outer}^2$. Gravitational acceleration is zero
for $r\leq R_\text{inner}$, so the potential there is constant and a
continuity requirement yields
```{math}
U(r)=2\pi G \rho_0 (R_\text{inner}^2-R_\text{outer}^2).
```

The gravity postprocessor in can be used to calculate the radial components of
gravity ($g_r$ and $U$) at arbitrary points using the sampling scheme
'*list of points*.' For this benchmark we calculate points along a
line from the center of the planet to a distant point, $r=0$ to
$r=10,000~\text{ km}$ ({numref}`fig:grav-thick-shell1`). Arbitrarily, the latitude and longitude
are both set to $13\text{ \degree}$ so as to avoid potential measurement
artifacts due to symmetry. The list of radii is defined as follows:

```{literalinclude} thick_shell.prm
```

The resulting measurements obtained for a mesh composed of 12 caps of $32^3$
cells (i.e., 393,216 total mesh cells) are shown in {numref}`fig:grav-thick-shell2` and are
in good agreement with the analytical profiles.


```{figure-md} fig:grav-thick-shell1
<img src="gravity_g.*" width="100%" />

Gravity benchmark for a thick shell. Gravitational potential computed on a line from the center of a constant density shell to a radius of 10,000 km. The gray area indicates the region $R_{inner} \leq r \leq R_{outer}$ inside the shell, where the density is not zero.
```

```{figure-md} fig:grav-thick-shell2
<img src="gravity_U.*" width="100%" />

Gravity benchmark for a thick shell. Gravitational acceleration computed on a line from the center of a constant density shell to a radius of 10,000 km. The gray area indicates the region $R_{inner} \leq r \leq R_{outer}$ inside the shell, where the density is not zero.
```

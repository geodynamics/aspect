# The Rayleigh-Taylor instability

*This section was contributed by Cedric Thieulot.*

This benchmark is carried out in {cite:t}`deubelbeiss:kaus:2008,gerya:2010,thieulot:2011`) and is based on the analytical solution by Ramberg {cite}`ramberg:1968`, which consists of a gravitationally unstable two-layer system. Free
slip are imposed on the sides while no-slip boundary conditions are imposed on
the top and the bottom of the box. Fluid 1 $(\rho_1,\eta_1)$ of thickness
$h_1$ overlays fluid 2 $(\rho_2,\eta_2)$ of thickness $h_2$ (with
$h_1+h_2=L_y$). An initial sinusoidal disturbance of the interface between
these layers is introduced and is characterized by an amplitude $\Delta$ and a
wavelength $\lambda=L_x/2$ as shown in {numref}`fig:RTi_setup`.

```{figure-md} fig:RTi_setup
<img src="setup.*" style="width:40.0%" />

 Setup of the Rayleigh-Taylor instability benchmark (taken from {cite:t}`thieulot:2011`)
```

Under this condition, the velocity of the diapiric growth $v_y$ is given by
the relation
```{math}
\frac{v_y}{\Delta} = - K \frac{\rho_1-\rho_2}{2 \eta_2} h_2 g
\qquad
\qquad
\text{with}
\qquad
\qquad
K=\frac{-d_{12}}{c_{11}j_{22}-d_{12}i_{21}}
```
where $K$ is the dimensionless
growth factor and
```{math}

c_{11} &= \frac{\eta_1 2 \phi_1^2}{\eta_2(\cosh 2\phi_1 - 1 - 2\phi_1^2)} - \frac{2\phi_2^2}{\cosh 2\phi_2 - 1 - 2 \phi_2^2}\\
d_{12} &= \frac{\eta_1(\sinh 2\phi_1 -2\phi_1)}{\eta_2(\cosh 2\phi_1 -1 -2\phi_1^2)} + \frac{\sinh 2\phi_2 - 2\phi_2}{\cosh 2\phi_2 -1 -2\phi_2^2} \\
i_{21} &= \frac{\eta_1\phi_2 (\sinh 2 \phi_1 + 2 \phi_1)}{\eta_2(\cosh 2\phi_1 -1 -2\phi_1^2)}
+ \frac{\phi_2 (\sinh 2\phi_2 + 2\phi_2)}{\cosh 2\phi_2 -1 -2\phi_2^2} \\
j_{22} &= \frac{\eta_1 2 \phi_1^2 \phi_2}{\eta_2(\cosh 2\phi_1 -1-2\phi_1^2)} - \frac{2\phi_2^3}{ \cosh 2\phi_2 -1 -2\phi_2^2}\\
\phi_1&=\frac{2\pi h_1}{\lambda} \\
\phi_2&=\frac{2\pi h_2}{\lambda}
```
We set
$L_x=L_y=512 \text{ km}$, $h_1=h_2=256 \text{ km}$,
$|\boldsymbol{g}|=10\text{m}/\text{s}^2$, $\Delta=3\text{ km}$,
$\rho_1=3300\text{ kg}/\text{m}^3$, $\rho_2=3000\text{kg}/\text{m}^3$,
$\eta_1=1e21\text{ Pa s}$. $\eta_2$ is varied between $10^{20}$ and $10^{23}$
and 3 values of $\lambda$ (64, 128, and 256km) are used. Adaptive mesh
refinement based on density is used to capture the interface between the two
fluids, as shown in {numref}`fig:RTi_grids_a` and {numref}`fig:RTi_grids_b`. This translates as follows in the input
file:

    subsection Mesh refinement
      set Initial global refinement = 4
      set Initial adaptive refinement = 6
      set Strategy = density
      set Refinement fraction = 0.6
    end


**[Description of benchmark files](../README.md)**

```{figure-md} fig:RTi_grids_a
<img src="grid.*" style="width:44.0%" />

Grid with initial global refinement 4 and adaptive refinement 6.
```

```{figure-md} fig:RTi_grids_b
<img src="grid2.*" style="width:52.0%" />

Density field with detail of the mesh.
```

The maximum vertical velocity is plotted against $\phi_1$ in {numref}`fig:RTi_vels`
and is found to match analytical results.

```{figure-md} fig:RTi_vels
<img src="plot.*" style="width:75.0%" />

 Maximum velocity for three values of the \phi_1 parameter.
```

:::{toctree}
../README.md
:::

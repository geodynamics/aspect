# The Rayleigh-Taylor instability

*This section was contributed by Cedric Thieulot.*

This benchmark is carried out in (Deubelbeiss and Kaus 2008; Gerya 2010;
Thieulot 2011) and is based on the analytical solution by Ramberg(Ramberg
1968), which consists of a gravitationally unstable two-layer system. Free
slip are imposed on the sides while no-slip boundary conditions are imposed on
the top and the bottom of the box. Fluid 1 $(\rho_1,\eta_1)$ of thickness
$h_1$ overlays fluid 2 $(\rho_2,\eta_2)$ of thickness $h_2$ (with
$h_1+h_2=L_y$). An initial sinusoidal disturbance of the interface between
these layers is introduced and is characterised by an amplitude $\Delta$ and a
wavelength $\lambda=L_x/2$ as shown in Figure&nbsp;[1].

```{figure-md} fig:RTi_setup
<img src="setup.*" style="width:40.0%" />

 Setup of the Rayleigh-Taylor instability benchmark (taken from (Thieulot 2011))
```

Under this condition, the velocity of the diapiric growth $v_y$ is given by
the relation $$\frac{v_y}{\Delta} = - K \frac{\rho_1-\rho_2}{2 \eta_2} h_2 g
\qquad
\qquad
\text{with}
\qquad
\qquad
K=\frac{-d_{12}}{c_{11}j_{22}-d_{12}i_{21}}$$ where $K$ is the dimensionless
growth factor and $$\begin{aligned}
c_{11} &=& \frac{\eta_1 2 \phi_1^2}{\eta_2(\cosh 2\phi_1 - 1 - 2\phi_1^2)} - \frac{2\phi_2^2}{\cosh 2\phi_2 - 1 - 2 \phi_2^2}\\
d_{12} &=& \frac{\eta_1(\sinh 2\phi_1 -2\phi_1)}{\eta_2(\cosh 2\phi_1 -1 -2\phi_1^2)} + \frac{\sinh 2\phi_2 - 2\phi_2}{\cosh 2\phi_2 -1 -2\phi_2^2} \\
i_{21} &=& \frac{\eta_1\phi_2 (\sinh 2 \phi_1 + 2 \phi_1)}{\eta_2(\cosh 2\phi_1 -1 -2\phi_1^2)}
+ \frac{\phi_2 (\sinh 2\phi_2 + 2\phi_2)}{\cosh 2\phi_2 -1 -2\phi_2^2} \\
j_{22} &=& \frac{\eta_1 2 \phi_1^2 \phi_2}{\eta_2(\cosh 2\phi_1 -1-2\phi_1^2)} - \frac{2\phi_2^3}{ \cosh 2\phi_2 -1 -2\phi_2^2}\\
\phi_1&=&\frac{2\pi h_1}{\lambda} \\
\phi_2&=&\frac{2\pi h_2}{\lambda}\end{aligned}$$ We set
$L_x=L_y=\SI{512}{km}$, $h_1=h_2=\SI{256}{km}$,
$|\boldsymbol{g}|=\SI{10}{m/s^2}$, $\Delta=\SI{3}{km}$,
$\rho_1=\SI{3300}{kg/m^3}$, $\rho_2=\SI{3000}{kg/m^3}$,
$\eta_1=\SI{1e21}{Pa.s}$. $\eta_2$ is varied between $10^{20}$ and $10^{23}$
and 3 values of $\lambda$ (64, 128, and 256km) are used. Adaptive mesh
refinement based on density is used to capture the interface between the two
fluids, as shown in Figure&nbsp;[3]. This translates as follows in the input
file:

    subsection Mesh refinement
      set Initial global refinement = 4
      set Initial adaptive refinement = 6
      set Strategy = density
      set Refinement fraction = 0.6
    end


```{figure-md} fig:RTi_grids
<img src="grid.*" style="width:44.0%" />

 Left: grid with initial global refinement 4 and adaptive refinement 6; Right: density field with detail of the mesh.
```

```{figure-md} fig:RTi_grids
<img src="grid2.*" style="width:52.0%" />

 Left: grid with initial global refinement 4 and adaptive refinement 6; Right: density field with detail of the mesh.
```

The maximum vertical velocity is plotted against $\phi_1$ in Figure&nbsp;[4]
and is found to match analytical results.

```{figure-md} fig:RTi_vels
<img src="plot.*" style="width:75.0%" />

 Maximum velocity for three values of the \phi_1 parameter.
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Deu08" class="csl-entry">

Deubelbeiss, Y., and B. J. P. Kaus. 2008. "Comparison of Eulerian and
Lagrangian Numerical Techniques for the Stokes Equations in the Presence of
Strongly Varying Viscosity." *Physics of the Earth and Planetary
Interiors* 171: 92--111.

</div>

<div id="ref-Ger10" class="csl-entry">

Gerya, T. 2010. *Introduction to Numerical Geodynamic Modelling*. Cambridge
University Press.

</div>

<div id="ref-ramb68" class="csl-entry">

Ramberg, Hans. 1968. "Instability of Layered Systems in the Field of
Gravity." *Phys. Earth Planet. Interiors* 1: 427--47.

</div>

<div id="ref-thie11" class="csl-entry">

Thieulot, C. 2011. "<span class="nocase">FANTOM: two- and
three-dimensional numerical modelling of creeping flows for the solution of
geological problems</span>."
*Phys.&nbsp;Earth.&nbsp;Planet.&nbsp;Inter.* 188: 47--68.

</div>

</div>

  [1]: #fig:RTi_setup
  [3]: #fig:RTi_grids
  [4]: #fig:RTi_vels

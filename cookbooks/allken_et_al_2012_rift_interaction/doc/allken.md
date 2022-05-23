(sec:cookbooks:allken)=
# rift interaction in brittle-ductile coupled systems

*This section was contributed by Cedric Thieulot*

The setup for this experiment is based of {cite:t}`alht11,alht12,alhf13`.
The Cartesian domain is $L_x\times L_y\times L_z=\SI{210}{\km}\times\SI{210}{\km}\times\SI{30}{\km}$
in size (see {numref}`fig:allken_setup`).
It is a simple 2-layer crustal model consisting of $h_{uc}=\SI{15}{\km}$ of upper crust (characterised by a Mohr-Coulomb rheology with cohesion
$C^0=\SI{20}{\mega\pascal}$ and angle of friction $\phi^0=\SI{16}{\degree}$) and $h_{lc}=\SI{15}{\km}$ of lower crust
(characterised by a Newtonian viscous rheology with a viscosity $\eta_{lc}=\SI{1e20}{\pascal\second}$).
The flow is assumed to be incompressible and isothermal.
Both material have the same density $\rho=\SI{2800}{\kg\per\cubic\meter}$ and the gravity is vertical with magnitude \SI{9.81}{\meter\per\square\second}.

```{figure-md} fig:allken-setup
<img src="allken_setup.png" style="width:80.0%" />

(a) Model setup showing box of dimension. (b) Rheological profile implemented: Mohr-Coulomb plasticity in the upper crust and a fixed viscosity. (c) Frictional plastic strain weakening behavior of the upper crust. Figure taken from {cite}`alht12`.
```

The boundary conditions are free slip on the bottom and on the $y=0$ and $y=L_y$ faces,
free surface at the top, and a velocity $(\pm v_{ext},0,0)$ is prescribed on the
$x=0$ and $x=L_x$ sides, with $v_{ext}=\SI{0.5}{\cm\per year}$, so as to generate extension in the domain.

Strain weakening is an essential part of this experiment: both cohesion and angle of friction see their values linearly decrease between strains $\varepsilon_1=0.25$ and $\varepsilon_2=1.25$ to reach a maximum weakening factor of $R=C^0/C^{sw}=\phi^0/\phi^{sw}=4$ (i.e. $C^{sw}=\SI{5}{\mega\pascal}$ and $\phi^{sw}=\SI{4}{\degree}$).
Two weak zones are placed at the base of the upper crust and are characterised by a prescribed strain value of $\varepsilon_2$. As shown in {numref}`fig:allken_setup` they are offset by a distance $\Delta$ which is a multiple (between 2 and 6) of $h_{uc}$.
The cross section of the two weak zones is $\SI{4}{\km}\times\SI{2}{\km}$ and they are \SI{50}{\km} long.

```{figure-md} fig:allken-results
<img src="vx.png" style="width:45%" />
<img src="vz.png" style="width:45%" />
<img src="visc.png" style="width:45%" />
<img src="sr.png" style="width:45%" />

Horizontal and vertical velocity components, effective viscosity and effective strain rate after 1 Myr.
```

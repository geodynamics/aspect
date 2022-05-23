(sec:cookbooks:allken)=
# Rift interaction in brittle-ductile coupled systems

*This section was contributed by Cedric Thieulot*

The setup for this experiment is based of {cite:t}`alht11,alht12,alhf13`.
The Cartesian domain is $L_x\times L_y\times L_z=210\text{ km}\times 210\text{ km}\times 30\text{ km}$
in size (see {numref}`fig:allken_setup`).
It is a simple 2-layer crustal model consisting of $h_{uc}=15\text{ km}$ of upper crust (characterised by a Mohr-Coulomb rheology with cohesion
$C^0=20\text{ MPa}$ and angle of friction $\phi^0=16^o$) and $h_{lc}=15\text{ km}$ of lower crust
(characterised by a Newtonian viscous rheology with a viscosity $\eta_{lc}=10^{20}\text{ Pa.s}$).
The flow is assumed to be incompressible and isothermal.
Both materials have the same density $\rho=2800\text{ kg}/\text{m}^3$ and the gravity is vertical with magnitude $9.81\text{ m}/\text{s}^2$.

```{figure-md} fig:allken-setup
<img src="allken_setup.png" style="width:80.0%" />

(a) Model setup showing box of dimension. (b) Rheological profile implemented: Mohr-Coulomb plasticity in the upper crust and a fixed viscosity. (c) Frictional plastic strain weakening behavior of the upper crust. Figure taken from {cite}`alht12`.
```

The boundary conditions are free slip on the bottom and on the $y=0$ and $y=L_y$ faces,
free surface at the top, and a velocity $(\pm v_{ext},0,0)$ is prescribed on the
$x=0$ and $x=L_x$ sides, with $v_{ext}=0.5\text{ cm}/\text{yr}$, so as to generate extension in the domain.

Strain weakening is an essential part of this experiment: both cohesion and angle of friction see their values linearly decrease between strains $\varepsilon_1=0.25$ and $\varepsilon_2=1.25$ to reach a maximum weakening factor of $R=C^0/C^{sw}=\phi^0/\phi^{sw}=4$ (i.e. $C^{sw}=5\text{ MPa}$ and $\phi^{sw}=4^o$).
Two weak zones are placed at the base of the upper crust and are characterised by a prescribed strain value of $\varepsilon_2$. As shown in {numref}`fig:allken_setup` they are offset by a distance $\Delta$ which is a multiple (between 2 and 6) of $h_{uc}$.
The cross section of the two weak zones is $4\text{ km}\times 2\text{ km}$ and they are $50\text{ km}$ long.
Results after $1\text{ Myr}$ of extension are shown in {numref}`fig:allken-result1`, {numref}`fig:allken-result2`, {numref}`fig:allken-result3` and {numref}`fig:allken-result4`.  

```{figure-md} fig:allken-result1
<img src="vx.png" style="width:60%" />

Velocity ($x$-component) field after 1 Myr.
```

```{figure-md} fig:allken-result2
<img src="vz.png" style="width:60%" />

Velocity ($z$-component) field after 1 Myr.
```

```{figure-md} fig:allken-result3
<img src="visc.png" style="width:60%" />

Effective viscosity field after 1 Myr.
```

```{figure-md} fig:allken-result4
<img src="sr.png" style="width:60%" />

Effective strain rate field after 1 Myr.
```

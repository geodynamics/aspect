# The slab detachment benchmark

*This section was contributed by Cedric Thieulot and Anne Glerum.*

Slab detachment (slab break-off) may occur in the final stages of subduction
as a consequence of the combination of a buoyant crust and strong slab pull.
It is often invoked to explain geophysical and geological observations such as
tomographic images of slab remnants and exhumed ultra-high-pressure rocks
{cite}`wortel:spakman:2000,vanhunen:allen:2011,garzanti:etal:2018`.

This benchmark is based on the setup by S. Schmalholtz {cite}`schmalholz:2011`,
which was subsequently run with by A. Glerum {cite}`glerum:etal:2018`. The
computational domain is a $1000 \text{ km}\times 660 \text{ km}$ box. No-slip
boundary conditions are imposed on the sides of the system, while free-slip
boundary conditions are imposed at the top and bottom.

```{figure-md} fig:slab_detachment_setup
<img src="drawing.png" alt="Initial geometry" width="100%"/>

 Slab detachment benchmark: Initial geometry
```

Two materials are present in the domain: the lithosphere and the mantle as
shown in {numref}`fig:slab_detachment_setup`. The gravity acceleration is Earth-like with
$g=9.81 \text{ m}\text{ s}^2$. The overriding plate is $80\text{ km}$ thick and is
placed at the top of the domain. The already subducted lithosphere extends
vertically into the mantle for $250 \text{ km}$. This slab has a density
$\rho_s=3300\text{ kg/m}^3$ and is characterized by a power-law flow law so
that its effective viscosity depends on the square root of the second
invariant of the strain rate $\dot\varepsilon$:
```{math}
\eta_{eff} = \eta_0 \, \dot\varepsilon^{1/n-1}
```
with $n=4$ and
$\eta_0=4.75\times 10^{11}\text{ Pa . s}$. The mantle occupies the rest of the domain and
has a constant viscosity $\eta_m=1\times 10^{21}\text{ Pa . s}$ and a density
$\rho_m=3150\text{ kg/m}^3$. Viscosity is capped between $1\times10^{21}\text{ Pa . s}$
and $1\times 10^{25} \text{ Pa . s}$. {numref}`fig:slab_detachment_evolution`
shows the various fields and their evolution through time. As shown in
{cite:t}`schmalholz:2011,glerum:etal:2018` the hanging slab necks, helped by the
localizing effect of the nonlinear rheology. Model results were shown to
compare favorably to the results of {cite:t}`schmalholz:2011` in
{cite:t}`glerum:etal:2018,hillebrand:etal:2014` and the effect of viscosity and material averaging was
explored in {cite}`glerum:etal:2018`.

:::{figure-md} fig:slab_detachment_evolution
<img src="results.*" alt="Screenshot"  width="85%"/>

Slab detachment benchmark: a,b) velocity and strain rate fields at $t=0$. c,d,e) and f,g,h) time evolution of the viscosity and slab composition fields at $t=0, 6, 12\text{ Myr}$.
:::

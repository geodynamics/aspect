# The Bunge et al. mantle convection experiments

*This section was contributed by Cedric Thieulot and Bob Myhill.*

Early mantle modeling studies of the 1970s, 1980s and 1990s were often
concerned with simple set-ups (Cartesian geometries, incompressible fluid,
free slip boundary conditions) and investigated the influence of the Rayleigh
number, the heating mode or the temperature dependence of the viscosity on
temperature, pressure and/or strain rate {cite}`BBC89,bunge:etal:1997,busse:1975,busse:1979,busse:etal:1993,vandenberg:etal:1993,young:1974`. In this cookbook, we use
the 'simple' material model to reproduce the set-up in {cite:t}`bunge:etal:1996`, which reported that even modest increases in
mantle viscosity with depth could have a marked effect on the style of mantle
convection. The prm file corresponding to this cookbook can be found at
[cookbooks/bunge_et_al_mantle_convection/bunge_et_al.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/bunge_et_al_mantle_convection/bunge_et_al.prm).

Although the original article showcases results obtained in a 3D hollow
sphere, we here run the models in an annular domain of inner radius
$R_{inner} = 3480~\text{ km}$ and outer radius
$R_{outer} = 6370~\text{ km}$. The surface temperature is set to
$T{_{surf}}$ = 1060 K and the bottom temperature to
$T{_{cmb}} = 3450$ K. The gravity vector is radial and its
magnitude is $g = 10 \text{m s}^{-2}$.

There is a single incompressible fluid in the domain, characterized by
$\rho_0 = 4500\text{ kg m}^{-3}$, $\alpha = 2.5\cdot10^{-5}\text{ K}$, $k = 4\text{ W m}^{-1}\text{K}^{-1}$, $C_p = 1000 \text{ J kg}^{-1}\text{K}^{-1}$ and its internal heating rate is
$Q{_{int}} = 1\cdot10^{-12}\text{ W kg}^{-1}$. The
interface between the upper mantle (viscosity $\eta_{um}$) and
the lower mantle (viscosity $\eta_{lm}$) is fixed at 670 km
depth. As in the article we consider four time-independent radial viscosity
profiles:

-   Isoviscous mantle:
    $\eta_{um}=\eta_{lm}=1.7\cdot 10^{24}$ Pa&nbsp;s

-   Mantle with step change in viscosity:
    $\eta_{um}=5.8\cdot 10^{22}$ Pa&nbsp;s,
    $\eta_{lm}=30\eta_{um}$

-   Isoviscous mantle:
    $\eta_{um}=\eta_{lm}=5.8\cdot 10^{22}$ Pa&nbsp;s

-   Mantle with step change in viscosity:
    $\eta_{um}=7\cdot 10^{21}$ Pa&nbsp;s,
    $\eta_{lm}=30\eta_{um}$

Separate ascii files `visc_depth_X.txt` with `X={a,b,c,d}` contain each of
these viscosity profiles. The resulting temperature fields after 5 billion
years of convection are shown in {numref}`fig:bunge_et_al`. Similar to the results
obtained by {cite:t}`bunge:etal:1996`, models in which the lower
mantle is more viscous than the upper mantle are distinctly colder than their
isoviscous equivalents, with more clearly defined upwellings. You can find a
movie of how the temperature evolves over this time period at
<https://youtu.be/5SPCU1sFGGc>.

```{figure-md} fig:bunge_et_al
<img src="temps.png" style="width:90.0%" />

 Bunge et al. benchmark. From left to right: temperature field at time $t=5\cdot 10^9$ years obtained with viscosity profiles a, b, c and d.
```

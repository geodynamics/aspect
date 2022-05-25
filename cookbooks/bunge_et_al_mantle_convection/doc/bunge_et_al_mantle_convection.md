# The Bunge et al.&nbsp;mantle convection experiments

*This section was contributed by Cedric Thieulot and Bob Myhill.*

Early mantle modeling studies of the 1970s, 1980s and 1990s were often
concerned with simple set-ups (Cartesian geometries, incompressible fluid,
free slip boundary conditions) and investigated the influence of the Rayleigh
number, the heating mode or the temperature dependence of the viscosity on
temperature, pressure and/or strain rate (Young 1974; F. H. Busse 1975; FH
Busse 1979; Blankenbach et al. 1989; F. Busse et al. 1993; Berg, Keken, and
Yuen 1993; Bunge, Richards, and Baumgardner 1997). In this cookbook, we use
the &lsquo;simple&rsquo; material model to reproduce the set-up in (Bunge,
Richards, and Baumgardner 1996), which reported that even modest increases in
mantle viscosity with depth could have a marked effect on the style of mantle
convection. The prm file corresponding to this cookbook can be found at
[cookbooks/bunge_et_al_mantle_convection/bunge_et_al.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/bunge_et_al_mantle_convection/bunge_et_al.prm).

Although the original article showcases results obtained in a 3D hollow
sphere, we here run the models in an annular domain of inner radius
$R\textsubscript{inner} = 3480~\si{\km}$ and outer radius
$R\textsubscript{outer} = 6370~\si{\km}$. The surface temperature is set to
$T{\textsubscript{surf}}$ = 1060 K and the bottom temperature to
$T{\textsubscript{cmb}} = 3450$ K. The gravity vector is radial and its
magnitude is $g = 10$ m&nbsp;s<sup>&minus;2</sup>.

There is a single incompressible fluid in the domain, characterized by
$\rho_0 = 4500$ kg&nbsp;m<sup>&minus;3</sup>, $\alpha = 2.5\cdot10^{-5}$
&nbsp;K, $k = 4$ W&nbsp;m&nbsp;K<sup>&minus;&minus;1</sup>, $C_p = 1000$
J&nbsp;kg&nbsp;K<sup>&minus;&minus;1</sup> and its internal heating rate is
$Q{\textsubscript{int}} = 1\cdot10^{-12}$ W&nbsp;kg<sup>&minus;1</sup>. The
interface between the upper mantle (viscosity $\eta\textsubscript{um}$) and
the lower mantle (viscosity $\eta\textsubscript{lm}$) is fixed at 670 km
depth. As in the article we consider four time-independent radial viscosity
profiles:

-   Isoviscous mantle:
    $\eta\textsubscript{um}=\eta\textsubscript{lm}=1.7\cdot 10^{24}$ Pa&nbsp;s

-   Mantle with step change in viscosity:
    $\eta\textsubscript{um}=5.8\cdot 10^{22}$ Pa&nbsp;s,
    $\eta\textsubscript{lm}=30\eta\textsubscript{um}$

-   Isoviscous mantle:
    $\eta\textsubscript{um}=\eta\textsubscript{lm}=5.8\cdot 10^{22}$ Pa&nbsp;s

-   Mantle with step change in viscosity:
    $\eta\textsubscript{um}=7\cdot 10^{21}$ Pa&nbsp;s,
    $\eta\textsubscript{lm}=30\eta\textsubscript{um}$

Separate ascii files `visc_depth_X.txt` with `X={a,b,c,d}` contain each of
these viscosity profiles. The resulting temperature fields after 5 billion
years of convection are shown in Fig.&nbsp;[1][]. Similar to the results
obtained by (Bunge, Richards, and Baumgardner 1996), models in which the lower
mantle is more viscous than the upper mantle are distinctly colder than their
isoviscous equivalents, with more clearly defined upwellings. You can find a
movie of how the temperature evolves over this time period at
<https://youtu.be/5SPCU1sFGGc>.

<figure>
<img src="cookbooks/bunge_et_al_mantle_convection/doc/temps.png" id="fig:bunge_et_al" style="width:90.0%" alt="Bunge et al.&#xA0;benchmark. From left to right: temperature field at time t=5\cdot 10^9 years obtained with viscosity profiles a, b, c and d." /><figcaption aria-hidden="true"><em>Bunge et al.&#xA0;benchmark. From left to right: temperature field at time <span class="math inline"><em>t</em>&#x2004;=&#x2004;5&#x2005;&#x22C5;&#x2005;10<sup>9</sup></span> years obtained with viscosity profiles a, b, c and d.</em></figcaption>
</figure>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-vavy93" class="csl-entry">

Berg, Arie P van den, Peter E van Keken, and David A Yuen. 1993. &ldquo;The
Effects of a Composite Non-Newtonian and Newtonian Rheology on Mantle
Convection.&rdquo; *Geophy.&nbsp;J.&nbsp;Int.* 115 (1): 62&ndash;78.

</div>

<div id="ref-BBC89" class="csl-entry">

Blankenbach, B., F. Busse, U. Christensen, L. Cserepes, D. Gunkel, U. Hansen,
H. Harder, et al. 1989. &ldquo;A Benchmark Comparison for Mantle Convection
Codes.&rdquo; *Geophys. J. Int.* 98: 23&ndash;38.

</div>

<div id="ref-burb96" class="csl-entry">

Bunge, H.-P., M. A. Richards, and J. R. Baumgardner. 1996. &ldquo;<span
class="nocase">Effect of depth-dependent viscosity on the planform of mantle
convection</span>.&rdquo; *Nature* 379: 436&ndash;38.

</div>

<div id="ref-burb97" class="csl-entry">

&mdash;&mdash;&mdash;. 1997. &ldquo;<span class="nocase">A sensitivity study
of three-dimensional spherical mantle convection at $10^8$ Rayleigh number:
Effects of depth-dependent viscosity, heating mode, and endothermic phase
change</span>.&rdquo; *J.&nbsp;Geophys.&nbsp;Res.* 102 (B6): 11, 991&ndash;12,
007. <https://doi.org/10.1029/96JB03806>.

</div>

<div id="ref-buss75" class="csl-entry">

Busse, F. H. 1975. &ldquo;<span class="nocase">Patterns of convection in
spherical shells</span>.&rdquo; *J. Fluid Mech.* 72 (1): 67&ndash;85.

</div>

<div id="ref-BC93" class="csl-entry">

Busse, F., U. Christensen, R. Clever, L. Cserepes, C. Gable, E. Giannandrea,
L. Guillou, et al. 1993. &ldquo;3d Convection at Infinite Prandtl Numbers in
Cartesian Geometry &mdash; a Benchmark Comparison.&rdquo; *Geophys. Astrophys.
Fluid Dynamics* 75: 39&ndash;59.

</div>

<div id="ref-buss79" class="csl-entry">

Busse, FH. 1979. &ldquo;High Prandtl Number Convection.&rdquo;
*Phys.&nbsp;Earth.&nbsp;Planet.&nbsp;Inter.* 19 (2): 149&ndash;57.

</div>

<div id="ref-youn74" class="csl-entry">

Young, Richard E. 1974. &ldquo;Finite-Amplitude Thermal Convection in a
Spherical Shell.&rdquo; *J. Fluid Mech.* 63 (4): 695&ndash;721.

</div>

</div>

  [cookbooks/bunge_et_al_mantle_convection/bunge_et_al.prm]: cookbooks/bunge_et_al_mantle_convection/bunge_et_al.prm
  [1]: #fig:bunge_et_al

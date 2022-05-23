#### The slab detachment benchmark

*This section was contributed by Cedric Thieulot and Anne Glerum.*

Slab detachment (slab break-off) may occur in the final stages of subduction
as a consequence of the combination of a buoyant crust and strong slab pull.
It is often invoked to explain geophysical and geological observations such as
tomographic images of slab remnants and exhumed ultra-high-pressure rocks
(Wortel and Spakman 2000; van Hunen and Allen 2011; Garzanti, Radeff, and
Malus&agrave; 2018).

This benchmark is based on the setup by S. Schmalholtz (Schmalholz 2011),
which was subsequently run with by A. Glerum (Glerum et al. 2018). The
computational domain is a $1000 \si{km}\times 660 \si{km}$ box. No-slip
boundary conditions are imposed on the sides of the system, while free-slip
boundary conditions are imposed at the top and bottom.

<figure>
<img src="cookbooks/benchmarks/slab_detachment/doc/drawing.png" id="fig:slab_detachment_setup" alt="Slab detachment benchmark: Initial geometry [fig:slab_detachment_setup]" /><figcaption aria-hidden="true"><em>Slab detachment benchmark: Initial geometry <span id="fig:slab_detachment_setup" label="fig:slab_detachment_setup">[fig:slab_detachment_setup]</span></em></figcaption>
</figure>

Two materials are present in the domain: the lithosphere and the mantle as
shown in Figure&nbsp;[1]. The gravity acceleration is Earth-like with
$g=9.81 \si{m}\si{s}^2$. The overriding plate is $80\si{km}$ thick and is
placed at the top of the domain. The already subducted lithosphere extends
vertically into the mantle for $250 \si{km}$. This slab has a density
$\rho_s=3300\si{kg}/\si{m}^3$ and is characterized by a power-law flow law so
that its effective viscosity depends on the square root of the second
invariant of the strainrate $\dot\varepsilon$:
$$\eta_{eff} = \eta_0 \, \dot\varepsilon^{1/n-1}$$ with $n=4$ and
$\eta_0=\SI{4.75e11}{Pa . s}$. The mantle occupies the rest of the domain and
has a constant viscosity $\eta_m=\SI{1e21}{Pa . s}$ and a density
$\rho_m=\SI{3150}{kg/m^3}$. Viscosity is capped between $\SI{1e21}{Pa . s}$
and $\SI{1e25}{Pa . s}$. Figure&nbsp;[\[fig:slab_detachment_evolution\]][2]
shows the various fields and their evolution through time. As shown in
(Schmalholz 2011; Glerum et al. 2018) the hanging slab necks, helped by the
localizing effect of the nonlinear rheology. Model results were shown to
compare favorably to the results of (Schmalholz 2011) in (Glerum et al. 2018;
Hillebrand et al. 2014) and the effect of viscosity and material averaging was
explored in (Glerum et al. 2018).

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-garm18" class="csl-entry">

Garzanti, E., G. Radeff, and M. G. Malus&agrave;. 2018. &ldquo;Slab Breakoff:
A Critical Appraisal of a Geological Theory as Applied in Space and
Time.&rdquo; *Earth-Science Reviews* 177: 303&ndash;19.
<https://doi.org/10.1016/j.earscirev.2017.11.012>.

</div>

<div id="ref-gltf18" class="csl-entry">

Glerum, A., C. Thieulot, M. Fraters, C. Blom, and W. Spakman. 2018.
&ldquo;Nonlinear Viscoplasticity in <span class="smallcaps">ASPECT</span>:
Benchmarking and Applications to Subduction.&rdquo; *Solid Earth* 9 (2):
267&ndash;94. <https://doi.org/10.5194/se-9-267-2018>.

</div>

<div id="ref-hitg14" class="csl-entry">

Hillebrand, B., C. Thieulot, T. Geenen, A. P. van den Berg, and W. Spakman.
2014. &ldquo;Using the Level Set Method in Geodynamical Modeling of
Multi-Material Flows and Earth&rsquo;s Free Surface.&rdquo; *Solid Earth* 5
(2): 1087&ndash;98. <https://doi.org/10.5194/se-5-1087-2014>.

</div>

<div id="ref-schm11" class="csl-entry">

Schmalholz, S. M. 2011. &ldquo;<span class="nocase">A simple analytical
solution for slab detachment</span>.&rdquo;
*Earth&nbsp;Planet.&nbsp;Sci.&nbsp;Lett.* 304: 45&ndash;54.

</div>

<div id="ref-vaal11" class="csl-entry">

van Hunen, J., and M. B. Allen. 2011. &ldquo;<span class="nocase">Continental
collision and slab break-off: A comparison of 3-D numerical models with
observations</span>.&rdquo; *Earth&nbsp;Planet.&nbsp;Sci.&nbsp;Lett.* 302:
27&ndash;37.

</div>

<div id="ref-wosp00" class="csl-entry">

Wortel, M. J. R., and W. Spakman. 2000. &ldquo;<span class="nocase">Subduction
and slab detachment in the Mediterranean-Carpathian region</span>.&rdquo;
*Science* 290: 1910&ndash;17.

</div>

</div>

  [1]: #fig:slab_detachment_setup
  [2]: #fig:slab_detachment_evolution

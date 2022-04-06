### Setups for teaching

Because <span class="smallcaps">ASPECT</span> is freely available, has an
extensive documentation and can be applied to a variety of problems in
geophysics, it can be a useful tool for teaching geophysics in general or
geodynamic modeling in particular. In the following section, we will present a
number of cookbooks that can be used for this purpose. Many of them are
modifications of existing cookbooks, but have been changed to run faster to be
more suitable for running them in the classroom, or they include additional
ideas for what parameters can be changed to learn more about the physical
behaviour that controls the model results.

##### Introduction to Geophysics

*This section was contributed by Juliane Dannberg, based on the course
&ldquo;Introduction to Geophysics&rdquo; at University of Florida.*

The course is designed to teach general concepts of geophysics, and it
includes the following cookbooks:

1.  (using the files in [cookbooks/convection-box-particles/][])

2.  (using the files in [cookbooks/heat_flow/][])

3.  (using the files in [cookbooks/onset_of_convection/][])

4.  (using the files in [cookbooks/magnetic_stripes/][])

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-FBTGS19" class="csl-entry">

Fraters, M. R. T., W. Bangerth, C. Thieulot, A. C. Glerum, and W. Spakman.
2019. &ldquo;Efficient and Practical Newton Solvers for Nonlinear Stokes
Systems in Geodynamics Problems.&rdquo; *Geophysics Journal International* 218
(2): 873&ndash;94. <https://doi.org/10.1093/gji/ggz183>.

</div>

<div id="ref-heister_aspect_methods2" class="csl-entry">

Heister, Timo, Juliane Dannberg, Rene Gassm&ouml;ller, and Wolfgang Bangerth.
2017. &ldquo;High Accuracy Mantle Convection Simulation Through Modern
Numerical Methods. II: Realistic Models and Problems.&rdquo; *Geophysical
Journal International* 210 (2): 833&ndash;51.
<https://doi.org/10.1093/gji/ggx195>.

</div>

<div id="ref-KHB12" class="csl-entry">

Kronbichler, M., T. Heister, and W. Bangerth. 2012. &ldquo;High Accuracy
Mantle Convection Simulation Through Modern Numerical Methods.&rdquo;
*Geophysical Journal International* 191: 12&ndash;29.
<https://doi.org/10.1111/j.1365-246X.2012.05609.x>.

</div>

<div id="ref-T15" class="csl-entry">

Tosi, N., C. Stein, L. Noack, C. H&uuml;ttig, P. Maierova, H. Samual, D. R.
Davies, et al. 2015. &ldquo;A Community Benchmark for Viscoplastic Thermal
Convection in a 2-d Square Box.&rdquo; *Geochem.&nbsp;Geophys.&nbsp;Geosyst.*
16: 2175&ndash;96.

</div>

</div>

[1] You can also extend <span class="smallcaps">ASPECT</span> using plugins
&ndash; i.e., pieces of code you compile separately and either link into the
<span class="smallcaps">ASPECT</span> executable itself, or reference from the
input file. This is discussed in Section&nbsp;[\[sec:extending\]][3].

[2] Internally, the geometry models <span class="smallcaps">ASPECT</span> uses
label every part of the boundary with what is called a *boundary indicator*
&ndash; a number that identifies pieces of the boundary. If you know which
number each piece has, you can list these numbers on the right hand sides of
the assignments of boundary types above. For example, the left boundary of the
box has boundary indicator zero (see
Section&nbsp;[\[parameters:Geometry_20model\]][6]), and using this number
instead of the `left` would have been equally valid. However, numbers are far
more difficult to remember than names, and consequently every geometry model
provides string aliases such as &ldquo;`left`&rdquo; for each boundary
indicator describing parts of the boundary. These symbolic aliases are
specific to the geometry &ndash; for the box, they are &ldquo;`left`,&rdquo;
&ldquo;`right`,&rdquo; &ldquo;`bottom`,&rdquo; etc., whereas for a spherical
shell they are &ldquo;`inner`&rdquo; and &ldquo;`outer`&rdquo; &ndash; but are
described in the documentation of every geometry model, see
Section&nbsp;[\[parameters:Geometry_20model\]][6].

[3] Verification is the first half of the *verification and validation* (V&V)
procedure: *verification* intends to ensure that the mathematical model is
solved correctly, while *validation* intends to ensure that the mathematical
model is correct. Obviously, much of the aim of computational geodynamics is
to validate the models that we have.

:::{toctree}
convection-box-particles.md
heat-flow.md
onsert_of_convection.md
magnetic_stripes.md
:::

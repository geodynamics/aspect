(sec:benchmarks:solitary_wave)=
# The solitary wave benchmark

*This section was contributed by Juliane Dannberg and is based on a section in
(Dannberg and Heister 2016) by Juliane Dannberg and Timo Heister.*

One of the most widely used benchmarks for codes that model the migration of
melt through a compacting and dilating matrix is the propagation of solitary
waves (e.g. (Simpson and Spiegelman 2011; Keller, May, and Kaus 2013;
Schmeling 2000)). The benchmark is intended to test the accuracy of the
solution of the two-phase flow equations as described in Section
[\[sec:melt_transport\]][1] and makes use of the fact that there is an
analytical solution for the shape of solitary waves that travel through a
partially molten rock with a constant background porosity without changing
their shape and with a constant wave speed. Here, we follow the setup of the
benchmark as it is described in (Barcilon and Richter 1986), which considers
one-dimensional solitary waves.

The model features a perturbation of higher porosity with the amplitude
$A \phi_0$ in a uniform low-porosity ($\phi=\phi_0$) background. Due to its
lower density, melt migrates upwards, dilating the solid matrix at its front
and compacting it at its end.

Assuming constant shear and compaction viscosities and using a permeability
law of the form $$\begin{aligned}
k_\phi &= k_0 \phi^3, && \text{ implying a Darcy coefficient }
K_D(\phi) = \frac{k_0}{\eta_f} \phi^3 , \\
\intertext{and the non-dimensionalization }
x &= \delta x'
  && \text{ with the compaction length } \delta = \sqrt{K_D(\phi_0)(\xi + \frac{4}{3}\eta)} , \\
\phi &= \phi_0 \phi '
  && \text{ with the background porosity } \phi_0 , \\
(\mathbf u_s, \mathbf u_f) &= u_0 (\mathbf u_s, \mathbf u_f)'
  && \text{ with the separation flux } \phi_0 u_0 = K_D(\phi_0) \Delta\rho g ,\end{aligned}$$
the analytical solution for the shape of the solitary wave can be written in
implicit form as: $$\begin{aligned}
x(\phi) &= \pm (A + 0.5)
\left[ -2 \sqrt{A-\phi} + \frac{1}{\sqrt{A-1}}
\ln \frac{\sqrt{A-1} - \sqrt{A-\phi}}{\sqrt{A-1} + \sqrt{A-\phi}} \right]\end{aligned}$$
and the phase speed $c$, scaled back to physical units, is $c = u_0 (2A+1)$.
This is only valid in the limit of small porosity $\phi_0 \ll 1$.
Figure&nbsp;[1][] illustrates the model setup.

<div class="center">

```{figure-md}
<embed src="cookbooks/benchmarks/solitary_wave/doc/setup.pdf" id="fig:setup-solitary-wave" style="width:65.0%" />

Setup of the solitary wave benchmark. The domain is <span class="math inline">400</span> m high and includes a low porosity (<span class="math inline"><em>&#x3D5;</em>&#x2004;=&#x2004;0.001</span>) background with an initial perturbation (<span class="math inline"><em>&#x3D5;</em>&#x2004;=&#x2004;0.1</span>). The solid density is <span class="math inline">3300&#x2006;<em>k</em><em>g</em>/<em>m</em><sup>3</sup></span> and the melt density is <span class="math inline">2500&#x2006;<em>k</em><em>g</em>/<em>m</em><sup>3</sup></span>. We apply the negative phase speed of the solitary wave <span class="math inline"><strong>u</strong><sub><em>s</em></sub>&#x2004;=&#x2004;&#x2005;&#x2212;&#x2005;<em>c</em>&#x2006;<strong>e</strong><sub><em>z</em></sub></span> as velocity boundary condition, so that the wave will stay at its original position while the background is moving.</em></figcaption>
```

</div>

The parameter file and material model for this setup can be found in
[benchmarks/solitary_wave/solitary_wave.prm][] and
[benchmarks/solitary_wave/solitary_wave.cc][]. The most relevant sections are
shown in the following paragraph.

``` prmfile
```

The benchmark uses a custom model to generate the initial condition for the
porosity field as specified by the analytical solution, and its own material
model, which includes the additional material properties needed by models with
melt migration, such as the permeability, melt density and compaction
viscosity. The solitary wave postprocessor compares the porosity and pressure
in the model to the analytical solution, and computes the errors for the shape
of the porosity, shape of the compaction pressure and the phase speed. We
apply the negative phase speed of the solitary wave as a boundary condition
for the solid velocity. This changes the reference frame, so that the solitary
wave stays in the center of the domain, while the solid moves downwards. The
temperature evolution does not play a role in this benchmark, so all
temperature and heating-related parameters are disabled or set to zero.

And extensive discussion of the results and convergence behavior can be found
in (Dannberg and Heister 2016).

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-BR86" class="csl-entry">

Barcilon, V., and F. M. Richter. 1986. &ldquo;Nonlinear Waves in Compacting
Media.&rdquo; *Journal of Fluid Mechanics* 164: 429&ndash;48.

</div>

<div id="ref-dannberg_melt" class="csl-entry">

Dannberg, J., and T. Heister. 2016. &ldquo;Compressible Magma/Mantle Dynamics:
3d, Adaptive Simulations in ASPECT.&rdquo; *Geophysical Journal International*
207 (3): 1343&ndash;66. <https://doi.org/10.1093/gji/ggw329>.

</div>

<div id="ref-KMK2013" class="csl-entry">

Keller, Tobias, Dave A. May, and Boris J. P. Kaus. 2013. &ldquo;Numerical
Modelling of Magma Dynamics Coupled to Tectonic Deformation of Lithosphere and
Crust.&rdquo; *Geophysical Journal International* 195 (3): 1406&ndash;42.
<https://doi.org/10.1093/gji/ggt306>.

</div>

<div id="ref-Schm00" class="csl-entry">

Schmeling, Harro. 2000. &ldquo;Partial Melting and Melt Segregation in a
Convecting Mantle.&rdquo; In *Physics and Chemistry of Partially Molten
Rocks*, 141&ndash;78. Springer.

</div>

<div id="ref-SS11" class="csl-entry">

Simpson, Gideon, and Marc Spiegelman. 2011. &ldquo;Solitary Wave Benchmarks in
Magma Dynamics.&rdquo; *Journal of Scientific Computing* 49 (3): 268&ndash;90.

</div>

</div>

  [1]: #sec:melt_transport
  [1]: #fig:setup-solitary-wave
  [benchmarks/solitary_wave/solitary_wave.prm]: benchmarks/solitary_wave/solitary_wave.prm
  [benchmarks/solitary_wave/solitary_wave.cc]: benchmarks/solitary_wave/solitary_wave.cc

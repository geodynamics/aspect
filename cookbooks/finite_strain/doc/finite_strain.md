(sec:cookbooks:finite_strain)=
# Tracking finite strain

*This section was contributed by Juliane Dannberg and Rene Gassm&ouml;ller*

In many geophysical settings, material properties, and in particular the
rheology, do not only depend on the current temperature, pressure and strain
rate, but also on the history of the system. This can be incorporated in
models by tracking history variables through compositional fields. In this
cookbook, we will show how to do this by tracking the strain that idealized
little grains of finite size accumulate over time at every (Lagrangian) point
in the model.

Here, we use a material model plugin that defines the compositional fields as
the components of the deformation gradient tensor $\mathbf F_{ij}$, and
modifies the right-hand side of the corresponding advection equations to
accumulate strain over time. This is done by adjusting the
`out.reaction_terms` variable:

``` c++
```

Let us denote the accumulated deformation at time step $n$ as $\mathbf F^n$.
We can calculate its time derivative as the product of two tensors, namely the
current velocity gradient $\mathbf G_{ij} = \frac{\partial u_i}{\partial x_j}$
and the deformation gradient $\mathbf F^{n-1}$ accumulated up to the previous
time step, in other words
$\frac{\partial \mathbf F}{\partial t} = \mathbf G \mathbf F$, and
$\mathbf F^0 = \mathbf I$, with $\mathbf I$ being the identity tensor. While
we refer to other studies (McKenzie and Jackson 1983; Dahlen and Tromp 1998;
Becker et al. 2003) for a derivation of this relationship, we can give an
intuitive example for the necessity to apply the velocity gradient to the
already accumulated deformation, instead of simply integrating the velocity
gradient over time. Consider a simple one-dimensional &ldquo;grain&rdquo; of
length $1.0$, in which case the deformation tensor only has one component, the
compression in $x$-direction. If one embeds this grain into a convergent flow
field for a compressible medium where the dimensionless velocity gradient is
$-0.5$ (e.g. a velocity of zero at its left end at $x=0.0$, and a velocity of
$-0.5$ at its right end at $x=1.0$), simply integrating the velocity gradient
would suggest that the grain reaches a length of zero after two units of time,
and would then &ldquo;flip&rdquo; its orientation, which is clearly
non-physical. What happens instead can be seen by solving the equation of
motion for the right end of the grain $\frac{dx}{dt} = v = -0.5 x$. Solving
this equation for $x$ leads to $x(t) = e^{-0.5t}$. This is therefore also the
solution for $\mathbf F$ since $\mathbf F x$ transforms the initial position
of $x(t=0)=1.0$ into the deformed position of $x(t=1) = e^{-0.5}$, which is
the definition of $\mathbf F$.

In more general cases a visualization of $\mathbf F$ is not intuitive, because
it contains rotational components that represent a rigid body rotation without
deformation. Following (Becker et al. 2003) we can polar-decompose the tensor
into a positive-definite and symmetric left stretching tensor $\mathbf L$, and
an orthogonal rotation tensor $\mathbf Q$, as
$\mathbf F = \mathbf L \mathbf Q$, therefore
$\mathbf L^2 = \mathbf L \mathbf L^T = \mathbf F \mathbf F^T$. The left
stretching tensor $\mathbf L$ (or finite strain tensor) then describes the
deformation we are interested in, and its eigenvalues $\lambda_i$ and
eigenvectors $\mathbf e_i$ describe the length and orientation of the
half-axes of the finite strain ellipsoid. Moreover, we will represent the
amount of relative stretching at every point by the ratio
$\ln(\lambda_1/\lambda_2)$, called the *natural strain* (Ribe 1992).

The full plugin implementing the integration of $\mathbf F$ can be found in
[cookbooks/finite_strain/finite_strain.cc][] and can be compiled with
`cmake . && make` in the [cookbooks/finite_strain][] directory. It can be
loaded in a parameter file as an &ldquo;Additional shared library,&rdquo; and
selected as material model. As it is derived from the &ldquo;simple&rdquo;
material model, all input parameters for the material properties are read in
from the subsection `Simple model`.

``` prmfile
```

```{figure-md}
<embed src="cookbooks/finite_strain/doc/finite_strain.pdf" id="fig:finite_strain" style="width:75.0%" />

Accumulated finite strain in an example convection model, as described in Section <a href="#sec:finite-strain" data-reference-type="ref" data-reference="sec:finite-strain">0.0.1</a> at a time of 67.6&#xA0;Ma. Top panel: Temperature distribution. Bottom panel: Natural strain distribution. Additional black crosses are the scaled eigenvectors of the stretching tensor <span class="math inline"><strong>L</strong></span>, showing the direction of stretching and compression.
```

The plugin was tested against analytical solutions for the deformation
gradient tensor in simple and pure shear as described in
[benchmarks/finite_strain/pure_shear.prm][] and
[benchmarks/finite_strain/simple_shear.prm][].

We will demonstrate its use at the example of a 2D Cartesian convection model
(Figure&nbsp;[1][]): Heating from the bottom leads to the ascent of plumes
from the boundary layer (top panel), and the amount of stretching is visible
in the distribution of natural strain (color in lower panel). Additionally,
the black crosses show the direction of stretching and compression (the
eigenvectors of $\mathbf L$). Material moves to the sides at the top of the
plume head, so that it is shortened in vertical direction (short vertical
lines) and stretched in horizontal direction (long horizontal lines). The
sides of the plume head show the opposite effect. Shear occurs mostly at the
edges of the plume head, in the plume tail, and in the bottom boundary layer
(black areas in the natural strain distribution).

The example used here shows how history variables can be integrated up over
the model evolution. While we do not use these variables actively in the
computation (in our example, there is no influence of the accumulated strain
on the rheology or any other material property), it would be trivial to extend
this material model in a way that material properties depend on the integrated
strain: Because the values of the compositional fields are part of what the
material model gets as inputs, they can easily be used for computing material
model outputs such as the viscosity.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Becker2003" class="csl-entry">

Becker, Thorsten W., James B. Kellogg, G&ouml;ran Ekstr&ouml;m, and Richard J.
O&rsquo;Connell. 2003. &ldquo;<span class="nocase">Comparison of azimuthal
seismic anisotropy from surface waves and finite strain from global
mantle-circulation models</span>.&rdquo; *Geophysical Journal International*
155 (2): 696&ndash;714. <https://doi.org/10.1046/j.1365-246X.2003.02085.x>.

</div>

<div id="ref-dahlen1998theoretical" class="csl-entry">

Dahlen, FA, and Jeroen Tromp. 1998. *Theoretical Global Seismology*. Princeton
University Press.

</div>

<div id="ref-McKenzie1983" class="csl-entry">

McKenzie, Dan, and James Jackson. 1983. &ldquo;<span class="nocase">The
relationship between strain rates, crustal thickening, palaeomagnetism, finite
strain and fault movements within a deforming zone</span>.&rdquo; *Earth and
Planetary Science Letters* 65 (1): 182&ndash;202.
<https://doi.org/10.1016/0012-821X(83)90198-X>.

</div>

<div id="ref-Ribe1992" class="csl-entry">

Ribe, Neil M. 1992. &ldquo;<span class="nocase">On the relation between
seismic anisotropy and finite strain</span>.&rdquo; *Journal of Geophysical
Research* 97 (B6): 8737. <https://doi.org/10.1029/92JB00551>.

</div>

</div>

  [cookbooks/finite_strain/finite_strain.cc]: cookbooks/finite_strain/finite_strain.cc
  [cookbooks/finite_strain]: cookbooks/finite_strain
  [benchmarks/finite_strain/pure_shear.prm]: benchmarks/finite_strain/pure_shear.prm
  [benchmarks/finite_strain/simple_shear.prm]: benchmarks/finite_strain/simple_shear.prm
  [1]: #fig:finite_strain

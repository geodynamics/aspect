#### The &ldquo;Stokes&rsquo; law&rdquo; benchmark

*This section was contributed by Juliane Dannberg.*

Stokes&rsquo; law was derived by George Gabriel Stokes in 1851 and describes
the frictional force a sphere with a density different than the surrounding
fluid experiences in a laminar flowing viscous medium. A setup for testing
this law is a sphere with the radius $r$ falling in a highly viscous fluid
with lower density. Due to its higher density the sphere is accelerated by the
gravitational force. While the frictional force increases with the velocity of
the falling particle, the buoyancy force remains constant. Thus, after some
time the forces will be balanced and the settling velocity of the sphere $v_s$
will remain constant:

$$\begin{aligned}
  \label{eq:stokes-law}
  \underbrace{6 \pi \, \eta \, r \, v_s}_{\text{frictional force}} =
  \underbrace{4/3 \pi \, r^3 \, \Delta\rho \, g,}_{\text{buoyancy force}}\end{aligned}$$
where $\eta$ is the dynamic viscosity of the fluid, $\Delta\rho$ is the
density difference between sphere and fluid and $g$ the gravitational
acceleration. The resulting settling velocity is then given by
$$\begin{aligned}
  \label{eq:stokes-velo}
  v_s = \frac{2}{9} \frac{\Delta\rho \, r^2 \, g}{\eta}.\end{aligned}$$
Because we do not take into account inertia in our numerical computation, the
falling particle will reach the constant settling velocity right after the
first timestep.

For the setup of this benchmark, we chose the following parameters:
$$\begin{aligned}
  \label{eq:stokes-parameters}
  r &= 200 \, \si{km}\\
  \Delta\rho &= 100 \, \si{kg}/\si{m}^3\\
  \eta &= 10^{22} \, \si{Pa.s}\\
  g &= 9.81 \, \si{m}/\si{s}^2.\end{aligned}$$ With these values, the exact
value of sinking velocity is $v_s =
\num{8.72e-10} \, \si{m}/\si{s}$.

To run this benchmark, we need to set up an input file that describes the
situation. In principle, what we need to do is to describe a spherical object
with a density that is larger than the surrounding material. There are
multiple ways of doing this. For example, we could simply set the initial
temperature of the material in the sphere to a lower value, yielding a higher
density with any of the common material models. Or, we could use &rsquo;s
facilities to advect along what are called &ldquo;compositional fields&rdquo;
and make the density dependent on these fields.

We will go with the second approach and tell to advect a single compositional
field. The initial conditions for this field will be zero outside the sphere
and one inside. We then need to also tell the material model to increase the
density by $\Delta\rho=100 kg\, m^{-3}$ times the concentration of the
compositional field. This can be done, like everything else, from the input
file.

All of this setup is then described by the following input file. (You can find
the input file to run this cookbook example in
[cookbooks/stokes/stokes.prm]. For your first runs you will probably want to
reduce the number of mesh refinement steps to make things run more quickly.)

``` prmfile
```

Using this input file, let us try to evaluate the results of the current
computations for the settling velocity of the sphere. You can visualize the
output in different ways, one of it being ParaView and shown in
Fig.&nbsp;[2] (an alternative is to use Visit as described in
Section&nbsp;[\[sec:viz\]][1]; 3d images of this simulation using Visit are
shown in Fig.&nbsp;[5]). Here, ParaView has the advantage that you can
calculate the average velocity of the sphere using the following filters:

1.  Threshold (Scalars: C_1, Lower Threshold 0.5, Upper Threshold 1),

2.  Integrate Variables,

3.  Cell Data to Point Data,

4.  Calculator (use the formula sqrt(velocity_x^2+
    velocity_y^2+velocity_z^2)/Volume).

If you then look at the Calculator object in the Spreadsheet View, you can see
the average sinking velocity of the sphere in the column &ldquo;Result&rdquo;
and compare it to the theoretical value
$v_s = \num{8.72e-10} \, \si{m}/\si{s}$. In this case, the numerical result is
$\num{8.865e-10} \,
\si{m}/\si{s}$ when you add a few more refinement steps to actually resolve
the 3d flow field adequately. The &ldquo;velocity statistics&rdquo;
postprocessor we have selected above also provides us with a maximal velocity
that is on the same order of magnitude. The difference between the analytical
and the numerical values can be explained by different at least the following
three points: (i) In our case the sphere is viscous and not rigid as assumed
in Stokes&rsquo; initial model, leading to a velocity field that varies inside
the sphere rather than being constant. (ii) Stokes&rsquo; law is derived using
an infinite domain but we have a finite box instead. (iii) The mesh may not
yet fine enough to provide a fully converges solution. Nevertheless, the fact
that we get a result that is accurate to less than 2% is a good indication
that implements the equations correctly.

<div class="center">

<img src="cookbooks/benchmarks/stokes/doc/stokes-velocity.png" title="fig:" id="fig:stokes-falling-sphere-2d" style="width:55.0%" alt="Stokes benchmark. Both figures show only a 2D slice of the three-dimensional model. Left: The compositional field and overlaid to it some velocity vectors. The composition is 1 inside a sphere with the radius of 200 km and 0 outside of this sphere. As the velocity vectors show, the sphere sinks in the viscous medium. Right: The density distribution of the model. The compositional density contrast of 100 kg/\si{m}^3 leads to a higher density inside of the sphere." />
<img src="cookbooks/benchmarks/stokes/doc/stokes-density.png" title="fig:" id="fig:stokes-falling-sphere-2d" style="width:44.0%" alt="Stokes benchmark. Both figures show only a 2D slice of the three-dimensional model. Left: The compositional field and overlaid to it some velocity vectors. The composition is 1 inside a sphere with the radius of 200 km and 0 outside of this sphere. As the velocity vectors show, the sphere sinks in the viscous medium. Right: The density distribution of the model. The compositional density contrast of 100 kg/\si{m}^3 leads to a higher density inside of the sphere." />

</div>

<div class="center">

<img src="cookbooks/benchmarks/stokes/doc/composition.png" title="fig:" id="fig:stokes-falling-sphere-3d" style="width:30.0%" alt="Stokes benchmark. Three-dimensional views of the compositional field (left), the adaptively refined mesh (center) and the resulting velocity field (right)." />
<img src="cookbooks/benchmarks/stokes/doc/mesh.png" title="fig:" id="fig:stokes-falling-sphere-3d" style="width:30.0%" alt="Stokes benchmark. Three-dimensional views of the compositional field (left), the adaptively refined mesh (center) and the resulting velocity field (right)." />
<img src="cookbooks/benchmarks/stokes/doc/velocity.png" title="fig:" id="fig:stokes-falling-sphere-3d" style="width:30.0%" alt="Stokes benchmark. Three-dimensional views of the compositional field (left), the adaptively refined mesh (center) and the resulting velocity field (right)." />

</div>

  [cookbooks/stokes/stokes.prm]: cookbooks/stokes/stokes.prm
  [2]: #fig:stokes-falling-sphere-2d
  [1]: #sec:viz
  [5]: #fig:stokes-falling-sphere-3d

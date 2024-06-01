(sec:cookbooks:artificial-viscosity-smoothing)=
# Artificial viscosity smoothing

*This section was contributed by Ryan Grove.*

Standard finite element discretizations of advection-diffusion equations
introduce unphysical oscillations around steep gradients. Therefore,
stabilization must be added to the discrete formulation to obtain correct
solutions. In ASPECT, we use the Entropy Viscosity scheme developed by
{cite:t}`guermond:etal:2011`.
In this scheme, an artificial viscosity is calculated on
every cell and used to try to combat these oscillations that cause unwanted
overshoot and undershoot. More information about how does this is located at
<https://dealii.org/developer/doxygen/deal.II/step_31.html>.

Instead of just looking at an individual cell&rsquo;s artificial viscosity,
improvements in the minimizing of the oscillations can be made by smoothing.
Smoothing is the act of finding the maximum artificial viscosity taken over a
cell $T$ and the neighboring cells across the faces of $T$, i.e.,
```{math}
\bar{v_h}(T) = \max_{K \in N(T)} v_h(K)
```
where $N(T)$ is the set containing $T$ and the neighbors across the faces of $T$.

This feature can be turned on by setting the `Use artificial viscosity smoothing`
flag inside the `Stabilization` subsection inside the `Discretization`
subsection in your parameter file.

To show how this can be used in practice, let us consider the simple
convection in a quarter of a 2d annulus cookbook in Section
{ref}`sec:cookbooks:shell_simple_2d`, a radial compositional field was added to help
show the advantages of using the artificial viscosity smoothing feature.

By applying the following changes shown below to the parameters of the already
existing file [cookbooks/shell_simple_2d/shell_simple_2d.prm](https://github.com/geodynamics/aspect/blob/main/cookbooks/shell_simple_2d/shell_simple_2d.prm),

```{literalinclude} shell_simple_2d_smoothing.part.prm
```

it is possible to produce pictures of the simple convection in a quarter of a
2d annulus such as the ones visualized in Figure&nbsp;{numref}`fig:smoothing`.

```{figure-md} fig:smoothing
<img src="artificial_viscosity_smoothing.*" style="width:100.0%" />

Example of the output of two similar runs.  The run on the left has the artificial viscosity smoothing turned on and the run on the right does not, as described in Section {ref}`sec:cookbooks:artificial-viscosity-smoothing`.
```

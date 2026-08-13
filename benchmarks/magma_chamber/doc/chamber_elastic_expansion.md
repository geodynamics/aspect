```{tags}
category:benchmark
feature:2d
feature:cartesian
feature:nonlinear-solver
feature:compositional-fields
feature:community-benchmark
```

(sec:benchmarks:magma_chamber)=
# Pressurized magma chamber

*This section was contributed by Annie Thompson and Cedric Thieulot.*

This is a 2D model of a pressurized magma chamber using a compressible equation of state and viscoelastic rheology. We create a pressurized chamber by including a non-zero thermal expansivity and adding heat to the chamber, causing the chamber to dilate. This benchmark tests the implementation of effectively elastic material in ASPECT. The elastic material in ASPECT matches the stress fields from an elastic finite element model. 

This is based on a commonly used 3D benchmark in the volcanology community for a pressurized magma chamber in elastic half-space, shown in {numref}`fig:chamber`. Examples can be found [here](http://www.driversofvolcanodeformation.org/) and in {cite:t}`crozier:etal:2023` under exercise 1.3.

```{figure-md} fig:chamber
<img src="chamber.png" width="50%" />

View of chamber and output transect locations from {cite}`crozier:etal:2023`. In ASPECT, we implement this chamber setup in the center of a 50x25 km box. Subsurface output lines are used for {numref}`fig:stress` and {numref}`fig:displacement`.
```

This implementation simulates a half-space using a 50x25 km domain that is significantly larger than the 1 km spherical magma chamber. The model has free-slip boundary conditions on the left, right, and bottom boundaries, and a free surface on the top boundary. It is calibrated to have a constant 10 MPa pressure inside the chamber by adjusting the heating rate and thermal expansivity of the material. To create elastic behavior, we set viscosity to 10<sup>20</sup> Pa s in the viscoelastic material model to make viscous behavior negligible. The material has a Poisson's ratio of 0.25, following the analytical solution.

As this is a 2D benchmark, it is not expected to match the 3D analytical solution from {cite:t}`zhong:dabrowski:jamtveit:2019` (below). Instead, we compare our results against the output of a purely elastic finite element model, Fieldstone from {cite:t}`fieldstone`. Fieldstone has been tested in 3D and matches the original 3D analytical solution. Below are comparisons between this benchmark's implementation, Fieldstone, and the analytical solution of {cite:t}`zhong:dabrowski:jamtveit:2019` ({numref}`fig:6-panel`, {numref}`fig:stress`, {numref}`fig:displacement`).


```{figure-md} fig:6-panel
<img src="6-panel.png" width="100%" />

Comparison of stress, displacement, and pressure fields around the chamber. Each inset is a combined plot of Fieldstone (left half of chamber) and ASPECT (right half of chamber) to help visualize the differences between model outputs.
```

```{figure-md} fig:stress
<img src="stress_plot.png" width="100%" />

Plots of stress along the transects shown in {numref}`fig:chamber` for this ASPECT benchmark, the 2D Fieldstone implementation, and the 3D analytical solution. The interior transect is 2km above the chamber and the surface transect is 4km above.
```
As seen in {numref}`fig:stress`, the stress fields from ASPECT match those of the elastic FEM. This confirms that ASPECT correctly solves the quasi-static equilibrium equation

```{math}
\nabla \cdot \boldsymbol{\sigma} = 0, \qquad
\boldsymbol{\sigma}\cdot\mathbf{n}\big|_{\partial\Omega_{chamber}} = -\Delta P\,\mathbf{n}
```

for the thermal eigenstrain loading used to pressurize the chamber: the divergence of the stress tensor correctly balances the internal overpressure {cite:t}`zhong:dabrowski:jamtveit:2019`.

It is important to note that the displacement fields are qualitatively similar, shown in {numref}`fig:displacement`, but there is a mismatch in magnitude. The primary difference between the two models is that in the elastic FEM, linear elastic material stress and strain are related algebraically through Hooke's law ($\boldsymbol\sigma = \lambda\,\mathrm{tr}(\boldsymbol\varepsilon)\,\mathbf I + 2\mu\,\boldsymbol\varepsilon$). ASPECT solves for velocity rather than displacement directly, and the displacements shown here are obtained as $\mathbf u \approx \mathbf v\,\Delta t^{e}$ from the velocity at a single elastic timestep.

```{figure-md} fig:displacement
<img src="displacement_plot.png" width="80%" />

Plots of displacement along the transects shown in {numref}`fig:chamber` for this ASPECT benchmark, the 2D Fieldstone implementation, and the 3D analytical solution. The interior transect is 2km above the chamber and the surface transect is 4km above. Displacement from ASPECT was calculated by multiplying velocity at 1 elastic timestep by the elastic timestep.
```
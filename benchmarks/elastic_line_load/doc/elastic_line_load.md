```{tags}
category:benchmark
feature:2d
feature:cartesian
feature:analytical-solution
feature:mesh-deformation
feature:visco-elastic-plastic
```

(sec:benchmarks:elastic-line-load)=
# 2D Elastic Line Loading Problem

*This section was contributed by Daniel Douglas, Cedric Thieulot, Prajakta Mohite and John Naliboff*

This example uses ASPECT to reproduce the modeling setup of {cite:t}`dase96`, which prescribes
a line "load" on the top of an elastic half-space. The analytic solution for $\sigma_{xx}$, $\sigma_{xy}$
and $\sigma_{yy}$ comes from Albert Flamant and assumes a homogeneous, isotropic elastic half-space. To
approximate this in ASPECT, we use the `Viscoelastic` material model with a high viscosity, a uniform shear
modulus, and a small time step. This approximates a homogeneous, isotropic elastic half-space because the
`Viscoelastic` material model defines deformation with a Maxwell rheology that combines the viscous part and
the elastic part of the deformation in serial such that:

```{math}
\frac{1}{E} \frac{d \sigma}{dt} + \frac{\sigma}{\eta} = \frac{d \epsilon}{dt}
```

where $E$ is Young's modulus, $\eta$ is the viscosity, $\sigma$ is the stress, and $\epsilon$ is the strain.
We do not define $E$ within ASPECT, but instead work with the shear modulus $G$ when using viscoelasticity.
$E$ is related to $G$ with:
```{math}
G = \frac{E}{2(1 + \nu)}
```
where $\nu$ is Poisson's ratio. For an incompressible material, which we assume is that case in this benchmark, $\nu$ = 0.5.
When $G$ is very large, the first term goes to 0 and $\sigma$ is related to the time derivative of the strain.
If $\eta$ is very large, then the second term goes to 0, the time derivatives cancel out, and $\sigma$ is
related to $\epsilon$. By making the viscosity large, we approximate a purely elastic response. The model setup,
and the determination of the analytical solution is outlined below. It is important to note that in the analytical solution for the stresses,
there is no dependence on the properties of the elastic material, only on the magnitude of the load and the distance of a given point from the load.
While developing this benchmark, we also investigated the effect of a compressible equation of state, which had no effect on the solution.
Additionally, since the stress applied to the top boundary never changes, varying the shear modulus has no effect on the solution.
This is consistent with the analytic solution, which is independent of material parameters, and is also intuitive when considering that
the rheology of an elastic material is linear. For the same applied traction, increasing the shear modulus then
decreases the strain, while decreasing the shear modulus increases the strain.
By playing around with the input file in this benchmark, you can confirm that this is the case by seeing how the
velocity magnitude changes with different values of the shear modulus.

The analytic solution which describes the stresses within the domain are outlined in Section 4.6 of {cite:t}`dase96`, but the expressions for each of the 3 independent components is:

```{math}
\sigma_{xx} = \frac{p_0}{\pi} \left[\theta_2 - \theta_1 - \frac{1}{2} \cdot (\text{sin}(2 \theta_2) - \text{sin}(2 \theta_1) ) \right]
```

```{math}
\sigma_{xy} = \frac{p_0}{\pi} \left[(\text{sin}^2(\theta_2) - \text{sin}^2(\theta_1) ) \right]
```

```{math}
\sigma_{yy} = \frac{p_0}{\pi} \left[\theta_2 - \theta_1 + \frac{1}{2} \cdot (\text{sin}(2 \theta_2) - \text{sin}(2 \theta_1) ) \right]
```

where $\theta_1$ is the angle from the right edge of the line load to the current point, and $\theta_2$ is the angle from the
left edge of the line load to the current point, and $p_0$ is the magnitude of the traction introduced by the load, in Pa (see {numref}`fig:theta1-theta2-ref`). Within this
benchmark, the plugin `flamant_solution.cc` creates a postprocessor that visualizes the analytical solution onto the ASPECT grid for
direct comparison to the ASPECT stresses.

```{figure-md} fig:theta1-theta2-ref
<img src="dase96b.*" alt="Figure showing theta1 theta2"/>

Figure from {cite:t}`dase96` showing how the analytic solution is computed using angles from the edges of the load, $\theta_1$ and $\theta_2$.
```

Before comparing the ASPECT output to the analytic solution, we first show what the ASPECT stress distribution looks like.
The first three figures show the 3 unique components of the stress tensor in 2D, $\sigma_{xx}$ {numref}`fig:sigma-xx-aspect`,
$\sigma_{xy}$ {numref}`fig:sigma-xy-aspect`, and $\sigma_{yy}$ {numref}`fig:sigma-yy-aspect`.

```{figure-md} fig:sigma-xx-aspect
<img src="sigma_xx_displacement.*" alt="The ASPECT output for $\sigma_{xx}$"/>

Calculated $\sigma_{xx}$ from ASPECT. Vectors show the displacement direction and are
scaled by the magnitude of the velocity.
```

```{figure-md} fig:sigma-xy-aspect
<img src="sigma_xy_displacement.*" alt="The ASPECT output for $\sigma_{xy}$"/>

Calculated $\sigma_{xy}$ from ASPECT. Vectors show the displacement direction and are
scaled by the magnitude of the velocity.
```

```{figure-md} fig:sigma-yy-aspect
<img src="sigma_yy_displacement.*" alt="The ASPECT output for $\sigma_{yy}$"/>

Calculated $\sigma_{yy}$ from ASPECT. Vectors show the displacement direction and are
scaled by the magnitude of the velocity.
```

The direction of displacement is intuitive given the uniform vertical load applied in the center of the model.
Displacement is directly down beneath the center of the load and angles away from the load with increasing
lateral distance from the load.

In the next three figures below, we compare the stresses output by ASPECT with the analytic
solution within a 150 km $\times$ 150 km subset of the domain, centered on the load. To quantify how well ASPECT reproduces the analytic
solution, we calculate a relative error with:

```{math}
\chi = \text{abs}\left( \frac{\sigma_i^{analytic} - \sigma_i^{ASPECT}}{\sigma_i^{analytic}} \right)
```

where $i$ = $xx$, $xy$, or $yy$ is a given stress component.

```{figure-md} fig:sigma-xx-error
<img src="sigma_xx_model_error.*" alt="The relative error for $\sigma_{xx}$"/>

The relative error for $\sigma_{xx}$.
```

```{figure-md} fig:sigma-xy-error
<img src="sigma_xy_model_error.*" alt="The relative error for $\sigma_{xy}$"/>

The relative error for $\sigma_{xy}$.
```

```{figure-md} fig:sigma-yy-error
<img src="sigma_yy_model_error.*" alt="The relative error for $\sigma_{yy}$"/>

The relative error for $\sigma_{yy}$.
```

ASPECT most accurately reproduces $\sigma_{yy}$, which is clearly visible in {numref}`fig:sigma-yy-error`.
Similarly, ASPECT accurately reproduces $\sigma_{xy}$ as well, the vertical line directly in the center of
{numref}`fig:sigma-xy-error` is a numerical artifact since the stress is 0 here, meaning that in our relative
error calculation we are dividing by 0. The highest deviation occurs for $\sigma_{xx}$, where the misfit increases
with depth, and is higher directly beneath the load ({numref}`fig:sigma-xx-error`). We can also visualize the stress
along the surface of the ASPECT model, as well as with depth directly beneath the center of the load ({numref}`fig:surface-depth-stresses`).
Along the surface, with the exception of the points at the left and right edges of the line load, there is good agreement
with the predicted stresses in ASPECT and with the analytic solution. The oscillations in the stresses at the edges of the
load are numerical artifacts introduced by the discontinuity in the traction at this point, from 0 to $p_0$. ASPECT also
predicts non-zero $\sigma_{xx}$ towards the edges of the model domain (top row of {numref}`fig:surface-depth-stresses`).
With depth beneath the middle of the load, ASPECT accurately reproduces all of the stress components at depths $\leq$40 km.
At depths greater than 40 km, $\sigma_{xx}$ begins to noticeably deviate from the analytic solution, which is likely due to the
increasing proximity to the boundary conditions. The analytic solution assumes an infinite elastic half-space, but in this
example we place a physical boundary at 300 km depth.

```{figure-md} fig:surface-depth-stresses
<img src="surface_and_depth_stresses.png" alt="Plotting stresses along the surface and beneath the center of the load."/>

Plotting stresses along the surface and beneath the center of the load. Note the very small magnitudes of stress in the
bottom middle panel.
```

Output from ASPECT's topography postprocessor compared to the analytic solution
is shown in {numref}`fig:flexure-comparison`. Both the flexural amplitude and the
flexural wavelength are accurately recovered.

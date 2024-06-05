# Particle integration scheme benchmark

*This section was contributed by Gabriel Johnston and Rene Gassmöller*

This benchmark is designed to test the convergence of various particle integration schemes (Forward-Euler, RK2, RK4) as detailed in {cite}`gassmoller:etal:2018`. The first benchmark involves 100 particles moving in circular motion around a ring in a 2D box. The error is measured as the distance between the final and initial positions, which should ideally be zero. The second benchmark involves particles advected in a flow field defined by {math:numref}`eq:flow`.

```{math}
:label: eq:flow
\mathbf{u}_x(t) = 1 + e^t \quad \text{and} \quad \mathbf{u}_y(t) = 1 + e^t.
```
The final analytic position is represented by {math:numref}`eq:fposition` and the distance between the initial and final position is defined by {math:numref}`eq:fidistance`.

```{math}
:label: eq:fposition
\mathbf{x}(1) = \mathbf{x}_0 + \int_{0}^{1} \mathbf{u}(t) \, dt
```

```{math}
:label: eq:fidistance
d = \| \mathbf{x}_1 - \mathbf{x}_0 \|_2 = \sqrt{2} e
```

```{figure-md} fig:particle-integration-scheme
<img src="convergence.*" width="100%" />

Figure shows the convergence order of particle advection schemes in spatially variable flow (a) and temporally changing flow (b). Additionally, the figure (right) depicts the RK4 scheme using the analytical velocity solution for the halfway time instead of interpolating the discrete time steps. Figure adapted from {cite}`gassmoller:etal:2018`
```

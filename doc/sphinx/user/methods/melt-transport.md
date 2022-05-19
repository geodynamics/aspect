(sec:methods:melt-transport)=
# Calculations with melt transport

The original formulation of the equations in {ref}`2.1][] describes
the movement of solid mantle material. These computations also allow for
taking into account how partially molten material changes the material
properties and the energy balance through the release of latent heat. However,
this will not consider melt extraction or any relative movement between melt
and solid and there might be problems where the transport of melt is of
interest. Thus, ASPECT allows for solving
additional equations describing the behavior of silicate melt percolating
through and interacting with a viscously deforming host rock. This requires
the advection of a compositional field representing the volume fraction of
melt present at any given time (the porosity $\phi$), and also a change of the
mechanical part of the system. The latter is implemented using the approach of
(Keller, May, and Kaus 2013) and changes the Stokes system to
```{math}
\begin{align}
  \label{eq:stokes-1-melt}
  -\nabla \cdot \left[2\eta \left(\varepsilon(\mathbf{u}_s)
                                  - \frac{1}{3}(\nabla \cdot \mathbf{u}_s)\mathbf 1\right)
                \right] + \nabla p_f + \nabla p_c  &=
  \rho \mathbf g
  & \qquad
  & \textrm{in $\Omega$},
  \\
  \label{eq:stokes-2-melt}
  \nabla \cdot \mathbf{u}_s - \nabla \cdot K_D \nabla p_f
  - K_D \nabla p_f \cdot \frac{\nabla \rho_f}{\rho_f}
  &=
  - \nabla \cdot K_D \rho_f \mathbf g
  \notag
  \\
  &\quad
  + \Gamma \left( \frac{1}{\rho_f} - \frac{1}{\rho_s} \right)
  \\
  &\quad
  - \frac{\phi }{\rho_f} \mathbf{u}_s \cdot \nabla\rho_f
  - \frac{1 - \phi }{\rho_s} \mathbf{u}_s \cdot \nabla\rho_s
  \notag
  \\
  &\quad
  - K_D \mathbf g \cdot \nabla \rho_f
  & \qquad
  & \textrm{in $\Omega$},
  \notag
  \\
  \label{eq:stokes-3-melt}
  \nabla \cdot \mathbf{u}_s + \frac{p_c}{\xi}
  &=
  0.
\end{align}
```

We use the indices $s$ to indicate properties of the solid and $f$ for the
properties of the fluid. The equations are solved for the solid velocity
$\mathbf{u}_s$, the fluid pressure $p_f$, and an additional variable, the
compaction pressure $p_c$, which is related to the fluid and solid pressure
through the relation $p_c = (1-\phi) (p_s-p_f)$. $K_D$ is the Darcy
coefficient, which is defined as the quotient of the permeability and the
fluid viscosity and $\Gamma$ is the melting rate. $\eta$ and $\xi$ are the
shear and compaction viscosities and can depend on the porosity, temperature,
pressure, strain rate and composition. However, there are various laws for
these quantities and so they are implemented in the material model. Common
formulations for the dependence on porosity are
$\eta = (1-\phi) \eta_0 e^{-\alpha_\phi \phi}$ with
$\alpha_\phi \approx 25...30$ and $\xi = \eta_0 \phi^{-n}$ with $n \approx 1$.

To avoid the density gradients in Equation{math:numref}`eq::stokes-2-melt`40],
which would have to be specified individually for each material model by the
user, we can use the same method as for the mass conservation (described in
{ref}`2.10.3][]) and assume the change in solid density is dominated
by the change in static pressure, which can be written as
$\nabla p_s \approx \nabla p_{\text{static}} \approx \rho_s \textbf{g}$. This
finally allows us to write
```{math}
\frac{1}{\rho_s} \nabla \rho_s
\approx \frac{1}{\rho_s} \frac{\partial \rho_s}{\partial p_s} \nabla p_s
\approx \frac{1}{\rho_s} \frac{\partial \rho_s}{\partial p_s} \nabla p_s
\approx \frac{1}{\rho_s} \frac{\partial \rho_s}{\partial p_s} \rho_s \textbf{g}
\approx \beta_s \rho_s \textbf{g}.
```
where $\beta_s$ is the compressibility of
the solid. In the paper that describes the implementation (Dannberg and
Heister 2016), $\kappa$ is used for the compressibility. We change the
variable here to be consistent throughout the manual.

For the fluid pressure, choosing a good approximation depends on the model
parameters and setup (see (Dannberg and Heister 2016)). Hence, we make
$\nabla \rho_{f}$ a model input parameter, which can be adapted based on the
forces that are expected to be dominant in the model. We can then replace the
second equation by
```{math}
\begin{align}
\nabla \cdot \mathbf{u}_s - \nabla \cdot K_D \nabla p_f
  - K_D \nabla p_f \cdot \frac{\nabla \rho_f}{\rho_f}
  &=
  - \nabla \cdot (K_D\rho_f \mathbf g)
  \\
  &\quad
  + \Gamma \left( \frac{1}{\rho_f} - \frac{1}{\rho_s} \right)
  \notag
  \\
  &\quad
  - \frac{\phi }{\rho_f} \mathbf{u}_s \cdot \nabla\rho_f
  - (\mathbf{u}_s \cdot \mathbf g ) (1 - \phi) \beta_s \rho_s
  \notag
  \\
  &\quad
  - K_D \mathbf g \cdot \nabla \rho_f .
  \notag
\end{align}
```
The melt velocity is computed as
```{math}
\mathbf{u}_f =  \mathbf{u}_s - \frac{K_D}{\phi} (\nabla p_f - \rho_f g),
```
but is only used for postprocessing purposes and for computing the time step
length.

<div class="center">

</div>

Moreover, melt transport requires an advection equation for the porosity field
$\phi$:
```{math}
\begin{align}
  \label{eq:porosity}
  \rho_s \frac{\partial (1 - \phi)}{\partial t} + \nabla \cdot \left[ \rho_s (1 - \phi) \mathbf{u}_s \right]
  &=
  - \Gamma
  & \quad
  & \textrm{in $\Omega$},
  i=1\ldots C
\end{align}
```

In order to solve this equation in the same way as the other advection
equations, we replace the second term of the equation by:

```{math}
\nabla \cdot \left[ \rho_s (1 - \phi) \mathbf{u}_s \right]
= \left( 1-\phi \right) \left( \rho_s \nabla \cdot \mathbf{u}_s
+ \nabla \rho_s \cdot \mathbf{u}_s \right)
- \nabla \phi \cdot \rho_s \mathbf{u}_s
```
Then we use the same method as
described above and assume again that the change in density is dominated by
the change in static pressure
```{math}
\frac{1}{\rho_s} \nabla \rho_s \cdot \mathbf{u}_s
\approx \beta_s \rho_s \textbf{g} \cdot \mathbf{u}_s``` so we get
```{math}
\frac{\partial \phi}{\partial t} + \mathbf{u}_s \cdot \nabla \phi
= \frac{\Gamma}{\rho_s}
+ (1 - \phi) (\nabla \cdot \mathbf{u}_s + \beta_s \rho_s \textbf{g} \cdot \mathbf{u}_s ).
```

More details on the implementation can be found in (Dannberg and Heister
2016). A benchmark case demonstrating the propagation of solitary waves can be
found in {ref}`sec:benchmark-solitary_wave`41].

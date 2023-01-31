(sec:methods:choosing-a-formulation:mass-conservation-approx)=
# Mass conservation approximation

First, we have to choose how to approximate the conservation of mass: $\nabla \cdot (\rho \mathbf u) = 0$ (see Equation {math:numref}`eq:stokes-2`).
We provide the following options, which can be selected in the parameter file in the subsection `Formulation/Mass conservation` (see also {ref}`parameters:Formulation/Mass_20conservation`):

-   "incompressible":
```{math}
\nabla \cdot \textbf{u} = 0,
```

-   "isothermal compression":
    ```{math}
    \nabla \cdot \textbf{u} = -\rho \beta \textbf{g} \cdot \textbf{u},
    ```
    where $\beta = \frac{1}{\rho} \frac{\partial \rho}{\partial p}$ is the compressibility, and is defined in the material model.
    Despite the name, this approximation can be used either for isothermal compression (where $\beta = \beta_T$) or isentropic compression (where $\beta = \beta_S$).
    The material model determines which compressibility is used.
    This is an explicit compressible mass equation where the velocity $\textbf{u}$ on the right-hand side is an extrapolated velocity from the last timesteps.

-   "hydrostatic compression":
    ```{math}
    \nabla \cdot \textbf{u}
    = - \left( \frac{1}{\rho} \left( \frac{\partial \rho}{\partial p} \right)_{T} \rho \textbf{g} + \frac{1}{\rho} \left( \frac{\partial \rho}{\partial T} \right)_{p} \nabla T \right) \cdot \textbf{u}
    = - \left( \beta_T \rho \textbf{g} - \alpha \nabla T \right) \cdot \textbf{u}
    ```
    where $\beta_T = \frac{1}{\rho} \left(\frac{\partial \rho}{\partial p} \right)_{T}$ is the isothermal compressibility, $\alpha = - \frac{1}{\rho} \left(\frac{\partial \rho}{\partial T} \right)_{p}$ is the thermal expansion coefficient, and both are defined in the material model.
    The approximation made here is that $\nabla p = \rho \textbf{g}$.

-   "reference density profile":
    ```{math}
    \nabla \cdot \textbf{u} = -\frac{1}{\bar{\rho}} \frac{\partial \bar{\rho}}{\partial z} \frac{\textbf{g}}{\|\textbf{g}\|} \cdot \textbf{u},
    ```
    where the reference profiles for the density $\bar{\rho}$ and the density gradient $\frac{\partial \bar{\rho}}{\partial z}$ provided by the adiabatic conditions model ({ref}`sec:methods:initial-conditions`) are used.
    Note that the gravity is assumed to point downwards in depth direction.
    This is the explicit mass equation where the velocity $\textbf{u}$ on the right-hand side is an extrapolated velocity from the last timesteps.

-   "implicit reference density profile":
    ```{math}
    \nabla \cdot \textbf{u} + \frac{1}{\bar{\rho}} \frac{\partial \bar{\rho}}{\partial z} \frac{\textbf{g}}{\|\textbf{g}\|} \cdot \textbf{u} = 0,
    ```
    which uses the same approximation for the density as "reference density profile," but implements this term on the left-hand side instead of the right-hand side of the mass conservation equation.
    This effectively uses the current velocity $\textbf{u}$ instead of an explicitly extrapolated velocity from the last timesteps.

-   "ask material model," which uses "isothermal compression" if the material model reports that it is compressible and "incompressible" otherwise.

:::{admonition} The stress tensor approximation
:class: note

If a medium is incompressible, that is, if the mass conservation equation reads $\nabla \cdot \textbf{u} = 0$, then the shear stress in the momentum and temperature equation simplifies from
```{math}
\tau = 2\eta\left(\varepsilon\left(u\right) - \frac{1}{3}\left(\nabla\cdot\textbf{u}\right)1\right)
```
to
```{math}
\tau = 2\eta\varepsilon\left(\textbf{u}\right)
```
:::

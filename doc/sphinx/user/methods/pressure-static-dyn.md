(sec:methods:static-v-dynamic)=
# Static or dynamic pressure?

One could reformulate equation {math:numref}`eq:stokes-1` somewhat.
To this end, let us say that we would want to represent the pressure $p$ as the sum of two parts that we will call static and dynamic, $p=p_s+p_d$.
If we assume that $p_s$ is already given, then we can replace {math:numref}`eq:stokes-1` by
```{math}
\begin{aligned}
  -\nabla \cdot 2\eta \nabla \mathbf u + \nabla p_d = \rho\mathbf g - \nabla p_s.
\end{aligned}
```
One typically chooses $p_s$ as the pressure one would get if the whole medium were at rest - i.e., as the hydrostatic pressure.
This pressure can be computed noting that {math:numref}`eq:stokes-1` reduces to
```{math}
\begin{aligned}
  \nabla p_s = \rho(p_s,T_s,\mathbf x)\mathbf g = \bar\rho \mathbf g
\end{aligned}
```
in the absence of any motion where $T_s$ is some static temperature field (see also {ref}`sec:methods:initial-conditions`).
This, our rewritten version of {math:numref}`eq:stokes-1` would look like this:
```{math}
\begin{aligned}
  -\nabla \cdot 2\eta   \nabla \mathbf u + \nabla p_d = \left[\rho(p,T,\mathbf x)-\rho(p_s,T_s,\mathbf x)\right]\mathbf g.
\end{aligned}
```
In this formulation, it is clear that the quantity that drives the fluid flow is in fact the *buoyancy* caused by the *variation* of densities, not the density itself.

This reformulation has a number of advantages and disadvantages:

-   One can notice that in many realistic cases, the dynamic component $p_d$ of the pressure is orders of magnitude smaller than the static component $p_s$.
  For example, in the Earth, the two are separated by around 6 orders of magnitude at the bottom of the Earth's mantle.
  Consequently, if one wants to solve the linear system that arises from discretization of the original equations, one has to solve it to a significant degree of accuracy (6-7 digits) to get the dynamic part of the pressure correct to even one digit.
  This entails a very significant numerical effort, and one that is not necessary if we can split the pressure in a way so that the pre-computed static pressure $p_s$ (or, rather, the density using the static pressure and temperature from which $p_s$ results) absorbs the dominant part and one only has to compute the remaining, dynamic pressure to 2 or 3 digits of accuracy, rather than the corresponding 7-8 for the total pressure.

-   On the other hand, the pressure $p_d$ one computes this way is not immediately comparable to quantities that we use to look up pressure-dependent quantities such as the density.
  Rather, one needs to first find the static pressure as well (see {ref}`sec:methods:initial-conditions`) and add the two together before they can be used to look up material properties or to compare them with experimental results.
  Consequently, if the pressure a program outputs (either for visualization, or in the internal interfaces to parts of the code where users can implement pressure- and temperature-dependent material properties) is only the dynamic component, then all of the consumers of this information need to convert it into the total pressure when comparing with physical experiments.
  Since any code implementing realistic material models has a great many of these places, there is a large potential for inadvertent errors and bugs.

-   Finally, the definition of a reference density $\rho(p_s,T_s,\mathbf x)$ derived from static pressures and temperatures is only simple if we have incompressible models and under the assumption that the temperature-induced density variations are small compared to the overall density.
In this case, we can choose $\rho(p_s,T_s,\mathbf x)=\rho_0$ with a constant reference density $\rho_0$.
On the other hand, for more complicated models, it is not a priori clear which density to choose since we first need to compute static pressures and temperatures - quantities that satisfy equations that introduce boundary layers, may include phase changes releasing latent heat, and where the density may have discontinuities at certain depths, see {ref}`sec:methods:initial-conditions`.

    Thus, if we compute adiabatic pressures and temperatures $\bar p_s,\bar T_s$ under the assumption of a thermal boundary layer worth 900 Kelvin at the top, and we get a corresponding density profile $\bar\rho=\rho(\bar p_s,\bar T_s, \mathbf x)$, but after running for a few million years the temperature turns out to be so that the top boundary layer has a jump of only 800 Kelvin with corresponding adiabatic pressures and temperatures $\hat p_s,\hat T_s$, then a more appropriate density profile would be $\hat\rho=\rho(\hat p_s,\hat T_s, \mathbf x)$.

    The problem is that it may well be that the erroneously computed density profile $\hat \rho$ does *not* lead to a separation where $|p_d|\ll|p_s|$ because, especially if the material undergoes phase changes, there will be entire areas of the computational domain in which $|\rho-\hat \rho_s|\ll |\rho|$ but $|\rho-\bar \rho_s|\not\ll |\rho|$.
    Consequently the benefits of lesser requirements on the iterative linear solver would not be realized.

We do note that most of the codes available today, and that we are aware of, split the pressure into static and dynamic parts either internally or require the user to specify the density profile as the difference between the true and the hydrostatic density.
This may, in part, be due to the fact that historically most codes were written to solve problems in which the medium was considered incompressible, i.e., where the definition of a static density was simple.

On the other hand, we intend ASPECT to be a code that can solve more general models for which this definition is not as simple.
As a consequence, we have chosen to solve the equations as stated originally - i.e., we solve for the *full* pressure rather than just its *dynamic* component.
With most traditional methods, this would lead to a catastrophic loss of accuracy in the dynamic pressure since it is many orders of magnitude smaller than the total pressure at the bottom of the earth mantle.
We avoid this problem in ASPECT by using a cleverly chosen iterative solver that ensures that the full pressure we compute is accurate enough so that the dynamic pressure can be extracted from it with the same accuracy one would get if one were to solve for only the dynamic component.
The methods that ensure this are described in detail in {cite:t}`kronbichler:etal:2012` and in particular in the appendix of that paper.

:::{note}
By default, ASPECT uses the full pressure in the equations, and only prescribing density deviations from a reference state on the right-hand side of {math:numref}`eq:stokes-1` would lead to negative densities in the energy equation {math:numref}`eq:temperature`.
However, when using one of the approximations described in {ref}`sec:methods:approximate-equations`, the energy balance uses the reference density $\bar\rho$ instead of the full density, which makes it possible to formulate the Stokes system in terms of the dynamic instead of the full pressure.
In order to do this, one would have to use a material model (see {ref}`sec:extending:plugin-types:material-models`) in which the density is in fact a density variation, and then the pressure solution variable would only be the dynamic pressure.
:::

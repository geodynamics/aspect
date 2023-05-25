(sec:methods:coefficients)=
# Coefficients

The equations above contain a significant number of coefficients that we will discuss in the following.
In the most general form, many of these coefficients depend nonlinearly on the solution variables pressure $p$, temperature $T$ and, in the case of the viscosity, on the strain rate $\varepsilon(\mathbf{u})$.
If compositional fields $\mathfrak c=\{c_1,\ldots,c_C\}$ are present (i.e., if $C>0$), coefficients may also depend on them.
Alternatively, they may be parameterized as a function of the spatial variable $\mathbf{x}$.
ASPECT allows both kinds of parameterizations.

Note that below we will discuss examples of the dependence of coefficients on other quantities; which dependence is actually implemented in the code is a different matter.
As we will discuss in {ref}`parameters` and {ref}`cha:extending`, some versions of these models are already implemented and can be selected from the input parameter file; others are easy to add to ASPECT by providing self-contained descriptions of a set of coefficients that the rest of the code can then use without a need for further modifications.

Concretely, we consider the following coefficients and dependencies:

-   *The viscosity $\eta=\eta(p,T,\varepsilon(\mathbf u),\mathfrak c,\mathbf x)$:* Units $\text{ Pa . s} = \text{ kg}\frac{1}{\text{ m . s}}$.

    The viscosity is the proportionality factor that relates total forces (external gravity minus pressure gradients) and fluid velocities $\mathbf u$.
    The simplest models assume that $\eta$ is constant, with the constant often chosen to be on the order of $10^{21} \text{ Pa . s}$.

    More complex (and more realistic) models assume that the viscosity depends on pressure, temperature and strain rate.
    Since this dependence is often difficult to quantify, one modeling approach is to make $\eta$ spatially dependent.

-   *The density $\rho=\rho(p,T,\mathfrak c,\mathbf x)$:* Units $\frac{\text{ kg}}{\text{ m}^3}$.

    In general, the density depends on pressure and temperature, both through pressure compression, thermal expansion, and phase changes the material   may undergo as it moves through the pressure-temperature phase diagram.

    The simplest parameterization for the density is to assume a linear dependence on temperature, yielding the form $\rho(T)=\rho_{\text{ref}}[1-\alpha (T-T_{\text{ref}})]$ where $\rho_{\text{ref}}$ is the reference density at temperature $T_{\text{ref}}$ and $\alpha$ is the linear thermal expansion coefficient.
    For the Earth's mantle, typical values for this parameterization would be $\rho_{\text{ref}}=3300\frac{\text{ kg}}{\text{ m}^3}$, $T_{\text{ref}}=293 \text{ K}$, $\alpha=2\times 10^{-5} \frac{1}{\mathrm{K}}$.

-   *The gravity vector $\mathbf g=\mathbf g(\mathbf x)$:* Units $\frac{\text{ m}}{\textrm{s}^2}$.

    Simple models assume a radially inward gravity vector of constant magnitude (e.g., the surface gravity of Earth, $9.81 \frac{\text{ m}}{\textrm{s}^2}$), or one that can be computed analytically assuming a homogeneous mantle density.

    A physically self-consistent model would compute the gravity vector as $\mathbf g = -\nabla \varphi$ with a gravity potential $\varphi$ that
    satisfies $-\Delta\varphi=4\pi G\rho$ with the density $\rho$ from above and $G$ the universal constant of gravity.
    This would provide a gravity vector that changes as a function of time.
    Such a model is not currently implemented.

-   *The specific isobaric heat capacity
    $C_p=C_p(p,T,\mathfrak c,\mathbf x)$:* Units J/kg/K =
    m<sup>2</sup>/s<sup>2</sup>/K.

    The specific heat capacity denotes the amount of energy needed to increase the temperature of one kilogram of material by one Kelvin at constant pressure.
    Wikipedia lists a value of $790 \text{ J/kg/K}$ for granite[^footnote1].
    For the Earth's mantle, a value of $1250 \text{ J/kg/K}$ is within the range suggested by the literature.

-   *The thermal conductivity $k=k(p,T,\mathfrak c,\mathbf x)$:* Units $\frac{\textrm{W}}{\text{ m}\cdot\text{ K}}=\frac{\text{ kg}\cdot\text{ m}}{\textrm{s}^3\cdot\text{ K}}$.

    The thermal conductivity denotes the amount of thermal energy flowing through a unit area for a given temperature gradient.
    It depends on the material and as such will from a physical perspective depend on pressure and temperature due to phase changes of the material as well as through different mechanisms for heat transport (see, for example, the partial transparency of perovskite, the most abundant material in the Earth's mantle, at pressures above around 120 GPa {cite}`badro:etal:2004`.

    As a rule of thumb for its order of magnitude, Wikipedia quotes values of $1.83$-$2.90\frac{\textrm{W}}{\text{ m}\cdot\text{ K}}$ for sandstone and $1.73$-$3.98\frac{\textrm{W}}{\text{ m}\cdot\text{ K}}$ for granite[^footnote2].
    The values in the mantle are almost certainly higher than this though probably not by much.
    The exact value is not really all that important: heat transport through convection is several orders of magnitude more important than through thermal conduction.

    The thermal conductivity $k$ is often expressed in terms of the *thermal diffusivity* $\kappa$ using the relation $k = \rho C_p \kappa$.

-   *The intrinsic specific heat production $H=H(\mathbf x)$:* Units $\frac{\textrm{W}}{\text{ kg}}=\frac{\text{ m}^2}{\textrm{s}^3}$.

    This term denotes the intrinsic heating of the material, for example due to the decay of radioactive material.
    As such, it depends not on pressure or temperature, but may depend on the location due to different chemical composition of material in the Earth's mantle.
    The literature suggests a value of $\gamma=7.4\times 10^{-12}\frac{\textrm{W}}{\text{ kg}}$.

-   *The thermal expansion coefficient $\alpha=\alpha(p,T,\mathfrak c ,\mathbf x)$:* Units $\frac{1}{\text{ K}}$.

    This term denotes by how much the material under consideration expands due to temperature increases at constant pressure.
    This coefficient is defined as $\alpha = -\frac{1}{\rho} \left(\frac{\partial \rho}{\partial T}\right)_{p}$, where the negative sign is due the fact that the density *decreases* as a function of temperature.
    Alternatively, if one considers the *volume* $V=V(T)$ a piece of material of mass $M$ occupies, $V=\frac{M}{\rho}$, then the thermal expansion coefficient is defined as the relative increase in volume, $\alpha=\frac{1}{V}\frac{\partial V(T)}{\partial T}$, because $\frac{\partial V(T)}{\partial T} = \frac{\partial \frac{M}{\rho}}{\partial T} = -\frac{M}{\rho^2} \frac{\partial \rho}{\partial T} = -\frac{V}{\rho} \frac{\partial \rho}{\partial T}$.

    The literature suggests that values of $\alpha=1\times 10^{-5}\frac{1}{\text{ K}}$ at the core-mantle boundary and $\alpha=4\times 10^{-5}\frac{1}{\text{ K}}$ are appropriate for Earth.

-   *The isothermal compressibility $\beta_T=\beta_T(p,T,\mathfrak c ,\mathbf x)$:* Units $\frac{1}{\textrm{Pa}}$.

    This term quantifies how much the material under consideration contracts due to pressure increases at constant temperature.
    This coefficient is defined as $\beta_T = \frac{1}{\rho} \left( \frac{\partial \rho}{\partial p} \right)_{T}$.
    Alternatively, if one considers the *volume* $V=V(p, T)$ a piece of material of mass $M$ occupies, $V=\frac{M}{\rho}$, then the isothermal compressibility is defined as the relative increase in volume, $\beta=\frac{1}{V}\left(\frac{\partial V(p, T)}{\partial p}\right)_{T}$, because $\frac{\partial V(p, T)}{\partial p} = \frac{\partial \frac{M}{\rho}}{\partial p} = -\frac{M}{\rho^2} \frac{\partial \rho}{\partial p} = -\frac{V}{\rho} \frac{\partial \rho}{\partial p}$.

    Values of $\beta=10^{-12}$ - $10^{-11} \frac{1}{\textrm{Pa}}$ are reasonable for Earth's mantle, with values decreasing by about a factor of 5 between the shallow lithosphere and core-mantle boundary.

-   *The isentropic/adiabatic compressibility* $\beta_S=\beta_S(p,T,\mathfrak c ,\mathbf x)$:* Units $\frac{1}{\textrm{Pa}}$.

    This term quantifies how much the material under consideration contracts due to pressure increases at constant entropy.
    This coefficient is defined as $\beta_S = \frac{1}{\rho} \left( \frac{\partial \rho}{\partial p} \right)_{S}$.
    Alternatively, if one considers the *volume* $V=V(p, T)$ a piece of material of mass $M$ occupies, $V=\frac{M}{\rho}$, then the isentropic compressibility is defined as the relative increase in volume, $\beta=\frac{1}{V}\left(\frac{\partial V(p, T)}{\partial p}\right)_{S}$, because $\frac{\partial V(p, T)}{\partial p} = \frac{\partial \frac{M}{\rho}}{\partial p} = -\frac{M}{\rho^2} \frac{\partial \rho}{\partial p} = -\frac{V}{\rho} \frac{\partial \rho}{\partial p}$.
    The isentropic and isothermal compressibility are related by the expression:
    %
    \begin{equation}
    \beta_S = \beta_T - \frac{\alpha^2 T}{\rho C_p}
    \end{equation}
    %
    The ratio of the compressibilities decreases with increasing temperature and increases with increasing pressure.
    In the Earth's convecting mantle, $\beta_S/\beta_T = 0.92$-$0.98$.
    Different mineral assemblages have different values of this ratio under the same conditions.
    For example, the upper-lower boundary may exhibit a 3-4% drop in $\beta_S / \beta_T$ as a result of a 40% lower $C_p$ of bridgmanite-periclase assemblages relative to the olivine polymorphs.

-   *The change in entropy $\Delta S$ at a phase transition together with the derivatives of the phase function $X=X(p,T,\mathfrak c,\mathbf x)$ with regard to temperature and pressure:* Units J/kg/K<sup>2</sup> ($-\Delta S \frac{\partial X}{\partial T}$) and m<sup>3</sup>/kg/K ($\Delta S \frac{\partial X}{\partial p}$).

    When material undergoes a phase transition, the entropy changes due to release or consumption of latent heat.
    However, phase transitions occur gradually and for a given chemical composition it depends on temperature and pressure which phase prevails.
    Thus, the latent heat release can be calculated from the change of entropy $\Delta S$ and the derivatives of the phase function $\frac{\partial X}{\partial T}$ and $\frac{\partial X}{\partial p}$.
    These values have to be provided by the material model, separately for the coefficient $-\Delta S \frac{\partial X}{\partial T}$ on the left-hand side and $\Delta S \frac{\partial X}{\partial p}$ on the right-hand side of the temperature equation.
    However, they may be either approximated with the help of an analytic phase function, employing data from a thermodynamic database or in any other way that seems appropriate to the user.

:::{toctree}
self-consistency.md
averaging.md
:::

[^footnote1]: See <http://en.wikipedia.org/wiki/Specific_heat>
[^footnote2]: See <http://en.wikipedia.org/wiki/Thermal_conductivity> and <http://en.wikipedia.org/wiki/List_of_thermal_conductivities>.

(sec:methods:dimensionalize)=
# Dimensional or non-dimensionalized equations?

Equations {math:numref}`eq:stokes-1` - {math:numref}`eq:temperature` are stated in their physically correct form.
One would usually interpret them in a way that the various coefficients such as the viscosity, density and thermal conductivity $\eta,\rho,\kappa$ are given in their correct physical units, typically expressed in a system such as the meter, kilogram, second (MKS) system that is part of the [SI](https://en.wikipedia.org/wiki/SI) system.
This is certainly how we envision ASPECT to be used: with geometries, material models, boundary conditions and initial values to be given in their correct physical units.
As a consequence, when ASPECT prints information about the simulation onto the screen, it typically does so by using a postfix such as `m/s` to indicate a velocity or `W/m^2` to indicate a heat flux.

:::{note}
For mantle convection simulations, it is often convenient to work with time units of *years* instead of *seconds*.
The flag "`Use years in output instead of seconds`" ({ref}`parameters:global`) in the input file determines how input and output parameters with units of time or velocity are interpreted.
For details, see {ref}`sec:methods:dimensionalize:years-or-seconds`.
:::

That said, in reality, ASPECT has no preferred system of units as long as every material constant, geometry, time, etc., are all expressed in the same system.
In other words, it is entirely legitimate to implement geometry and material models in which the dimension of the domain is one, density and viscosity are one, and the density variation as a function of temperature is scaled by the Rayleigh number - i.e., to use the usual non-dimensionalization of the equations {math:numref}`eq:stokes-1` - {math:numref}`eq:temperature`.
Some of the cookbooks in {ref}`cha:cookbooks` use this non-dimensional form; for example, the simplest cookbook in {ref}`sec:cookbooks:convection-box` as well as the SolCx, SolKz and inclusion benchmarks in {ref}`sec:benchmarks:solcx`, are such cases.
Whenever this is the case, output showing units `m/s` or `W/m^2` clearly no longer have a literal meaning.
Rather, the unit postfix must in this case simply be interpreted to mean that the number that precedes the first is a velocity and a heat flux in the second case.

In other words, whether a computation uses physical or non-dimensional units really depends on the geometry, material, initial and boundary condition description of the particular case under consideration - ASPECT will simply use whatever it is given.
Whether one or the other is the more appropriate description is a decision we purposefully leave to the user.
There are of course good reasons to use non-dimensional descriptions of realistic problems, rather than to use the original form in which all coefficients remain in their physical units.
On the other hand, there are also downsides:

-   Non-dimensional descriptions, such as when using the [Rayleigh](https://en.wikipedia.org/wiki/Rayleigh_number) number to indicate the relative strength of convective to diffusive thermal transport, have the advantage that they allow one to reduce a system to its essence.
For example, it is clear that we get the same behavior if one increases both the viscosity and the thermal expansion coefficient by a factor of two because the resulting Rayleigh number; similarly, if we were to increase the size of the domain by a factor of 2 and thermal diffusion coefficient by a factor of 8.
In both of these cases, the non-dimensional equations are exactly the same. On the other hand, the equations in their physical unit form are different and one may not see that the result of this variations in coefficients will be exactly the same as before.
Using non-dimensional variables therefore reduces the space of independent parameters one may have to consider when doing parameter studies.

-   From a practical perspective, equations {math:numref}`eq:stokes-1` - {math:numref}`eq:temperature` are often ill-conditioned in their original form: the two sides of each equation have physical units different from those of the other equations, and their numerical values are often vastly different[^footnote1].
Of course, these values can not be compared: they have different physical units, and the ratios between these values depends on whether we choose to measure lengths in meters or kilometers, for example.
Nevertheless, when implementing these equations in software, at one point or another, we have to work with numbers and at this point the physical units are lost.
If one does not take care at this point, it is easy to get software in which all accuracy is lost due to round-off errors.
On the other hand, non-dimensionalization typically avoids this since it normalizes all quantities so that values that appear in computations are typically on the order of one.

-   On the downside, the numbers non-dimensionalized equations produce are not immediately comparable to ones we know from physical experiments.
This is of little concern if all we have to do is convert every output number of our program back to physical units.
On the other hand, it is more difficult and a source of many errors if this has to be done inside the program, for example, when looking up the viscosity as a pressure-, temperature- and strain-rate-dependent function: one first has to convert pressure, temperature and strain rate from non-dimensional to physical units, look up the corresponding viscosity in a table, and then convert the viscosity back to non-dimensional quantities.
Getting this right at every one of the dozens or hundreds of places inside a program and using the correct (but distinct) conversion factors for each of these quantities is both a challenge and a possible source of errors.

-   From a mathematical viewpoint, it is typically clear how an equation needs to be non-dimensionalized if all coefficients are constant.
However, how is one to normalize the equations if, as is the case in the Earth's mantle, the viscosity varies by several orders of magnitude?
In cases like these, one has to choose a reference viscosity, density, etc.
While the resulting non-dimensionalization retains the universality of parameters in the equations, as discussed above, it is not entirely clear that this would also retain the numerical stability if the reference values are poorly chosen.

As a consequence of such considerations, most codes in the past have used non-dimensionalized models.
This was aided by the fact that until recently and with notable exceptions, many models had constant coefficients and the difficulties associated with variable coefficients were not a concern.
On the other hand, our goal with ASPECT is for it to be a code that solves realistic problems using complex models and that is easy to use.
Thus, we allow users to input models in physical or non-dimensional units, at their discretion. We believe that this makes the description of realistic models simpler.
On the other hand, ensuring numerical stability is not something users should have to be concerned about, and is taken care of in the implementation of ASPECT's core (see the corresponding section in {cite:t}`kronbichler:etal:2012`).

:::{toctree}
years-or-seconds.md
:::


[^footnote1]: To illustrate this, consider convection in the Earth as a back-of-the-envelope example.
With the length scale of the mantle $L=3\times 10^{6}\;\text{ m}$, viscosity $\eta=10^{24} \; \text{ kg}/\text{ m}/\text{ s}$, density $\rho=3\times 10^{3} \; \text{ kg}/\text{ m}^3$ and a typical velocity of $U=0.1\;\text{ m}/\text{year}=3\times 10^{-9}\; \text{ m}/\text{ s}$, we get that the friction term in {math:numref}`eq:stokes-1` has size $\eta U/L^2 \approx 3\times 10^{2} \; \text{ kg}/\text{ m}^2/\text{ s}^2$.
On the other hand, the term $\nabla\cdot(\rho u)$ in the continuity equation {math:numref}`eq:stokes-2` has size $\rho U/L\approx 3\times 10^{-12} \; \text{ kg}/\text{ s}/\text{ m}^3$.
In other words, their *numerical values* are 14 orders of magnitude apart.

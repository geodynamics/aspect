(parameters:Discretization)=
# Discretization


## **Subsection:** Discretization


::::{dropdown} __Parameter:__ {ref}`Composition polynomial degree<parameters:Discretization/Composition_20polynomial_20degree>`
:name: parameters:Discretization/Composition_20polynomial_20degree
**Default value:** 2

**Pattern:** [List of <[Integer range 0...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The polynomial degree to use for the composition variable(s). As an example, a value of 2 for this parameter will yield either the element $Q_2$ or $DGQ_2$ for the compositional field(s), depending on whether we use continuous or discontinuous field(s).

For continuous elements, the value needs to be 1 or larger as $Q_1$ is the lowest order element, while $DGQ_0$ is a valid choice. Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`Stokes velocity polynomial degree<parameters:Discretization/Stokes_20velocity_20polynomial_20degree>`
:name: parameters:Discretization/Stokes_20velocity_20polynomial_20degree
**Default value:** 2

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The polynomial degree to use for the velocity variables in the Stokes system. The polynomial degree for the pressure variable will then be one less in order to make the velocity/pressure pair conform with the usual LBB (Babu{\v s}ka-Brezzi) condition. In other words, we are using a Taylor-Hood element for the Stokes equations and this parameter indicates the polynomial degree of it. As an example, a value of 2 for this parameter will yield the element $Q_2^d \times Q_1$ for the $d$ velocity components and the pressure, respectively (unless the &lsquo;Use locally conservative discretization&rsquo; parameter is set, which modifies the pressure element).

Be careful if you choose 1 as the degree. The resulting element is not stable and it may lead to artifacts in the solution. Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`Temperature polynomial degree<parameters:Discretization/Temperature_20polynomial_20degree>`
:name: parameters:Discretization/Temperature_20polynomial_20degree
**Default value:** 2

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The polynomial degree to use for the temperature variable. As an example, a value of 2 for this parameter will yield either the element $Q_2$ or $DGQ_2$ for the temperature field, depending on whether we use a continuous or discontinuous field. Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`Use discontinuous composition discretization<parameters:Discretization/Use_20discontinuous_20composition_20discretization>`
:name: parameters:Discretization/Use_20discontinuous_20composition_20discretization
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Whether to use a composition discretization that is discontinuous as opposed to continuous. This then requires the assembly of face terms between cells, and weak imposition of boundary terms for the composition field via the discontinuous Galerkin method.
::::

::::{dropdown} __Parameter:__ {ref}`Use discontinuous temperature discretization<parameters:Discretization/Use_20discontinuous_20temperature_20discretization>`
:name: parameters:Discretization/Use_20discontinuous_20temperature_20discretization
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a temperature discretization that is discontinuous as opposed to continuous. This then requires the assembly of face terms between cells, and weak imposition of boundary terms for the temperature field via the interior-penalty discontinuous Galerkin method.
::::

::::{dropdown} __Parameter:__ {ref}`Use equal order interpolation for Stokes<parameters:Discretization/Use_20equal_20order_20interpolation_20for_20Stokes>`
:name: parameters:Discretization/Use_20equal_20order_20interpolation_20for_20Stokes
**Default value:** false

**Pattern:** [Bool]

**Documentation:** By default (i.e., when this parameter is set to its default value &lsquo;false&rsquo;) ASPECT uses finite element combinations in which the pressure shape functions are polynomials one degree lower than the shape functions for the velocity. An example is the Taylor-Hood element that uses $Q_k$ elements for the velocity and $Q_{k-1}$ for the pressure. This is because using the *same* polynomial degree for both the velocity and the pressure turns out to violate some mathematical properties necessary to make the problem solvable. (In particular, thecondition in question goes by the name &ldquo;inf-sup&rdquo; or Babu{\v s}ka-Brezzi or LBB condition.) A consequence of violating this condition is that the pressure may show oscillations and not converge to the correct pressure.

That said, people have often used $Q_1$ elements for both the velocity and pressure anyway. This is commonly referred to as using the $Q_1-Q_1$ method. It is, by default, not stable as mentioned above, but it can be made stable by adding a small amount of compressibility to the model. There are numerous ways to do that. Today, the way that is generally considered to be the best approach is the one by Dohrmann and Bochev {cite}`DohrmannBochev2004`.

When this parameter is set to &ldquo;true&rdquo;, then ASPECT will use this method by using $Q_k\times Q_k$ elements for velocity and pressure, respectively, where $k$ is the value provided for the parameter &ldquo;Stokes velocity polynomial degree&rdquo;.

:::{note}
While ASPECT *allows* you to use this method, it is generally understood that this is not a great idea as it leads to rather low accuracy in general as documented in {cite}`thba22`. It also leads to substantial problems when using free surfaces. As a consequence, the presence of this parameter should not be seen as an endorsement of the method, or a suggestion to actually use it. It simply makes the method available.
:::
::::

::::{dropdown} __Parameter:__ {ref}`Use locally conservative discretization<parameters:Discretization/Use_20locally_20conservative_20discretization>`
:name: parameters:Discretization/Use_20locally_20conservative_20discretization
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a Stokes discretization that is locally conservative at the expense of a larger number of degrees of freedom (true), or to go with a cheaper discretization that does not locally conserve mass, although it is globally conservative (false).

When using a locally conservative discretization, the finite element space for the pressure is discontinuous between cells and is the polynomial space $P_{-(k-1)}$ of polynomials of degree $k-1$ in each variable separately. Here, $k$ is the value given in the parameter &ldquo;Stokes velocity polynomial degree&rdquo;, and consequently the polynomial degree for the pressure, $k-1$, is one lower than that for the velocity.

As a consequence of choosing this element for the pressure rather than the more commonly used $Q_{k-1}$ element that is continuous, it can be shown that if the medium is considered incompressible then the computed discrete velocity field $\mathbf u_h$ satisfies the property $\int_ {\partial K} \mathbf u_h \cdot \mathbf n = 0$ for every cell $K$, i.e., for each cell inflow and outflow exactly balance each other as one would expect for an incompressible medium. In other words, the velocity field is *locally conservative*.

On the other hand, if this parameter is set to &ldquo;false&rdquo;(the default), then the finite element space is chosen as $Q_{k-1}$. This choice does not yield the local conservation property but has the advantage of requiring fewer degrees of freedom. Furthermore, the error is generally smaller with this choice.

For an in-depth discussion of these issues and a quantitative evaluation of the different choices, see {cite}`kronbichler:etal:2012`.
::::

(parameters:Discretization/Stabilization_20parameters)=
## **Subsection:** Discretization / Stabilization parameters
::::{dropdown} __Parameter:__ {ref}`Discontinuous penalty<parameters:Discretization/Stabilization_20parameters/Discontinuous_20penalty>`
:name: parameters:Discretization/Stabilization_20parameters/Discontinuous_20penalty
**Default value:** 10.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value used to penalize discontinuities in the discontinuous Galerkin method. This is used only for the temperature field, and not for the composition field, as pure advection does not use the interior penalty method. This is largely empirically decided -- it must be large enough to ensure the bilinear form is coercive, but not so large as to penalize discontinuity at all costs.
::::

::::{dropdown} __Parameter:__ {ref}`Global composition maximum<parameters:Discretization/Stabilization_20parameters/Global_20composition_20maximum>`
:name: parameters:Discretization/Stabilization_20parameters/Global_20composition_20maximum
**Default value:** 1.7976931348623157e+308

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The maximum global composition values that will be used in the bound preserving limiter for the discontinuous solutions from composition advection fields. The number of the input &rsquo;Global composition maximum&rsquo; values separated by &rsquo;,&rsquo; has to be one or the same as the number of the compositional fields. When only one value is supplied, this same value is assumed for all compositional fields.
::::

::::{dropdown} __Parameter:__ {ref}`Global composition minimum<parameters:Discretization/Stabilization_20parameters/Global_20composition_20minimum>`
:name: parameters:Discretization/Stabilization_20parameters/Global_20composition_20minimum
**Default value:** -1.7976931348623157e+308

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The minimum global composition value that will be used in the bound preserving limiter for the discontinuous solutions from composition advection fields. The number of the input &rsquo;Global composition minimum&rsquo; values separated by &rsquo;,&rsquo; has to be one or the same as the number of the compositional fields. When only one value is supplied, this same value is assumed for all compositional fields.
::::

::::{dropdown} __Parameter:__ {ref}`Global temperature maximum<parameters:Discretization/Stabilization_20parameters/Global_20temperature_20maximum>`
:name: parameters:Discretization/Stabilization_20parameters/Global_20temperature_20maximum
**Default value:** 1.7976931348623157e+308

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum global temperature value that will be used in the bound preserving limiter for the discontinuous solutions from temperature advection fields.
::::

::::{dropdown} __Parameter:__ {ref}`Global temperature minimum<parameters:Discretization/Stabilization_20parameters/Global_20temperature_20minimum>`
:name: parameters:Discretization/Stabilization_20parameters/Global_20temperature_20minimum
**Default value:** -1.7976931348623157e+308

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The minimum global temperature value that will be used in the bound preserving limiter for the discontinuous solutions from temperature advection fields.
::::

::::{dropdown} __Parameter:__ {ref}`List of compositional fields with disabled boundary entropy viscosity<parameters:Discretization/Stabilization_20parameters/List_20of_20compositional_20fields_20with_20disabled_20boundary_20entropy_20viscosity>`
:name: parameters:Discretization/Stabilization_20parameters/List_20of_20compositional_20fields_20with_20disabled_20boundary_20entropy_20viscosity
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** Select for which compositional fields to skip the entropy viscosity stabilization at dirichlet boundaries. This is only advisable for compositional fields that have intrinsic physical diffusion terms, otherwise oscillations may develop. The parameter should contain a list of compositional field names.
::::

::::{dropdown} __Parameter:__ {ref}`Stabilization method<parameters:Discretization/Stabilization_20parameters/Stabilization_20method>`
:name: parameters:Discretization/Stabilization_20parameters/Stabilization_20method
**Default value:** entropy viscosity

**Pattern:** [Selection entropy viscosity|SUPG ]

**Documentation:** Select the method for stabilizing the advection equation. The original method implemented is &rsquo;entropy viscosity&rsquo; as described in {cite}`kronbichler:etal:2012`. SUPG is currently experimental.
::::

::::{dropdown} __Parameter:__ {ref}`Use artificial viscosity smoothing<parameters:Discretization/Stabilization_20parameters/Use_20artificial_20viscosity_20smoothing>`
:name: parameters:Discretization/Stabilization_20parameters/Use_20artificial_20viscosity_20smoothing
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If set to false, the artificial viscosity of a cell is computed and is computed on every cell separately as discussed in {cite}`kronbichler:etal:2012`. If set to true, the maximum of the artificial viscosity in the cell as well as the neighbors of the cell is computed and used instead.
::::

::::{dropdown} __Parameter:__ {ref}`Use limiter for discontinuous composition solution<parameters:Discretization/Stabilization_20parameters/Use_20limiter_20for_20discontinuous_20composition_20solution>`
:name: parameters:Discretization/Stabilization_20parameters/Use_20limiter_20for_20discontinuous_20composition_20solution
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Whether to apply the bound preserving limiter as a correction after having the discontinuous composition solution. The limiter will only have an effect if the &rsquo;Global composition maximum&rsquo; and &rsquo;Global composition minimum&rsquo; parameters are defined in the .prm file. This limiter keeps the discontinuous solution in the range given by Global composition maximum&rsquo; and &rsquo;Global composition minimum&rsquo;. The number of input values in this parameter separated by &rsquo;,&rsquo; has to be one or the number of the compositional fields. When only one value is supplied, this same value is assumed for all compositional fields, otherwise each value represents if the limiter should be applied to the respective compositional field. Because this limiter modifies the solution it no longer satisfies the assembled equation. Therefore, the nonlinear residual for this field is meaningless, and in nonlinear solvers we will ignore the residual for this field to evaluate if the nonlinear solver has converged.
::::

::::{dropdown} __Parameter:__ {ref}`Use limiter for discontinuous temperature solution<parameters:Discretization/Stabilization_20parameters/Use_20limiter_20for_20discontinuous_20temperature_20solution>`
:name: parameters:Discretization/Stabilization_20parameters/Use_20limiter_20for_20discontinuous_20temperature_20solution
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to apply the bound preserving limiter as a correction after computing the discontinuous temperature solution. The limiter will only have an effect if the &rsquo;Global temperature maximum&rsquo; and &rsquo;Global temperature minimum&rsquo; parameters are defined in the .prm file. This limiter keeps the discontinuous solution in the range given by &rsquo;Global temperature maximum&rsquo; and &rsquo;Global temperature minimum&rsquo;. Because this limiter modifies the solution it no longer satisfies the assembled equation. Therefore, the nonlinear residual for this field is meaningless, and in nonlinear solvers we will ignore the residual for this field to evaluate if the nonlinear solver has converged.
::::

::::{dropdown} __Parameter:__ {ref}`alpha<parameters:Discretization/Stabilization_20parameters/alpha>`
:name: parameters:Discretization/Stabilization_20parameters/alpha
**Default value:** 2

**Pattern:** [Integer range 1...2 (inclusive)]

**Documentation:** The exponent $\alpha$ in the entropy viscosity stabilization. Valid options are 1 or 2. The recommended setting is 2. (This parameter does not correspond to any variable in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see {cite}`kronbichler:etal:2012`. Rather, the paper always uses 2 as the exponent in the definition of the entropy, following equation (15) of the paper. The full approach is discussed in {cite}`guermond:etal:2011`.) Note that this is not the thermal expansion coefficient, also commonly referred to as $\alpha$.Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`beta<parameters:Discretization/Stabilization_20parameters/beta>`
:name: parameters:Discretization/Stabilization_20parameters/beta
**Default value:** 0.052

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The $\beta$ factor in the artificial viscosity stabilization. This parameter controls the maximum dissipation of the entropy viscosity, which is the part that only scales with the cell diameter and the maximum velocity in the cell, but does not depend on the solution field itself or its residual. An appropriate value for 2d is 0.052 and 0.078 for 3d. (For historical reasons, the name used here is different from the one used in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see {cite}`kronbichler:etal:2012`. This parameter can be given as a single value or as a list with as many entries as one plus the number of compositional fields. In the former case all advection fields use the same stabilization parameters, in the latter case each field (temperature first, then all compositions) use individual parameters. This can be useful to reduce the stabilization for the temperature, which already has some physical diffusion. This parameter corresponds to the factor $\alpha_{\text{max}}$ in the formulas following equation (15) of the paper.) Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`cR<parameters:Discretization/Stabilization_20parameters/cR>`
:name: parameters:Discretization/Stabilization_20parameters/cR
**Default value:** 0.11

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The $c_R$ factor in the entropy viscosity stabilization. This parameter controls the part of the entropy viscosity that depends on the solution field itself and its residual in addition to the cell diameter and the maximum velocity in the cell. This parameter can be given as a single value or as a list with as many entries as one plus the number of compositional fields. In the former case all advection fields use the same stabilization parameters, in the latter case each field (temperature first, then all compositions) use individual parameters. This can be useful to reduce the stabilization for the temperature, which already has some physical diffusion. (For historical reasons, the name used here is different from the one used in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see {cite}`kronbichler:etal:2012`. This parameter corresponds to the factor $\alpha_E$ in the formulas following equation (15) of the paper.) Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`gamma<parameters:Discretization/Stabilization_20parameters/gamma>`
:name: parameters:Discretization/Stabilization_20parameters/gamma
**Default value:** 0.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The strain rate scaling factor in the artificial viscosity stabilization. This parameter determines how much the strain rate (in addition to the velocity) should influence the stabilization. (This parameter does not correspond to any variable in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see {cite}`kronbichler:etal:2012`. Rather, the paper always uses 0, i.e. they specify the maximum dissipation $\nu_h^\text{max}$ as $\nu_h^\text{max}\vert_K = \alpha_{\text{max}} h_K \|\mathbf u\|_{\infty,K}$. Here, we use $\|\lvert\mathbf u\rvert + \gamma h_K \lvert\varepsilon (\mathbf u)\rvert\|_{\infty,K}$ instead of $\|\mathbf u\|_{\infty,K}$. Units: None.
::::

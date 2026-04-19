(parameters:Heating_20model)=
# Heating model


## **Subsection:** Heating model


::::{dropdown} __Parameter:__ {ref}`List of model names<parameters:Heating_20model/List_20of_20model_20names>`
:name: parameters:Heating_20model/List_20of_20model_20names
**Default value:**

**Pattern:** [MultipleSelection adiabatic heating|adiabatic heating of melt|anisotropic shear heating|compositional heating|constant heating|function|latent heat|latent heat melt|radioactive decay|shear heating|shear heating with melt|tidal heating ]

**Documentation:** A comma separated list of heating models that will be used to calculate the heating terms in the energy equation. The results of each of these criteria, i.e., the heating source terms and the latent heat terms for the left hand side will be added.

The following heating models are available:

&lsquo;adiabatic heating&rsquo;: Implementation of a standard and a simplified model of adiabatic heating.

&lsquo;adiabatic heating of melt&rsquo;: Implementation of a standard and a simplified model of adiabatic heating of melt. The full model implements the heating term
$\alpha T (-\phi \mathbf u_s \cdot \nabla p) + \alpha T (\phi \mathbf u_f \cdot \nabla p)$.
For full adiabatic heating, this has to be used in combination with the heating model &lsquo;adiabatic heating&rsquo; to also include adiabatic heating for the solid part, and the full heating term is then $\alpha T ((1-\phi) \mathbf u_s \cdot \nabla p) + \alpha T (\phi \mathbf u_f \cdot \nabla p)$.

&lsquo;anisotropic shear heating&rsquo;: Implementation of a standard model for shear heating extended for an anisotropic viscosity tensor. If the material model provides a stress-strain director tensor, then the strain-rate is multiplied with this tensor to compute the stress that is used when computing the shear heating.

&lsquo;compositional heating&rsquo;: Implementation of a model in which magnitude of internal heat production is determined from fixed values assigned to each compositional field. These values are interpreted as having units \si{\watt\per\meter\cubed}.

&lsquo;constant heating&rsquo;: Implementation of a model in which the heating rate is constant.

&lsquo;function&rsquo;: Implementation of a model in which the heating rate is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Heating model|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

The formula is interpreted as having units W/kg.

Since the symbol $t$ indicating time may appear in the formulas for the heating rate, it is interpreted as having units seconds unless the global parameter &ldquo;Use years instead of seconds&rdquo; is set.

&lsquo;latent heat&rsquo;: Implementation of a standard model for latent heat.

&lsquo;latent heat melt&rsquo;: Implementation of a standard model for latent heat of melting. This assumes that there is a compositional field called porosity, and it uses the reaction term of this field (the fraction of material that melted in the current time step) multiplied by a constant entropy change for melting all of the material as source term of the heating model.
If there is no field called porosity, the heating terms are 0.

&lsquo;radioactive decay&rsquo;: Implementation of a model in which the internal heating rate is radioactive decaying in the following rule:
\[(\text{initial concentration})\cdot 0.5^{\text{time}/(\text{half life})}\]
The crust and mantle can have different concentrations, and the crust can be defined either by depth or by a certain compositional field.
The formula is interpreted as having units W/kg.

&lsquo;shear heating&rsquo;: Implementation of a standard model for shear heating. Adds the term: $  2 \eta \left( \varepsilon - \frac{1}{3} \text{tr} \varepsilon \mathbf 1 \right) : \left( \varepsilon - \frac{1}{3} \text{tr} \varepsilon \mathbf 1 \right)$ to the right-hand side of the temperature equation.

&lsquo;shear heating with melt&rsquo;: Implementation of a standard model for shear heating of migrating melt, including bulk (compression) heating $\xi \left( \nabla \cdot \mathbf u_s \right)^2 $ and heating due to melt segregation $\frac{\eta_f \phi^2}{k} \left( \mathbf u_f - \mathbf u_s \right)^2 $. For full shear heating, this has to be used in combination with the heating model shear heating to also include shear heating for the solid part.

&lsquo;tidal heating&rsquo;: A tidal heating implementation related to diurnal tides. The default equation ignores regional (radial/lateral) changes. This equation is the Eq.12 from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099). Selecting &rsquo;latitudinal variation&rsquo; from &rsquo;Custom distribution of tidal strain rate&rsquo; allows simplified latitudinal variation with cosine function, &rsquo;Maximum tidal strain rate&rsquo; and &rsquo;Minimum tidal strain rate. Latitudinal variation of tidal strain rate is shown in Fig.3 from Nimmo et al. (2007) (https://doi.org/10.1016/j.icarus.2007.04.021). Unit: W/m^3.
::::

(parameters:Heating_20model/Adiabatic_20heating)=
## **Subsection:** Heating model / Adiabatic heating
::::{dropdown} __Parameter:__ {ref}`Use simplified adiabatic heating<parameters:Heating_20model/Adiabatic_20heating/Use_20simplified_20adiabatic_20heating>`
:name: parameters:Heating_20model/Adiabatic_20heating/Use_20simplified_20adiabatic_20heating
**Default value:** false

**Pattern:** [Bool]

**Documentation:** A flag indicating whether the adiabatic heating should be simplified from $\alpha T (\mathbf u \cdot \nabla p)$ to $ \alpha \rho T (\mathbf u \cdot \mathbf g) $.
::::

(parameters:Heating_20model/Adiabatic_20heating_20of_20melt)=
## **Subsection:** Heating model / Adiabatic heating of melt
::::{dropdown} __Parameter:__ {ref}`Use simplified adiabatic heating<parameters:Heating_20model/Adiabatic_20heating_20of_20melt/Use_20simplified_20adiabatic_20heating>`
:name: parameters:Heating_20model/Adiabatic_20heating_20of_20melt/Use_20simplified_20adiabatic_20heating
**Default value:** false

**Pattern:** [Bool]

**Documentation:** A flag indicating whether the adiabatic heating should be simplified from $\alpha T (\mathbf u \cdot \nabla p)$ to $ \alpha \rho T (\mathbf u \cdot \mathbf g) $.
::::

(parameters:Heating_20model/Compositional_20heating)=
## **Subsection:** Heating model / Compositional heating
::::{dropdown} __Parameter:__ {ref}`Compositional heating values<parameters:Heating_20model/Compositional_20heating/Compositional_20heating_20values>`
:name: parameters:Heating_20model/Compositional_20heating/Compositional_20heating_20values
**Default value:** 0.

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** List of heat production per unit volume values for background and compositional fields, for a total of N+1 values, where the first value corresponds to the background material, and N is the number of compositional fields. Units: \si{\watt\per\meter\cubed}.
::::

::::{dropdown} __Parameter:__ {ref}`Use compositional field for heat production averaging<parameters:Heating_20model/Compositional_20heating/Use_20compositional_20field_20for_20heat_20production_20averaging>`
:name: parameters:Heating_20model/Compositional_20heating/Use_20compositional_20field_20for_20heat_20production_20averaging
**Default value:** 1

**Pattern:** [List of <[Integer range 0...1 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of integers with as many entries as compositional fields plus one. The first entry corresponds to the background material, each following entry corresponds to a particular compositional field. If the entry for a field is &rsquo;1&rsquo; this field is considered during the computation of volume fractions, if it is &rsquo;0&rsquo; the field is ignored. This is useful if some compositional fields are used to track properties like finite strain that should not contribute to heat production. The first entry determines whether the background field contributes to heat production or not (essentially similar to setting its &rsquo;Compositional heating values&rsquo; to zero, but included for consistency in the length of the input lists).
::::

(parameters:Heating_20model/Constant_20heating)=
## **Subsection:** Heating model / Constant heating
::::{dropdown} __Parameter:__ {ref}`Radiogenic heating rate<parameters:Heating_20model/Constant_20heating/Radiogenic_20heating_20rate>`
:name: parameters:Heating_20model/Constant_20heating/Radiogenic_20heating_20rate
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The specific rate of heating due to radioactive decay (or other bulk sources you may want to describe). This parameter corresponds to the variable $H$ in the temperature equation stated in the manual, and the heating term is $\rho H$. Units: W/kg.
::::

(parameters:Heating_20model/Function)=
## **Subsection:** Heating model / Function
::::{dropdown} __Parameter:__ {ref}`Coordinate system<parameters:Heating_20model/Function/Coordinate_20system>`
:name: parameters:Heating_20model/Function/Coordinate_20system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.
::::

::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Heating_20model/Function/Function_20constants>`
:name: parameters:Heating_20model/Function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Heating_20model/Function/Function_20expression>`
:name: parameters:Heating_20model/Function/Function_20expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Heating_20model/Function/Variable_20names>`
:name: parameters:Heating_20model/Function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Heating_20model/Latent_20heat_20melt)=
## **Subsection:** Heating model / Latent heat melt
::::{dropdown} __Parameter:__ {ref}`Melting entropy change<parameters:Heating_20model/Latent_20heat_20melt/Melting_20entropy_20change>`
:name: parameters:Heating_20model/Latent_20heat_20melt/Melting_20entropy_20change
**Default value:** -300.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The entropy change for the phase transition from solid to melt. Units: $\frac{\text{J}}{\text{K}\text{kg}}$.
::::

::::{dropdown} __Parameter:__ {ref}`Retrieve entropy change from material model<parameters:Heating_20model/Latent_20heat_20melt/Retrieve_20entropy_20change_20from_20material_20model>`
:name: parameters:Heating_20model/Latent_20heat_20melt/Retrieve_20entropy_20change_20from_20material_20model
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Instead of using the entropy change given in the &rsquo;Melting entropy change&rsquo; query the EnthalpyAdditionalOutputs in the material model to compute the entropy change for the phase transition from solid to melt.Units: $J/(kg K)$.
::::

(parameters:Heating_20model/Radioactive_20decay)=
## **Subsection:** Heating model / Radioactive decay
::::{dropdown} __Parameter:__ {ref}`Crust composition number<parameters:Heating_20model/Radioactive_20decay/Crust_20composition_20number>`
:name: parameters:Heating_20model/Radioactive_20decay/Crust_20composition_20number
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Which composition field should be treated as crust
::::

::::{dropdown} __Parameter:__ {ref}`Crust defined by composition<parameters:Heating_20model/Radioactive_20decay/Crust_20defined_20by_20composition>`
:name: parameters:Heating_20model/Radioactive_20decay/Crust_20defined_20by_20composition
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether crust defined by composition or depth
::::

::::{dropdown} __Parameter:__ {ref}`Crust depth<parameters:Heating_20model/Radioactive_20decay/Crust_20depth>`
:name: parameters:Heating_20model/Radioactive_20decay/Crust_20depth
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Depth of the crust when crust if defined by depth. Units: \si{\meter}.
::::

::::{dropdown} __Parameter:__ {ref}`Half decay times<parameters:Heating_20model/Radioactive_20decay/Half_20decay_20times>`
:name: parameters:Heating_20model/Radioactive_20decay/Half_20decay_20times
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Half decay times. Units: (Seconds), or (Years) if set &lsquo;use years instead of seconds&rsquo;.
::::

::::{dropdown} __Parameter:__ {ref}`Heating rates<parameters:Heating_20model/Radioactive_20decay/Heating_20rates>`
:name: parameters:Heating_20model/Radioactive_20decay/Heating_20rates
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Heating rates of different elements (W/kg)
::::

::::{dropdown} __Parameter:__ {ref}`Initial concentrations crust<parameters:Heating_20model/Radioactive_20decay/Initial_20concentrations_20crust>`
:name: parameters:Heating_20model/Radioactive_20decay/Initial_20concentrations_20crust
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Initial concentrations of different elements (ppm)
::::

::::{dropdown} __Parameter:__ {ref}`Initial concentrations mantle<parameters:Heating_20model/Radioactive_20decay/Initial_20concentrations_20mantle>`
:name: parameters:Heating_20model/Radioactive_20decay/Initial_20concentrations_20mantle
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Initial concentrations of different elements (ppm)
::::

::::{dropdown} __Parameter:__ {ref}`Number of elements<parameters:Heating_20model/Radioactive_20decay/Number_20of_20elements>`
:name: parameters:Heating_20model/Radioactive_20decay/Number_20of_20elements
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Number of radioactive elements
::::

(parameters:Heating_20model/Shear_20heating)=
## **Subsection:** Heating model / Shear heating
::::{dropdown} __Parameter:__ {ref}`Cohesion for maximum shear stress<parameters:Heating_20model/Shear_20heating/Cohesion_20for_20maximum_20shear_20stress>`
:name: parameters:Heating_20model/Shear_20heating/Cohesion_20for_20maximum_20shear_20stress
**Default value:** 2e7

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Cohesion for maximum shear stress that should be used for the computation of shear heating. It can be useful to limit the shear stress in models where velocities are prescribed, and actual stresses in the Earth would be lower than the stresses introduced by the boundary conditions. Only used if &rsquo;Limit stress contribution to shear heating&rsquo; is true. Units: Pa.
::::

::::{dropdown} __Parameter:__ {ref}`Friction angle for maximum shear stress<parameters:Heating_20model/Shear_20heating/Friction_20angle_20for_20maximum_20shear_20stress>`
:name: parameters:Heating_20model/Shear_20heating/Friction_20angle_20for_20maximum_20shear_20stress
**Default value:** 0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Friction angle for maximum shear stress that should be used for the computation of shear heating. It can be useful to limit the shear stress in models where velocities are prescribed, and actual stresses in the Earth would be lower than the stresses introduced by the boundary conditions. Only used if &rsquo;Limit stress contribution to shear heating&rsquo; is true. Units: none.
::::

::::{dropdown} __Parameter:__ {ref}`Limit stress contribution to shear heating<parameters:Heating_20model/Shear_20heating/Limit_20stress_20contribution_20to_20shear_20heating>`
:name: parameters:Heating_20model/Shear_20heating/Limit_20stress_20contribution_20to_20shear_20heating
**Default value:** false

**Pattern:** [Bool]

**Documentation:** In models with prescribed boundary velocities, stresses can become unrealistically large. Using these large stresses when calculating the amount of shear heating would then lead to an unreasonable increase in temperature. This parameter indicates if the stress being used to compute the amount of shear heating should be limited based on a Drucker-Prager yield criterion with the cohesion given by the &rsquo;Cohesion for maximum shear stress&rsquo; parameter and the friction angle given by the &rsquo;Friction angle for maximum shear stress&rsquo; parameter.
::::

(parameters:Heating_20model/Tidal_20heating)=
## **Subsection:** Heating model / Tidal heating
::::{dropdown} __Parameter:__ {ref}`Constant tidal strain rate<parameters:Heating_20model/Tidal_20heating/Constant_20tidal_20strain_20rate>`
:name: parameters:Heating_20model/Tidal_20heating/Constant_20tidal_20strain_20rate
**Default value:** 2e-10

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Time-averaged strain rate by a lunar diurnal tide that is simplified to be constant regardless of location. Default value is for the diurnal tide of Europa from Tobie et al. (2003). Units: 1/s.
::::

::::{dropdown} __Parameter:__ {ref}`Custom distribution of tidal strain rate<parameters:Heating_20model/Tidal_20heating/Custom_20distribution_20of_20tidal_20strain_20rate>`
:name: parameters:Heating_20model/Tidal_20heating/Custom_20distribution_20of_20tidal_20strain_20rate
**Default value:** constant

**Pattern:** [Selection constant|latitudinal variation ]

**Documentation:** Choose how the time-averaged tidal strain rate is distributed. If &rsquo;constant&rsquo;, the tidal strain rate is fixed to &rsquo;Constant tidal strain rate&rsquo;. If &rsquo;latitudinal variation&rsquo;, &rsquo;Maximum tidal strain rate&rsquo; and &rsquo;Minimum tidal strain rate&rsquo; are used.
::::

::::{dropdown} __Parameter:__ {ref}`Elastic shear modulus<parameters:Heating_20model/Tidal_20heating/Elastic_20shear_20modulus>`
:name: parameters:Heating_20model/Tidal_20heating/Elastic_20shear_20modulus
**Default value:** 3.3e9

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Elastic shear modulus of the material. For simplicity, this parameter will be used even if elasticity is set in the material model. Default value is for Europa&rsquo;s icy shell. Units: Pa.
::::

::::{dropdown} __Parameter:__ {ref}`Maximum tidal strain rate<parameters:Heating_20model/Tidal_20heating/Maximum_20tidal_20strain_20rate>`
:name: parameters:Heating_20model/Tidal_20heating/Maximum_20tidal_20strain_20rate
**Default value:** 2.81e-10

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum time-averaged tidal strain rate by lunar diurnal tide at the poles. This parameter will be used when &rsquo;Custom distribution of tidal strain rate&rsquo; is &rsquo;latitudinal variation&rsquo;. Default value is for Europa at pole from Nimmo et al. (2007). Units: 1/s.
::::

::::{dropdown} __Parameter:__ {ref}`Minimum tidal strain rate<parameters:Heating_20model/Tidal_20heating/Minimum_20tidal_20strain_20rate>`
:name: parameters:Heating_20model/Tidal_20heating/Minimum_20tidal_20strain_20rate
**Default value:** 1.67e-10

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum time-averaged tidal strain rate by lunar diurnal tide at the equator. This parameter will be used when &rsquo;Custom distribution of tidal strain rate&rsquo; is &rsquo;latitudinal variation&rsquo;. Default value is for Europa at Equator from Nimmo et al. (2007). Units: 1/s.
::::

::::{dropdown} __Parameter:__ {ref}`Tidal frequency<parameters:Heating_20model/Tidal_20heating/Tidal_20frequency>`
:name: parameters:Heating_20model/Tidal_20heating/Tidal_20frequency
**Default value:** 2.048e-5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The orbital/tidal frequency that produces the heating. Default value is the diurnal tidal frequency of Europa, ~3.551 days. Units: 1/s.
::::

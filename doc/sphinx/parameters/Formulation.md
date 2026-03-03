(parameters:Formulation)=
# Formulation


## **Subsection:** Formulation


::::{dropdown} __Parameter:__ {ref}`Enable additional Stokes RHS<parameters:Formulation/Enable_20additional_20Stokes_20RHS>`
:name: parameters:Formulation/Enable_20additional_20Stokes_20RHS
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to ask the material model for additional terms for the right-hand side of the Stokes equation. This feature is likely only used when implementing force vectors for manufactured solution problems and requires filling additional outputs of type AdditionalMaterialOutputsStokesRHS.
::::

::::{dropdown} __Parameter:__ {ref}`Enable elasticity<parameters:Formulation/Enable_20elasticity>`
:name: parameters:Formulation/Enable_20elasticity
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to include the additional elastic terms on the right-hand side of the Stokes equation.
::::

::::{dropdown} __Parameter:__ {ref}`Enable prescribed dilation<parameters:Formulation/Enable_20prescribed_20dilation>`
:name: parameters:Formulation/Enable_20prescribed_20dilation
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to include additional terms on the right-hand side of the Stokes equation to set a given compression term specified in the MaterialModel output PrescribedPlasticDilation.
::::

::::{dropdown} __Parameter:__ {ref}`Formulation<parameters:Formulation/Formulation>`
:name: parameters:Formulation/Formulation
**Default value:** custom

**Pattern:** [Selection isentropic compression|custom|anelastic liquid approximation|Boussinesq approximation ]

**Documentation:** Select a formulation for the basic equations. Different published formulations are available in ASPECT (see the list of possible values for this parameter in the manual for available options). Two ASPECT specific options are
  * &lsquo;isentropic compression&rsquo;: ASPECT&rsquo;s original formulation, using the explicit compressible mass equation, and the full density for the temperature equation.
  * &lsquo;custom&rsquo;: A custom selection of &lsquo;Mass conservation&rsquo; and &lsquo;Temperature equation&rsquo;.
:::{warning}
The &lsquo;custom&rsquo; option is implemented for advanced users that want full control over the equations solved. It is possible to choose inconsistent formulations and no error checking is performed on the consistency of the resulting equations.
:::

:::{note}
The &lsquo;anelastic liquid approximation&rsquo; option here can also be used to set up the &lsquo;truncated anelastic liquid approximation&rsquo; as long as this option is chosen together with a material model that defines a density that depends on temperature and depth and not on the pressure.
:::
::::

::::{dropdown} __Parameter:__ {ref}`Mass conservation<parameters:Formulation/Mass_20conservation>`
:name: parameters:Formulation/Mass_20conservation
**Default value:** ask material model

**Pattern:** [Selection incompressible|isentropic compression|hydrostatic compression|reference density profile|implicit reference density profile|projected density field|ask material model ]

**Documentation:** Possible approximations for the density derivatives in the mass conservation equation. Note that this parameter is only evaluated if &lsquo;Formulation&rsquo; is set to &lsquo;custom&rsquo;. Other formulations ignore the value of this parameter.
::::

::::{dropdown} __Parameter:__ {ref}`Temperature equation<parameters:Formulation/Temperature_20equation>`
:name: parameters:Formulation/Temperature_20equation
**Default value:** real density

**Pattern:** [Selection real density|reference density profile ]

**Documentation:** Possible approximations for the density in the temperature equation. Possible approximations are &lsquo;real density&rsquo; and &lsquo;reference density profile&rsquo;. Note that this parameter is only evaluated if &lsquo;Formulation&rsquo; is set to &lsquo;custom&rsquo;. Other formulations ignore the value of this parameter.
::::

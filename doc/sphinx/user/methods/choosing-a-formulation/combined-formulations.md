(sec:methods:combined-formulations)=
# Combined formulations

Not all combinations of the different approximations discussed above are physically reasonable, and to help users choose between these options, we provide a number of combined "Formulations" that are equivalent to the approximate equations discussed previously ({ref}`sec:methods:approximate-equations`).
They can be selected in the subsection `Formulation/Formulation` (see also {ref}`parameters:Formulation/Formulation`):

-   "anelastic liquid approximation": This formulation sets the mass conservation approximation to "reference density profile," the temperature equation approximation to "reference density profile" and checks that both adiabatic and shear heating are included in the list of heating plugins used in the model, using the simplified version of the adiabatic heating term (see {ref}`parameters:Heating_20model/Adiabatic_20heating`).
The default setting for the adiabatic conditions is an adiabatic temperature profile, and hydrostatic pressure and density profiles.
This option should be chosen together with a material model that defines a density that depends on temperature and pressure (and potentially depth), which would be equivalent to the anelastic liquid approximation ({ref}`sec:methods:approximate-equations:ala`), or with a material model that defines a density that depends on temperature and depth (and not on the pressure), which would be equivalent to the truncated anelastic liquid approximation ({ref}`sec:methods:approximate-equations:tala`).

-   "Boussinesq approximation": This formulation sets the mass conservation approximation to "incompressible," the temperature equation approximation to "reference density profile" and checks that neither adiabatic nor shear heating are included in the list of heating plugins used in the model.
    The default setting for the adiabatic conditions is a constant temperature, and hydrostatic pressure and density profiles.
    This option should be chosen together with a material model that defines a density that only depends on temperature and depth (and not on the pressure).
    This is equivalent to the Boussinesq approximation ({ref}`sec:methods:approximate-equations:ba`).

-   "isothermal compression": This formulation sets the mass conservation approximation to "isothermal compression," the temperature equation approximation to "real density" and checks that both adiabatic and shear heating are included in the list of heating plugins used in the model.
    The default setting for the adiabatic conditions is an adiabatic temperature profile, and hydrostatic pressure and density profiles.
    The density can depend on any of the solution variables.
    This is equivalent to the isothermal compression approximation ({ref}`sec:methods:approximate-equations:ica`).

-   "custom": By default, this formulation sets the mass conservation approximation to "ask material model" and the temperature equation approximation to "real density."
    The adiabatic conditions model uses an adiabatic temperature profile if adiabatic heating is included in the model, and a constant temperature if adiabatic heating is not included. Pressure and density profiles are hydrostatic.
    The density can depend on any of the solution variables.
    However, this option can also be used to arbitrarily combine the different approximations described in this section.
    Users should be careful when using this option, as some combinations may lead to unphysical model behavior.

An example cookbook that shows a comparison between different approximations
is discussed in {ref}`sec:cookbooks:burnman`.

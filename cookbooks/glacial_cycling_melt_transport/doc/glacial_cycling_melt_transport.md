(sec:cookbooks:glacial-cycling-melt-trasport)=
# Glacial Unloading and Melt Migration with Viscoplastic Rheologies

*This section was contributed by Prajakta Mohite & John Naliboff.*

This cookbook demonstrates interactions between time and spatially-dependent surface loading, lithospheric deformation, and reactive melt transport. It is motivated by observations and prior modeling investigations highlighting the role of glacial loading and unloading on volcanic activity, and designed to serve as a potential template for investigating the effects of glacial cycles on tectonic and magmatic processes, such as in Antarctica and Iceland.
For example, the model in this cookbook can be used or extended to investigate the following topics:

1. Quantify how changes in surface loading drive stress evolution in the crust and upper mantle.

2. Examine the onset, distribution, and temporal evolution of melt during and after rapid unloading.

3. Compare the relative effects of different rheological laws and material properties on crustal and mantle dynamics.

## The input file

The model in this cookbook is designed to investigate the mechanical and magmatic response of the continental lithosphere and asthenosphere to rapid deglaciation, using a physically motivated and time-dependent parabolic ice load. The model  domain is 400 (X) x 200 (Y) km wide, and contains layers representing the continental crust, mantle lithosphere, and asthenosphere. The surface load is imposed as a 200 km wide parabolic ice sheet with a maximum height of 1 km, which is centered along the horizontal extent of the model and applied through a traction boundary condition. The load remains constant for a prescribed period of 1 Myr, after which it decreases linearly to zero over an interval of 20 kyr, thereby simulating rapid deglaciation.

This is defined in the parameter file as:
```{literalinclude} traction_boundary.part.prm
```
where ```{math}/rho``` , ``` /g```, and thickness represent ice density, gravitational acceleration, and ice thickness, respectively, and the spatial function defines the parabolic load geometry. Respectively, the model side and bottom boundaries contain no-slip and free-slip velocity boundary conditions.

The highly nonlinear system of equations is solved with the nonlinear solver scheme 'iterated Advection and Stokes', following prior investigations and cookbooks {ref}`sec:cookbooks:global-melt` modeling two-phase reactive melt transport in ASPECT {cite:t}`dannberg:heister:2016`. Respectively, the linear and nonlinear solver tolerance are set to, although we note that the model dynamics may slightly or moderately differ when using stricter values. Here, the values are selected as a  compromise between accuracy and model run times. We note that the introduction of a free surface, as compared to a free-slip top boundary, introduces significant challenges for both the linear and nonlinear solvers. Likewise, we note that using an open traction boundary at the model base with the traction magnitude equal to the adiabatic pressure may also introduce significant nonlinear solver issues.

The initial temperature field is defined using a depth-dependent conductive temperature profile through the lithosphere in a similar fashion to the continental extension cookbook, and an adiabatic gradient of 0.5 C/km from the base of the lithosphere (80 km depth) to the model base (200 km depth). Although temperatures are specified for the model sides, these values are not used, as the boundaries are insulating (zero net heat flux). Temperature evolves in the model domain through flow induced by spatiotemporal variations in surface loading, and includes the effects of adiabatic and shear heating introduced through the use of the extended Boussinesq approximation [cookbooks/continental_extension/](https://github.com/geodynamics/aspect/tree/main/cookbooks/continental_extension).
The model rheology follows nonlinear dislocation creep similar to the setup in the continental extension cookbook, with distinct flow parameters for the upper crust, lower crust, and mantle. Melt generation is parameterized using the {cite:t}`katz:etal:2003` model, with a prescribed melt freezing rate of 0.5 and maximum melt migration depth of 30 km.

The configuration below allows us to combine two-phase flow transport (Katz 2003) with viscoplastic rheologies[reaction-model:katz2003_mantle_melting.cc](https://github.com/geodynamics/aspect/blob/main/source/material_model/reaction_model/katz2003_mantle_melting.cc).

This is defined in the parameter file as:
```{literalinclude} reactive_fluid_transport_model.part.prm
```

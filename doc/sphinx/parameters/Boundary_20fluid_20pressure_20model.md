(parameters:Boundary_20fluid_20pressure_20model)=
# Boundary fluid pressure model


## **Subsection:** Boundary fluid pressure model


(parameters:Boundary_20fluid_20pressure_20model/Plugin_20name)=
### __Parameter name:__ Plugin name
**Default value:** density

**Pattern:** [Selection density ]

**Documentation:** Select one of the following plugins:

&lsquo;density&rsquo;: A plugin that prescribes the fluid pressure gradient at the boundary based on fluid/solid density from the material model.

(parameters:Boundary_20fluid_20pressure_20model/Density)=
## **Subsection:** Boundary fluid pressure model / Density
(parameters:Boundary_20fluid_20pressure_20model/Density/Density_20formulation)=
### __Parameter name:__ Density formulation
**Default value:** solid density

**Pattern:** [Selection solid density|fluid density|average density ]

**Documentation:** The density formulation used to compute the fluid pressure gradient at the model boundary.

&lsquo;solid density&rsquo; prescribes the gradient of the fluid pressure as solid density times gravity (which is the lithostatic pressure) and leads to approximately the same pressure in the melt as in the solid, so that fluid is only flowing in or out due to differences in dynamic pressure.

&lsquo;fluid density&rsquo; prescribes the gradient of the fluid pressure as fluid density times gravity and causes melt to flow in with the same velocity as inflowing solid material, or no melt flowing in or out if the solid velocity normal to the boundary is zero.

&rsquo;average density&rsquo; prescribes the gradient of the fluid pressure as the averaged fluid and solid density times gravity (which is a better approximation for the lithostatic pressure than just the solid density) and leads to approximately the same pressure in the melt as in the solid, so that fluid is only flowing in or out due to differences in dynamic pressure.

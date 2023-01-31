(parameters:Melt_20settings)=
# Melt settings


## **Subsection:** Melt settings


(parameters:Melt_20settings/Average_20melt_20velocity)=
### __Parameter name:__ Average melt velocity
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to cell-wise average the material properties that are used to compute the melt velocity or not. The melt velocity is computed as the sum of the solid velocity and the phase separation flux $ - K_D / \phi (\nabla p_f - \rho_f \mathbf g)$. If this parameter is set to true, $K_D$ and $\phi$ will be averaged cell-wise in the computation of the phase separation flux. This is useful because in some models the melt velocity can have spikes close to the interface between regions of melt and no melt, as both $K_D$ and $\phi$ go to zero for vanishing melt fraction. As the melt velocity is used for computing the time step size, and in models that use heat transport by melt or shear heating of melt, setting this parameter to true can speed up the model and make it mode stable. In computations where accuracy and convergence behavior of the melt velocity is important (like in benchmark cases with an analytical solution), this parameter should probably be set to &rsquo;false&rsquo;.

(parameters:Melt_20settings/Heat_20advection_20by_20melt)=
### __Parameter name:__ Heat advection by melt
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a porosity weighted average of the melt and solid velocity to advect heat in the temperature equation or not. If this is set to true, additional terms are assembled on the left-hand side of the temperature advection equation. Only used if Include melt transport is true. If this is set to false, only the solid velocity is used (as in models without melt migration).

(parameters:Melt_20settings/Include_20melt_20transport)=
### __Parameter name:__ Include melt transport
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to include the transport of melt into the model or not. If this is set to true, two additional pressures (the fluid pressure and the compaction pressure) will be added to the finite element. Including melt transport in the simulation also requires that there is one compositional field that has the name &lsquo;porosity&rsquo;. This field will be used for computing the additional pressures and the melt velocity, and has a different advection equation than other compositional fields, as it is effectively advected with the melt velocity.

(parameters:Melt_20settings/Melt_20scaling_20factor_20threshold)=
### __Parameter name:__ Melt scaling factor threshold
**Default value:** 1e-7

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The factor by how much the Darcy coefficient K\_D in a cell can be smaller than the reference Darcy coefficient for this cell still to be considered a melt cell (for which the melt transport equations are solved). For smaller Darcy coefficients, the Stokes equations (without melt) are solved instead. Only used if &ldquo;Include melt transport&rdquo; is true.

(parameters:Melt_20settings/Use_20discontinuous_20compaction_20pressure)=
### __Parameter name:__ Use discontinuous compaction pressure
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to use a discontinuous element for the compaction pressure or not. From our preliminary tests, continuous elements seem to work better in models where the porosity is > 0 everywhere in the domain, and discontinuous elements work better in models where in parts of the domain the porosity = 0.

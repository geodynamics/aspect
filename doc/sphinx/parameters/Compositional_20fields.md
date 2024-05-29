(parameters:Compositional_20fields)=
# Compositional fields


## **Subsection:** Compositional fields


(parameters:Compositional_20fields/Compositional_20field_20methods)=
### __Parameter name:__ Compositional field methods
**Default value:**

**Pattern:** [List of <[Selection field|particles|volume of fluid|static|melt field|darcy field|prescribed field|prescribed field with diffusion ]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting the solution method of each compositional field. Each entry of the list must be one of the currently implemented field methods.

These choices correspond to the following methods by which compositional fields gain their values:\begin{itemize}\item &ldquo;field&rdquo;: If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the velocity field, and applying reaction rates to it. In other words, this corresponds to the usual notion of a composition field as mentioned in {ref}`sec:methods:compositional-fields`.
\item &ldquo;particles&rdquo;: If a compositional field is marked with this method, then its values are obtained in each time step by interpolating the corresponding properties from the particles located on each cell. The time evolution therefore happens because particles move along with the velocity field, and particle properties can react with each other as well. See {ref}`sec:methods:particles` for more information about how particles behave.
\item &ldquo;volume of fluid&ldquo;: If a compositional field is marked with this method, then its values are obtained in each timestep by reconstructing a polynomial finite element approximation on each cell from a volume of fluid interface tracking method, which is used to compute the advection updates.
\item &ldquo;static&rdquo;: If a compositional field is marked this way, then it does not evolve at all. Its values are simply set to the initial conditions, and will then never change.
\item &ldquo;melt field&rdquo;: If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the melt velocity, and applying reaction rates to it. In other words, this corresponds to the usual notion of a composition field as mentioned in {ref}`sec:methods:compositional-fields`, except that it is advected with the melt velocity instead of the solid velocity. This method can only be chosen if melt transport is active in the model.
\item &ldquo;darcy field&rdquo;: If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the fluid velocity prescribed by Darcy&rsquo;s Law, and applying reaction rates to it. In other words this corresponds to the usual notion of a composition field as mentioned in {ref}`sec:methods:compositional-fields`, except that it is advected with the Darcy velocity instead of the solid velocity. This method requires there to be a compositional field named porosity that is advected the darcy field method. We calculate the fluid velocity $u_f$ using an approximation of Darcy&rsquo;s Law: $u_f = u_s - K_D / \phi * (rho_s * g - rho_f * g)$.
\item &ldquo;prescribed field&rdquo;: The value of these fields is determined in each time step from the material model. If a compositional field is marked with this method, then the value of a specific additional material model output, called the &lsquo;PrescribedFieldOutputs&rsquo; is interpolated onto the field. This field does not change otherwise, it is not advected with the flow.
\item &ldquo;prescribed field with diffusion&rdquo;: If a compositional field is marked this way, the value of a specific additional material model output, called the &lsquo;PrescribedFieldOutputs&rsquo; is interpolated onto the field, as in the &ldquo;prescribed field&rdquo; method. Afterwards, the field is diffused based on a solver parameter, the diffusion length scale, smoothing the field. Specifically, the field is updated by solving the equation $(I-l^2 \Delta) C_\text{smoothed} = C_\text{prescribed}$, where $l$ is the diffusion length scale. Note that this means that the amount of diffusion is independent of the time step size, and that the field is not advected with the flow.\end{itemize}

(parameters:Compositional_20fields/List_20of_20normalized_20fields)=
### __Parameter name:__ List of normalized fields
**Default value:**

**Pattern:** [List of <[Integer range 0...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of integers smaller than or equal to the number of compositional fields. All compositional fields in this list will be normalized before the first timestep. The normalization is implemented in the following way: First, the sum of the fields to be normalized is calculated at every point and the global maximum is determined. Second, the compositional fields to be normalized are divided by this maximum.

(parameters:Compositional_20fields/Mapped_20particle_20properties)=
### __Parameter name:__ Mapped particle properties
**Default value:**

**Pattern:** [Map of <[Anything]>:<[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting the particle properties that will be projected to those compositional fields that are of the &ldquo;particles&rdquo; field type.

The format of valid entries for this parameter is that of a map given as &ldquo;key1: value1, key2: value2 [component2], key3: value3 [component4], ...&rdquo; where each key must be a valid field name of the &ldquo;particles&rdquo; type, and each value must be one of the currently selected particle properties. Component is a component index of the particle property that is 0 by default, but can be set up to n-1, where n is the number of vector components of this particle property. The component indicator only needs to be set if not the first component of the particle property should be mapped (e.g. the $y$-component of the velocity at the particle positions).

(parameters:Compositional_20fields/Names_20of_20fields)=
### __Parameter name:__ Names of fields
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A user-defined name for each of the compositional fields requested.

(parameters:Compositional_20fields/Number_20of_20fields)=
### __Parameter name:__ Number of fields
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of fields that will be advected along with the flow field, excluding velocity, pressure and temperature.

(parameters:Compositional_20fields/Types_20of_20fields)=
### __Parameter name:__ Types of fields
**Default value:** unspecified

**Pattern:** [List of <[Selection chemical composition|stress|strain|grain size|porosity|density|entropy|generic|unspecified ]> of length 0...4294967295 (inclusive)]

**Documentation:** A type for each of the compositional fields requested. Each entry of the list must be one of several recognized types: chemical composition, stress, strain, grain size, porosity, density, entropy, general and unspecified. The generic type is intended to be a placeholder type that has no effect on the running of any material model, while the unspecified type is intended to tell ASPECT that the user has not explicitly indicated the type of field (facilitating parameter file checking). Plugins such as material models can use these types to affect how that plugin functions.

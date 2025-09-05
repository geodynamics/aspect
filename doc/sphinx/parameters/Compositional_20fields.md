(parameters:Compositional_20fields)=
# Compositional fields


## **Subsection:** Compositional fields


(parameters:Compositional_20fields/Compositional_20field_20methods)=
### __Parameter name:__ Compositional field methods
**Default value:**

**Pattern:** [List of <[Selection field|particles|volume of fluid|static|melt field|darcy field|prescribed field|prescribed field with diffusion ]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting the solution method of each compositional field. Each entry of the list must be one of the currently implemented field methods.

These choices correspond to the following methods by which compositional fields gain their values:* &ldquo;field&rdquo;: If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the velocity field, and applying reaction rates to it. In other words, this corresponds to the usual notion of a composition field as mentioned in {ref}`sec:methods:compositional-fields`.
* &ldquo;particles&rdquo;: If a compositional field is marked with this method, then its values are obtained in each time step by interpolating the corresponding properties from the particles located on each cell. The time evolution therefore happens because particles move along with the velocity field, and particle properties can react with each other as well. See {ref}`sec:methods:particles` for more information about how particles behave.
* &ldquo;volume of fluid&ldquo;: If a compositional field is marked with this method, then its values are obtained in each timestep by reconstructing a polynomial finite element approximation on each cell from a volume of fluid interface tracking method, which is used to compute the advection updates.
* &ldquo;static&rdquo;: If a compositional field is marked this way, then it does not evolve at all. Its values are simply set to the initial conditions, and will then never change.
* &ldquo;melt field&rdquo;: If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the melt velocity, and applying reaction rates to it. In other words, this corresponds to the usual notion of a composition field as mentioned in {ref}`sec:methods:compositional-fields`, except that it is advected with the melt velocity instead of the solid velocity. This method can only be chosen if melt transport is active in the model.
* &ldquo;darcy field&rdquo;: If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the fluid velocity prescribed by Darcy&rsquo;s Law, and applying reaction rates to it. In other words this corresponds to the usual notion of a composition field as mentioned in {ref}`sec:methods:compositional-fields`, except that it is advected with the Darcy velocity instead of the solid velocity. This method requires there to be a compositional field named porosity that is advected the darcy field method. We calculate the fluid velocity $u_f$ using an approximation of Darcy&rsquo;s Law: $u_f = u_s - K_D / \phi * (rho_s * g - rho_f * g)$.
* &ldquo;prescribed field&rdquo;: The value of these fields is determined in each time step from the material model. If a compositional field is marked with this method, then the value of a specific additional material model output, called the &lsquo;PrescribedFieldOutputs&rsquo; is interpolated onto the field. This field does not change otherwise, it is not advected with the flow.
* &ldquo;prescribed field with diffusion&rdquo;: If a compositional field is marked this way, the value of a specific additional material model output, called the &lsquo;PrescribedFieldOutputs&rsquo; is interpolated onto the field, as in the &ldquo;prescribed field&rdquo; method. Afterwards, the field is diffused based on a solver parameter, the diffusion length scale, smoothing the field. Specifically, the field is updated by solving the equation $(I-l^2 \Delta) C_\text{smoothed} = C_\text{prescribed}$, where $l$ is the diffusion length scale. Note that this means that the amount of diffusion is independent of the time step size, and that the field is not advected with the flow.

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

**Documentation:** A comma separated list denoting a &ldquo;type&rdquo; for each of the compositional fields requested. ASPECT uses these types to determine how fields are handled when evaluating the material model and when solving the equations as described below.

Each entry of the list must be one of several recognized types: * &ldquo;chemical composition&rdquo;: This type of field represents the bulk composition of the modeled material. This type of field is generally considered by the material model to determine the equation of state, rheology, and reactions.
* &ldquo;stress&rdquo;: This type of field represents stress in the material. Whether the fields represents a scalar stress invariant or tensor components, and which type of stress is represented depends on the interpretation of the material model.
* &ldquo;strain&rdquo;: This type of field represents accumulated strain. It behaves similar to the type &ldquo;stress&rdquo; discussed above except tracking the accumulated strain.
* &ldquo;grain size&rdquo;: This type of field represents an average mineral grain size of the material. It will only be considered in material models that include models for grain size evolution.
* &ldquo;porosity&rdquo;: This type of field represents porosity in a two-phase flow or Darcy flow system. Note that setting the type of a compositional field to &ldquo;porosity&rdquo; does not automatically enable melt transport, which is done with the parameter &ldquo;Melt settings/Include melt transport&rdquo;.
* &ldquo;density&rdquo;: This type of field is a finite-element field representation of the density in the model. This field type is not usually used except for the projected density approximation of the compressible Stokes equations, which uses this field type to compute gradients and time-derivatives of the density.
* &ldquo;entropy&rdquo;: This type of field represents entropy. If one or more entropy fields are found in a model, they automatically replace temperature as the main thermodynamic state variable in the model. The temperature equation is then automatically changed to a pure diffusion equation, which is coupled to the entropy advection equation as described in the paper {cite}`dannberg:etal:2022`.
* &ldquo;generic&rdquo;: The generic type is intended to be a placeholder type that is not used by any component of ASPECT unless in user-provided source code.
* &ldquo;unspecified&rdquo;: The unspecified type is intended to tell ASPECT that the user has not explicitly indicated the type of this field. ASPECT will then try to detect the type automatically based on the name, but will default to &ldquo;chemical composition&rdquo; if the name does not correspond to a known type.

Note that while ASPECT&rsquo;s functionality can make use of the field types, not all of the code will make use of it. It is the user&rsquo;s responsibility to check that the chosen material model and other plugins interpret the compositional fields as intended.

(parameters:Mesh_20refinement)=
# Mesh refinement


## **Subsection:** Mesh refinement


(parameters:Mesh_20refinement/Adapt_20by_20fraction_20of_20cells)=
### __Parameter name:__ Adapt by fraction of cells
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Use fraction of the total number of cells instead of fraction of the total error as the limit for refinement and coarsening.

(parameters:Mesh_20refinement/Additional_20refinement_20times)=
### __Parameter name:__ Additional refinement times
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of times so that if the end time of a time step is beyond this time, an additional round of mesh refinement is triggered. This is mostly useful to make sure we can get through the initial transient phase of a simulation on a relatively coarse mesh, and then refine again when we are in a time range that we are interested in and where we would like to use a finer mesh. Units: Each element of the list has units years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Mesh_20refinement/Coarsening_20fraction)=
### __Parameter name:__ Coarsening fraction
**Default value:** 0.05

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** Cells are sorted from largest to smallest by their total error (determined by the Strategy). Then the cells with the smallest error (bottom of this sorted list) that account for the given fraction of the error are coarsened.

(parameters:Mesh_20refinement/Initial_20adaptive_20refinement)=
### __Parameter name:__ Initial adaptive refinement
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of adaptive refinement steps performed after initial global refinement but while still within the first time step. These refinement steps (n) are added to the value for initial global refinement (m) so that the final mesh has cells that are at most on refinement level $n+m$.

(parameters:Mesh_20refinement/Initial_20global_20refinement)=
### __Parameter name:__ Initial global refinement
**Default value:** 2

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of global refinement steps performed on the initial coarse mesh, before the problem is first solved there.

Note that it is possible to supply conflicting refinement and coarsening settings, such as an &rsquo;Initial global refinement&rsquo; of 4 and a &rsquo;Maximum refinement function&rsquo; strategy that limits the refinement locally to 2. In this case, the tagging strategies such as the &rsquo;Maximum refinement function&rsquo; will remove refinement flags in each initial global refinement step, such that the resulting mesh is not necessarily uniform or of the level given by the &rsquo;Initial global refinement&rsquo; parameter.

(parameters:Mesh_20refinement/Minimum_20refinement_20level)=
### __Parameter name:__ Minimum refinement level
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The minimum refinement level each cell should have, and that can not be exceeded by coarsening. Should not be higher than the &rsquo;Initial global refinement&rsquo; parameter.

(parameters:Mesh_20refinement/Normalize_20individual_20refinement_20criteria)=
### __Parameter name:__ Normalize individual refinement criteria
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If multiple refinement criteria are specified in the &ldquo;Strategy&rdquo; parameter, then they need to be combined somehow to form the final refinement indicators. This is done using the method described by the &ldquo;Refinement criteria merge operation&rdquo; parameter which can either operate on the raw refinement indicators returned by each strategy (i.e., dimensional quantities) or using normalized values where the indicators of each strategy are first normalized to the interval $[0,1]$ (which also makes them non-dimensional). This parameter determines whether this normalization will happen.

(parameters:Mesh_20refinement/Refinement_20criteria_20merge_20operation)=
### __Parameter name:__ Refinement criteria merge operation
**Default value:** max

**Pattern:** [Selection plus|max ]

**Documentation:** If multiple mesh refinement criteria are computed for each cell (by passing a list of more than element to the `Strategy` parameter in this section of the input file) then one will have to decide which criteria should win when deciding which cells to refine. The operation that determines how to combine these competing criteria is the one that is selected here. The options are:

\begin{itemize}
\item `plus`: Add the various error indicators together and refine those cells on which the sum of indicators is largest.
\item `max`: Take the maximum of the various error indicators and refine those cells on which the maximal indicators is largest.
\end{itemize}The refinement indicators computed by each strategy are modified by the &ldquo;Normalize individual refinement criteria&rdquo; and &ldquo;Refinement criteria scale factors&rdquo; parameters.

(parameters:Mesh_20refinement/Refinement_20criteria_20scaling_20factors)=
### __Parameter name:__ Refinement criteria scaling factors
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of scaling factors by which every individual refinement criterion will be multiplied by. If only a single refinement criterion is selected (using the &ldquo;Strategy&rdquo; parameter, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other.

If &ldquo;Normalize individual refinement criteria&rdquo; is set to true, then the criteria will first be normalized to the interval $[0,1]$ and then multiplied by the factors specified here. You will likely want to choose the factors to be not too far from 1 in that case, say between 1 and 10, to avoid essentially disabling those criteria with small weights. On the other hand, if the criteria are not normalized to $[0,1]$ using the parameter mentioned above, then the factors you specify here need to take into account the relative numerical size of refinement indicators (which in that case carry physical units).

You can experimentally play with these scaling factors by choosing to output the refinement indicators into the graphical output of a run.

If the list of indicators given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are indicators chosen in the &ldquo;Strategy&rdquo; parameter.

(parameters:Mesh_20refinement/Refinement_20fraction)=
### __Parameter name:__ Refinement fraction
**Default value:** 0.3

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** Cells are sorted from largest to smallest by their total error (determined by the Strategy). Then the cells with the largest error (top of this sorted list) that account for given fraction of the error are refined.

(parameters:Mesh_20refinement/Run_20postprocessors_20on_20initial_20refinement)=
### __Parameter name:__ Run postprocessors on initial refinement
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not the postprocessors should be executed after each of the initial adaptive refinement cycles that are run at the start of the simulation. This is useful for plotting/analyzing how the mesh refinement parameters are working for a particular model.

(parameters:Mesh_20refinement/Skip_20setup_20initial_20conditions_20on_20initial_20refinement)=
### __Parameter name:__ Skip setup initial conditions on initial refinement
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not the initial conditions should be set up during the adaptive refinement cycles that are run at the start of the simulation.

(parameters:Mesh_20refinement/Skip_20solvers_20on_20initial_20refinement)=
### __Parameter name:__ Skip solvers on initial refinement
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not solvers should be executed during the initial adaptive refinement cycles that are run at the start of the simulation.

(parameters:Mesh_20refinement/Strategy)=
### __Parameter name:__ Strategy
**Default value:** thermal energy density

**Pattern:** [MultipleSelection artificial viscosity|boundary|compaction length|composition|composition approximate gradient|composition gradient|composition threshold|density|isosurfaces|maximum refinement function|minimum refinement function|nonadiabatic temperature|nonadiabatic temperature threshold|particle density|slope|strain rate|temperature|thermal energy density|topography|velocity|viscosity|volume of fluid interface ]

**Documentation:** A comma separated list of mesh refinement criteria that will be run whenever mesh refinement is required. The results of each of these criteria, i.e., the refinement indicators they produce for all the cells of the mesh will then be normalized to a range between zero and one and the results of different criteria will then be merged through the operation selected in this section.

The following criteria are available:

&lsquo;artificial viscosity&rsquo;: A mesh refinement criterion that computes refinement indicators from the artificial viscosity of the temperature or compositional fields based on user specified weights.

&lsquo;boundary&rsquo;: A class that implements a mesh refinement criterion which always flags all cells on specified boundaries for refinement. This is useful to provide high accuracy for processes at or close to the edge of the model domain.

To use this refinement criterion, you may want to combine it with other refinement criteria, setting the &rsquo;Normalize individual refinement criteria&rsquo; flag and using the &lsquo;max&rsquo; setting for &rsquo;Refinement criteria merge operation&rsquo;.

&lsquo;compaction length&rsquo;: A mesh refinement criterion for models with melt transport that computes refinement indicators based on the compaction length, defined as $\delta = \sqrt{\frac{(\xi + 4 \eta/3) k}{\eta_f}}$. $\xi$ is the bulk viscosity, $\eta$ is the shear viscosity, $k$ is the permeability and $\eta_f$ is the melt viscosity. If the cell width or height exceeds a multiple (which is specified as an input parameter) of this compaction length, the cell is marked for refinement.

&lsquo;composition&rsquo;: A mesh refinement criterion that computes refinement indicators from the compositional fields. If there is more than one compositional field, then it simply takes the sum of the indicators computed from each of the compositional field.

The way these indicators are computed is by evaluating the &lsquo;Kelly error indicator&rsquo; on each compositional field. This error indicator takes the finite element approximation of the compositional field and uses it to compute an approximation of the second derivatives of the composition for each cell. This approximation is then multiplied by an appropriate power of the cell&rsquo;s diameter to yield an indicator for how large the error is likely going to be on this cell. This construction rests on the observation that for many partial differential equations, the error on each cell is proportional to some power of the cell&rsquo;s diameter times the second derivatives of the solution on that cell.

For complex equations such as those we solve here, this observation may not be strictly true in the mathematical sense, but it often yields meshes that are surprisingly good.

&lsquo;composition approximate gradient&rsquo;: A mesh refinement criterion that computes refinement indicators from the gradients of compositional fields. If there is more than one compositional field, then it simply takes the sum of the indicators times a user-specified weight for each field.

In contrast to the &lsquo;composition gradient&rsquo; refinement criterion, the current criterion does not compute the gradient at quadrature points on each cell, but by a finite difference approximation between the centers of cells. Consequently, it also works if the compositional fields are computed using discontinuous finite elements.

&lsquo;composition gradient&rsquo;: A mesh refinement criterion that computes refinement indicators from the gradients of compositional fields. If there is more than one compositional field, then it simply takes the sum of the indicators times a user-specified weight for each field.

This refinement criterion computes the gradient of the compositional field at quadrature points on each cell, and then averages them in some way to obtain a refinement indicator for each cell. This will give a reasonable approximation of the true gradient of the compositional field if you are using a continuous finite element.

On the other hand, for discontinuous finite elements (see the &lsquo;Use discontinuous composition discretization&rsquo; parameter in the &lsquo;Discretization&rsquo; section), the gradient at quadrature points does not include the contribution of jumps in the compositional field between cells, and consequently will not be an accurate approximation of the true gradient. As an extreme example, consider the case of using piecewise constant finite elements for compositional fields; in that case, the gradient of the solution at quadrature points inside each cell will always be exactly zero, even if the finite element solution is different from each cell to the next. Consequently, the current refinement criterion will likely not be useful in this situation. That said, the &lsquo;composition approximate gradient&rsquo; refinement criterion exists for exactly this purpose.

&lsquo;composition threshold&rsquo;: A mesh refinement criterion that computes refinement indicators from the compositional fields. One threshold per compositional is given in the input file, and if any field exceeds its threshold, the cell is marked for refinement.

&lsquo;density&rsquo;: A mesh refinement criterion that computes refinement indicators from a field that describes the spatial variability of the density, $\rho$. Because this quantity may not be a continuous function ($\rho$ and $C_p$ may be discontinuous functions along discontinuities in the medium, for example due to phase changes), we approximate the gradient of this quantity to refine the mesh. The error indicator defined here takes the magnitude of the approximate gradient and scales it by $h_K^{1+d/2}$ where $h_K$ is the diameter of each cell and $d$ is the dimension. This scaling ensures that the error indicators converge to zero as $h_K\rightarrow 0$ even if the energy density is discontinuous, since the gradient of a discontinuous function grows like $1/h_K$.

&lsquo;isosurfaces&rsquo;: A mesh refinement criterion that computes coarsening and refinement indicators between two isosurfaces of specific field entries (e.g. temperature, composition).

The way these indicators are derived between pairs of isosurfaces is by checking whether the solutions of specific fields are within the ranges of the isosurface values given. If these conditions hold, then coarsening and refinement indicators are set such that the mesh refinement levels lies within the range of levels given. Usage of this plugin allows the user to put a conditional minimum and maximum refinement function onto fields that they are interested in.

For now, only temperature and compositional fields are allowed as field entries. The key words could be &rsquo;Temperature&rsquo; or one of the names of the compositional fields which are either specified by user or set up as C\_0, C\_1, etc.

Usage: A list of isosurfaces separated by semi-colons (;). Each isosurface entry consists of multiple entries separated by a comma. The first two entries indicate the minimum and maximum refinement levels respectively. The entries after the first two describe the fields the isosurface applies to, followed by a colon (:), which again is followed by the minimum and maximum field values separated by a bar (|). An example for two isosurface entries is &rsquo;0, 2, Temperature: 300 | 600; 2, 2, C\_1: 0.5 | 1&rsquo;. If both isoterm entries are triggered at the same location and the current refinement level is 1, it means that the first isoline will not set any flag and the second isoline will set a refinement flag. This means the cell will be refined. If both the coarsening and refinement flags are set, preference is given to refinement.

The minimum and maximum refinement levels per isosurface can be provided in absolute values relative to the global minimum and maximum refinement. This is done with the &rsquo;min&rsquo; and &rsquo;max&rsquo; key words. For example: &rsquo;set Isosurfaces = max-2,  max,    Temperature: 0 | 600 ; min + 1,min+2, Temperature: 1600 | 3000,   C\_2 : 0.0 | 0.5&rsquo;.

&lsquo;maximum refinement function&rsquo;: A mesh refinement criterion that ensures a maximum refinement level described by an explicit formula with the depth or position as argument. Which coordinate representation is used is determined by an input parameter. Whatever the coordinate system chosen, the function you provide in the input file will by default depend on variables &lsquo;x&rsquo;, &lsquo;y&rsquo; and &lsquo;z&rsquo; (if in 3d). However, the meaning of these symbols depends on the coordinate system. In the Cartesian coordinate system, they simply refer to their natural meaning. If you have selected &lsquo;depth&rsquo; for the coordinate system, then &lsquo;x&rsquo; refers to the depth variable and &lsquo;y&rsquo; and &lsquo;z&rsquo; will simply always be zero. If you have selected a spherical coordinate system, then &lsquo;x&rsquo; will refer to the radial distance of the point to the origin, &lsquo;y&rsquo; to the azimuth angle and &lsquo;z&rsquo; to the polar angle measured positive from the north pole. Note that the order of spherical coordinates is r,phi,theta and not r,theta,phi, since this allows for dimension independent expressions. Each coordinate system also includes a final &lsquo;t&rsquo; variable which represents the model time, evaluated in years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set, otherwise evaluated in seconds. After evaluating the function, its values are rounded to the nearest integer.

The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;minimum refinement function&rsquo;: A mesh refinement criterion that ensures a minimum refinement level described by an explicit formula with the depth or position as argument. Which coordinate representation is used is determined by an input parameter. Whatever the coordinate system chosen, the function you provide in the input file will by default depend on variables &lsquo;x&rsquo;, &lsquo;y&rsquo; and &lsquo;z&rsquo; (if in 3d). However, the meaning of these symbols depends on the coordinate system. In the Cartesian coordinate system, they simply refer to their natural meaning. If you have selected &lsquo;depth&rsquo; for the coordinate system, then &lsquo;x&rsquo; refers to the depth variable and &lsquo;y&rsquo; and &lsquo;z&rsquo; will simply always be zero. If you have selected a spherical coordinate system, then &lsquo;x&rsquo; will refer to the radial distance of the point to the origin, &lsquo;y&rsquo; to the azimuth angle and &lsquo;z&rsquo; to the polar angle measured positive from the north pole. Note that the order of spherical coordinates is r,phi,theta and not r,theta,phi, since this allows for dimension independent expressions. Each coordinate system also includes a final &lsquo;t&rsquo; variable which represents the model time, evaluated in years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set, otherwise evaluated in seconds. After evaluating the function, its values are rounded to the nearest integer.

The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;nonadiabatic temperature&rsquo;: A mesh refinement criterion that computes refinement indicators from the excess temperature(difference between temperature and adiabatic temperature.

&lsquo;nonadiabatic temperature threshold&rsquo;: A mesh refinement criterion that computes refinement indicators from the temperature difference between the actual temperature and the adiabatic conditions (the nonadiabatic temperature). If the temperature anomaly exceeds the threshold given in the input file, the cell is marked for refinement.

&lsquo;particle density&rsquo;: A mesh refinement criterion that computes refinement indicators based on the density of particles. In practice this plugin equilibrates the number of particles per cell, leading to fine cells in high particle density regions and coarse cells in low particle density regions. This plugin is mostly useful for models with inhomogeneous particle density, e.g. when tracking an initial interface with a high particle density, or when the spatial particle density denotes the region of interest. Additionally, this plugin tends to balance the computational load between processes in parallel computations, because the particle and mesh density is more aligned.

&lsquo;slope&rsquo;: A class that implements a mesh refinement criterion intended for use with deforming mesh boundaries, like the free surface. It calculates a local slope based on the angle between the surface normal and the local gravity vector. Cells with larger angles are marked for refinement.

To use this refinement criterion, you may want to combine it with other refinement criteria, setting the &rsquo;Normalize individual refinement criteria&rsquo; flag and using the &lsquo;max&rsquo; setting for &rsquo;Refinement criteria merge operation&rsquo;.

&lsquo;strain rate&rsquo;: A mesh refinement criterion that computes the refinement indicators equal to the strain rate norm computed at the center of the elements.

&lsquo;temperature&rsquo;: A mesh refinement criterion that computes refinement indicators from the temperature field.

The way these indicators are computed is by evaluating the &lsquo;Kelly error indicator&rsquo; on the temperature field. This error indicator takes the finite element approximation of the temperature field and uses it to compute an approximation of the second derivatives of the temperature for each cell. This approximation is then multiplied by an appropriate power of the cell&rsquo;s diameter to yield an indicator for how large the error is likely going to be on this cell. This construction rests on the observation that for many partial differential equations, the error on each cell is proportional to some power of the cell&rsquo;s diameter times the second derivatives of the solution on that cell.

For complex equations such as those we solve here, this observation may not be strictly true in the mathematical sense, but it often yields meshes that are surprisingly good.

&lsquo;thermal energy density&rsquo;: A mesh refinement criterion that computes refinement indicators from a field that describes the spatial variability of the thermal energy density, $\rho C_p T$. Because this quantity may not be a continuous function ($\rho$ and $C_p$ may be discontinuous functions along discontinuities in the medium, for example due to phase changes), we approximate the gradient of this quantity to refine the mesh. The error indicator defined here takes the magnitude of the approximate gradient and scales it by $h_K^{1.5}$ where $h_K$ is the diameter of each cell. This scaling ensures that the error indicators converge to zero as $h_K\rightarrow 0$ even if the energy density is discontinuous, since the gradient of a discontinuous function grows like $1/h_K$.

&lsquo;topography&rsquo;: A class that implements a mesh refinement criterion, which always flags all cells in the uppermost layer for refinement. This is useful to provide high accuracy for processes at or close to the surface.

To use this refinement criterion, you may want to combine it with other refinement criteria, setting the &rsquo;Normalize individual refinement criteria&rsquo; flag and using the &lsquo;max&rsquo; setting for &rsquo;Refinement criteria merge operation&rsquo;.

&lsquo;velocity&rsquo;: A mesh refinement criterion that computes refinement indicators from the velocity field.

The way these indicators are computed is by evaluating the &lsquo;Kelly error indicator&rsquo; on the velocity field. This error indicator takes the finite element approximation of the velocity field and uses it to compute an approximation of the second derivatives of the velocity for each cell. This approximation is then multiplied by an appropriate power of the cell&rsquo;s diameter to yield an indicator for how large the error is likely going to be on this cell. This construction rests on the observation that for many partial differential equations, the error on each cell is proportional to some power of the cell&rsquo;s diameter times the second derivatives of the solution on that cell.

For complex equations such as those we solve here, this observation may not be strictly true in the mathematical sense, but it often yields meshes that are surprisingly good.

&lsquo;viscosity&rsquo;: A mesh refinement criterion that computes refinement indicators from a field that describes the spatial variability of the logarithm of the viscosity, $\log\eta$. (We choose the logarithm of the viscosity because it can vary by orders of magnitude.)Because this quantity may not be a continuous function ($\eta$ may be a discontinuous function along discontinuities in the medium, for example due to phase changes), we approximate the gradient of this quantity to refine the mesh. The error indicator defined here takes the magnitude of the approximate gradient and scales it by $h_K^{1+d/2}$ where $h_K$ is the diameter of each cell and $d$ is the dimension. This scaling ensures that the error indicators converge to zero as $h_K\rightarrow 0$ even if the viscosity is discontinuous, since the gradient of a discontinuous function grows like $1/h_K$.

&lsquo;volume of fluid interface&rsquo;: A class that implements a mesh refinement criterion, which ensures a minimum level of refinement near the volume of fluid interface boundary.

(parameters:Mesh_20refinement/Time_20steps_20between_20mesh_20refinement)=
### __Parameter name:__ Time steps between mesh refinement
**Default value:** 10

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of time steps after which the mesh is to be adapted again based on computed error indicators. If 0 then the mesh will never be changed.

(parameters:Mesh_20refinement/Artificial_20viscosity)=
## **Subsection:** Mesh refinement / Artificial viscosity
(parameters:Mesh_20refinement/Artificial_20viscosity/Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of scaling factors by which every individual compositional field will be multiplied. These factors are used to weigh the various indicators relative to each other and to the temperature.

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to 0. If the list is not empty then it needs to have as many entries as there are compositional fields.

(parameters:Mesh_20refinement/Artificial_20viscosity/Temperature_20scaling_20factor)=
### __Parameter name:__ Temperature scaling factor
**Default value:** 0.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** A scaling factor for the artificial viscosity  of the temperature equation. Use 0.0 to disable.

(parameters:Mesh_20refinement/Boundary)=
## **Subsection:** Mesh refinement / Boundary
(parameters:Mesh_20refinement/Boundary/Boundary_20refinement_20indicators)=
### __Parameter name:__ Boundary refinement indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries where there should be mesh refinement.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

(parameters:Mesh_20refinement/Compaction_20length)=
## **Subsection:** Mesh refinement / Compaction length
(parameters:Mesh_20refinement/Compaction_20length/Mesh_20cells_20per_20compaction_20length)=
### __Parameter name:__ Mesh cells per compaction length
**Default value:** 1.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The desired ratio between compaction length and size of the mesh cells, or, in other words, how many cells the mesh should (at least) have per compaction length. Every cell where this ratio is smaller than the value specified by this parameter (in places with fewer mesh cells per compaction length) is marked for refinement.

(parameters:Mesh_20refinement/Composition)=
## **Subsection:** Mesh refinement / Composition
(parameters:Mesh_20refinement/Composition/Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of scaling factors by which every individual compositional field will be multiplied. If only a single compositional field exists, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other.

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are compositional fields.

(parameters:Mesh_20refinement/Composition_20approximate_20gradient)=
## **Subsection:** Mesh refinement / Composition approximate gradient
(parameters:Mesh_20refinement/Composition_20approximate_20gradient/Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of scaling factors by which every individual compositional field gradient will be multiplied. If only a single compositional field exists, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other.

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are compositional fields.

(parameters:Mesh_20refinement/Composition_20gradient)=
## **Subsection:** Mesh refinement / Composition gradient
(parameters:Mesh_20refinement/Composition_20gradient/Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of scaling factors by which every individual compositional field gradient will be multiplied. If only a single compositional field exists, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other.

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are compositional fields.

(parameters:Mesh_20refinement/Composition_20threshold)=
## **Subsection:** Mesh refinement / Composition threshold
(parameters:Mesh_20refinement/Composition_20threshold/Compositional_20field_20thresholds)=
### __Parameter name:__ Compositional field thresholds
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of thresholds, one for each compositional field to be evaluated against.

(parameters:Mesh_20refinement/Isosurfaces)=
## **Subsection:** Mesh refinement / Isosurfaces
(parameters:Mesh_20refinement/Isosurfaces/Isosurfaces)=
### __Parameter name:__ Isosurfaces
**Default value:**

**Pattern:** [Anything]

**Documentation:** A list of isosurfaces separated by semi-colons (;). Each isosurface entry consists of multiple entries separated by a comma. The first two entries indicate the minimum and maximum refinement levels respectively. The entries after the first two describe the fields the isosurface applies to, followed by a colon (:), which is again followed by the minimum and maximum property values separated by bar (|). An example for an isosurface is &rsquo;0, 2, Temperature: 300 | 600; 2, 2, C\_1: 0.5 | 1&rsquo;. In this example the mesh refinement is kept between level 0 and level 2 if the temperature is between 300 and 600 and at level 2 when the compositional field C\_1 is between 0.5 and 1. If both happen at the same location and the current refinement level is 1, it means that the first isoline will not set any flag and the second isoline will set a refinement flag. This means the cell will be refined. If both the coarsening and refinement flags are set, preference is given to refinement.

The first two entries for each isosurface, describing the minimum and maximum grid levels, can be two numbers or contain one of the key values &rsquo;min&rsquo; and &rsquo;max&rsquo;. This indicates the key will be replaced with the global minimum and maximum refinement levels. The &rsquo;min&rsquo; and &rsquo;max&rsquo; keys also accept adding values to be added or subtracted from them respectively. This is done by adding a &rsquo;+&rsquo; or &rsquo;-&rsquo; and a number behind them (e.g. min+2 or max-1). Note that you can&rsquo;t subtract a value from a minimum value or add a value to the maximum value. If, for example, &lsquo;max-4&lsquo; drops below the minimum or &lsquo;min+4&lsquo; goes above the maximum, it will simply use the global minimum and maximum values respectively. The same holds for any mesh refinement level below the global minimum or above the global maximum.

(parameters:Mesh_20refinement/Maximum_20refinement_20function)=
## **Subsection:** Mesh refinement / Maximum refinement function
(parameters:Mesh_20refinement/Maximum_20refinement_20function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** depth

**Pattern:** [Selection depth|cartesian|spherical ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;depth&rsquo;, &lsquo;cartesian&rsquo; and &lsquo;spherical&rsquo;. &lsquo;depth&rsquo; will create a function, in which only the first variable is non-zero, which is interpreted to be the depth of the point. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle.

(parameters:Mesh_20refinement/Maximum_20refinement_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Mesh_20refinement/Maximum_20refinement_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Mesh_20refinement/Maximum_20refinement_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Mesh_20refinement/Minimum_20refinement_20function)=
## **Subsection:** Mesh refinement / Minimum refinement function
(parameters:Mesh_20refinement/Minimum_20refinement_20function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** depth

**Pattern:** [Selection depth|cartesian|spherical ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;depth&rsquo;, &lsquo;cartesian&rsquo; and &lsquo;spherical&rsquo;. &lsquo;depth&rsquo; will create a function, in which only the first variable is non-zero, which is interpreted to be the depth of the point. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle.

(parameters:Mesh_20refinement/Minimum_20refinement_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Mesh_20refinement/Minimum_20refinement_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Mesh_20refinement/Minimum_20refinement_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Mesh_20refinement/Nonadiabatic_20temperature_20threshold)=
## **Subsection:** Mesh refinement / Nonadiabatic temperature threshold
(parameters:Mesh_20refinement/Nonadiabatic_20temperature_20threshold/Temperature_20anomaly_20type)=
### __Parameter name:__ Temperature anomaly type
**Default value:** absolute value

**Pattern:** [Selection negative only|positive only|absolute value ]

**Documentation:** What type of temperature anomaly should be considered when evaluating against the threshold: Only negative anomalies (negative only), only positive anomalies (positive only) or the absolute value of the nonadiabatic temperature.

(parameters:Mesh_20refinement/Nonadiabatic_20temperature_20threshold/Threshold)=
### __Parameter name:__ Threshold
**Default value:** 100

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** A threshold that the nonadiabatic temperature will be evaluated against. Units: \si{\kelvin}

(parameters:Mesh_20refinement/Volume_20of_20fluid_20interface)=
## **Subsection:** Mesh refinement / Volume of fluid interface
(parameters:Mesh_20refinement/Volume_20of_20fluid_20interface/Strict_20coarsening)=
### __Parameter name:__ Strict coarsening
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If true, then explicitly coarsen any cells not neighboring the VolumeOfFluid interface.

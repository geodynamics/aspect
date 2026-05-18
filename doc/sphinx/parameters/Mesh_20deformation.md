(parameters:Mesh_20deformation)=
# Mesh deformation


## **Subsection:** Mesh deformation


::::{dropdown} __Parameter:__ {ref}`Additional tangential mesh velocity boundary indicators<parameters:Mesh_20deformation/Additional_20tangential_20mesh_20velocity_20boundary_20indicators>`
:name: parameters:Mesh_20deformation/Additional_20tangential_20mesh_20velocity_20boundary_20indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries where there the mesh is allowed to move tangential to the boundary. All tangential mesh movements along those boundaries that have tangential material velocity boundary conditions are allowed by default, this parameters allows to generate mesh movements along other boundaries that are open, or have prescribed material velocities or tractions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.
::::

::::{dropdown} __Parameter:__ {ref}`Mesh deformation boundary indicators<parameters:Mesh_20deformation/Mesh_20deformation_20boundary_20indicators>`
:name: parameters:Mesh_20deformation/Mesh_20deformation_20boundary_20indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries where there the mesh is allowed to move according to the specified mesh deformation objects.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

The format is id1: object1 \& object2, id2: object3 \& object2, where objects are one of &lsquo;ascii data&rsquo;: Implementation of a model in which the initial mesh deformation (initial topography) is derived from a file containing data in ascii format. The following geometry models are currently supported: box, chunk, spherical. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;Topography [m]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Topography [m]&rsquo; in a 3d model, which means that there has to be a single column containing the topography. Note that the data in the input file needs to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle in radians  and &lsquo;y&rsquo; by the polar angle in radians measured positive from the north pole. The grid will be assumed to be a longitude-colatitude grid. Note that the order of spherical coordinates is &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;boundary function&rsquo;: A plugin, which prescribes the surface mesh to deform according to an analytically prescribed function. Note that the function prescribes a deformation velocity, i.e. the return value of this plugin is later multiplied by the time step length to compute the displacement increment in this time step. Although the function&rsquo;s time variable is interpreted as years when Use years instead of seconds is set to true, the boundary deformation velocity should still be given in m/s. The format of the functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;diffusion&rsquo;: A plugin that computes the deformation of surface vertices according to the solution of the hillslope diffusion problem. Specifically, at the end of each timestep, or after a specific number of timesteps, this plugin solves the following equation: \begin{align*}  \frac{\partial h}{\partial t} = \kappa \left( \frac{\partial^{2} h}{\partial x^{2}} + \frac{\partial^{2} h}{\partial y^{2}} \right), \end{align*} where $\kappa$ is the hillslope diffusion coefficient (diffusivity), and $h(x,y)$ the height of a point along the top boundary with respect to the surface of the unperturbed domain.

Using this definition, the plugin then solves for one time step, i.e., using as initial condition $h(t_{n-1})$ the current surface elevation, and computing $h(t_n)$ from it by solving the equation above over the time interval $t_n-t_{n-1}$. From this, one can then compute a surface velocity $v = \frac{h(t_n)-h(t_{n-1})}{t_n-t_{n-1}}$.

This surface velocity is used to deform the surface and as a boundary condition for solving the Laplace equation to determine the mesh velocity in the domain interior. Diffusion can be applied every timestep, mimicking surface processes of erosion and deposition, or at a user-defined timestep interval to purely smooth the surface topography to avoid too great a distortion of mesh elements when a free surface is also used.

&lsquo;fastscape&rsquo;: A plugin that uses the program FastScape to compute the deformation of the mesh surface. FastScape is a surface processes code that computes the erosion, transport and deposition of sediments both on land and in the marine domain. These surface processes include river incision (through the stream power law), hillslope diffusion and marine diffusion, as described in Braun and Willett 2013; Yuan et al. 2019; Yuan et al. 2019b.
Upon initialization, FastScape requires the initial topography of the surface boundary of ASPECT&rsquo;s model domain and several user-specified erosional and depositional parameters. In each ASPECT timestep, FastScape is then fed ASPECT&rsquo;s material velocity at the surface boundary. The z-component of this velocity is used to uplift the FastScape surface, while the horizontal components are used to advect the topography in the x-y plane.
After solving its governing equations (this can be done in several timesteps that are smaller than the ASPECT timestep), FastScape returns a new topography of the surface. The difference in topography before and after the call to FastScape divided by the ASPECT timestep provides the mesh velocity at the domain&rsquo;s surface that is used to displace the surface and internal mesh nodes.
FastScape can be used in both 2D and 3D ASPECT simulations. In 2D, one can think of the coupled model as a T-model.The ASPECT domain spans the x - z plane, while FastScape acts on the horizontal x-y plane. This means that to communicate ASPECT&rsquo;s material velocities to FastScape, FastScape mesh nodes with the same x-coordinate (so lying along the y-direction) get the same velocities. In turn, the FastScape topography is collapsed back onto the line of the ASPECT surface boundary by averaging the topography over the y-direction. In 3D no such actions are necessary.
The FastScape manual (https://fastscape.org/fastscapelib-fortran/) provides more information on the input parameters.

&lsquo;free surface&rsquo;: A plugin that computes the deformation of surface vertices according to the solution of the flow problem. In particular this means if the surface of the domain is left open to flow, this flow will carry the mesh with it. The implementation was described in {cite}`rose_freesurface`, with the stabilization of the free surface originally described in {cite}`kaus:etal:2010`.
::::

::::{dropdown} __Parameter:__ {ref}`Mesh deformation mapping order<parameters:Mesh_20deformation/Mesh_20deformation_20mapping_20order>`
:name: parameters:Mesh_20deformation/Mesh_20deformation_20mapping_20order
**Default value:** auto

**Pattern:** [Anything]

**Documentation:** Polynomial degree used by the MappingQEulerian object for mesh deformation. Set to &lsquo;auto&lsquo; to choose a robust default. For geometries with curved elements, the mapping order is the larger of 4 and the Stokes velocity polynomial degree. For geometries without curved elements, the mapping order is 1. Set this parameter to an integer >= 1 to explicitly enforce a mapping order. In most cases, &lsquo;auto&lsquo; is recommended. Explicit values are mainly useful for stability investigations or reproducibility studies.
::::

(parameters:Mesh_20deformation/Ascii_20data_20model)=
## **Subsection:** Mesh deformation / Ascii data model
::::{dropdown} __Parameter:__ {ref}`Data directory<parameters:Mesh_20deformation/Ascii_20data_20model/Data_20directory>`
:name: parameters:Mesh_20deformation/Ascii_20data_20model/Data_20directory
**Default value:** $ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT. A trailing slash at the end of the directory path is optional; the plugin will automatically append a &rsquo;/&rsquo; when the parameters are parsed if it is missing.
::::

::::{dropdown} __Parameter:__ {ref}`Data file name<parameters:Mesh_20deformation/Ascii_20data_20model/Data_20file_20name>`
:name: parameters:Mesh_20deformation/Ascii_20data_20model/Data_20file_20name
**Default value:** box_3d_%s.0.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.
::::

::::{dropdown} __Parameter:__ {ref}`Scale factor<parameters:Mesh_20deformation/Ascii_20data_20model/Scale_20factor>`
:name: parameters:Mesh_20deformation/Ascii_20data_20model/Scale_20factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.
::::

(parameters:Mesh_20deformation/Boundary_20function)=
## **Subsection:** Mesh deformation / Boundary function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Mesh_20deformation/Boundary_20function/Function_20constants>`
:name: parameters:Mesh_20deformation/Boundary_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Mesh_20deformation/Boundary_20function/Function_20expression>`
:name: parameters:Mesh_20deformation/Boundary_20function/Function_20expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Mesh_20deformation/Boundary_20function/Variable_20names>`
:name: parameters:Mesh_20deformation/Boundary_20function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Mesh_20deformation/Diffusion)=
## **Subsection:** Mesh deformation / Diffusion
::::{dropdown} __Parameter:__ {ref}`Hillslope transport coefficient<parameters:Mesh_20deformation/Diffusion/Hillslope_20transport_20coefficient>`
:name: parameters:Mesh_20deformation/Diffusion/Hillslope_20transport_20coefficient
**Default value:** 1e-6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The hillslope transport coefficient $\kappa$ used to diffuse the free surface, either as a  stabilization step or to mimic erosional and depositional processes. Units: $\si{m^2/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Time steps between diffusion<parameters:Mesh_20deformation/Diffusion/Time_20steps_20between_20diffusion>`
:name: parameters:Mesh_20deformation/Diffusion/Time_20steps_20between_20diffusion
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of time steps between each application of diffusion.
::::

(parameters:Mesh_20deformation/Fastscape)=
## **Subsection:** Mesh deformation / Fastscape
::::{dropdown} __Parameter:__ {ref}`Additional fastscape refinement<parameters:Mesh_20deformation/Fastscape/Additional_20fastscape_20refinement>`
:name: parameters:Mesh_20deformation/Fastscape/Additional_20fastscape_20refinement
**Default value:** 0

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** How many levels above the ASPECT mesh the FastScape mesh should be refined.
::::

::::{dropdown} __Parameter:__ {ref}`Additional output variables<parameters:Mesh_20deformation/Fastscape/Additional_20output_20variables>`
:name: parameters:Mesh_20deformation/Fastscape/Additional_20output_20variables
**Default value:** river incision rate

**Pattern:** [Selection river incision rate|deposition coefficient|uplift rate ]

**Documentation:** Select one additional Fastscape variable to output in the Fastcape vtk. Output are in units of per year.
::::

::::{dropdown} __Parameter:__ {ref}`Average out of plane surface topography in 2d<parameters:Mesh_20deformation/Fastscape/Average_20out_20of_20plane_20surface_20topography_20in_202d>`
:name: parameters:Mesh_20deformation/Fastscape/Average_20out_20of_20plane_20surface_20topography_20in_202d
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If this is set to false, then a 2D model will only consider the center slice FastScape gives. If set to true, then ASPECT will average the mesh along Y excluding the ghost nodes.
::::

::::{dropdown} __Parameter:__ {ref}`Fastscape seed<parameters:Mesh_20deformation/Fastscape/Fastscape_20seed>`
:name: parameters:Mesh_20deformation/Fastscape/Fastscape_20seed
**Default value:** 1000

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Seed used for adding an initial noise to FastScape topography based on the initial noise magnitude.
::::

::::{dropdown} __Parameter:__ {ref}`Initial noise magnitude<parameters:Mesh_20deformation/Fastscape/Initial_20noise_20magnitude>`
:name: parameters:Mesh_20deformation/Fastscape/Initial_20noise_20magnitude
**Default value:** 5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum topography change from the initial noise. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Maximum surface refinement level<parameters:Mesh_20deformation/Fastscape/Maximum_20surface_20refinement_20level>`
:name: parameters:Mesh_20deformation/Fastscape/Maximum_20surface_20refinement_20level
**Default value:** 1

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** This should be set to the highest ASPECT refinement level expected at the surface.
::::

::::{dropdown} __Parameter:__ {ref}`Maximum timestep length<parameters:Mesh_20deformation/Fastscape/Maximum_20timestep_20length>`
:name: parameters:Mesh_20deformation/Fastscape/Maximum_20timestep_20length
**Default value:** 10e3

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum timestep for FastScape. Units: ${yrs}$
::::

::::{dropdown} __Parameter:__ {ref}`Node tolerance<parameters:Mesh_20deformation/Fastscape/Node_20tolerance>`
:name: parameters:Mesh_20deformation/Fastscape/Node_20tolerance
**Default value:** 0.001

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Node tolerance for how close an ASPECT node must be to a FastScape node for the value to be transferred.
::::

::::{dropdown} __Parameter:__ {ref}`Number of fastscape timesteps per aspect timestep<parameters:Mesh_20deformation/Fastscape/Number_20of_20fastscape_20timesteps_20per_20aspect_20timestep>`
:name: parameters:Mesh_20deformation/Fastscape/Number_20of_20fastscape_20timesteps_20per_20aspect_20timestep
**Default value:** 5

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Initial number of fastscape time steps per ASPECT timestep, this value will double if the FastScape timestep is above the maximum FastScape timestep.
::::

::::{dropdown} __Parameter:__ {ref}`Sediment rain rates<parameters:Mesh_20deformation/Fastscape/Sediment_20rain_20rates>`
:name: parameters:Mesh_20deformation/Fastscape/Sediment_20rain_20rates
**Default value:** 0,0

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Sediment rain rates given as a list 1 greater than the number of sediment rain time intervals. E.g, If the time interval is given at 5 Myr, there will be one value for 0-5 Myr model time and a second value for 5+ Myr. Units: ${m/yr}$
::::

::::{dropdown} __Parameter:__ {ref}`Sediment rain time intervals<parameters:Mesh_20deformation/Fastscape/Sediment_20rain_20time_20intervals>`
:name: parameters:Mesh_20deformation/Fastscape/Sediment_20rain_20time_20intervals
**Default value:** 0

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of times to change the sediment rain rate. Units: ${yrs}$
::::

::::{dropdown} __Parameter:__ {ref}`Surface refinement difference<parameters:Mesh_20deformation/Fastscape/Surface_20refinement_20difference>`
:name: parameters:Mesh_20deformation/Fastscape/Surface_20refinement_20difference
**Default value:** 0

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** The difference between the lowest and highest refinement level at the surface. E.g., if three resolution levels are expected, this would be set to two.
::::

::::{dropdown} __Parameter:__ {ref}`Uplift and advect with fastscape<parameters:Mesh_20deformation/Fastscape/Uplift_20and_20advect_20with_20fastscape>`
:name: parameters:Mesh_20deformation/Fastscape/Uplift_20and_20advect_20with_20fastscape
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Flag to use FastScape advection and uplift.
::::

::::{dropdown} __Parameter:__ {ref}`Use ghost nodes<parameters:Mesh_20deformation/Fastscape/Use_20ghost_20nodes>`
:name: parameters:Mesh_20deformation/Fastscape/Use_20ghost_20nodes
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Flag to use ghost nodes.
::::

::::{dropdown} __Parameter:__ {ref}`Use marine component<parameters:Mesh_20deformation/Fastscape/Use_20marine_20component>`
:name: parameters:Mesh_20deformation/Fastscape/Use_20marine_20component
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Flag to use the marine component of FastScape.
::::

::::{dropdown} __Parameter:__ {ref}`Vertical exaggeration<parameters:Mesh_20deformation/Fastscape/Vertical_20exaggeration>`
:name: parameters:Mesh_20deformation/Fastscape/Vertical_20exaggeration
**Default value:** -1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Vertical exaggeration for FastScape&rsquo;s VTK file. -1 outputs topography, basement, and sealevel.
::::

::::{dropdown} __Parameter:__ {ref}`Y extent in 2d<parameters:Mesh_20deformation/Fastscape/Y_20extent_20in_202d>`
:name: parameters:Mesh_20deformation/Fastscape/Y_20extent_20in_202d
**Default value:** 100000

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** FastScape Y extent when using a 2D ASPECT model. Units: ${m}$
::::

(parameters:Mesh_20deformation/Fastscape/Boundary_20conditions)=
## **Subsection:** Mesh deformation / Fastscape / Boundary conditions
::::{dropdown} __Parameter:__ {ref}`Back<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Back>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Back
**Default value:** 1

**Pattern:** [Integer range 0...1 (inclusive)]

**Documentation:** Back (top) boundary condition, where 1 is fixed and 0 is reflective.
::::

::::{dropdown} __Parameter:__ {ref}`Back front ghost nodes periodic<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Back_20front_20ghost_20nodes_20periodic>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Back_20front_20ghost_20nodes_20periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to set the ghost nodes at the FastScape back and front boundary to periodic even if &rsquo;Back&rsquo; and &rsquo;Front&rsquo; are set to fixed boundary.
::::

::::{dropdown} __Parameter:__ {ref}`Back mass flux<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Back_20mass_20flux>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Back_20mass_20flux
**Default value:** 0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Flux per unit length through the back boundary. Units: ${m^2/yr}$
::::

::::{dropdown} __Parameter:__ {ref}`Front<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Front>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Front
**Default value:** 1

**Pattern:** [Integer range 0...1 (inclusive)]

**Documentation:** Front (bottom) boundary condition, where 1 is fixed and 0 is reflective.
::::

::::{dropdown} __Parameter:__ {ref}`Front mass flux<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Front_20mass_20flux>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Front_20mass_20flux
**Default value:** 0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Flux per unit length through the front boundary. Units: ${m^2/yr}$
::::

::::{dropdown} __Parameter:__ {ref}`Left<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Left>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Left
**Default value:** 1

**Pattern:** [Integer range 0...1 (inclusive)]

**Documentation:** Left boundary condition, where 1 is fixed and 0 is reflective.
::::

::::{dropdown} __Parameter:__ {ref}`Left mass flux<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Left_20mass_20flux>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Left_20mass_20flux
**Default value:** 0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Flux per unit length through the left boundary. Units: ${m^2/yr}$
::::

::::{dropdown} __Parameter:__ {ref}`Left right ghost nodes periodic<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Left_20right_20ghost_20nodes_20periodic>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Left_20right_20ghost_20nodes_20periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to set the ghost nodes at the FastScape left and right boundary to periodic even if &rsquo;Left&rsquo; and &rsquo;Right&rsquo; are set to fixed boundary.
::::

::::{dropdown} __Parameter:__ {ref}`Right<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Right>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Right
**Default value:** 1

**Pattern:** [Integer range 0...1 (inclusive)]

**Documentation:** Right boundary condition, where 1 is fixed and 0 is reflective.
::::

::::{dropdown} __Parameter:__ {ref}`Right mass flux<parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Right_20mass_20flux>`
:name: parameters:Mesh_20deformation/Fastscape/Boundary_20conditions/Right_20mass_20flux
**Default value:** 0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Flux per unit length through the right boundary. Units: ${m^2/yr}$
::::

(parameters:Mesh_20deformation/Fastscape/Erosional_20parameters)=
## **Subsection:** Mesh deformation / Fastscape / Erosional parameters
::::{dropdown} __Parameter:__ {ref}`Bedrock deposition coefficient<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Bedrock_20deposition_20coefficient>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Bedrock_20deposition_20coefficient
**Default value:** 1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Deposition coefficient for bedrock.
::::

::::{dropdown} __Parameter:__ {ref}`Bedrock diffusivity<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Bedrock_20diffusivity>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Bedrock_20diffusivity
**Default value:** 1e-2

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Transport coefficient (diffusivity) for bedrock. Units: ${m^2/yr}$ if &ldquo;Use years instead of seconds&rdquo; is true; otherwise, the units are ${m^2/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Bedrock river incision rate<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Bedrock_20river_20incision_20rate>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Bedrock_20river_20incision_20rate
**Default value:** 1e-5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** River incision rate for bedrock in the Stream Power Law. Units: ${m^(1-2drainage_area_exponent)/yr}$ if &ldquo;Use years instead of seconds&rdquo; is true; otherwise, the units are ${m^(1-2drainage_area_exponent)/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Drainage area exponent<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Drainage_20area_20exponent>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Drainage_20area_20exponent
**Default value:** 0.4

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The drainage area exponent for the Stream Power Law (m).
::::

::::{dropdown} __Parameter:__ {ref}`Elevation factor<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Elevation_20factor>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Elevation_20factor
**Default value:** 1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Amount to multiply the bedrock river incision rate and transport coefficient by past the given orographic elevation control.
::::

::::{dropdown} __Parameter:__ {ref}`Erosional base level<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Erosional_20base_20level>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Erosional_20base_20level
**Default value:** 0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** When &rsquo;Use a fixed erosional base level&rsquo; is set to true, all ghost nodes of fixed FastScape boundaries where no mass flux is specified by the user (FastScape boundary condition set to 1 and &rsquo;Left/Right/Bottom/Top mass flux&rsquo; set to 0) will be fixed to this elevation. The reflecting boundaries (FastScape boundary condition set to 0) will not be affected, nor are the boundaries where a mass flux is specified.
Units: m
::::

::::{dropdown} __Parameter:__ {ref}`Flag to use orographic controls<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Flag_20to_20use_20orographic_20controls>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Flag_20to_20use_20orographic_20controls
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not to apply orographic controls.
::::

::::{dropdown} __Parameter:__ {ref}`Multi-direction slope exponent<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Multi_2ddirection_20slope_20exponent>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Multi_2ddirection_20slope_20exponent
**Default value:** 1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Exponent to determine the distribution from the SPL to neighbor nodes, with 10 being steepest decent and 1 being more varied.
::::

::::{dropdown} __Parameter:__ {ref}`Orographic elevation control<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Orographic_20elevation_20control>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Orographic_20elevation_20control
**Default value:** 2000

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Above this height, the elevation factor is applied. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Orographic wind barrier height<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Orographic_20wind_20barrier_20height>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Orographic_20wind_20barrier_20height
**Default value:** 500

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** When terrain reaches this height the wind barrier factor is applied. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Sediment deposition coefficient<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Sediment_20deposition_20coefficient>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Sediment_20deposition_20coefficient
**Default value:** -1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Deposition coefficient for sediment. A value smaller than 0 sets this to the same as the bedrock deposition coefficient.
::::

::::{dropdown} __Parameter:__ {ref}`Sediment diffusivity<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Sediment_20diffusivity>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Sediment_20diffusivity
**Default value:** -1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Transport coefficient (diffusivity) for sediment. -1 sets this to the bedrock diffusivity. Units: ${m^2/yr}$ if &ldquo;Use years instead of seconds&rdquo; is true; otherwise, the units are ${m^2/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Sediment river incision rate<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Sediment_20river_20incision_20rate>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Sediment_20river_20incision_20rate
**Default value:** -1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** River incision rate for sediment in the Stream Power Law. A value smaller than 0 sets this to the bedrock river incision rate. Units: $m^(1-2drainage_area_exponent)/yr}$ if &ldquo;Use years instead of seconds&rdquo; is true; otherwise, the units are $m^(1-2drainage_area_exponent)/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Slope exponent<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Slope_20exponent>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Slope_20exponent
**Default value:** 1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The slope exponent for the Stream Power Law (n). Generally m/n should equal approximately 0.4
::::

::::{dropdown} __Parameter:__ {ref}`Stack orographic controls<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Stack_20orographic_20controls>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Stack_20orographic_20controls
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether or not to apply both controls to a point, or only a maximum of one set as the wind barrier.
::::

::::{dropdown} __Parameter:__ {ref}`Use a fixed erosional base level<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Use_20a_20fixed_20erosional_20base_20level>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Use_20a_20fixed_20erosional_20base_20level
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not to use an erosional base level that differs from sea level. Setting this parameter to true will set all ghost nodes of fixed FastScape boundaries to the height you specify in &rsquo;set Erosional base level&rsquo;.
This can make sense for a continental model where the model surrounding topography is assumed above sea level, e.g. highlands. If the sea level would be used as an erosional base level in this case, all topography erodes away with lots of &rsquo;sediment volume&rsquo; lost through the sides of the model. This is mostly important, when there are mountains in the middle of the model, while it is less important when there is lower relief in the middle of the model.
In the FastScape  visualization files, setting the extra base level may show up as a strong slope at the fixed boundaries of the model. However, in the ASPECT visualization files it will not show up, as the ghost nodes only exist in FastScape.
::::

::::{dropdown} __Parameter:__ {ref}`Use kd distribution function<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Use_20kd_20distribution_20function>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Use_20kd_20distribution_20function
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to define Bedrock transport coefficient (diffusivity) using a distribution function. If false, a constant kd value will be used, which can be specified by setting the parameter &ldquo;Bedrock diffusivity&rdquo;. Units: ${m^2/yr}$ if &ldquo;Use years instead of seconds&rdquo; is true; otherwise, the units are ${m^2/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Use kf distribution function<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Use_20kf_20distribution_20function>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Use_20kf_20distribution_20function
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to define bedrock river incision rate using a distribution function. If false, a constant kf value will be used, which can be specified by setting the parameter &ldquo;Bedrock river incision rate&rdquo;. Units: ${m^(1-2drainage_area_exponent)/yr}$ if &ldquo;Use years instead of seconds&rdquo; is true; otherwise, the units are ${m^(1-2drainage_area_exponent)/s}$.
::::

::::{dropdown} __Parameter:__ {ref}`Wind barrier factor<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Wind_20barrier_20factor>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Wind_20barrier_20factor
**Default value:** 1

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Amount to multiply the bedrock river incision rate and transport coefficient by past given wind barrier height.
::::

::::{dropdown} __Parameter:__ {ref}`Wind direction<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Wind_20direction>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/Wind_20direction
**Default value:** west

**Pattern:** [Selection east|west|south|north ]

**Documentation:** This parameter assumes a wind direction, deciding which side is reduced from the wind barrier.
::::

(parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function)=
## **Subsection:** Mesh deformation / Fastscape / Erosional parameters / kd distribution function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function/Function_20constants>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function/Function_20expression>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function/Function_20expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function/Variable_20names>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kd_20distribution_20function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function)=
## **Subsection:** Mesh deformation / Fastscape / Erosional parameters / kf distribution function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function/Function_20constants>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function/Function_20expression>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function/Function_20expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function/Variable_20names>`
:name: parameters:Mesh_20deformation/Fastscape/Erosional_20parameters/kf_20distribution_20function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Mesh_20deformation/Fastscape/Marine_20parameters)=
## **Subsection:** Mesh deformation / Fastscape / Marine parameters
::::{dropdown} __Parameter:__ {ref}`Depth averaging thickness<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Depth_20averaging_20thickness>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Depth_20averaging_20thickness
**Default value:** 1e2

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Depth averaging for the sand-silt equation. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Sand e-folding depth<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sand_20e_2dfolding_20depth>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sand_20e_2dfolding_20depth
**Default value:** 1e3

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** E-folding depth for the exponential of the sand porosity law. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Sand porosity<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sand_20porosity>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sand_20porosity
**Default value:** 0.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Porosity of sand.
::::

::::{dropdown} __Parameter:__ {ref}`Sand transport coefficient<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sand_20transport_20coefficient>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sand_20transport_20coefficient
**Default value:** 5e2

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Transport coefficient (diffusivity) for sand. Units: ${m^2/yr}$
::::

::::{dropdown} __Parameter:__ {ref}`Sea level<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level
**Default value:** 0.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Constant sea level relative to the ASPECT surface, where the maximum Z or Y extent in ASPECT is a sea level of zero. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Silt e-folding depth<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20e_2dfolding_20depth>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20e_2dfolding_20depth
**Default value:** 1e3

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** E-folding depth for the exponential of the silt porosity law. Units: ${m}$
::::

::::{dropdown} __Parameter:__ {ref}`Silt fraction<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20fraction>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20fraction
**Default value:** 0.5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Fraction of silt for material leaving continent. Formerly called Sand-silt ratio.
::::

::::{dropdown} __Parameter:__ {ref}`Silt porosity<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20porosity>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20porosity
**Default value:** 0.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Porosity of silt.
::::

::::{dropdown} __Parameter:__ {ref}`Silt transport coefficient<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20transport_20coefficient>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Silt_20transport_20coefficient
**Default value:** 2.5e2

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Transport coefficient (diffusivity) for silt. Units: ${m^2/yr}$
::::

::::{dropdown} __Parameter:__ {ref}`Use sea level function<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Use_20sea_20level_20function>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Use_20sea_20level_20function
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to define sea level using a time-dependent function. If false, a constant value will be used.
::::

(parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function)=
## **Subsection:** Mesh deformation / Fastscape / Marine parameters / Sea level function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function/Function_20constants>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function/Function_20expression>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function/Function_20expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function/Variable_20names>`
:name: parameters:Mesh_20deformation/Fastscape/Marine_20parameters/Sea_20level_20function/Variable_20names
**Default value:** x,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Mesh_20deformation/Free_20surface)=
## **Subsection:** Mesh deformation / Free surface
::::{dropdown} __Parameter:__ {ref}`Free surface stabilization theta<parameters:Mesh_20deformation/Free_20surface/Free_20surface_20stabilization_20theta>`
:name: parameters:Mesh_20deformation/Free_20surface/Free_20surface_20stabilization_20theta
**Default value:** 0.5

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** Theta parameter described in {cite}`kaus:etal:2010`. An unstabilized free surface can overshoot its equilibrium position quite easily and generate unphysical results.  One solution is to use a quasi-implicit correction term to the forces near the free surface.  This parameter describes how much the free surface is stabilized with this term, where zero is no stabilization, and one is fully implicit.
::::

::::{dropdown} __Parameter:__ {ref}`Surface velocity projection<parameters:Mesh_20deformation/Free_20surface/Surface_20velocity_20projection>`
:name: parameters:Mesh_20deformation/Free_20surface/Surface_20velocity_20projection
**Default value:** normal

**Pattern:** [Selection normal|vertical ]

**Documentation:** After each time step the free surface must be advected in the direction of the velocity field. Mass conservation requires that the mesh velocity is in the normal direction of the surface. However, for steep topography or large curvature, advection in the normal direction can become ill-conditioned, and instabilities in the mesh can form. Projection of the mesh velocity onto the local vertical direction can preserve the mesh quality better, but at the cost of slightly poorer mass conservation of the domain.
::::

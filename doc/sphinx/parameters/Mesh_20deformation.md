(parameters:Mesh_20deformation)=
# Mesh deformation


## **Subsection:** Mesh deformation


(parameters:Mesh_20deformation/Additional_20tangential_20mesh_20velocity_20boundary_20indicators)=
### __Parameter name:__ Additional tangential mesh velocity boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries where there the mesh is allowed to move tangential to the boundary. All tangential mesh movements along those boundaries that have tangential material velocity boundary conditions are allowed by default, this parameters allows to generate mesh movements along other boundaries that are open, or have prescribed material velocities or tractions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

(parameters:Mesh_20deformation/Mesh_20deformation_20boundary_20indicators)=
### __Parameter name:__ Mesh deformation boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries where there the mesh is allowed to move according to the specified mesh deformation objects.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

The format is id1: object1 \& object2, id2: object3 \& object2, where objects are one of &lsquo;ascii data&rsquo;: Implementation of a model in which the initial mesh deformation (initial topography) is derived from a file containing data in ascii format. The following geometry models are currently supported: box, chunk, spherical. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;Topography [m]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Topography [m]&rsquo; in a 3d model, which means that there has to be a single column containing the topography. Note that the data in the input file needs to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle in radians  and &lsquo;y&rsquo; by the polar angle in radians measured positive from the north pole. The grid will be assumed to be a longitude-colatitude grid. Note that the order of spherical coordinates is &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;boundary function&rsquo;: A plugin, which prescribes the surface mesh to deform according to an analytically prescribed function. Note that the function prescribes a deformation velocity, i.e. the return value of this plugin is later multiplied by the time step length to compute the displacement increment in this time step. Although the function&rsquo;s time variable is interpreted as years when Use years in output instead of seconds is set to true, the boundary deformation velocity should still be given in m/s. The format of the functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;diffusion&rsquo;: A plugin that computes the deformation of surface vertices according to the solution of the hillslope diffusion problem. Specifically, at the end of each timestep, or after a specific number of timesteps, this plugin solves the following equation: \begin{align*}  \frac{\partial h}{\partial t} = \kappa \left( \frac{\partial^{2} h}{\partial x^{2}} + \frac{\partial^{2} h}{\partial y^{2}} \right), \end{align*} where $\kappa$ is the hillslope diffusion coefficient (diffusivity), and $h(x,y)$ the height of a point along the top boundary with respect to the surface of the unperturbed domain.

Using this definition, the plugin then solves for one time step, i.e., using as initial condition $h(t_{n-1})$ the current surface elevation, and computing $h(t_n)$ from it by solving the equation above over the time interval $t_n-t_{n-1}$. From this, one can then compute a surface velocity $v = \frac{h(t_n)-h(t_{n-1})}{t_n-t_{n-1}}$.

This surface velocity is used to deform the surface and as a boundary condition for solving the Laplace equation to determine the mesh velocity in the domain interior. Diffusion can be applied every timestep, mimicking surface processes of erosion and deposition, or at a user-defined timestep interval to purely smooth the surface topography to avoid too great a distortion of mesh elements when a free surface is also used.

&lsquo;free surface&rsquo;: A plugin that computes the deformation of surface vertices according to the solution of the flow problem. In particular this means if the surface of the domain is left open to flow, this flow will carry the mesh with it. The implementation was described in {cite}`rose_freesurface`, with the stabilization of the free surface originally described in {cite}`kaus:etal:2010`.

(parameters:Mesh_20deformation/Ascii_20data_20model)=
## **Subsection:** Mesh deformation / Ascii data model
(parameters:Mesh_20deformation/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Mesh_20deformation/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_3d_%s.0.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Mesh_20deformation/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Mesh_20deformation/Boundary_20function)=
## **Subsection:** Mesh deformation / Boundary function
(parameters:Mesh_20deformation/Boundary_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Mesh_20deformation/Boundary_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Mesh_20deformation/Boundary_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Mesh_20deformation/Diffusion)=
## **Subsection:** Mesh deformation / Diffusion
(parameters:Mesh_20deformation/Diffusion/Hillslope_20transport_20coefficient)=
### __Parameter name:__ Hillslope transport coefficient
**Default value:** 1e-6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The hillslope transport coefficient $\kappa$ used to diffuse the free surface, either as a  stabilization step or to mimic erosional and depositional processes. Units: $\si{m^2/s}$.

(parameters:Mesh_20deformation/Diffusion/Time_20steps_20between_20diffusion)=
### __Parameter name:__ Time steps between diffusion
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of time steps between each application of diffusion.

(parameters:Mesh_20deformation/Free_20surface)=
## **Subsection:** Mesh deformation / Free surface
(parameters:Mesh_20deformation/Free_20surface/Free_20surface_20stabilization_20theta)=
### __Parameter name:__ Free surface stabilization theta
**Default value:** 0.5

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** Theta parameter described in {cite}`kaus:etal:2010`. An unstabilized free surface can overshoot its equilibrium position quite easily and generate unphysical results.  One solution is to use a quasi-implicit correction term to the forces near the free surface.  This parameter describes how much the free surface is stabilized with this term, where zero is no stabilization, and one is fully implicit.

(parameters:Mesh_20deformation/Free_20surface/Surface_20velocity_20projection)=
### __Parameter name:__ Surface velocity projection
**Default value:** normal

**Pattern:** [Selection normal|vertical ]

**Documentation:** After each time step the free surface must be advected in the direction of the velocity field. Mass conservation requires that the mesh velocity is in the normal direction of the surface. However, for steep topography or large curvature, advection in the normal direction can become ill-conditioned, and instabilities in the mesh can form. Projection of the mesh velocity onto the local vertical direction can preserve the mesh quality better, but at the cost of slightly poorer mass conservation of the domain.

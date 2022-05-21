# Parameters


## **Parameters in section** 


(parameters:Additional_20shared_20libraries)=
### __Parameter name:__ Additional shared libraries
**Default value:**  

**Pattern:** [List of <[FileName (Type: input)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of names of additional shared libraries that should be loaded upon starting up the program. The names of these files can contain absolute or relative paths (relative to the directory in which you call ASPECT). In fact, file names that do not contain any directory information (i.e., only the name of a file such as <myplugin.so> will not be found if they are not located in one of the directories listed in the \texttt{LD_LIBRARY_PATH} environment variable. In order to load a library in the current directory, use <./myplugin.so> instead.

The typical use of this parameter is so that you can implement additional plugins in your own directories, rather than in the ASPECT source directories. You can then simply compile these plugins into a shared library without having to re-compile all of ASPECT. See the section of the manual discussing writing extensions for more information on how to compile additional files into a shared library. 

(parameters:Adiabatic_20surface_20temperature)=
### __Parameter name:__ Adiabatic surface temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** In order to make the problem in the first time step easier to solve, we need a reasonable guess for the temperature and pressure. To obtain it, we use an adiabatic pressure and temperature field. This parameter describes what the `adiabatic' temperature would be at the surface of the domain (i.e. at depth zero). Note that this value need not coincide with the boundary condition posed at this point. Rather, the boundary condition may differ significantly from the adiabatic value, and then typically induce a thermal boundary layer.

For more information, see the section in the manual that discusses the general mathematical model. 

(parameters:CFL_20number)=
### __Parameter name:__ CFL number
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** In computations, the time step $k$ is chosen according to $k = c \min_K \frac {h_K} {\|u\|_{\infty,K} p_T}$ where $h_K$ is the diameter of cell $K$, and the denominator is the maximal magnitude of the velocity on cell $K$ times the polynomial degree $p_T$ of the temperature discretization. The dimensionless constant $c$ is called the CFL number in this program. For time discretizations that have explicit components, $c$ must be less than a constant that depends on the details of the time discretization and that is no larger than one. On the other hand, for implicit discretizations such as the one chosen here, one can choose the time step as large as one wants (in particular, one can choose $c>1$) though a CFL number significantly larger than one will yield rather diffusive solutions. Units: None. 

(parameters:Dimension)=
### __Parameter name:__ Dimension
**Default value:** 2 

**Pattern:** [Integer range 2...3 (inclusive)] 

**Documentation:** The number of space dimensions you want to run this program in. ASPECT can run in 2 and 3 space dimensions. 

(parameters:End_20time)=
### __Parameter name:__ End time
**Default value:** 5.69e+300 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The end time of the simulation. The default value is a number so that when converted from years to seconds it is approximately equal to the largest number representable in floating point arithmetic. For all practical purposes, this equals infinity. Units: Years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Max_20nonlinear_20iterations)=
### __Parameter name:__ Max nonlinear iterations
**Default value:** 10 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The maximal number of nonlinear iterations to be performed. 

(parameters:Max_20nonlinear_20iterations_20in_20pre_2drefinement)=
### __Parameter name:__ Max nonlinear iterations in pre_2drefinement
**Default value:** 2147483647 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The maximal number of nonlinear iterations to be performed in the pre-refinement steps. This does not include the last refinement step before moving to timestep 1. When this parameter has a larger value than max nonlinear iterations, the latter is used. 

(parameters:Maximum_20first_20time_20step)=
### __Parameter name:__ Maximum first time step
**Default value:** 5.69e+300 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Set a maximum time step size for only the first timestep. Generally the time step based on the CFL number should be sufficient, but for complicated models or benchmarking it may be useful to limit the first time step to some value, especially when using the free surface, which needs to settle to prevent instabilities. This should in that case be combined with a value set for ``Maximum relative increase in time step''. The default value is a value so that when converted from years into seconds it equals the largest number representable by a floating point number, implying an unlimited time step. Units: Years or seconds, depending on the ``Use years in output instead of seconds'' parameter. 

(parameters:Maximum_20relative_20increase_20in_20time_20step)=
### __Parameter name:__ Maximum relative increase in time step
**Default value:** 2147483647 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Set a percentage with which the time step is limited to increase. Generally the time step based on the CFL number should be sufficient, but for complicated models which may suddenly drastically change behavior, it may be useful to limit the increase in the time step, without limiting the time step size of the whole simulation to a particular number. For example, if this parameter is set to $50$, then that means that the time step can at most increase by 50\% from one time step to the next, or by a factor of 1.5. Units: \%. 

(parameters:Maximum_20time_20step)=
### __Parameter name:__ Maximum time step
**Default value:** 5.69e+300 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Set a maximum time step size for the solver to use. Generally the time step based on the CFL number should be sufficient, but for complicated models or benchmarking it may be useful to limit the time step to some value. The default value is a value so that when converted from years into seconds it equals the largest number representable by a floating point number, implying an unlimited time step.Units: Years or seconds, depending on the ``Use years in output instead of seconds'' parameter. 

(parameters:Nonlinear_20solver_20scheme)=
### __Parameter name:__ Nonlinear solver scheme
**Default value:** single Advection, single Stokes 

**Pattern:** [Selection single Advection, single Stokes|iterated Advection and Stokes|single Advection, iterated Stokes|no Advection, iterated Stokes|no Advection, single Stokes|no Advection, iterated defect correction Stokes|single Advection, iterated defect correction Stokes|iterated Advection and defect correction Stokes|iterated Advection and Newton Stokes|single Advection, iterated Newton Stokes|single Advection, no Stokes|IMPES|iterated IMPES|iterated Stokes|Newton Stokes|Stokes only|Advection only|first timestep only, single Stokes|no Advection, no Stokes ] 

**Documentation:** The kind of scheme used to resolve the nonlinearity in the system. `single Advection, single Stokes' means that no nonlinear iterations are done, and the temperature, compositional fields and Stokes equations are solved exactly once per time step, one after the other. The `iterated Advection and Stokes' scheme iterates this decoupled approach by alternating the solution of the temperature, composition and Stokes systems. The `single Advection, iterated Stokes' scheme solves the temperature and composition equation once at the beginning of each time step and then iterates out the solution of the Stokes equation. The `no Advection, iterated Stokes' scheme only solves the Stokes system, iterating out the solution, and ignores compositions and the temperature equation (careful, the material model must not depend on the temperature or composition; this is mostly useful for Stokes benchmarks).  The `no Advection, single Stokes' scheme only solves the Stokes system once per timestep. This is also mostly useful for Stokes benchmarks. The `single Advection, no Stokes' scheme only solves the temperature and other advection systems once, and instead of solving for the Stokes system, a prescribed velocity and pressure is used. The `iterated Advection and Newton Stokes' scheme iterates by alternating the solution of the temperature, composition and Stokes equations, using Picard iterations for the temperature and composition, and Newton iterations for the Stokes system. The `single Advection, iterated Newton Stokes' scheme solves the temperature and composition equations once at the beginning of each time step and then iterates out the solution of the Stokes equation, using Newton iterations for the Stokes system. The `iterated Advection and defect correction Stokes' scheme iterates by alternating the solution of the temperature, composition and Stokes equations, using Picard iterations for the temperature and composition, and defect correction Picard iterations for the Stokes system. The `single Advection, iterated defect correction Stokes' scheme solves the temperature and composition equations once at the beginning of each time step and then iterates out the solution of the Stokes equation, using defect correction Picard iterations for the Stokes system. The `no Advection, iterated defect correction Stokes' scheme solves the temperature and composition equations once at the beginning of each time step and then iterates out the solution of the Stokes equation, using defect correction Picard iterations for the Stokes system. The `first timestep only, single Stokes' scheme solves the Stokes equations exactly once, at the first time step. No nonlinear iterations are done, and the temperature and composition systems are not solved. 

The `IMPES' scheme is deprecated and only allowed for reasons of backwards compatibility. It is the same as `single Advection, single Stokes' .The `iterated IMPES' scheme is deprecated and only allowed for reasons of backwards compatibility. It is the same as `iterated Advection and Stokes'. The `iterated Stokes' scheme is deprecated and only allowed for reasons of backwards compatibility. It is the same as `single Advection, iterated Stokes'. The `Stokes only' scheme is deprecated and only allowed for reasons of backwards compatibility. It is the same as `no Advection, iterated Stokes'. The `Advection only' scheme is deprecated and only allowed for reasons of backwards compatibility. It is the same as `single Advection, no Stokes'. The `Newton Stokes' scheme is deprecated and only allowed for reasons of backwards compatibility. It is the same as `iterated Advection and Newton Stokes'. 

(parameters:Nonlinear_20solver_20tolerance)=
### __Parameter name:__ Nonlinear solver tolerance
**Default value:** 1e-5 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** A relative tolerance up to which the nonlinear solver will iterate. This parameter is only relevant if the `Nonlinear solver scheme' does nonlinear iterations, in other words, if it is set to something other than `single Advection, single Stokes' or `single Advection, no Stokes'. 

(parameters:Output_20directory)=
### __Parameter name:__ Output directory
**Default value:** output 

**Pattern:** [DirectoryName] 

**Documentation:** The name of the directory into which all output files should be placed. This may be an absolute or a relative path. 

(parameters:Pressure_20normalization)=
### __Parameter name:__ Pressure normalization
**Default value:** surface 

**Pattern:** [Selection surface|volume|no ] 

**Documentation:** If and how to normalize the pressure after the solution step. This is necessary because depending on boundary conditions, in many cases the pressure is only determined by the model up to a constant. On the other hand, we often would like to have a well-determined pressure, for example for table lookups of material properties in models or for comparing solutions. If the given value is `surface', then normalization at the end of each time steps adds a constant value to the pressure in such a way that the average pressure at the surface of the domain is what is set in the `Surface pressure' parameter; the surface of the domain is determined by asking the geometry model whether a particular face of the geometry has a zero or small `depth'. If the value of this parameter is `volume' then the pressure is normalized so that the domain average is zero. If `no' is given, the no pressure normalization is performed. 

(parameters:Resume_20computation)=
### __Parameter name:__ Resume computation
**Default value:** false 

**Pattern:** [Selection true|false|auto ] 

**Documentation:** A flag indicating whether the computation should be resumed from a previously saved state (if true) or start from scratch (if false). If auto is selected, models will be resumed if there is an existing checkpoint file, otherwise started from scratch. 

(parameters:Start_20time)=
### __Parameter name:__ Start time
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The start time of the simulation. Units: Years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Surface_20pressure)=
### __Parameter name:__ Surface pressure
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The value the pressure is normalized to in each time step when `Pressure normalization' is set to `surface' with default value 0. This setting is ignored in all other cases.

The mathematical equations that describe thermal convection only determine the pressure up to an arbitrary constant. On the other hand, for comparison and for looking up material parameters it is important that the pressure be normalized somehow. We do this by enforcing a particular average pressure value at the surface of the domain, where the geometry model determines where the surface is. This parameter describes what this average surface pressure value is supposed to be. By default, it is set to zero, but one may want to choose a different value for example for simulating only the volume of the mantle below the lithosphere, in which case the surface pressure should be the lithostatic pressure at the bottom of the lithosphere.

For more information, see the section in the manual that discusses the general mathematical model. 

(parameters:Timing_20output_20frequency)=
### __Parameter name:__ Timing output frequency
**Default value:** 100 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** How frequently in timesteps to output timing information. This is generally adjusted only for debugging and timing purposes. If the value is set to zero it will also output timing information at the initiation timesteps. 

(parameters:Use_20conduction_20timestep)=
### __Parameter name:__ Use conduction timestep
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Mantle convection simulations are often focused on convection dominated systems. However, these codes can also be used to investigate systems where heat conduction plays a dominant role. This parameter indicates whether the simulator should also use heat conduction in determining the length of each time step. 

(parameters:Use_20operator_20splitting)=
### __Parameter name:__ Use operator splitting
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If set to true, the advection and reactions of compositional fields and temperature are solved separately, and can use different time steps. Note that this will only work if the material/heating model fills the reaction\_rates/heating\_reaction\_rates structures. Operator splitting can be used with any existing solver schemes that solve the temperature/composition equations. 

(parameters:Use_20years_20in_20output_20instead_20of_20seconds)=
### __Parameter name:__ Use years in output instead of seconds
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** When computing results for mantle convection simulations, it is often difficult to judge the order of magnitude of results when they are stated in MKS units involving seconds. Rather, some kinds of results such as velocities are often stated in terms of meters per year (or, sometimes, centimeters per year). On the other hand, for non-dimensional computations, one wants results in their natural unit system as used inside the code. If this flag is set to `true' conversion to years happens; if it is `false', no such conversion happens.

Contrary to the word ``output'' in the name of this parameter, a number of plugins also use this parameter to determine how to interpret their \textit{inputs}. For example, when `true', several of the boundary velocity models described in Section~\ref{parameters:Boundary_20velocity_20model} interpret both specific times in years instead of seconds, and velocities in meters per year instead of meters per second. 

(parameters:World_20builder_20file)=
### __Parameter name:__ World builder file
**Default value:**  

**Pattern:** [FileName (Type: input)] 

**Documentation:** Name of the world builder file. If empty, the world builder is not initialized. 

(parameters:Adiabatic_20conditions_20model)=
## **Parameters in section** Adiabatic conditions model
(parameters:Adiabatic_20conditions_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** compute profile 

**Pattern:** [Selection ascii data|compute profile|function ] 

**Documentation:** Select one of the following models:

`ascii data': A model in which the adiabatic profile is read from a file that describes the reference state. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of points in the reference state as for example `# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide columns named `temperature', `pressure', and `density'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

`compute profile': A model in which the adiabatic profile is calculated by solving the hydrostatic equations for pressure and temperature in depth. The gravity is assumed to be in depth direction and the composition is either given by the initial composition at reference points or computed as a reference depth-function. All material parameters are computed by the material model plugin. The surface conditions are either constant or changing over time as prescribed by a user-provided function.

`function': A model in which the adiabatic profile is specified by a user defined function. The supplied function has to contain temperature, pressure, and density as a function of depth in this order. 

(parameters:Adiabatic_20conditions_20model:Ascii_20data_20model)=
## **Parameters in section** Adiabatic conditions model/Ascii data model
(parameters:Adiabatic_20conditions_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/tests/adiabatic-conditions/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Adiabatic_20conditions_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Adiabatic_20conditions_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile)=
## **Parameters in section** Adiabatic conditions model/Compute profile
(parameters:Adiabatic_20conditions_20model:Compute_20profile:Composition_20reference_20profile)=
### __Parameter name:__ Composition reference profile
**Default value:** initial composition 

**Pattern:** [Selection initial composition|function ] 

**Documentation:** Select how the reference profile for composition is computed. This profile is used to evaluate the material model, when computing the pressure and temperature profile. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Number_20of_20points)=
### __Parameter name:__ Number of points
**Default value:** 1000 

**Pattern:** [Integer range 5...2147483647 (inclusive)] 

**Documentation:** The number of points we use to compute the adiabatic profile. The higher the number of points, the more accurate the downward integration from the adiabatic surface temperature will be. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Use_20surface_20condition_20function)=
### __Parameter name:__ Use surface condition function
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use the 'Surface condition function' to determine surface conditions, or the 'Adiabatic surface temperature' and 'Surface pressure' parameters. If this is set to true the reference profile is updated every timestep. The function expression of the function should be independent of space, but can depend on time 't'. The function must return two components, the first one being reference surface pressure, the second one being reference surface temperature. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Surface_20condition_20function)=
## **Parameters in section** Adiabatic conditions model/Compute profile/Surface condition function
(parameters:Adiabatic_20conditions_20model:Compute_20profile:Surface_20condition_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Surface_20condition_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Adiabatic_20conditions_20model:Compute_20profile:Surface_20condition_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Adiabatic_20conditions_20model:Function)=
## **Parameters in section** Adiabatic conditions model/Function
(parameters:Adiabatic_20conditions_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Adiabatic_20conditions_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0.0; 0.0; 1.0 

**Pattern:** [Anything] 

**Documentation:** Expression for the adiabatic temperature, pressure, and density separated by semicolons as a function of `depth'. 

(parameters:Adiabatic_20conditions_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** depth 

**Pattern:** [Anything] 

**Documentation:**  

(parameters:Boundary_20composition_20model)=
## **Parameters in section** Boundary composition model
(parameters:Boundary_20composition_20model:Allow_20fixed_20composition_20on_20outflow_20boundaries)=
### __Parameter name:__ Allow fixed composition on outflow boundaries
**Default value:** false for models without melt 

**Pattern:** [Selection true|false|false for models without melt ] 

**Documentation:** When the composition is fixed on a given boundary as determined by the list of 'Fixed composition boundary indicators', there might be parts of the boundary where material flows out and one may want to prescribe the composition only on those parts of the boundary where there is inflow. This parameter determines if compositions are only prescribed at these inflow parts of the boundary (if false) or everywhere on a given boundary, independent of the flow direction (if true). By default, this parameter is set to false, except in models with melt transport (see below). Note that in this context, `fixed' refers to the fact that these are the boundary indicators where Dirichlet boundary conditions are applied, and does not imply that the boundary composition is time-independent. 

Mathematically speaking, the compositional fields satisfy an advection equation that has no diffusion. For this equation, one can only impose Dirichlet boundary conditions (i.e., prescribe a fixed compositional field value at the boundary) at those boundaries where material flows in. This would correspond to the ``false'' setting of this parameter, which is correspondingly the default. On the other hand, on a finite dimensional discretization such as the one one obtains from the finite element method, it is possible to also prescribe values on outflow boundaries, even though this may make no physical sense. This would then correspond to the ``true'' setting of this parameter.

A warning for models with melt transport: In models with fluid flow, some compositional fields (in particular the porosity) might be transported with the fluid velocity, and would need to set the constraints based on the fluid velocity. However, this is currently not possible, because we reuse the same matrix for all compositional fields, and therefore can not use different constraints for different fields. Consequently, we set this parameter to true by default in models where melt transport is enabled. Be aware that if you change this default setting, you will not use the melt velocity, but the solid velocity to determine on which parts of the boundaries there is outflow. 

(parameters:Boundary_20composition_20model:Fixed_20composition_20boundary_20indicators)=
### __Parameter name:__ Fixed composition boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries on which the composition is fixed and described by the boundary composition object selected in its own section of this input file. All boundary indicators used by the geometry but not explicitly listed here will end up with no-flux (insulating) boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed composition, but not what composition should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryComposition group, unless an existing implementation in this group already provides what you want. 

(parameters:Boundary_20composition_20model:List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**  

**Pattern:** [MultipleSelection ascii data|box|box with lithosphere boundary indicators|function|initial composition|spherical constant ] 

**Documentation:** A comma-separated list of boundary composition models that will be used to initialize the composition. These plugins are loaded in the order given, and modify the existing composition field via the operators listed in 'List of model operators'.

The following boundary composition models are available:

`ascii data': Implementation of a model in which the boundary composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `composition1', `composition2', etc. in a 2d model and `x', `y', `composition1', `composition2', etc., in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates.If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`box': A model in which the composition is chosen constant on the sides of a box which are selected by the parameters Left/Right/Top/Bottom/Front/Back composition

`box with lithosphere boundary indicators': A model in which the composition is chosen constant on all the sides of a box. Additional boundary indicators are added to the lithospheric parts of the vertical boundaries. This model is to be used with the 'Two Merged Boxes' Geometry Model.

`function': Implementation of a model in which the boundary composition is given in terms of an explicit formula that is elaborated in the parameters in section ``Boundary composition model|Function''. 

Since the symbol $t$ indicating time may appear in the formulas for the prescribed composition, it is interpreted as having units seconds unless the global input parameter ``Use years in output instead of seconds'' is set, in which case we interpret the formula expressions as having units year.

The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`initial composition': A model in which the composition at the boundary is chosen to be the same as given in the initial conditions.

Because this class simply takes what the initial composition had described, this class can not know certain pieces of information such as the minimal and maximal composition on the boundary. For operations that require this, for example in post-processing, this boundary composition model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section ``Boundary composition model/Initial composition''.

`spherical constant': A model in which the composition is chosen constant on the inner and outer boundaries of a surface, spherical shell, chunk or ellipsoidal chunk. Parameters are read from subsection 'Spherical constant'. 

(parameters:Boundary_20composition_20model:List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add 

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ] 

**Documentation:** A comma-separated list of operators that will be used to append the listed composition models onto the previous models. If only one operator is given, the same operator is applied to all models. 

(parameters:Boundary_20composition_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection ascii data|box|box with lithosphere boundary indicators|function|initial composition|spherical constant|unspecified ] 

**Documentation:** Select one of the following models:

`ascii data': Implementation of a model in which the boundary composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `composition1', `composition2', etc. in a 2d model and `x', `y', `composition1', `composition2', etc., in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates.If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`box': A model in which the composition is chosen constant on the sides of a box which are selected by the parameters Left/Right/Top/Bottom/Front/Back composition

`box with lithosphere boundary indicators': A model in which the composition is chosen constant on all the sides of a box. Additional boundary indicators are added to the lithospheric parts of the vertical boundaries. This model is to be used with the 'Two Merged Boxes' Geometry Model.

`function': Implementation of a model in which the boundary composition is given in terms of an explicit formula that is elaborated in the parameters in section ``Boundary composition model|Function''. 

Since the symbol $t$ indicating time may appear in the formulas for the prescribed composition, it is interpreted as having units seconds unless the global input parameter ``Use years in output instead of seconds'' is set, in which case we interpret the formula expressions as having units year.

The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`initial composition': A model in which the composition at the boundary is chosen to be the same as given in the initial conditions.

Because this class simply takes what the initial composition had described, this class can not know certain pieces of information such as the minimal and maximal composition on the boundary. For operations that require this, for example in post-processing, this boundary composition model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section ``Boundary composition model/Initial composition''.

`spherical constant': A model in which the composition is chosen constant on the inner and outer boundaries of a surface, spherical shell, chunk or ellipsoidal chunk. Parameters are read from subsection 'Spherical constant'.

\textbf{Warning}: This parameter provides an old and deprecated way of specifying boundary composition models and shouldn't be used. Please use 'List of model names' instead. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model)=
## **Parameters in section** Boundary composition model/Ascii data model
(parameters:Boundary_20composition_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-composition/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model:Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Time step between following data files. Depending on the setting of the global `Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model:Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. `Ma BP'). If this flag is set to `True' the plugin will first load the file with the number `First data file number' and decrease the file number during the model run. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model:First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The `First data file model time' parameter has been deactivated and will be removed in a future release. Do not use this paramter and instead provide data files starting from the model start time. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model:First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than `First velocity file model time'. 

(parameters:Boundary_20composition_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Boundary_20composition_20model:Box)=
## **Parameters in section** Boundary composition model/Box
(parameters:Boundary_20composition_20model:Box:Bottom_20composition)=
### __Parameter name:__ Bottom composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the bottom boundary (at minimal $y$-value in 2d, or minimal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box:Left_20composition)=
### __Parameter name:__ Left composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the left boundary (at minimal $x$-value). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box:Right_20composition)=
### __Parameter name:__ Right composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the right boundary (at maximal $x$-value). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box:Top_20composition)=
### __Parameter name:__ Top composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the top boundary (at maximal $y$-value in 2d, or maximal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators)=
## **Parameters in section** Boundary composition model/Box with lithosphere boundary indicators
(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators:Bottom_20composition)=
### __Parameter name:__ Bottom composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the bottom boundary (at minimal $y$-value in 2d, or minimal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators:Left_20composition)=
### __Parameter name:__ Left composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the left boundary (at minimal $x$-value). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators:Left_20composition_20lithosphere)=
### __Parameter name:__ Left composition lithosphere
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the left boundary (at minimal $x$-value). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators:Right_20composition)=
### __Parameter name:__ Right composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the right boundary (at maximal $x$-value). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators:Right_20composition_20lithosphere)=
### __Parameter name:__ Right composition lithosphere
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the right boundary (at maximal $x$-value). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Box_20with_20lithosphere_20boundary_20indicators:Top_20composition)=
### __Parameter name:__ Top composition
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of composition boundary values at the top boundary (at maximal $y$-value in 2d, or maximal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none. 

(parameters:Boundary_20composition_20model:Function)=
## **Parameters in section** Boundary composition model/Function
(parameters:Boundary_20composition_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are 'cartesian', 'spherical', and 'depth'. 'spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. 'depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Boundary_20composition_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Boundary_20composition_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Boundary_20composition_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Boundary_20composition_20model:Initial_20composition)=
## **Parameters in section** Boundary composition model/Initial composition
(parameters:Boundary_20composition_20model:Initial_20composition:Maximal_20composition)=
### __Parameter name:__ Maximal composition
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximal composition. Units: none. 

(parameters:Boundary_20composition_20model:Initial_20composition:Minimal_20composition)=
### __Parameter name:__ Minimal composition
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimal composition. Units: none. 

(parameters:Boundary_20composition_20model:Spherical_20constant)=
## **Parameters in section** Boundary composition model/Spherical constant
(parameters:Boundary_20composition_20model:Spherical_20constant:Inner_20composition)=
### __Parameter name:__ Inner composition
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Composition at the inner boundary (core mantle boundary). Units: none. 

(parameters:Boundary_20composition_20model:Spherical_20constant:Outer_20composition)=
### __Parameter name:__ Outer composition
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Composition at the outer boundary (lithosphere water/air). For a spherical geometry model, this is the only boundary. Units: none. 

(parameters:Boundary_20fluid_20pressure_20model)=
## **Parameters in section** Boundary fluid pressure model
(parameters:Boundary_20fluid_20pressure_20model:Plugin_20name)=
### __Parameter name:__ Plugin name
**Default value:** density 

**Pattern:** [Selection density ] 

**Documentation:** Select one of the following plugins:

`density': A plugin that prescribes the fluid pressure gradient at the boundary based on fluid/solid density from the material model. 

(parameters:Boundary_20fluid_20pressure_20model:Density)=
## **Parameters in section** Boundary fluid pressure model/Density
(parameters:Boundary_20fluid_20pressure_20model:Density:Density_20formulation)=
### __Parameter name:__ Density formulation
**Default value:** solid density 

**Pattern:** [Selection solid density|fluid density|average density ] 

**Documentation:** The density formulation used to compute the fluid pressure gradient at the model boundary.

`solid density' prescribes the gradient of the fluid pressure as solid density times gravity (which is the lithostatic pressure) and leads to approximately the same pressure in the melt as in the solid, so that fluid is only flowing in or out due to differences in dynamic pressure.

`fluid density' prescribes the gradient of the fluid pressure as fluid density times gravity and causes melt to flow in with the same velocity as inflowing solid material, or no melt flowing in or out if the solid velocity normal to the boundary is zero.

'average density' prescribes the gradient of the fluid pressure as the averaged fluid and solid density times gravity (which is a better approximation for the lithostatic pressure than just the solid density) and leads to approximately the same pressure in the melt as in the solid, so that fluid is only flowing in or out due to differences in dynamic pressure. 

(parameters:Boundary_20heat_20flux_20model)=
## **Parameters in section** Boundary heat flux model
(parameters:Boundary_20heat_20flux_20model:Fixed_20heat_20flux_20boundary_20indicators)=
### __Parameter name:__ Fixed heat flux boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries on which the heat flux is fixed and described by the boundary heat flux object selected in the 'Model name' parameter. All boundary indicators used by the geometry but not explicitly listed here or in the list of 'Fixed temperature boundary indicators' in the 'Boundary temperature model' will end up with no-flux (insulating) boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed heat flux, but not what heat flux should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryHeatFlux group, unless an existing implementation in this group already provides what you want. 

(parameters:Boundary_20heat_20flux_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** function 

**Pattern:** [Selection function ] 

**Documentation:** Select one of the following plugins:

`function': Implementation of a model in which the boundary heat flux is given in terms of an explicit formula that is elaborated in the parameters in section ``Boundary heat flux model|Function''. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

The formula you describe in the mentioned section is a scalar value for the heat flux that is assumed to be the flux normal to the boundary, and that has the unit W/(m$^2$) (in 3d) or W/m (in 2d). Negative fluxes are interpreted as the flow of heat into the domain, and positive fluxes are interpreted as heat flowing out of the domain.

The symbol $t$ indicating time that may appear in the formulas for the prescribed heat flux is interpreted as having units seconds unless the global parameter ``Use years in output instead of seconds'' has been set. 

(parameters:Boundary_20heat_20flux_20model:Function)=
## **Parameters in section** Boundary heat flux model/Function
(parameters:Boundary_20heat_20flux_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Boundary_20heat_20flux_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Boundary_20heat_20flux_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Boundary_20heat_20flux_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Boundary_20temperature_20model)=
## **Parameters in section** Boundary temperature model
(parameters:Boundary_20temperature_20model:Allow_20fixed_20temperature_20on_20outflow_20boundaries)=
### __Parameter name:__ Allow fixed temperature on outflow boundaries
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** When the temperature is fixed on a given boundary as determined by the list of 'Fixed temperature boundary indicators', there might be parts of the boundary where material flows out and one may want to prescribe the temperature only on the parts of the boundary where there is inflow. This parameter determines if temperatures are only prescribed at these inflow parts of the boundary (if false) or everywhere on a given boundary, independent of the flow direction (if true).Note that in this context, `fixed' refers to the fact that these are the boundary indicators where Dirichlet boundary conditions are applied, and does not imply that the boundary temperature is time-independent. 

Mathematically speaking, the temperature satisfies an advection-diffusion equation. For this type of equation, one can prescribe the temperature even on outflow boundaries as long as the diffusion coefficient is nonzero. This would correspond to the ``true'' setting of this parameter, which is correspondingly the default. In practice, however, this would only make physical sense if the diffusion coefficient is actually quite large to prevent the creation of a boundary layer. In addition, if there is no diffusion, one can only impose Dirichlet boundary conditions (i.e., prescribe a fixed temperature value at the boundary) at those boundaries where material flows in. This would correspond to the ``false'' setting of this parameter. 

(parameters:Boundary_20temperature_20model:Fixed_20temperature_20boundary_20indicators)=
### __Parameter name:__ Fixed temperature boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries on which the temperature is fixed and described by the boundary temperature object selected in the 'List of model names' parameter. All boundary indicators used by the geometry but not explicitly listed here will end up with no-flux (insulating) boundary conditions, or, if they are listed in the 'Fixed heat flux boundary indicators', with Neumann boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed temperature, but not what temperature should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryTemperature group, unless an existing implementation in this group already provides what you want. 

(parameters:Boundary_20temperature_20model:List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**  

**Pattern:** [MultipleSelection ascii data|box|box with lithosphere boundary indicators|constant|dynamic core|function|initial temperature|spherical constant ] 

**Documentation:** A comma-separated list of boundary temperature models that will be used to initialize the temperature. These plugins are loaded in the order given, and modify the existing temperature field via the operators listed in 'List of model operators'.

The following boundary temperature models are available:

`ascii data': Implementation of a model in which the boundary data is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `Temperature [K]' in a 2d model and  `x', `y', `Temperature [K]' in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`box': A model in which the temperature is chosen constant on the sides of a box which are selected by the parameters Left/Right/Top/Bottom/Front/Back temperature

`box with lithosphere boundary indicators': A model in which the temperature is chosen constant on all the sides of a box. Additional boundary indicators are added to the lithospheric parts of the vertical boundaries. This model is to be used with the 'Two Merged Boxes' Geometry Model.

`constant': A model in which the temperature is chosen constant on a given boundary indicator.  Parameters are read from the subsection 'Constant'.

`dynamic core': This is a boundary temperature model working only with spherical shell geometry and core statistics postprocessor. The temperature at the top is constant, and the core mantle boundary temperature is dynamically evolving through time by calculating the heat flux into the core and solving the core energy balance. The formulation is mainly following \cite{NPB+04}, and the plugin is used in Zhang et al. [2016]. The energy of core cooling and freeing of the inner core is included in the plugin. However, current plugin can not deal with the energy balance if the core is in the `snowing core' regime (i.e., the core solidifies from the top instead of bottom).

`function': Implementation of a model in which the boundary temperature is given in terms of an explicit formula that is elaborated in the parameters in section ``Boundary temperature model|Function''. 

Since the symbol $t$ indicating time may appear in the formulas for the prescribed temperatures, it is interpreted as having units seconds unless the global input parameter ``Use years in output instead of seconds'' is set, in which case we interpret the formula expressions as having units year.

Because this class simply takes what the function calculates, this class can not know certain pieces of information such as the minimal and maximal temperature on the boundary. For operations that require this, for example in post-processing, this boundary temperature model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section ``Boundary temperature model/Initial temperature''.

The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`initial temperature': A model in which the temperature at the boundary is chosen to be the same as given in the initial conditions.

Because this class simply takes what the initial temperature had described, this class can not know certain pieces of information such as the minimal and maximal temperature on the boundary. For operations that require this, for example in post-processing, this boundary temperature model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section ``Boundary temperature model/Initial temperature''.

`spherical constant': A model in which the temperature is chosen constant on the inner and outer boundaries of a spherical shell, ellipsoidal chunk or chunk. Parameters are read from subsection 'Spherical constant'. 

(parameters:Boundary_20temperature_20model:List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add 

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ] 

**Documentation:** A comma-separated list of operators that will be used to append the listed temperature models onto the previous models. If only one operator is given, the same operator is applied to all models. 

(parameters:Boundary_20temperature_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection ascii data|box|box with lithosphere boundary indicators|constant|dynamic core|function|initial temperature|spherical constant|unspecified ] 

**Documentation:** Select one of the following models:

`ascii data': Implementation of a model in which the boundary data is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `Temperature [K]' in a 2d model and  `x', `y', `Temperature [K]' in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`box': A model in which the temperature is chosen constant on the sides of a box which are selected by the parameters Left/Right/Top/Bottom/Front/Back temperature

`box with lithosphere boundary indicators': A model in which the temperature is chosen constant on all the sides of a box. Additional boundary indicators are added to the lithospheric parts of the vertical boundaries. This model is to be used with the 'Two Merged Boxes' Geometry Model.

`constant': A model in which the temperature is chosen constant on a given boundary indicator.  Parameters are read from the subsection 'Constant'.

`dynamic core': This is a boundary temperature model working only with spherical shell geometry and core statistics postprocessor. The temperature at the top is constant, and the core mantle boundary temperature is dynamically evolving through time by calculating the heat flux into the core and solving the core energy balance. The formulation is mainly following \cite{NPB+04}, and the plugin is used in Zhang et al. [2016]. The energy of core cooling and freeing of the inner core is included in the plugin. However, current plugin can not deal with the energy balance if the core is in the `snowing core' regime (i.e., the core solidifies from the top instead of bottom).

`function': Implementation of a model in which the boundary temperature is given in terms of an explicit formula that is elaborated in the parameters in section ``Boundary temperature model|Function''. 

Since the symbol $t$ indicating time may appear in the formulas for the prescribed temperatures, it is interpreted as having units seconds unless the global input parameter ``Use years in output instead of seconds'' is set, in which case we interpret the formula expressions as having units year.

Because this class simply takes what the function calculates, this class can not know certain pieces of information such as the minimal and maximal temperature on the boundary. For operations that require this, for example in post-processing, this boundary temperature model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section ``Boundary temperature model/Initial temperature''.

The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`initial temperature': A model in which the temperature at the boundary is chosen to be the same as given in the initial conditions.

Because this class simply takes what the initial temperature had described, this class can not know certain pieces of information such as the minimal and maximal temperature on the boundary. For operations that require this, for example in post-processing, this boundary temperature model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section ``Boundary temperature model/Initial temperature''.

`spherical constant': A model in which the temperature is chosen constant on the inner and outer boundaries of a spherical shell, ellipsoidal chunk or chunk. Parameters are read from subsection 'Spherical constant'.

\textbf{Warning}: This parameter provides an old and deprecated way of specifying boundary temperature models and shouldn't be used. Please use 'List of model names' instead. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model)=
## **Parameters in section** Boundary temperature model/Ascii data model
(parameters:Boundary_20temperature_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-temperature/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model:Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Time step between following data files. Depending on the setting of the global `Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model:Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. `Ma BP'). If this flag is set to `True' the plugin will first load the file with the number `First data file number' and decrease the file number during the model run. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model:First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The `First data file model time' parameter has been deactivated and will be removed in a future release. Do not use this paramter and instead provide data files starting from the model start time. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model:First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than `First velocity file model time'. 

(parameters:Boundary_20temperature_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Boundary_20temperature_20model:Box)=
## **Parameters in section** Boundary temperature model/Box
(parameters:Boundary_20temperature_20model:Box:Bottom_20temperature)=
### __Parameter name:__ Bottom temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the bottom boundary (at minimal $z$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box:Left_20temperature)=
### __Parameter name:__ Left temperature
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the left boundary (at minimal $x$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box:Right_20temperature)=
### __Parameter name:__ Right temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the right boundary (at maximal $x$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box:Top_20temperature)=
### __Parameter name:__ Top temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the top boundary (at maximal $x$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators)=
## **Parameters in section** Boundary temperature model/Box with lithosphere boundary indicators
(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators:Bottom_20temperature)=
### __Parameter name:__ Bottom temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the bottom boundary (at minimal $z$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators:Left_20temperature)=
### __Parameter name:__ Left temperature
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the left boundary (at minimal $x$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators:Left_20temperature_20lithosphere)=
### __Parameter name:__ Left temperature lithosphere
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the additional left lithosphere boundary (specified by user in Geometry Model). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators:Right_20temperature)=
### __Parameter name:__ Right temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the right boundary (at maximal $x$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators:Right_20temperature_20lithosphere)=
### __Parameter name:__ Right temperature lithosphere
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the additional right lithosphere boundary (specified by user in Geometry Model). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Box_20with_20lithosphere_20boundary_20indicators:Top_20temperature)=
### __Parameter name:__ Top temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the top boundary (at maximal $x$-value). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Constant)=
## **Parameters in section** Boundary temperature model/Constant
(parameters:Boundary_20temperature_20model:Constant:Boundary_20indicator_20to_20temperature_20mappings)=
### __Parameter name:__ Boundary indicator to temperature mappings
**Default value:**  

**Pattern:** [Map of <[Anything]>:<[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of mappings between boundary indicators and the temperature associated with the boundary indicators. The format for this list is ``indicator1 : value1, indicator2 : value2, ...'', where each indicator is a valid boundary indicator (either a number or the symbolic name of a boundary as provided by the geometry model) and each value is the temperature of that boundary. 

(parameters:Boundary_20temperature_20model:Dynamic_20core)=
## **Parameters in section** Boundary temperature model/Dynamic core
(parameters:Boundary_20temperature_20model:Dynamic_20core:Alpha)=
### __Parameter name:__ Alpha
**Default value:** 1.35e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Core thermal expansivity. Units: \si{\per\kelvin}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Beta_20composition)=
### __Parameter name:__ Beta composition
**Default value:** 1.1 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Compositional expansion coefficient $Beta_c$. See \cite{NPB+04} for more details. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:CMB_20pressure)=
### __Parameter name:__ CMB pressure
**Default value:** 0.14e12 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Pressure at CMB. Units: \si{\pascal}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Core_20conductivity)=
### __Parameter name:__ Core conductivity
**Default value:** 60. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Core heat conductivity $k_c$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Core_20density)=
### __Parameter name:__ Core density
**Default value:** 12.5e3 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Density of the core. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Core_20heat_20capacity)=
### __Parameter name:__ Core heat capacity
**Default value:** 840. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Heat capacity of the core. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Delta)=
### __Parameter name:__ Delta
**Default value:** 0.5 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** Partition coefficient of the light element. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Gravity_20acceleration)=
### __Parameter name:__ Gravity acceleration
**Default value:** 9.8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Gravitation acceleration at CMB. Units: \si{\meter\per\second\squared}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Initial_20light_20composition)=
### __Parameter name:__ Initial light composition
**Default value:** 0.01 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Initial light composition (eg. S,O) concentration in weight fraction. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Inner_20temperature)=
### __Parameter name:__ Inner temperature
**Default value:** 6000. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the inner boundary (core mantle boundary) at the beginning. Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:K0)=
### __Parameter name:__ K0
**Default value:** 4.111e11 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Core compressibility at zero pressure. See \cite{NPB+04} for more details. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Lh)=
### __Parameter name:__ Lh
**Default value:** 750e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The latent heat of core freeze. Units: \si{\joule\per\kilogram}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Max_20iteration)=
### __Parameter name:__ Max iteration
**Default value:** 30000 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The max iterations for nonlinear core energy solver. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Outer_20temperature)=
### __Parameter name:__ Outer temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the outer boundary (lithosphere water/air). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Rh)=
### __Parameter name:__ Rh
**Default value:** -27.7e6 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The heat of reaction. Units: \si{\joule\per\kilogram}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Rho0)=
### __Parameter name:__ Rho0
**Default value:** 7.019e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Core density at zero pressure. Units: \si{\kilogram\per\meter\cubed}. See \cite{NPB+04} for more details. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:dR_20over_20dt)=
### __Parameter name:__ dR over dt
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Initial inner core radius changing rate. Units: \si{\kilo\meter}/year. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:dT_20over_20dt)=
### __Parameter name:__ dT over dt
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Initial CMB temperature changing rate. Units: \si{\kelvin}/year. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:dX_20over_20dt)=
### __Parameter name:__ dX over dt
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Initial light composition changing rate. Units: 1/year. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters)=
## **Parameters in section** Boundary temperature model/Dynamic core/Geotherm parameters
(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters:Composition_20dependency)=
### __Parameter name:__ Composition dependency
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** If melting curve dependent on composition. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters:Theta)=
### __Parameter name:__ Theta
**Default value:** 0.11 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Melting curve (\cite{NPB+04} eq. (40)) parameter Theta. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters:Tm0)=
### __Parameter name:__ Tm0
**Default value:** 1695. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Melting curve (\cite{NPB+04} eq. (40)) parameter Tm0. Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters:Tm1)=
### __Parameter name:__ Tm1
**Default value:** 10.9 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Melting curve (\cite{NPB+04} eq. (40)) parameter Tm1. Units: \si{\per\tera\pascal}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters:Tm2)=
### __Parameter name:__ Tm2
**Default value:** -8.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Melting curve (\cite{NPB+04} eq. (40)) parameter Tm2. Units: \si{\per\tera\pascal\squared}. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Geotherm_20parameters:Use_20BW11)=
### __Parameter name:__ Use BW11
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If using the Fe-FeS system solidus from Buono \& Walker (2011) instead. 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Other_20energy_20source)=
## **Parameters in section** Boundary temperature model/Dynamic core/Other energy source
(parameters:Boundary_20temperature_20model:Dynamic_20core:Other_20energy_20source:File_20name)=
### __Parameter name:__ File name
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Data file name for other energy source into the core. The 'other energy source' is used for external core energy source.For example if someone want to test the early lunar core powered by precession (Dwyer, C. A., et al. (2011). A long-lived lunar dynamo driven by continuous mechanical stirring. Nature 479(7372): 212-214.)Format [Time(Gyr)   Energy rate(W)] 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Radioactive_20heat_20source)=
## **Parameters in section** Boundary temperature model/Dynamic core/Radioactive heat source
(parameters:Boundary_20temperature_20model:Dynamic_20core:Radioactive_20heat_20source:Half_20life_20times)=
### __Parameter name:__ Half life times
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Half decay times of different elements (Ga) 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Radioactive_20heat_20source:Heating_20rates)=
### __Parameter name:__ Heating rates
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Heating rates of different elements (W/kg) 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Radioactive_20heat_20source:Initial_20concentrations)=
### __Parameter name:__ Initial concentrations
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Initial concentrations of different elements (ppm) 

(parameters:Boundary_20temperature_20model:Dynamic_20core:Radioactive_20heat_20source:Number_20of_20radioactive_20heating_20elements)=
### __Parameter name:__ Number of radioactive heating elements
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Number of different radioactive heating elements in core 

(parameters:Boundary_20temperature_20model:Function)=
## **Parameters in section** Boundary temperature model/Function
(parameters:Boundary_20temperature_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Boundary_20temperature_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Boundary_20temperature_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Boundary_20temperature_20model:Function:Maximal_20temperature)=
### __Parameter name:__ Maximal temperature
**Default value:** 3773. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximal temperature. Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Function:Minimal_20temperature)=
### __Parameter name:__ Minimal temperature
**Default value:** 273. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimal temperature. Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Boundary_20temperature_20model:Initial_20temperature)=
## **Parameters in section** Boundary temperature model/Initial temperature
(parameters:Boundary_20temperature_20model:Initial_20temperature:Maximal_20temperature)=
### __Parameter name:__ Maximal temperature
**Default value:** 3773. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximal temperature. Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Initial_20temperature:Minimal_20temperature)=
### __Parameter name:__ Minimal temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimal temperature. Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Spherical_20constant)=
## **Parameters in section** Boundary temperature model/Spherical constant
(parameters:Boundary_20temperature_20model:Spherical_20constant:Inner_20temperature)=
### __Parameter name:__ Inner temperature
**Default value:** 6000. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the inner boundary (core mantle boundary). Units: \si{\kelvin}. 

(parameters:Boundary_20temperature_20model:Spherical_20constant:Outer_20temperature)=
### __Parameter name:__ Outer temperature
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Temperature at the outer boundary (lithosphere water/air). Units: \si{\kelvin}. 

(parameters:Boundary_20traction_20model)=
## **Parameters in section** Boundary traction model
(parameters:Boundary_20traction_20model:Prescribed_20traction_20boundary_20indicators)=
### __Parameter name:__ Prescribed traction boundary indicators
**Default value:**  

**Pattern:** [Map of <[Anything]>:<[Selection ascii data|function|initial lithostatic pressure|zero traction ]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list denoting those boundaries on which a traction force is prescribed, i.e., where known external forces act, resulting in an unknown velocity. This is often used to model ``open'' boundaries where we only know the pressure. This pressure then produces a force that is normal to the boundary and proportional to the pressure.

The format of valid entries for this parameter is that of a map given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where each key must be a valid boundary indicator (which is either an integer or the symbolic name the geometry model in use may have provided for this part of the boundary) and each value must be one of the currently implemented boundary traction models. ``selector'' is an optional string given as a subset of the letters `xyz' that allows you to apply the boundary conditions only to the components listed. As an example, '1 y: function' applies the type `function' to the y component on boundary 1. Without a selector it will affect all components of the traction. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model)=
## **Parameters in section** Boundary traction model/Ascii data model
(parameters:Boundary_20traction_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-traction/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model:Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Time step between following data files. Depending on the setting of the global `Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model:Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. `Ma BP'). If this flag is set to `True' the plugin will first load the file with the number `First data file number' and decrease the file number during the model run. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model:First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The `First data file model time' parameter has been deactivated and will be removed in a future release. Do not use this paramter and instead provide data files starting from the model start time. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model:First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than `First velocity file model time'. 

(parameters:Boundary_20traction_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Boundary_20traction_20model:Function)=
## **Parameters in section** Boundary traction model/Function
(parameters:Boundary_20traction_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Boundary_20traction_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Boundary_20traction_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Boundary_20traction_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Boundary_20traction_20model:Initial_20lithostatic_20pressure)=
## **Parameters in section** Boundary traction model/Initial lithostatic pressure
(parameters:Boundary_20traction_20model:Initial_20lithostatic_20pressure:Number_20of_20integration_20points)=
### __Parameter name:__ Number of integration points
**Default value:** 1000 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of integration points over which we integrate the lithostatic pressure downwards. 

(parameters:Boundary_20traction_20model:Initial_20lithostatic_20pressure:Representative_20point)=
### __Parameter name:__ Representative point
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The point where the pressure profile will be calculated. Cartesian coordinates $(x,y,z)$ when geometry is a box, otherwise enter radius, longitude, and in 3D latitude. Note that the coordinate related to the depth ($y$ in 2D cartesian, $z$ in 3D cartesian and radius in spherical coordinates) is not used. Units: \si{\meter} or degrees. 

(parameters:Boundary_20velocity_20model)=
## **Parameters in section** Boundary velocity model
(parameters:Boundary_20velocity_20model:Prescribed_20velocity_20boundary_20indicators)=
### __Parameter name:__ Prescribed velocity boundary indicators
**Default value:**  

**Pattern:** [Map of <[Anything]>:<[Selection ascii data|function|gplates|zero velocity ]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list denoting those boundaries on which the velocity is prescribed, i.e., where unknown external forces act to prescribe a particular velocity. This is often used to prescribe a velocity that equals that of overlying plates.

The format of valid entries for this parameter is that of a map given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where each key must be a valid boundary indicator (which is either an integer or the symbolic name the geometry model in use may have provided for this part of the boundary) and each value must be one of the currently implemented boundary velocity models. ``selector'' is an optional string given as a subset of the letters `xyz' that allows you to apply the boundary conditions only to the components listed. As an example, '1 y: function' applies the type `function' to the y component on boundary 1. Without a selector it will affect all components of the velocity.

Note that the no-slip boundary condition is a special case of the current one where the prescribed velocity happens to be zero. It can thus be implemented by indicating that a particular boundary is part of the ones selected using the current parameter and using ``zero velocity'' as the boundary values. Alternatively, you can simply list the part of the boundary on which the velocity is to be zero with the parameter ``Zero velocity boundary indicator'' in the current parameter section.

Note that when ``Use years in output instead of seconds'' is set to true, velocity should be given in m/yr. The following boundary velocity models are available:

`ascii data': Implementation of a model in which the boundary velocity is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `velocity${}_x$', `velocity${}_y$' in a 2d model or `x', `y', `velocity${}_x$', `velocity${}_y$', `velocity${}_z$' in a 3d model. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates.If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions. Velocities can be specified using cartesian (by default) or spherical unit vectors. No matter which geometry model is chosen, the unit of the velocities is assumed to be m/s or m/yr depending on the `Use years in output instead of seconds' flag. If you provide velocities in cm/yr, set the `Scale factor' option to 0.01.

`function': Implementation of a model in which the boundary velocity is given in terms of an explicit formula that is elaborated in the parameters in section ``Boundary velocity model|Function''. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

The formula you describe in the mentioned section is a semicolon separated list of velocities for each of the $d$ components of the velocity vector. These $d$ formulas are interpreted as having units m/s, unless the global input parameter ``Use years in output instead of seconds'' is set, in which case we interpret the formula expressions as having units m/year.

Likewise, since the symbol $t$ indicating time may appear in the formulas for the prescribed velocities, it is interpreted as having units seconds unless the global parameter above has been set.

`gplates': Implementation of a model in which the boundary velocity is derived from files that are generated by the GPlates program.

`zero velocity': Implementation of a model in which the boundary velocity is zero. This is commonly referred to as a ``stick boundary condition'', indicating that the material ``sticks'' to the material on the other side of the boundary. 

(parameters:Boundary_20velocity_20model:Tangential_20velocity_20boundary_20indicators)=
### __Parameter name:__ Tangential velocity boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries on which the velocity is tangential and unrestrained, i.e., free-slip where no external forces act to prescribe a particular tangential velocity (although there is a force that requires the flow to be tangential).

The names of the boundaries listed here can either by numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model. 

(parameters:Boundary_20velocity_20model:Zero_20velocity_20boundary_20indicators)=
### __Parameter name:__ Zero velocity boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries on which the velocity is zero.

The names of the boundaries listed here can either by numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model)=
## **Parameters in section** Boundary velocity model/Ascii data model
(parameters:Boundary_20velocity_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-velocity/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Time step between following data files. Depending on the setting of the global `Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. `Ma BP'). If this flag is set to `True' the plugin will first load the file with the number `First data file number' and decrease the file number during the model run. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The `First data file model time' parameter has been deactivated and will be removed in a future release. Do not use this paramter and instead provide data files starting from the model start time. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than `First velocity file model time'. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Boundary_20velocity_20model:Ascii_20data_20model:Use_20spherical_20unit_20vectors)=
### __Parameter name:__ Use spherical unit vectors
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Specify velocity as r, phi, and theta components instead of x, y, and z. Positive velocities point up, east, and north (in 3D) or out and clockwise (in 2D). This setting only makes sense for spherical geometries. 

(parameters:Boundary_20velocity_20model:Function)=
## **Parameters in section** Boundary velocity model/Function
(parameters:Boundary_20velocity_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Boundary_20velocity_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Boundary_20velocity_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Boundary_20velocity_20model:Function:Use_20spherical_20unit_20vectors)=
### __Parameter name:__ Use spherical unit vectors
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Specify velocity as $r$, $\phi$, and $\theta$ components instead of $x$, $y$, and $z$. Positive velocities point up, east, and north (in 3D) or out and clockwise (in 2D). This setting only makes sense for spherical geometries. 

(parameters:Boundary_20velocity_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Boundary_20velocity_20model:GPlates_20model)=
## **Parameters in section** Boundary velocity model/GPlates model
(parameters:Boundary_20velocity_20model:GPlates_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a '/') or relative to the current directory. The path may also include the special text '$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.  

(parameters:Boundary_20velocity_20model:GPlates_20model:Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Time step between following velocity files. Depending on the setting of the global 'Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years. 

(parameters:Boundary_20velocity_20model:GPlates_20model:Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. 'Ma BP'). If this flag is set to 'True' the plugin will first load the file with the number 'First velocity file number' and decrease the file number during the model run. 

(parameters:Boundary_20velocity_20model:GPlates_20model:First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Time from which on the velocity file with number 'First velocity file number' is used as boundary condition. Previous to this time, a no-slip boundary condition is assumed. Depending on the setting of the global 'Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. 

(parameters:Boundary_20velocity_20model:GPlates_20model:First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than 'First velocity file model time'. 

(parameters:Boundary_20velocity_20model:GPlates_20model:Lithosphere_20thickness)=
### __Parameter name:__ Lithosphere thickness
**Default value:** 100000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Determines the depth of the lithosphere, so that the GPlates velocities can be applied at the sides of the model as well as at the surface. 

(parameters:Boundary_20velocity_20model:GPlates_20model:Point_20one)=
### __Parameter name:__ Point one
**Default value:** 1.570796,0.0 

**Pattern:** [Anything] 

**Documentation:** Point that determines the plane in which a 2D model lies in. Has to be in the format `a,b' where a and b are theta (polar angle) and phi in radians. This value is not utilized in 3D geometries, and can therefore be set to the default or any user-defined quantity. 

(parameters:Boundary_20velocity_20model:GPlates_20model:Point_20two)=
### __Parameter name:__ Point two
**Default value:** 1.570796,1.570796 

**Pattern:** [Anything] 

**Documentation:** Point that determines the plane in which a 2D model lies in. Has to be in the format `a,b' where a and b are theta (polar angle) and phi in radians. This value is not utilized in 3D geometries, and can therefore be set to the default or any user-defined quantity. 

(parameters:Boundary_20velocity_20model:GPlates_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the boundary velocity. You might want to use this to scale the velocities to a reference model (e.g. with free-slip boundary) or another plate reconstruction. 

(parameters:Boundary_20velocity_20model:GPlates_20model:Velocity_20file_20name)=
### __Parameter name:__ Velocity file name
**Default value:** phi.%d 

**Pattern:** [Anything] 

**Documentation:** The file name of the material data. Provide file in format: (Velocity file name).\%d.gpml where \%d is any sprintf integer qualifier, specifying the format of the current file number. 

(parameters:Checkpointing)=
## **Parameters in section** Checkpointing
(parameters:Checkpointing:Steps_20between_20checkpoint)=
### __Parameter name:__ Steps between checkpoint
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of timesteps between performing checkpoints. If 0 and time between checkpoint is not specified, checkpointing will not be performed. Units: None. 

(parameters:Checkpointing:Time_20between_20checkpoint)=
### __Parameter name:__ Time between checkpoint
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The wall time between performing checkpoints. If 0, will use the checkpoint step frequency instead. Units: Seconds. 

(parameters:Compositional_20fields)=
## **Parameters in section** Compositional fields
(parameters:Compositional_20fields:Compositional_20field_20methods)=
### __Parameter name:__ Compositional field methods
**Default value:**  

**Pattern:** [List of <[Selection field|particles|volume of fluid|static|melt field|darcy field|prescribed field|prescribed field with diffusion ]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list denoting the solution method of each compositional field. Each entry of the list must be one of the currently implemented field methods.

These choices correspond to the following methods by which compositional fields gain their values:\begin{itemize}\item ``field'': If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the velocity field, and applying reaction rates to it. In other words, this corresponds to the usual notion of a composition field as mentioned in Section~\ref{sec:compositional}. 
\item ``particles'': If a compositional field is marked with this method, then its values are obtained in each time step by interpolating the corresponding properties from the particles located on each cell. The time evolution therefore happens because particles move along with the velocity field, and particle properties can react with each other as well. See Section~\ref{sec:particles} for more information about how particles behave.
\item ``volume of fluid``: If a compositional field is marked with this method, then its values are obtained in each timestep by reconstructing a polynomial finite element approximation on each cell from a volume of fluid interface tracking method, which is used to compute the advection updates.
\item ``static'': If a compositional field is marked this way, then it does not evolve at all. Its values are simply set to the initial conditions, and will then never change.
\item ``melt field'': If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the melt velocity, and applying reaction rates to it. In other words, this corresponds to the usual notion of a composition field as mentioned in Section~\ref{sec:compositional}, except that it is advected with the melt velocity instead of the solid velocity. This method can only be chosen if melt transport is active in the model.
\item ``darcy field'': If a compositional field is marked with this method, then its values are computed in each time step by advecting along the values of the previous time step using the fluid velocity prescribed by Darcy's Law, and applying reaction rates to it. In other words this corresponds to the usual notion of a composition field as mentioned in Section~\ref{sec:compositional}, except that it is advected with the Darcy velocity instead of the solid velocity. This method requires there to be a compositional field named porosity that is advected the darcy field method. We calculate the fluid velocity $u_f$ using an approximation of Darcy's Law: $u_f = u_s - K_D / \phi * (rho_s * g - rho_f * g)$.
\item ``prescribed field'': The value of these fields is determined in each time step from the material model. If a compositional field is marked with this method, then the value of a specific additional material model output, called the `PrescribedFieldOutputs' is interpolated onto the field. This field does not change otherwise, it is not advected with the flow.
\item ``prescribed field with diffusion'': If a compositional field is marked this way, the value of a specific additional material model output, called the `PrescribedFieldOutputs' is interpolated onto the field, as in the ``prescribed field'' method. Afterwards, the field is diffused based on a solver parameter, the diffusion length scale, smoothing the field. Specifically, the field is updated by solving the equation $(I-l^2 \Delta) C_\text{smoothed} = C_\text{prescribed}$, where $l$ is the diffusion length scale. Note that this means that the amount of diffusion is independent of the time step size, and that the field is not advected with the flow.\end{itemize} 

(parameters:Compositional_20fields:List_20of_20normalized_20fields)=
### __Parameter name:__ List of normalized fields
**Default value:**  

**Pattern:** [List of <[Integer range 0...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of integers smaller than or equal to the number of compositional fields. All compositional fields in this list will be normalized before the first timestep. The normalization is implemented in the following way: First, the sum of the fields to be normalized is calculated at every point and the global maximum is determined. Second, the compositional fields to be normalized are divided by this maximum. 

(parameters:Compositional_20fields:Mapped_20particle_20properties)=
### __Parameter name:__ Mapped particle properties
**Default value:**  

**Pattern:** [Map of <[Anything]>:<[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list denoting the particle properties that will be projected to those compositional fields that are of the ``particles'' field type.

The format of valid entries for this parameter is that of a map given as ``key1: value1, key2: value2 [component2], key3: value3 [component4], ...'' where each key must be a valid field name of the ``particles'' type, and each value must be one of the currently selected particle properties. Component is a component index of the particle property that is 0 by default, but can be set up to n-1, where n is the number of vector components of this particle property. The component indicator only needs to be set if not the first component of the particle property should be mapped (e.g. the $y$-component of the velocity at the particle positions). 

(parameters:Compositional_20fields:Names_20of_20fields)=
### __Parameter name:__ Names of fields
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A user-defined name for each of the compositional fields requested. 

(parameters:Compositional_20fields:Number_20of_20fields)=
### __Parameter name:__ Number of fields
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of fields that will be advected along with the flow field, excluding velocity, pressure and temperature. 

(parameters:Compositional_20fields:Types_20of_20fields)=
### __Parameter name:__ Types of fields
**Default value:** unspecified 

**Pattern:** [List of <[Selection chemical composition|stress|grain size|porosity|density|generic|unspecified ]> of length 0...4294967295 (inclusive)] 

**Documentation:** A type for each of the compositional fields requested. Each entry of the list must be one of several recognized types: chemical composition, stress, grain size, porosity, general and unspecified. The generic type is intended to be a placeholder type that has no effect on the running of any material model, while the unspecified type is intended to tell ASPECT that the user has not explicitly indicated the type of field (facilitating parameter file checking). If a plugin such as a material model uses these types, the choice of type will affect how that module functions. 

(parameters:Discretization)=
## **Parameters in section** Discretization
(parameters:Discretization:Composition_20polynomial_20degree)=
### __Parameter name:__ Composition polynomial degree
**Default value:** 2 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The polynomial degree to use for the composition variable(s). As an example, a value of 2 for this parameter will yield either the element $Q_2$ or $DGQ_2$ for the compositional field(s), depending on whether we use continuous or discontinuous field(s). 

For continuous elements, the value needs to be 1 or larger as $Q_1$ is the lowest order element, while $DGQ_0$ is a valid choice. Units: None. 

(parameters:Discretization:Stokes_20velocity_20polynomial_20degree)=
### __Parameter name:__ Stokes velocity polynomial degree
**Default value:** 2 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The polynomial degree to use for the velocity variables in the Stokes system. The polynomial degree for the pressure variable will then be one less in order to make the velocity/pressure pair conform with the usual LBB (Babu{\v s}ka-Brezzi) condition. In other words, we are using a Taylor-Hood element for the Stokes equations and this parameter indicates the polynomial degree of it. As an example, a value of 2 for this parameter will yield the element $Q_2^d \times Q_1$ for the $d$ velocity components and the pressure, respectively (unless the `Use locally conservative discretization' parameter is set, which modifies the pressure element). 

Be careful if you choose 1 as the degree. The resulting element is not stable and it may lead to artifacts in the solution. Units: None. 

(parameters:Discretization:Temperature_20polynomial_20degree)=
### __Parameter name:__ Temperature polynomial degree
**Default value:** 2 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The polynomial degree to use for the temperature variable. As an example, a value of 2 for this parameter will yield either the element $Q_2$ or $DGQ_2$ for the temperature field, depending on whether we use a continuous or discontinuous field. Units: None. 

(parameters:Discretization:Use_20discontinuous_20composition_20discretization)=
### __Parameter name:__ Use discontinuous composition discretization
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a composition discretization that is discontinuous as opposed to continuous. This then requires the assembly of face terms between cells, and weak imposition of boundary terms for the composition field via the discontinuous Galerkin method. 

(parameters:Discretization:Use_20discontinuous_20temperature_20discretization)=
### __Parameter name:__ Use discontinuous temperature discretization
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a temperature discretization that is discontinuous as opposed to continuous. This then requires the assembly of face terms between cells, and weak imposition of boundary terms for the temperature field via the interior-penalty discontinuous Galerkin method. 

(parameters:Discretization:Use_20equal_20order_20interpolation_20for_20Stokes)=
### __Parameter name:__ Use equal order interpolation for Stokes
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** By default (i.e., when this parameter is set to its default value `false') \aspect{} uses finite element combinations in which the pressure shape functions are polynomials one degree lower than the shape functions for the velocity. An example is the Taylor-Hood element that uses $Q_k$ elements for the velocity and $Q_{k-1}$ for the pressure. This is because using the \textit{same} polynomial degree for both the velocity and the pressure turns out to violate some mathematical properties necessary to make the problem solvable. (In particular, thecondition in question goes by the name ``inf-sup'' or Babu{\v s}ka-Brezzi or LBB condition.) A consequence of violating this condition is that the pressure may show oscillations and not converge to the correct pressure.

That said, people have often used $Q_1$ elements for both the velocity and pressure anyway. This is commonly referred to as using the $Q_1-Q_1$ method. It is, by default, not stable as mentioned above, but it can be made stable by adding a small amount of compressibility to the model. There are numerous ways to do that. Today, the way that is generally considered to be the best approach is the one by Dohrmann and Bochev \cite{DohrmannBochev2004}.

When this parameter is set to ``true'', then \aspect{} will use this method by using $Q_k\times Q_k$ elements for velocity and pressure, respectively, where $k$ is the value provided for the parameter ``Stokes velocity polynomial degree''.

\note{While \aspect{} \textit{allows} you to use this   method, it is generally understood that this is not a   great idea as it leads to rather low accuracy in   general as documented in \cite{thba22}.   It also leads to substantial problems when   using free surfaces. As a consequence, the presence   of this parameter should not be seen as an   endorsement of the method, or a suggestion to   actually use it. It simply makes the method available.} 

(parameters:Discretization:Use_20locally_20conservative_20discretization)=
### __Parameter name:__ Use locally conservative discretization
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a Stokes discretization that is locally conservative at the expense of a larger number of degrees of freedom (true), or to go with a cheaper discretization that does not locally conserve mass, although it is globally conservative (false).

When using a locally conservative discretization, the finite element space for the pressure is discontinuous between cells and is the polynomial space $P_{-(k-1)}$ of polynomials of degree $k-1$ in each variable separately. Here, $k$ is the value given in the parameter ``Stokes velocity polynomial degree'', and consequently the polynomial degree for the pressure, $k-1$, is one lower than that for the velocity.

As a consequence of choosing this element for the pressure rather than the more commonly used $Q_{k-1}$ element that is continuous, it can be shown that if the medium is considered incompressible then the computed discrete velocity field $\mathbf u_h$ satisfies the property $\int_ {\partial K} \mathbf u_h \cdot \mathbf n = 0$ for every cell $K$, i.e., for each cell inflow and outflow exactly balance each other as one would expect for an incompressible medium. In other words, the velocity field is \textit{locally conservative}.

On the other hand, if this parameter is set to ``false''(the default), then the finite element space is chosen as $Q_{k-1}$. This choice does not yield the local conservation property but has the advantage of requiring fewer degrees of freedom. Furthermore, the error is generally smaller with this choice.

For an in-depth discussion of these issues and a quantitative evaluation of the different choices, see \cite{KHB12}. 

(parameters:Discretization:Stabilization_20parameters)=
## **Parameters in section** Discretization/Stabilization parameters
(parameters:Discretization:Stabilization_20parameters:Discontinuous_20penalty)=
### __Parameter name:__ Discontinuous penalty
**Default value:** 10. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value used to penalize discontinuities in the discontinuous Galerkin method. This is used only for the temperature field, and not for the composition field, as pure advection does not use the interior penalty method. This is largely empirically decided -- it must be large enough to ensure the bilinear form is coercive, but not so large as to penalize discontinuity at all costs. 

(parameters:Discretization:Stabilization_20parameters:Global_20composition_20maximum)=
### __Parameter name:__ Global composition maximum
**Default value:** 1.7976931348623157e+308 

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The maximum global composition values that will be used in the bound preserving limiter for the discontinuous solutions from composition advection fields. The number of the input 'Global composition maximum' values separated by ',' has to be one or the same as the number of the compositional fields. When only one value is supplied, this same value is assumed for all compositional fields. 

(parameters:Discretization:Stabilization_20parameters:Global_20composition_20minimum)=
### __Parameter name:__ Global composition minimum
**Default value:** -1.7976931348623157e+308 

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The minimum global composition value that will be used in the bound preserving limiter for the discontinuous solutions from composition advection fields. The number of the input 'Global composition minimum' values separated by ',' has to be one or the same as the number of the compositional fields. When only one value is supplied, this same value is assumed for all compositional fields. 

(parameters:Discretization:Stabilization_20parameters:Global_20temperature_20maximum)=
### __Parameter name:__ Global temperature maximum
**Default value:** 1.7976931348623157e+308 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum global temperature value that will be used in the bound preserving limiter for the discontinuous solutions from temperature advection fields. 

(parameters:Discretization:Stabilization_20parameters:Global_20temperature_20minimum)=
### __Parameter name:__ Global temperature minimum
**Default value:** -1.7976931348623157e+308 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum global temperature value that will be used in the bound preserving limiter for the discontinuous solutions from temperature advection fields. 

(parameters:Discretization:Stabilization_20parameters:List_20of_20compositional_20fields_20with_20disabled_20boundary_20entropy_20viscosity)=
### __Parameter name:__ List of compositional fields with disabled boundary entropy viscosity
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** Select for which compositional fields to skip the entropy viscosity stabilization at dirichlet boundaries. This is only advisable for compositional fieldsthat have intrinsic physical diffusion terms, otherwise oscillations may develop. The parameter should contain a list of compositional field names. 

(parameters:Discretization:Stabilization_20parameters:Stabilization_20method)=
### __Parameter name:__ Stabilization method
**Default value:** entropy viscosity 

**Pattern:** [Selection entropy viscosity|SUPG ] 

**Documentation:** Select the method for stabilizing the advection equation. The original method implemented is 'entropy viscosity' as described in \cite {KHB12}. SUPG is currently experimental. 

(parameters:Discretization:Stabilization_20parameters:Use_20artificial_20viscosity_20smoothing)=
### __Parameter name:__ Use artificial viscosity smoothing
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If set to false, the artificial viscosity of a cell is computed and is computed on every cell separately as discussed in \cite{KHB12}. If set to true, the maximum of the artificial viscosity in the cell as well as the neighbors of the cell is computed and used instead. 

(parameters:Discretization:Stabilization_20parameters:Use_20limiter_20for_20discontinuous_20composition_20solution)=
### __Parameter name:__ Use limiter for discontinuous composition solution
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to apply the bound preserving limiter as a correction after having the discontinuous composition solution. Currently we apply this only to the compositional solution if the 'Global composition maximum' and 'Global composition minimum' are already defined in the .prm file. This limiter keeps the discontinuous solution in the range given by Global composition maximum' and 'Global composition minimum'. 

(parameters:Discretization:Stabilization_20parameters:Use_20limiter_20for_20discontinuous_20temperature_20solution)=
### __Parameter name:__ Use limiter for discontinuous temperature solution
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to apply the bound preserving limiter as a correction after computing the discontinuous temperature solution. Currently we apply this only to the temperature solution if the 'Global temperature maximum' and 'Global temperature minimum' are already defined in the .prm file. This limiter keeps the discontinuous solution in the range given by 'Global temperature maximum' and 'Global temperature minimum'. 

(parameters:Discretization:Stabilization_20parameters:alpha)=
### __Parameter name:__ alpha
**Default value:** 2 

**Pattern:** [Integer range 1...2 (inclusive)] 

**Documentation:** The exponent $\alpha$ in the entropy viscosity stabilization. Valid options are 1 or 2. The recommended setting is 2. (This parameter does not correspond to any variable in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see \cite{KHB12}. Rather, the paper always uses 2 as the exponent in the definition of the entropy, following equation (15) of the paper. The full approach is discussed in \cite{GPP11}.) Note that this is not the thermal expansion coefficient, also commonly referred to as $\alpha$.Units: None. 

(parameters:Discretization:Stabilization_20parameters:beta)=
### __Parameter name:__ beta
**Default value:** 0.052 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The $\beta$ factor in the artificial viscosity stabilization. This parameter controls the maximum dissipation of the entropy viscosity, which is the part that only scales with the cell diameter and the maximum velocity in the cell, but does not depend on the solution field itself or its residual. An appropriate value for 2d is 0.052 and 0.78 for 3d. (For historical reasons, the name used here is different from the one used in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see \cite{KHB12}. This parameter can be given as a single value or as a list with as many entries as one plus the number of compositional fields. In the former case all advection fields use the same stabilization parameters, in the latter case each field (temperature first, then all compositions) use individual parameters. This can be useful to reduce the stabilization for the temperature, which already has some physical diffusion. This parameter corresponds to the factor $\alpha_{\text{max}}$ in the formulas following equation (15) of the paper.) Units: None. 

(parameters:Discretization:Stabilization_20parameters:cR)=
### __Parameter name:__ cR
**Default value:** 0.11 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The $c_R$ factor in the entropy viscosity stabilization. This parameter controls the part of the entropy viscosity that depends on the solution field itself and its residual in addition to the cell diameter and the maximum velocity in the cell. This parameter can be given as a single value or as a list with as many entries as one plus the number of compositional fields. In the former case all advection fields use the same stabilization parameters, in the latter case each field (temperature first, then all compositions) use individual parameters. This can be useful to reduce the stabilization for the temperature, which already has some physical diffusion. (For historical reasons, the name used here is different from the one used in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see \cite{KHB12}. This parameter corresponds to the factor $\alpha_E$ in the formulas following equation (15) of the paper.) Units: None. 

(parameters:Discretization:Stabilization_20parameters:gamma)=
### __Parameter name:__ gamma
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The strain rate scaling factor in the artificial viscosity stabilization. This parameter determines how much the strain rate (in addition to the velocity) should influence the stabilization. (This parameter does not correspond to any variable in the 2012 paper by Kronbichler, Heister and Bangerth that describes ASPECT, see \cite{KHB12}. Rather, the paper always uses 0, i.e. they specify the maximum dissipation $\nu_h^\text{max}$ as $\nu_h^\text{max}\vert_K = \alpha_{\text{max}} h_K \|\mathbf u\|_{\infty,K}$. Here, we use $\|\lvert\mathbf u\rvert + \gamma h_K \lvert\varepsilon (\mathbf u)\rvert\|_{\infty,K}$ instead of $\|\mathbf u\|_{\infty,K}$. Units: None. 

(parameters:Formulation)=
## **Parameters in section** Formulation
(parameters:Formulation:Enable_20additional_20Stokes_20RHS)=
### __Parameter name:__ Enable additional Stokes RHS
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to ask the material model for additional terms for the right-hand side of the Stokes equation. This feature is likely only used when implementing force vectors for manufactured solution problems and requires filling additional outputs of type AdditionalMaterialOutputsStokesRHS. 

(parameters:Formulation:Enable_20elasticity)=
### __Parameter name:__ Enable elasticity
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to include the additional elastic terms on the right-hand side of the Stokes equation. 

(parameters:Formulation:Enable_20prescribed_20dilation)=
### __Parameter name:__ Enable prescribed dilation
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to include additional terms on the right-hand side of the Stokes equation to set a given compression term specified in the MaterialModel output PrescribedPlasticDilation. 

(parameters:Formulation:Formulation)=
### __Parameter name:__ Formulation
**Default value:** custom 

**Pattern:** [Selection isentropic compression|custom|anelastic liquid approximation|Boussinesq approximation ] 

**Documentation:** Select a formulation for the basic equations. Different published formulations are available in ASPECT (see the list of possible values for this parameter in the manual for available options). Two ASPECT specific options are
\begin{enumerate}
  \item `isentropic compression': ASPECT's original formulation, using the explicit compressible mass equation, and the full density for the temperature equation.
  \item `custom': A custom selection of `Mass conservation' and `Temperature equation'.
\end{enumerate}

\note{Warning: The `custom' option is implemented for advanced users that want full control over the equations solved. It is possible to choose inconsistent formulations and no error checking is performed on the consistency of the resulting equations.}

\note{The `anelastic liquid approximation' option here can also be used to set up the `truncated anelastic liquid approximation' as long as this option is chosen together with a material model that defines a density that depends on temperature and depth and not on the pressure.} 

(parameters:Formulation:Mass_20conservation)=
### __Parameter name:__ Mass conservation
**Default value:** ask material model 

**Pattern:** [Selection incompressible|isentropic compression|hydrostatic compression|reference density profile|implicit reference density profile|projected density field|ask material model ] 

**Documentation:** Possible approximations for the density derivatives in the mass conservation equation. Note that this parameter is only evaluated if `Formulation' is set to `custom'. Other formulations ignore the value of this parameter. 

(parameters:Formulation:Temperature_20equation)=
### __Parameter name:__ Temperature equation
**Default value:** real density 

**Pattern:** [Selection real density|reference density profile ] 

**Documentation:** Possible approximations for the density in the temperature equation. Possible approximations are `real density' and `reference density profile'. Note that this parameter is only evaluated if `Formulation' is set to `custom'. Other formulations ignore the value of this parameter. 

(parameters:Geometry_20model)=
## **Parameters in section** Geometry model
(parameters:Geometry_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection box|box with lithosphere boundary indicators|chunk|chunk with lithosphere boundary indicators|ellipsoidal chunk|sphere|spherical shell|unspecified ] 

**Documentation:** Select one of the following models:

`box': A box geometry parallel to the coordinate directions. The extent of the box in each coordinate direction is set in the parameter file. The box geometry labels its 2*dim sides as follows: in 2d, boundary indicators 0 through 3 denote the left, right, bottom and top boundaries; in 3d, boundary indicators 0 through 5 indicate left, right, front, back, bottom and top boundaries (see also the documentation of the deal.II class ``ReferenceCell''). You can also use symbolic names ``left'', ``right'', etc., to refer to these boundaries in input files. It is also possible to add initial topography to the box model. Note however that this is done after the last initial adaptive refinement cycle. Also, initial topography is supposed to be small, as it is not taken into account when depth or a representative point is computed. 

`box with lithosphere boundary indicators': A box geometry parallel to the coordinate directions. The extent of the box in each coordinate direction is set in the parameter file. This geometry model labels its sides with 2*dim+2*(dim-1) boundary indicators: in 2d, boundary indicators 0 through 3 denote the left, right, bottom and top boundaries, while indicators4 and 5 denote the upper part of the left and right vertical boundary, respectively. In 3d, boundary indicators 0 through 5 indicate left, right, front, back, bottom and top boundaries (see also the documentation of the deal.II class ``ReferenceCell''), while indicators 6, 7, 8 and 9 denote the left, right, front and back upper parts of the vertical boundaries, respectively. You can also use symbolic names ``left'', ``right'', ``left lithosphere'', etc., to refer to these boundaries in input files.

Note that for a given ``Global refinement level'' and no user-specified ``Repetitions'', the lithosphere part of the mesh will be more refined. 

The additional boundary indicators for the lithosphere allow for selecting boundary conditions for the lithosphere different from those for the underlying mantle. An example application of this geometry is to prescribe a velocity on the lithospheric plates, but use open boundary conditions underneath. 

`chunk': A geometry which can be described as a chunk of a spherical shell, bounded by lines of longitude, latitude and radius. The minimum and maximum longitude, latitude (if in 3d) and depth of the chunk is set in the parameter file. The chunk geometry labels its 2*dim sides as follows: ``west'' and ``east'': minimum and maximum longitude, ``south'' and ``north'': minimum and maximum latitude, ``inner'' and ``outer'': minimum and maximum radii. 

The dimensions of the model are specified by parameters of the following form: Chunk (minimum || maximum) (longitude || latitude): edges of geographical quadrangle (in degrees)Chunk (inner || outer) radius: Radii at bottom and top of chunk(Longitude || Latitude || Radius) repetitions: number of cells in each coordinate direction.

When used in 2d, this geometry does not imply the use of a spherical coordinate system. Indeed, in 2d the geometry is simply a sector of an annulus in a Cartesian coordinate system and consequently would correspond to a sector of a cross section of the fluid filled space between two infinite cylinders where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in Section~\ref{sec:meaning-of-2d}. It is also possible to add initial topography to the chunk geometry, based on an ascii data file. 

`chunk with lithosphere boundary indicators': A geometry which can be described as a chunk of a spherical shell, bounded by lines of longitude, latitude and radius. The side boundaries have two boundary indicators, so the user can prescribe different boundary conditions on these boundaries. The minimum and maximum longitude, (latitude) and depth of the chunk are set in the parameter file. The chunk geometry labels its 2*dim+2*(dim-1) sides as follows: ``lower west'' and ``lower east'': minimum and maximum longitude of the lower part of the east and west side boundaries, ``upper west and upper east'': minimum and maximum longitude of the upper part of the east and west side boundaries, ``lower south'' and ``lower north'': minimum and maximum latitude of the lower part of the south and north side boundaries, ``upper south'' and ``upper north'': minimum and maximum latitude of the upper part of the south and north side boundaries, 

The dimensions of the model are specified by parameters of the following form: Chunk (minimum || maximum) (longitude || latitude): edges of geographical quadrangle (in degrees). Chunk (inner || outer || middle boundary) radius: Radii at bottom and top of chunk and the radius at which the lower boundary indicator along a side boundary transitions into the upper boundary indicator. (Longitude || Latitude) repetitions: number of cells in each coordinate direction.(Inner || Outer) chunk radius repetitions: number of cells in the radial coordinate direction for the lower part of the domain (up to the Middle boundary radius) and for the upper part of the domain. 

When used in 2d, this geometry does not imply the use of a spherical coordinate system. Indeed, in 2d the geometry is simply a sector of an annulus in a Cartesian coordinate system and consequently would correspond to a sector of a cross section of the fluid filled space between two infinite cylinders where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in Section~\ref{sec:meaning-of-2d}. It is also possible to add initial topography to the chunk geometry, based on an ascii data file. 

`ellipsoidal chunk': A 3D chunk geometry that accounts for Earth's ellipticity (default assuming the WGS84 ellipsoid definition) which can be defined in non-coordinate directions. In the description of the ellipsoidal chunk, two of the ellipsoidal axes have the same length so that there is only a semi-major axis and a semi-minor axis. The user has two options for creating an ellipsoidal chunk geometry: 1) by defining two opposing points (SW and NE or NW and SE) a coordinate parallel ellipsoidal chunk geometry will be created. 2) by defining three points a non-coordinate parallel ellipsoidal chunk will be created. The points are defined in the input file by longitude:latitude. It is also possible to define additional subdivisions of the mesh in each direction. The boundary of the domain is formed by linear interpolation in longitude-latitude space between adjacent points (i.e. [lon, lat](f) = [lon1*f + lon2*(1-f), lat1*f + lat2*(1-f)], where f is a value between 0 and 1). Faces of the model are defined as 0, west; 1,east; 2, south; 3, north; 4, inner; 5, outer.

This geometry model supports initial topography for deforming the initial mesh.

`sphere': A geometry model for a sphere with a user specified radius. This geometry has only a single boundary, so the only valid boundary indicator to specify in input files is ``0''. It can also be referenced by the symbolic name ``surface'' in input files.

Despite the name, this geometry does not imply the use of a spherical coordinate system when used in 2d. Indeed, in 2d the geometry is simply a circle in a Cartesian coordinate system and consequently would correspond to a cross section of the fluid filled interior of an infinite cylinder where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in Section~\ref{sec:meaning-of-2d}.

`spherical shell': A geometry representing a spherical shell or a piece of it. Inner and outer radii are read from the parameter file in subsection 'Spherical shell'.

The spherical shell may be generated as per the original code (with respect to the inner and outer radius, and an initial number of cells along circumference) or following a custom mesh scheme: list of radial values or number of slices. A surface mesh is first generated and refined as desired, before it is extruded radially. A list of radial values subdivides the spherical shell at specified radii. The number of slices subdivides the spherical shell into N slices of equal thickness. The custom spherical shell only works with an opening angle of 360 degrees.

Despite the name, this geometry does not imply the use of a spherical coordinate system when used in 2d. Indeed, in 2d the geometry is simply an annulus in a Cartesian coordinate system and consequently would correspond to a cross section of the fluid filled space between two infinite cylinders where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in Section~\ref{sec:meaning-of-2d}.

The model assigns boundary indicators as follows: In 2d, inner and outer boundaries get boundary indicators zero and one, and if the opening angle set in the input file is less than 360, then left and right boundaries are assigned indicators two and three. These boundaries can also be referenced using the symbolic names `inner', `outer' and (if applicable) `left', `right'.

In 3d, inner and outer indicators are treated as in 2d. If the opening angle is chosen as 90 degrees, i.e., the domain is the intersection of a spherical shell and the first octant, then indicator 2 is at the face $x=0$, 3 at $y=0$, and 4 at $z=0$. These last three boundaries can then also be referred to as `east', `west' and `south' symbolically in input files. 

(parameters:Geometry_20model:Box)=
## **Parameters in section** Geometry model/Box
(parameters:Geometry_20model:Box:Box_20origin_20X_20coordinate)=
### __Parameter name:__ Box origin X coordinate
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** X coordinate of box origin. Units: \si{\meter}. 

(parameters:Geometry_20model:Box:Box_20origin_20Y_20coordinate)=
### __Parameter name:__ Box origin Y coordinate
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Y coordinate of box origin. Units: \si{\meter}. 

(parameters:Geometry_20model:Box:Box_20origin_20Z_20coordinate)=
### __Parameter name:__ Box origin Z coordinate
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Z coordinate of box origin. This value is ignored if the simulation is in 2d. Units: \si{\meter}. 

(parameters:Geometry_20model:Box:X_20extent)=
### __Parameter name:__ X extent
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Extent of the box in x-direction. Units: \si{\meter}. 

(parameters:Geometry_20model:Box:X_20periodic)=
### __Parameter name:__ X periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in X direction 

(parameters:Geometry_20model:Box:X_20repetitions)=
### __Parameter name:__ X repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in X direction. 

(parameters:Geometry_20model:Box:Y_20extent)=
### __Parameter name:__ Y extent
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Extent of the box in y-direction. Units: \si{\meter}. 

(parameters:Geometry_20model:Box:Y_20periodic)=
### __Parameter name:__ Y periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in Y direction 

(parameters:Geometry_20model:Box:Y_20repetitions)=
### __Parameter name:__ Y repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in Y direction. 

(parameters:Geometry_20model:Box:Z_20extent)=
### __Parameter name:__ Z extent
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Extent of the box in z-direction. This value is ignored if the simulation is in 2d. Units: \si{\meter}. 

(parameters:Geometry_20model:Box:Z_20periodic)=
### __Parameter name:__ Z periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in Z direction 

(parameters:Geometry_20model:Box:Z_20repetitions)=
### __Parameter name:__ Z repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in Z direction. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators)=
## **Parameters in section** Geometry model/Box with lithosphere boundary indicators
(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Box_20origin_20X_20coordinate)=
### __Parameter name:__ Box origin X coordinate
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** X coordinate of box origin. Units: \si{\meter}. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Box_20origin_20Y_20coordinate)=
### __Parameter name:__ Box origin Y coordinate
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Y coordinate of box origin. Units: \si{\meter}. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Box_20origin_20Z_20coordinate)=
### __Parameter name:__ Box origin Z coordinate
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Z coordinate of box origin. This value is ignored if the simulation is in 2d. Units: \si{\meter}. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Lithospheric_20thickness)=
### __Parameter name:__ Lithospheric thickness
**Default value:** 0.2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The thickness of the lithosphere used to create additional boundary indicators to set specific boundary conditions for the lithosphere.  

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:X_20extent)=
### __Parameter name:__ X extent
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Extent of the box in x-direction. Units: \si{\meter}. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:X_20periodic)=
### __Parameter name:__ X periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in X direction. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:X_20periodic_20lithosphere)=
### __Parameter name:__ X periodic lithosphere
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in X direction in the lithosphere. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:X_20repetitions)=
### __Parameter name:__ X repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in X direction of the lower box. The same number of repetitions will be used in the upper box. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Y_20extent)=
### __Parameter name:__ Y extent
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Extent of the box in y-direction. Units: \si{\meter}. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Y_20periodic)=
### __Parameter name:__ Y periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in Y direction. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Y_20periodic_20lithosphere)=
### __Parameter name:__ Y periodic lithosphere
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in Y direction in the lithosphere. This value is ignored if the simulation is in 2d.  

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Y_20repetitions)=
### __Parameter name:__ Y repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in Y direction of the lower box. If the simulation is in 3d, the same number of repetitions will be used in the upper box. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Y_20repetitions_20lithosphere)=
### __Parameter name:__ Y repetitions lithosphere
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in Y direction in the lithosphere. This value is ignored if the simulation is in 3d. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Z_20extent)=
### __Parameter name:__ Z extent
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Extent of the box in z-direction. This value is ignored if the simulation is in 2d. Units: \si{\meter}. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Z_20periodic)=
### __Parameter name:__ Z periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the box should be periodic in Z direction. This value is ignored if the simulation is in 2d. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Z_20repetitions)=
### __Parameter name:__ Z repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in Z direction of the lower box. This value is ignored if the simulation is in 2d. 

(parameters:Geometry_20model:Box_20with_20lithosphere_20boundary_20indicators:Z_20repetitions_20lithosphere)=
### __Parameter name:__ Z repetitions lithosphere
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in Z direction in the lithosphere. This value is ignored if the simulation is in 2d. 

(parameters:Geometry_20model:Chunk)=
## **Parameters in section** Geometry model/Chunk
(parameters:Geometry_20model:Chunk:Chunk_20inner_20radius)=
### __Parameter name:__ Chunk inner radius
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Radius at the bottom surface of the chunk. Units: \si{\meter}. 

(parameters:Geometry_20model:Chunk:Chunk_20maximum_20latitude)=
### __Parameter name:__ Chunk maximum latitude
**Default value:** 1. 

**Pattern:** [Double -90...90 (inclusive)] 

**Documentation:** Maximum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees. 

(parameters:Geometry_20model:Chunk:Chunk_20maximum_20longitude)=
### __Parameter name:__ Chunk maximum longitude
**Default value:** 1. 

**Pattern:** [Double -180...360 (inclusive)] 

**Documentation:** Maximum longitude of the chunk. Units: degrees. 

(parameters:Geometry_20model:Chunk:Chunk_20minimum_20latitude)=
### __Parameter name:__ Chunk minimum latitude
**Default value:** 0. 

**Pattern:** [Double -90...90 (inclusive)] 

**Documentation:** Minimum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees. 

(parameters:Geometry_20model:Chunk:Chunk_20minimum_20longitude)=
### __Parameter name:__ Chunk minimum longitude
**Default value:** 0. 

**Pattern:** [Double -180...360 (inclusive)] 

**Documentation:** Minimum longitude of the chunk. Units: degrees. 

(parameters:Geometry_20model:Chunk:Chunk_20outer_20radius)=
### __Parameter name:__ Chunk outer radius
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Radius at the top surface of the chunk. Units: \si{\meter}. 

(parameters:Geometry_20model:Chunk:Latitude_20repetitions)=
### __Parameter name:__ Latitude repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in latitude. This value is ignored if the simulation is in 2d 

(parameters:Geometry_20model:Chunk:Longitude_20repetitions)=
### __Parameter name:__ Longitude repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in longitude. 

(parameters:Geometry_20model:Chunk:Radius_20repetitions)=
### __Parameter name:__ Radius repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in radius. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators)=
## **Parameters in section** Geometry model/Chunk with lithosphere boundary indicators
(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20inner_20radius)=
### __Parameter name:__ Chunk inner radius
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Radius at the bottom surface of the chunk. Units: \si{\meter}. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20maximum_20latitude)=
### __Parameter name:__ Chunk maximum latitude
**Default value:** 1. 

**Pattern:** [Double -90...90 (inclusive)] 

**Documentation:** Maximum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20maximum_20longitude)=
### __Parameter name:__ Chunk maximum longitude
**Default value:** 1. 

**Pattern:** [Double -180...360 (inclusive)] 

**Documentation:** Maximum longitude of the chunk. Units: degrees. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20middle_20boundary_20radius)=
### __Parameter name:__ Chunk middle boundary radius
**Default value:** 1 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Radius at the top surface of the lower chunk, where it merges with the upper chunk. Units: \si{\meter}. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20minimum_20latitude)=
### __Parameter name:__ Chunk minimum latitude
**Default value:** 0. 

**Pattern:** [Double -90...90 (inclusive)] 

**Documentation:** Minimum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20minimum_20longitude)=
### __Parameter name:__ Chunk minimum longitude
**Default value:** 0. 

**Pattern:** [Double -180...360 (inclusive)] 

**Documentation:** Minimum longitude of the chunk. Units: degrees. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Chunk_20outer_20radius)=
### __Parameter name:__ Chunk outer radius
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Radius at the top surface of the chunk. Units: \si{\meter}. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Inner_20chunk_20radius_20repetitions)=
### __Parameter name:__ Inner chunk radius repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in radial direction for the lower chunk. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Latitude_20repetitions)=
### __Parameter name:__ Latitude repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in latitude. This value is ignored if the simulation is in 2d 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Longitude_20repetitions)=
### __Parameter name:__ Longitude repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in longitude. 

(parameters:Geometry_20model:Chunk_20with_20lithosphere_20boundary_20indicators:Outer_20chunk_20radius_20repetitions)=
### __Parameter name:__ Outer chunk radius repetitions
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of cells in radial direction for the upper chunk. 

(parameters:Geometry_20model:Ellipsoidal_20chunk)=
## **Parameters in section** Geometry model/Ellipsoidal chunk
(parameters:Geometry_20model:Ellipsoidal_20chunk:Depth)=
### __Parameter name:__ Depth
**Default value:** 500000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Bottom depth of model region. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:Depth_20subdivisions)=
### __Parameter name:__ Depth subdivisions
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of subdivisions of the coarse (initial) mesh in depth. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:East_2dWest_20subdivisions)=
### __Parameter name:__ East_2dWest subdivisions
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of subdivisions of the coarse (initial) mesh in the East-West direction. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:Eccentricity)=
### __Parameter name:__ Eccentricity
**Default value:** 8.1819190842622e-2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Eccentricity of the ellipsoid. Zero is a perfect sphere, default (8.1819190842622e-2) is WGS84. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:NE_20corner)=
### __Parameter name:__ NE corner
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Longitude:latitude in degrees of the North-East corner point of model region.The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:NW_20corner)=
### __Parameter name:__ NW corner
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Longitude:latitude in degrees of the North-West corner point of model region. The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:North_2dSouth_20subdivisions)=
### __Parameter name:__ North_2dSouth subdivisions
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of subdivisions of the coarse (initial) mesh in the North-South direction. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:SE_20corner)=
### __Parameter name:__ SE corner
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Longitude:latitude in degrees of the South-East corner point of model region. The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:SW_20corner)=
### __Parameter name:__ SW corner
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Longitude:latitude in degrees of the South-West corner point of model region. The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel. 

(parameters:Geometry_20model:Ellipsoidal_20chunk:Semi_2dmajor_20axis)=
### __Parameter name:__ Semi_2dmajor axis
**Default value:** 6378137.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The semi-major axis (a) of an ellipsoid. This is the radius for a sphere (eccentricity=0). Default WGS84 semi-major axis. 

(parameters:Geometry_20model:Initial_20topography_20model)=
## **Parameters in section** Geometry model/Initial topography model
(parameters:Geometry_20model:Initial_20topography_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** zero topography 

**Pattern:** [Selection ascii data|function|prm polygon|zero topography ] 

**Documentation:** Select one of the following models:

`ascii data': Implementation of a model in which the surface topography is derived from a file containing data in ascii format. The following geometry models are currently supported: box, chunk. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `Topography [m]' in a 2d model and  `x', `y', `Topography [m]' in a 3d model, which means that there has to be a single column containing the topography. Note that the data in the input file needs to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle in radians  and `y' by the polar angle in radians measured positive from the north pole. The grid will be assumed to be a longitude-colatitude grid. Note that the order of spherical coordinates is `phi', `theta' and not `theta', `phi', since this allows for dimension independent expressions.

`function': Implementation of a model in which the initial topography is described by a function in cartesian or spherical coordinates. 

`prm polygon': An initial topography model that defines the initial topography as constant inside each of a set of polygonal parts of the surface. The polygons, and their associated surface elevation, are defined in the `Geometry model/Initial topography/Prm polygon' section.

`zero topography': Implementation of a model in which the initial topography is zero.  

(parameters:Geometry_20model:Initial_20topography_20model:Ascii_20data_20model)=
## **Parameters in section** Geometry model/Initial topography model/Ascii data model
(parameters:Geometry_20model:Initial_20topography_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Geometry_20model:Initial_20topography_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.0.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Geometry_20model:Initial_20topography_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Geometry_20model:Initial_20topography_20model:Function)=
## **Parameters in section** Geometry model/Initial topography model/Function
(parameters:Geometry_20model:Initial_20topography_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian' and `spherical'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle.  

(parameters:Geometry_20model:Initial_20topography_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Geometry_20model:Initial_20topography_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Geometry_20model:Initial_20topography_20model:Function:Maximum_20topography_20value)=
### __Parameter name:__ Maximum topography value
**Default value:** 2000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum value the topography given by the function can take.  

(parameters:Geometry_20model:Initial_20topography_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Geometry_20model:Initial_20topography_20model:Prm_20polygon)=
## **Parameters in section** Geometry model/Initial topography model/Prm polygon
(parameters:Geometry_20model:Initial_20topography_20model:Prm_20polygon:Topography_20parameters)=
### __Parameter name:__ Topography parameters
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Set the topography height and the polygon which should be set to that height. The format is : "The topography height 	extgreater The point list describing a polygon \& The next topography height 	extgreater the next point list describing a polygon." The format for the point list describing the polygon is "x1,y1;x2,y2". For example for two triangular areas of 100 and -100 meters high set: '100 	extgreater 0,0;5,5;0,10 \& -100 	extgreater 10,10;10,15;20,15'. Units of the height are always in meters. The units of the coordinates are dependent on the geometry model. In the box model they are in meters, in the chunks they are in degrees, etc. Please refer to the manual of the individual geometry model to so see how the topography is implemented. 

(parameters:Geometry_20model:Sphere)=
## **Parameters in section** Geometry model/Sphere
(parameters:Geometry_20model:Sphere:Radius)=
### __Parameter name:__ Radius
**Default value:** 6371000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Radius of the sphere. Units: \si{\meter}. 

(parameters:Geometry_20model:Spherical_20shell)=
## **Parameters in section** Geometry model/Spherical shell
(parameters:Geometry_20model:Spherical_20shell:Cells_20along_20circumference)=
### __Parameter name:__ Cells along circumference
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of cells in circumferential direction that are created in the coarse mesh in 2d. If zero, this number is chosen automatically in a way that produces meshes in which cells have a reasonable aspect ratio for models in which the depth of the mantle is roughly that of the Earth. For planets with much shallower mantles and larger cores, you may want to chose a larger number to avoid cells that are elongated in tangential and compressed in radial direction.

In 3d, the number of cells is computed differently and does not have an easy interpretation. Valid values for this parameter in 3d are 0 (let this class choose), 6, 12 and 96. Other possible values may be discussed in the documentation of the deal.II function GridGenerator::hyper_shell. The parameter is best left at its default in 3d.

In either case, this parameter is ignored unless the opening angle of the domain is 360 degrees. This parameter is also ignored when using a custom mesh subdivision scheme. 

(parameters:Geometry_20model:Spherical_20shell:Custom_20mesh_20subdivision)=
### __Parameter name:__ Custom mesh subdivision
**Default value:** none 

**Pattern:** [Selection none|list of radial values|number of slices ] 

**Documentation:** Choose how the spherical shell mesh is generated. By default, a coarse mesh is generated with respect to the inner and outer radius, and an initial number of cells along circumference. In the other cases, a surface mesh is first generated and refined as desired, before it is extruded radially following the specified subdivision scheme. 

(parameters:Geometry_20model:Spherical_20shell:Initial_20lateral_20refinement)=
### __Parameter name:__ Initial lateral refinement
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Initial lateral refinement for the custom mesh subdivision schemes.The number of refinement steps performed on the initial coarse surface mesh, before the surface is extruded radially. This parameter allows the user more control over the ratio between radial and lateral refinement of the mesh. 

(parameters:Geometry_20model:Spherical_20shell:Inner_20radius)=
### __Parameter name:__ Inner radius
**Default value:** 3481000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Inner radius of the spherical shell. Units: \si{\meter}.

\note{The default value of 3,481,000 m equals the radius of a sphere with equal volume as Earth (i.e., 6371 km) minus the average depth of the core-mantle boundary (i.e., 2890 km).} 

(parameters:Geometry_20model:Spherical_20shell:List_20of_20radial_20values)=
### __Parameter name:__ List of radial values
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of radial values for the custom mesh scheme. Units: $\si{m}$. A list of radial values subdivides the spherical shell at specified radii. The list must be strictly ascending, and the first value must be greater than the inner radius while the last must be less than the outer radius. 

(parameters:Geometry_20model:Spherical_20shell:Number_20of_20slices)=
### __Parameter name:__ Number of slices
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of slices for the custom mesh subdivision scheme. The number of slices subdivides the spherical shell into N slices of equal thickness. Must be greater than 0. 

(parameters:Geometry_20model:Spherical_20shell:Opening_20angle)=
### __Parameter name:__ Opening angle
**Default value:** 360. 

**Pattern:** [Double 0...360 (inclusive)] 

**Documentation:** Opening angle in degrees of the section of the shell that we want to build. The only opening angles that are allowed for this geometry are 90, 180, and 360 in 2d; and 90 and 360 in 3d. Units: degrees. 

(parameters:Geometry_20model:Spherical_20shell:Outer_20radius)=
### __Parameter name:__ Outer radius
**Default value:** 6336000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Outer radius of the spherical shell. Units: \si{\meter}.

\note{The default value of 6,336,000 m equals the radius of a sphere with equal volume as Earth (i.e., 6371 km) minus the average depth of the mantle-crust interface (i.e., 35 km).} 

(parameters:Geometry_20model:Spherical_20shell:Phi_20periodic)=
### __Parameter name:__ Phi periodic
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether the shell should be periodic in the phi direction. 

(parameters:Gravity_20model)=
## **Parameters in section** Gravity model
(parameters:Gravity_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection ascii data|function|radial constant|radial earth-like|radial linear|vertical|unspecified ] 

**Documentation:** Select one of the following models:

`ascii data': Gravity is read from a file that describes the reference state. The default profile follows the preliminary reference Earth model (PREM, Dziewonski and Anderson, 1981). Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of points in the reference state as for example `# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide a column named `gravity'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

`function': Gravity is given in terms of an explicit formula that is elaborated in the parameters in section ``Gravity model|Function''. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`radial constant': A gravity model in which the gravity has a constant magnitude and the direction is radial (pointing inward if the value is positive). The magnitude is read from the parameter file in subsection 'Radial constant'.

`radial earth-like': This plugin has been removed due to its misleading name. The included profile was hard-coded and was less earth-like than the `ascii data' plugin, which uses the profile of the Preliminary Reference Earth Model (PREM). Use `ascii data' instead of `radial earth-like'.

`radial linear': A gravity model which is radial (pointing inward if the gravity is positive) and the magnitude changes linearly with depth. The magnitude of gravity at the surface and bottom is read from the input file in a section ``Gravity model/Radial linear''.

`vertical': A gravity model in which the gravity direction is vertical (pointing downward for positive values) and at a constant magnitude by default equal to one. 

(parameters:Gravity_20model:Ascii_20data_20model)=
## **Parameters in section** Gravity model/Ascii data model
(parameters:Gravity_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/gravity-model/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Gravity_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** prem.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Gravity_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Gravity_20model:Function)=
## **Parameters in section** Gravity model/Function
(parameters:Gravity_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Gravity_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Gravity_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Gravity_20model:Radial_20constant)=
## **Parameters in section** Gravity model/Radial constant
(parameters:Gravity_20model:Radial_20constant:Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 9.81 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Magnitude of the gravity vector in $m/s^2$. For positive values the direction is radially inward towards the center of the earth. 

(parameters:Gravity_20model:Radial_20linear)=
## **Parameters in section** Gravity model/Radial linear
(parameters:Gravity_20model:Radial_20linear:Magnitude_20at_20bottom)=
### __Parameter name:__ Magnitude at bottom
**Default value:** 10.7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Magnitude of the radial gravity vector at the bottom of the domain. `Bottom' means themaximum depth in the chosen geometry, and for example represents the core-mantle boundary in the case of the `spherical shell' geometry model, and the center in the case of the `sphere' geometry model. Units: \si{\meter\per\second\squared}. 

(parameters:Gravity_20model:Radial_20linear:Magnitude_20at_20surface)=
### __Parameter name:__ Magnitude at surface
**Default value:** 9.8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Magnitude of the radial gravity vector at the surface of the domain. Units: \si{\meter\per\second\squared}. 

(parameters:Gravity_20model:Vertical)=
## **Parameters in section** Gravity model/Vertical
(parameters:Gravity_20model:Vertical:Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Value of the gravity vector in $m/s^2$ directed along negative y (2D) or z (3D) axis (if the magnitude is positive. 

(parameters:Heating_20model)=
## **Parameters in section** Heating model
(parameters:Heating_20model:List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**  

**Pattern:** [MultipleSelection adiabatic heating|adiabatic heating of melt|compositional heating|constant heating|function|latent heat|latent heat melt|radioactive decay|shear heating|shear heating with melt ] 

**Documentation:** A comma separated list of heating models that will be used to calculate the heating terms in the energy equation. The results of each of these criteria, i.e., the heating source terms and the latent heat terms for the left hand side will be added.

The following heating models are available:

`adiabatic heating': Implementation of a standard and a simplified model of adiabatic heating.

`adiabatic heating of melt': Implementation of a standard and a simplified model of adiabatic heating of melt. The full model implements the heating term 
$\alpha T (-\phi \mathbf u_s \cdot \nabla p) + \alpha T (\phi \mathbf u_f \cdot \nabla p)$.
For full adiabatic heating, this has to be used in combination with the heating model `adiabatic heating' to also include adiabatic heating for the solid part, and the full heating term is then $\alpha T ((1-\phi) \mathbf u_s \cdot \nabla p) + \alpha T (\phi \mathbf u_f \cdot \nabla p)$.

`compositional heating': Implementation of a model in which magnitude of internal heat production is determined from fixed values assigned to each compositional field. These values are interpreted as having units \si{\watt\per\meter\cubed}.

`constant heating': Implementation of a model in which the heating rate is constant.

`function': Implementation of a model in which the heating rate is given in terms of an explicit formula that is elaborated in the parameters in section ``Heating model|Function''. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

The formula is interpreted as having units W/kg.

Since the symbol $t$ indicating time may appear in the formulas for the heating rate, it is interpreted as having units seconds unless the global parameter ``Use years in output instead of seconds'' is set.

`latent heat': Implementation of a standard model for latent heat.

`latent heat melt': Implementation of a standard model for latent heat of melting. This assumes that there is a compositional field called porosity, and it uses the reaction term of this field (the fraction of material that melted in the current time step) multiplied by a constant entropy change for melting all of the material as source term of the heating model.
If there is no field called porosity, the heating terms are 0.

`radioactive decay': Implementation of a model in which the internal heating rate is radioactive decaying in the following rule:
\[(\text{initial concentration})\cdot 0.5^{\text{time}/(\text{half life})}\]
The crust and mantle can have different concentrations, and the crust can be defined either by depth or by a certain compositional field.
The formula is interpreted as having units W/kg.

`shear heating': Implementation of a standard model for shear heating. Adds the term: $  2 \eta \left( \varepsilon - \frac{1}{3} \text{tr} \varepsilon \mathbf 1 \right) : \left( \varepsilon - \frac{1}{3} \text{tr} \varepsilon \mathbf 1 \right)$ to the right-hand side of the temperature equation.

`shear heating with melt': Implementation of a standard model for shear heating of migrating melt, including bulk (compression) heating $\xi \left( \nabla \cdot \mathbf u_s \right)^2 $ and heating due to melt segregation $\frac{\eta_f \phi^2}{k} \left( \mathbf u_f - \mathbf u_s \right)^2 $. For full shear heating, this has to be used in combination with the heating model shear heating to also include shear heating for the solid part. 

(parameters:Heating_20model:Adiabatic_20heating)=
## **Parameters in section** Heating model/Adiabatic heating
(parameters:Heating_20model:Adiabatic_20heating:Use_20simplified_20adiabatic_20heating)=
### __Parameter name:__ Use simplified adiabatic heating
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** A flag indicating whether the adiabatic heating should be simplified from $\alpha T (\mathbf u \cdot \nabla p)$ to $ \alpha \rho T (\mathbf u \cdot \mathbf g) $. 

(parameters:Heating_20model:Adiabatic_20heating_20of_20melt)=
## **Parameters in section** Heating model/Adiabatic heating of melt
(parameters:Heating_20model:Adiabatic_20heating_20of_20melt:Use_20simplified_20adiabatic_20heating)=
### __Parameter name:__ Use simplified adiabatic heating
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** A flag indicating whether the adiabatic heating should be simplified from $\alpha T (\mathbf u \cdot \nabla p)$ to $ \alpha \rho T (\mathbf u \cdot \mathbf g) $. 

(parameters:Heating_20model:Compositional_20heating)=
## **Parameters in section** Heating model/Compositional heating
(parameters:Heating_20model:Compositional_20heating:Compositional_20heating_20values)=
### __Parameter name:__ Compositional heating values
**Default value:** 0. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of heat production per unit volume values for background and compositional fields, for a total of N+1 values, where the first value corresponds to the background material, and N is the number of compositional fields. Units: \si{\watt\per\meter\cubed}. 

(parameters:Heating_20model:Compositional_20heating:Use_20compositional_20field_20for_20heat_20production_20averaging)=
### __Parameter name:__ Use compositional field for heat production averaging
**Default value:** 1 

**Pattern:** [List of <[Integer range 0...1 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of integers with as many entries as compositional fields plus one. The first entry corresponds to the background material, each following entry corresponds to a particular compositional field. If the entry for a field is '1' this field is considered during the computation of volume fractions, if it is '0' the field is ignored. This is useful if some compositional fields are used to track properties like finite strain that should not contribute to heat production. The first entry determines whether the background field contributes to heat production or not (essentially similar to setting its 'Compositional heating values' to zero, but included for consistency in the length of the input lists). 

(parameters:Heating_20model:Constant_20heating)=
## **Parameters in section** Heating model/Constant heating
(parameters:Heating_20model:Constant_20heating:Radiogenic_20heating_20rate)=
### __Parameter name:__ Radiogenic heating rate
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The specific rate of heating due to radioactive decay (or other bulk sources you may want to describe). This parameter corresponds to the variable $H$ in the temperature equation stated in the manual, and the heating term is $ho H$. Units: W/kg. 

(parameters:Heating_20model:Function)=
## **Parameters in section** Heating model/Function
(parameters:Heating_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Heating_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Heating_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Heating_20model:Latent_20heat_20melt)=
## **Parameters in section** Heating model/Latent heat melt
(parameters:Heating_20model:Latent_20heat_20melt:Melting_20entropy_20change)=
### __Parameter name:__ Melting entropy change
**Default value:** -300. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The entropy change for the phase transition from solid to melt. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Heating_20model:Radioactive_20decay)=
## **Parameters in section** Heating model/Radioactive decay
(parameters:Heating_20model:Radioactive_20decay:Crust_20composition_20number)=
### __Parameter name:__ Crust composition number
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Which composition field should be treated as crust 

(parameters:Heating_20model:Radioactive_20decay:Crust_20defined_20by_20composition)=
### __Parameter name:__ Crust defined by composition
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether crust defined by composition or depth 

(parameters:Heating_20model:Radioactive_20decay:Crust_20depth)=
### __Parameter name:__ Crust depth
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Depth of the crust when crust if defined by depth. Units: \si{\meter}. 

(parameters:Heating_20model:Radioactive_20decay:Half_20decay_20times)=
### __Parameter name:__ Half decay times
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Half decay times. Units: (Seconds), or (Years) if set `use years instead of seconds'. 

(parameters:Heating_20model:Radioactive_20decay:Heating_20rates)=
### __Parameter name:__ Heating rates
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Heating rates of different elements (W/kg) 

(parameters:Heating_20model:Radioactive_20decay:Initial_20concentrations_20crust)=
### __Parameter name:__ Initial concentrations crust
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Initial concentrations of different elements (ppm) 

(parameters:Heating_20model:Radioactive_20decay:Initial_20concentrations_20mantle)=
### __Parameter name:__ Initial concentrations mantle
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Initial concentrations of different elements (ppm) 

(parameters:Heating_20model:Radioactive_20decay:Number_20of_20elements)=
### __Parameter name:__ Number of elements
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Number of radioactive elements 

(parameters:Initial_20composition_20model)=
## **Parameters in section** Initial composition model
(parameters:Initial_20composition_20model:List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**  

**Pattern:** [MultipleSelection adiabatic density|ascii data|ascii data layered|function|porosity|world builder ] 

**Documentation:** A comma-separated list of initial composition models that together describe the initial composition field. These plugins are loaded in the order given, and modify the existing composition field via the operators listed in 'List of model operators'.

The following composition models are available:

`adiabatic density': Specify the initial composition as the adiabatic reference density at each position. Note that only the field of the type 'density' will be filled. For all other fields this plugin returns 0.0.

`ascii data': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2', etc. in a 2d model and `x', `y', `z', `composition1', `composition2', etc. in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model.Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`ascii data layered': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Each file defines a surface on which compositional fields are defined. Between the surfaces, the fields can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2' etc. in a 2d model and `x', `y', `z', `composition1', `composition2' etc. in a 3d model; i.e. the columns before the compositional field always contains the position of the surface along the vertical direction. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle and `y' (if 3D) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3D is `phi', `theta', `r', `T'and not `theta', `phi', `r', `T' as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

`function': Specify the composition in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`porosity': A class that implements initial conditions for the porosity field by computing the equilibrium melt fraction for the given initial condition and reference pressure profile. Note that this plugin only works if there is a compositional field called `porosity', and the used material model implements the 'MeltFractionModel' interface. For all compositional fields except porosity this plugin returns 0.0, and they are therefore not changed as long as the default `add' operator is selected for this plugin.

`world builder': Specify the initial composition through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter 'World builder file'. It is possible to use the World Builder only for selected compositional fields by specifying the parameter 'List of relevant compositions'. 

(parameters:Initial_20composition_20model:List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add 

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ] 

**Documentation:** A comma-separated list of operators that will be used to append the listed composition models onto the previous models. If only one operator is given, the same operator is applied to all models. 

(parameters:Initial_20composition_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection adiabatic density|ascii data|ascii data layered|function|porosity|world builder|unspecified ] 

**Documentation:** Select one of the following models:

`adiabatic density': Specify the initial composition as the adiabatic reference density at each position. Note that only the field of the type 'density' will be filled. For all other fields this plugin returns 0.0.

`ascii data': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2', etc. in a 2d model and `x', `y', `z', `composition1', `composition2', etc. in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model.Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`ascii data layered': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Each file defines a surface on which compositional fields are defined. Between the surfaces, the fields can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2' etc. in a 2d model and `x', `y', `z', `composition1', `composition2' etc. in a 3d model; i.e. the columns before the compositional field always contains the position of the surface along the vertical direction. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle and `y' (if 3D) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3D is `phi', `theta', `r', `T'and not `theta', `phi', `r', `T' as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

`function': Specify the composition in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`porosity': A class that implements initial conditions for the porosity field by computing the equilibrium melt fraction for the given initial condition and reference pressure profile. Note that this plugin only works if there is a compositional field called `porosity', and the used material model implements the 'MeltFractionModel' interface. For all compositional fields except porosity this plugin returns 0.0, and they are therefore not changed as long as the default `add' operator is selected for this plugin.

`world builder': Specify the initial composition through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter 'World builder file'. It is possible to use the World Builder only for selected compositional fields by specifying the parameter 'List of relevant compositions'.

\textbf{Warning}: This parameter provides an old and deprecated way of specifying initial composition models and shouldn't be used. Please use 'List of model names' instead. 

(parameters:Initial_20composition_20model:Volume_20of_20fluid_20initialization_20type)=
### __Parameter name:__ Volume of fluid initialization type
**Default value:**  

**Pattern:** [Map of <[Anything]>:<[Selection composition|level set ]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list denoting the method to be used to initialize a composition field specified to be advected using the volume of fluid method.

The format of valid entries for this parameter is that of a map given as ``key1:value1, key2:value2`` where each key must be the name of a compositional field using the volume of fluid advection method, and the value is one of ``composition`` or ``level set``. ``composition`` is the default

When ``composition is specified, the initial model is treated as a standard composition field with bounds between 0 and 1 assumed, The initial fluid fractions are then based on an iterated midpoint quadrature. Resultant volume fractions outside of the bounds will be coerced to the nearest valid value (ie 0 or 1). If ``level set`` is specified, the intial data will be assumed to be in the form of a signed distance level set function (i.e. a function which is positive when in the fluid, negative outside, and zero on the interface and the magnitude is always the distance to the interface so the gradient is one everywhere). 

(parameters:Initial_20composition_20model:Ascii_20data_20model)=
## **Parameters in section** Initial composition model/Ascii data model
(parameters:Initial_20composition_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.  

(parameters:Initial_20composition_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** initial_composition_top_mantle_box_3d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20composition_20model:Ascii_20data_20model:Data_20file_20names)=
### __Parameter name:__ Data file names
**Default value:** initial_composition_top_mantle_box_3d.txt 

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** The file names of the model data (comma separated).  

(parameters:Initial_20composition_20model:Ascii_20data_20model:First_20point_20on_20slice)=
### __Parameter name:__ First point on slice
**Default value:** 0.0,1.0,0.0 

**Pattern:** [Anything] 

**Documentation:** Point that determines the plane in which the 2D slice lies in. This variable is only used if 'Slice dataset in 2D plane' is true. The slice will go through this point, the point defined by the parameter 'Second point on slice', and the center of the model domain. After the rotation, this first point will lie along the (0,1,0) axis of the coordinate system. The coordinates of the point have to be given in Cartesian coordinates. 

(parameters:Initial_20composition_20model:Ascii_20data_20model:Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** linear 

**Pattern:** [Selection piecewise constant|linear ] 

**Documentation:** Method to interpolate between layer boundaries. Select from piecewise constant or linear. Piecewise constant takes the value from the nearest layer boundary above the data point. The linear option interpolates linearly between layer boundaries. Above and below the domain given by the layer boundaries, the values aregiven by the top and bottom layer boundary. 

(parameters:Initial_20composition_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20composition_20model:Ascii_20data_20model:Second_20point_20on_20slice)=
### __Parameter name:__ Second point on slice
**Default value:** 1.0,0.0,0.0 

**Pattern:** [Anything] 

**Documentation:** Second point that determines the plane in which the 2D slice lies in. This variable is only used if 'Slice dataset in 2D plane' is true. The slice will go through this point, the point defined by the parameter 'First point on slice', and the center of the model domain. The coordinates of the point have to be given in Cartesian coordinates. 

(parameters:Initial_20composition_20model:Ascii_20data_20model:Slice_20dataset_20in_202D_20plane)=
### __Parameter name:__ Slice dataset in 2D plane
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a 2D data slice of a 3D data file or the entire data file. Slicing a 3D dataset is only supported for 2D models. 

(parameters:Initial_20composition_20model:Function)=
## **Parameters in section** Initial composition model/Function
(parameters:Initial_20composition_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Initial_20composition_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Initial_20composition_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Initial_20composition_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Initial_20composition_20model:World_20builder)=
## **Parameters in section** Initial composition model/World builder
(parameters:Initial_20composition_20model:World_20builder:List_20of_20relevant_20compositions)=
### __Parameter name:__ List of relevant compositions
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of names of compositional fields for which to determine the initial composition using the World Builder. As World Builder evaluations can be expensive, this parameter allows to only evaluate the fields that are relevant. This plugin returns 0.0 for all compositions that are not selected in the list. By default the list is empty and the world builder is evaluated for all compositional fields. 

(parameters:Initial_20temperature_20model)=
## **Parameters in section** Initial temperature model
(parameters:Initial_20temperature_20model:List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**  

**Pattern:** [MultipleSelection S40RTS perturbation|SAVANI perturbation|adiabatic|adiabatic boundary|ascii data|ascii data layered|ascii profile|continental geotherm|function|harmonic perturbation|inclusion shape perturbation|lithosphere mask|mandelbox|patch on S40RTS|perturbed box|polar box|spherical gaussian perturbation|spherical hexagonal perturbation|world builder ] 

**Documentation:** A comma-separated list of initial temperature models that will be used to initialize the temperature. These plugins are loaded in the order given, and modify the existing temperature field via the operators listed in 'List of model operators'.

The following initial temperature models are available:

`S40RTS perturbation': An initial temperature field in which the temperature is perturbed following the S20RTS or S40RTS shear wave velocity model by Ritsema and others, which can be downloaded here \url{http://www.earth.lsa.umich.edu/~jritsema/research.html}. Information on the vs model can be found in Ritsema, J., Deuss, A., van Heijst, H.J. \& Woodhouse, J.H., 2011. S40RTS: a degree-40 shear-velocity model for the mantle from new Rayleigh wave dispersion, teleseismic traveltime and normal-mode splitting function measurements, Geophys. J. Int. 184, 1223-1236. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the 'Vs to density scaling' parameter or depth-dependent and read in from a file. To convert density the user can specify the 'Thermal expansion coefficient in initial temperature scaling' parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the `vs to density scaling' parameter and $\alpha$ is the 'Thermal expansion coefficient in initial temperature scaling' parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model). If a depth is specified in 'Remove temperature heterogeneity down to specified depth', there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of points in the reference state as for example '# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named `depth' and `vs\_to\_density'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.
If the plugin is used in 2D it will use an equatorial slice of the seismic tomography model.

`SAVANI perturbation': An initial temperature field in which the temperature is perturbed following the SAVANI shear wave velocity model by Auer and others, which can be downloaded here \url{http://n.ethz.ch/~auerl/savani.tar.bz2}. Information on the vs model can be found in Auer, L., Boschi, L., Becker, T.W., Nissen-Meyer, T. \& Giardini, D., 2014. Savani: A variable resolution whole-mantle model of anisotropic shear velocity variations based on multiple data sets. Journal of Geophysical Research: Solid Earth 119.4 (2014): 3006-3034. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the 'Vs to density scaling' parameter or depth-dependent and read in from a file. To convert density the user can specify the 'Thermal expansion coefficient in initial temperature scaling' parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the `vs to density scaling' parameter and $\alpha$ is the 'Thermal expansion coefficient in initial temperature scaling' parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).If a depth is specified in 'Remove temperature heterogeneity down to specified depth', there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of points in the reference state as for example '# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named `depth' and `vs\_to\_density'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

`adiabatic': Temperature is prescribed as an adiabatic profile with upper and lower thermal boundary layers, whose ages are given as input parameters. Note that this plugin uses the 'Adiabatic conditions model' to compute the adiabat. Thus, the results depend on variables defined outside of this specific subsection; e.g. the globally defined 'Adiabatic surface temperature', and the variables defined in the 'Material model' section including densities, heat capacities and thermal expansivities.

`adiabatic boundary': An initial temperature condition that allows for discretizing the model domain into two layers separated by a user-defined isothermal boundary. The user includes an input ascii data file that is formatted as 3 columns of `longitude(radians)', `colatitude(radians)', and `isotherm depth(meters)', where `isotherm depth' represents the depth of an initial temperature of 1673.15 K (by default). The first lines in the data file may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 69 121'. Note that the coordinates need to be sorted in a specific order: the `longitude' coordinate needs to ascend first, followed by the `colatitude' coordinate in order to assign the correct data (isotherm depth) to the prescribed coordinates. The temperature is defined from the surface (273.15 K) to the isotherm depth (1673.15 K) as a linear gradient. Below the isotherm depth the temperature increases approximately adiabatically (0.0005 K per meter). This plugin should work for all geometry models, but is currently only tested for spherical models.

`ascii data': Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `Temperature [K]' in a 2d model and  `x', `y', `z', `Temperature [K]' in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`ascii data layered': Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Each file defines a surface on which temperature is defined. Between the surfaces, the temperatures can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `Temperature [K]' in a 2d model and `x', `y', `z', `Temperature [K]' in a 3d model; i.e. the last two columns always contain the position of the isotherm along the vertical direction, and the temperature at that point. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle and `y' (if 3D) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3D is `phi', `theta', `r', `T'and not `theta', `phi', `r', `T' as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

`ascii profile': Implementation of a model in which the initial temperature is read from a file that provides these values as a function of depth. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of points in the temperature profile, for example `# POINTS: 10'. Following the comment lines, there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide columns named `depth' and`temperature'.Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

`continental geotherm': This is a temperature initial condition that computes a continental geotherm based on the solution of the steady-state conductive equation $k\frac{d^2 T}{dy^2}+\rho H = 0$ as described in e.g. Turcotte and Schubert, Ch. 4.6, or Chapman (1986). As boundary conditions, we take the surface temperature and the temperature of the Lithosphere-Asthenosphere Boundary (LAB). 
The geotherm is computed for a homogeneous lithosphere composed of an upper crust, lower crust and mantle layer. The crustal layers are assumed to have a constant radioactive heating, and all layers are assumed to have a constant thermal conductivity. Layer thicknesses, surface temperature and LAB temperature should be specified by the user. For consistency, the density, heat production and thermal conductivity of each layer are read from the visco plastic material model and the compositional heating model. 
For any depths below the depth of the LAB, a unrealistically high temperature is returned, such that this plugin can be combined with another temperature plugin through the 'minimum' operator. 
Note that the current implementation only works for a 3-layer lithosphere, even though in principle the heat conduction equation can be solved for any number of layers. The naming of the compositional fields that represent the layers is also very specific, namely `upper\_crust', `lower\_crust', and `lithospheric\_mantle'. 
Make sure the top and bottom temperatures of the lithosphere agree with temperatures set in for example the temperature boundary conditions.

`function': Specify the initial temperature in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`harmonic perturbation': An initial temperature field in which the temperature is perturbed following a harmonic function (spherical harmonic or sine depending on geometry and dimension) in lateral and radial direction from an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).

`inclusion shape perturbation': An initial temperature field in which there is an inclusion in a constant-temperature box field. The size, shape, gradient, position, and temperature of the inclusion are defined by parameters.

`lithosphere mask': Implementation of a model in which the initial temperature is set to a specified lithosphere temperature above the lithosphere-asthenosphere boundary (specified by an ascii file or maximum lithosphere depth value). Below this the initial temperature is set as NaN.  Note the required format of the input data file: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of grid points in each dimension as for example '# POINTS: 3 3'. For a spherical model, the order of the data columns has to be 'phi', 'theta', 'depth (m)', where phi is the azimuth angle and theta is the polar angle measured positive from the north pole. This plug-in can be combined with another using the 'replace if valid' operator. 

`mandelbox': Fractal-shaped temperature field.

`patch on S40RTS': Implementation of a model in which the initial temperature is derived from a file containing shear wave velocity perturbations in ascii format (e.g. a high resolution upper mantle tomography) combined with S40RTS. Note the required format of the input ascii input data: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of grid points in each dimension as for example '# POINTS: 3 3 3'. The order of the data columns has to be  `x', `y', `z', 'Vs Perturbation' in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. In the spherical model data will be handled as Cartesian, however, `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions. See S40RTS documentation for details on input parameters in the S40RTS perturbation subsection. The boundary between the two tomography models is smoothed using a depth weighted combination of Vs values within the region of smoothing. 

`perturbed box': An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is chosen in such a way that the initial temperature is constant to one along the entire boundary.

`polar box': An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is such that there are two poles on opposing corners of the box. 

`spherical gaussian perturbation': An initial temperature field in which the temperature is perturbed by a single Gaussian added to an otherwise spherically symmetric state. Additional parameters are read from the parameter file in subsection 'Spherical gaussian perturbation'.

`spherical hexagonal perturbation': An initial temperature field in which the temperature is perturbed following an $N$-fold pattern in a specified direction from an otherwise spherically symmetric state. The class's name comes from previous versions when the only option was $N=6$.

`world builder': Specify the initial temperature through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter 'World builder file'. 

(parameters:Initial_20temperature_20model:List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add 

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ] 

**Documentation:** A comma-separated list of operators that will be used to append the listed temperature models onto the previous models. If only one operator is given, the same operator is applied to all models. 

(parameters:Initial_20temperature_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection S40RTS perturbation|SAVANI perturbation|adiabatic|adiabatic boundary|ascii data|ascii data layered|ascii profile|continental geotherm|function|harmonic perturbation|inclusion shape perturbation|lithosphere mask|mandelbox|patch on S40RTS|perturbed box|polar box|spherical gaussian perturbation|spherical hexagonal perturbation|world builder|unspecified ] 

**Documentation:** Select one of the following models:

`S40RTS perturbation': An initial temperature field in which the temperature is perturbed following the S20RTS or S40RTS shear wave velocity model by Ritsema and others, which can be downloaded here \url{http://www.earth.lsa.umich.edu/~jritsema/research.html}. Information on the vs model can be found in Ritsema, J., Deuss, A., van Heijst, H.J. \& Woodhouse, J.H., 2011. S40RTS: a degree-40 shear-velocity model for the mantle from new Rayleigh wave dispersion, teleseismic traveltime and normal-mode splitting function measurements, Geophys. J. Int. 184, 1223-1236. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the 'Vs to density scaling' parameter or depth-dependent and read in from a file. To convert density the user can specify the 'Thermal expansion coefficient in initial temperature scaling' parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the `vs to density scaling' parameter and $\alpha$ is the 'Thermal expansion coefficient in initial temperature scaling' parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model). If a depth is specified in 'Remove temperature heterogeneity down to specified depth', there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of points in the reference state as for example '# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named `depth' and `vs\_to\_density'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.
If the plugin is used in 2D it will use an equatorial slice of the seismic tomography model.

`SAVANI perturbation': An initial temperature field in which the temperature is perturbed following the SAVANI shear wave velocity model by Auer and others, which can be downloaded here \url{http://n.ethz.ch/~auerl/savani.tar.bz2}. Information on the vs model can be found in Auer, L., Boschi, L., Becker, T.W., Nissen-Meyer, T. \& Giardini, D., 2014. Savani: A variable resolution whole-mantle model of anisotropic shear velocity variations based on multiple data sets. Journal of Geophysical Research: Solid Earth 119.4 (2014): 3006-3034. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the 'Vs to density scaling' parameter or depth-dependent and read in from a file. To convert density the user can specify the 'Thermal expansion coefficient in initial temperature scaling' parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the `vs to density scaling' parameter and $\alpha$ is the 'Thermal expansion coefficient in initial temperature scaling' parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).If a depth is specified in 'Remove temperature heterogeneity down to specified depth', there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of points in the reference state as for example '# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named `depth' and `vs\_to\_density'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

`adiabatic': Temperature is prescribed as an adiabatic profile with upper and lower thermal boundary layers, whose ages are given as input parameters. Note that this plugin uses the 'Adiabatic conditions model' to compute the adiabat. Thus, the results depend on variables defined outside of this specific subsection; e.g. the globally defined 'Adiabatic surface temperature', and the variables defined in the 'Material model' section including densities, heat capacities and thermal expansivities.

`adiabatic boundary': An initial temperature condition that allows for discretizing the model domain into two layers separated by a user-defined isothermal boundary. The user includes an input ascii data file that is formatted as 3 columns of `longitude(radians)', `colatitude(radians)', and `isotherm depth(meters)', where `isotherm depth' represents the depth of an initial temperature of 1673.15 K (by default). The first lines in the data file may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 69 121'. Note that the coordinates need to be sorted in a specific order: the `longitude' coordinate needs to ascend first, followed by the `colatitude' coordinate in order to assign the correct data (isotherm depth) to the prescribed coordinates. The temperature is defined from the surface (273.15 K) to the isotherm depth (1673.15 K) as a linear gradient. Below the isotherm depth the temperature increases approximately adiabatically (0.0005 K per meter). This plugin should work for all geometry models, but is currently only tested for spherical models.

`ascii data': Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `Temperature [K]' in a 2d model and  `x', `y', `z', `Temperature [K]' in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`ascii data layered': Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Each file defines a surface on which temperature is defined. Between the surfaces, the temperatures can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `Temperature [K]' in a 2d model and `x', `y', `z', `Temperature [K]' in a 3d model; i.e. the last two columns always contain the position of the isotherm along the vertical direction, and the temperature at that point. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle and `y' (if 3D) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3D is `phi', `theta', `r', `T'and not `theta', `phi', `r', `T' as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

`ascii profile': Implementation of a model in which the initial temperature is read from a file that provides these values as a function of depth. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of points in the temperature profile, for example `# POINTS: 10'. Following the comment lines, there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide columns named `depth' and`temperature'.Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

`continental geotherm': This is a temperature initial condition that computes a continental geotherm based on the solution of the steady-state conductive equation $k\frac{d^2 T}{dy^2}+\rho H = 0$ as described in e.g. Turcotte and Schubert, Ch. 4.6, or Chapman (1986). As boundary conditions, we take the surface temperature and the temperature of the Lithosphere-Asthenosphere Boundary (LAB). 
The geotherm is computed for a homogeneous lithosphere composed of an upper crust, lower crust and mantle layer. The crustal layers are assumed to have a constant radioactive heating, and all layers are assumed to have a constant thermal conductivity. Layer thicknesses, surface temperature and LAB temperature should be specified by the user. For consistency, the density, heat production and thermal conductivity of each layer are read from the visco plastic material model and the compositional heating model. 
For any depths below the depth of the LAB, a unrealistically high temperature is returned, such that this plugin can be combined with another temperature plugin through the 'minimum' operator. 
Note that the current implementation only works for a 3-layer lithosphere, even though in principle the heat conduction equation can be solved for any number of layers. The naming of the compositional fields that represent the layers is also very specific, namely `upper\_crust', `lower\_crust', and `lithospheric\_mantle'. 
Make sure the top and bottom temperatures of the lithosphere agree with temperatures set in for example the temperature boundary conditions.

`function': Specify the initial temperature in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`harmonic perturbation': An initial temperature field in which the temperature is perturbed following a harmonic function (spherical harmonic or sine depending on geometry and dimension) in lateral and radial direction from an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).

`inclusion shape perturbation': An initial temperature field in which there is an inclusion in a constant-temperature box field. The size, shape, gradient, position, and temperature of the inclusion are defined by parameters.

`lithosphere mask': Implementation of a model in which the initial temperature is set to a specified lithosphere temperature above the lithosphere-asthenosphere boundary (specified by an ascii file or maximum lithosphere depth value). Below this the initial temperature is set as NaN.  Note the required format of the input data file: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of grid points in each dimension as for example '# POINTS: 3 3'. For a spherical model, the order of the data columns has to be 'phi', 'theta', 'depth (m)', where phi is the azimuth angle and theta is the polar angle measured positive from the north pole. This plug-in can be combined with another using the 'replace if valid' operator. 

`mandelbox': Fractal-shaped temperature field.

`patch on S40RTS': Implementation of a model in which the initial temperature is derived from a file containing shear wave velocity perturbations in ascii format (e.g. a high resolution upper mantle tomography) combined with S40RTS. Note the required format of the input ascii input data: The first lines may contain any number of comments if they begin with '#', but one of these lines needs to contain the number of grid points in each dimension as for example '# POINTS: 3 3 3'. The order of the data columns has to be  `x', `y', `z', 'Vs Perturbation' in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. In the spherical model data will be handled as Cartesian, however, `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions. See S40RTS documentation for details on input parameters in the S40RTS perturbation subsection. The boundary between the two tomography models is smoothed using a depth weighted combination of Vs values within the region of smoothing. 

`perturbed box': An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is chosen in such a way that the initial temperature is constant to one along the entire boundary.

`polar box': An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is such that there are two poles on opposing corners of the box. 

`spherical gaussian perturbation': An initial temperature field in which the temperature is perturbed by a single Gaussian added to an otherwise spherically symmetric state. Additional parameters are read from the parameter file in subsection 'Spherical gaussian perturbation'.

`spherical hexagonal perturbation': An initial temperature field in which the temperature is perturbed following an $N$-fold pattern in a specified direction from an otherwise spherically symmetric state. The class's name comes from previous versions when the only option was $N=6$.

`world builder': Specify the initial temperature through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter 'World builder file'.

\textbf{Warning}: This parameter provides an old and deprecated way of specifying initial temperature models and shouldn't be used. Please use 'List of model names' instead. 

(parameters:Initial_20temperature_20model:Adiabatic)=
## **Parameters in section** Initial temperature model/Adiabatic
(parameters:Initial_20temperature_20model:Adiabatic:Age_20bottom_20boundary_20layer)=
### __Parameter name:__ Age bottom boundary layer
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The age of the lower thermal boundary layer, used for the calculation of the half-space cooling model temperature. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Initial_20temperature_20model:Adiabatic:Age_20top_20boundary_20layer)=
### __Parameter name:__ Age top boundary layer
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The age of the upper thermal boundary layer, used for the calculation of the half-space cooling model temperature. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Initial_20temperature_20model:Adiabatic:Amplitude)=
### __Parameter name:__ Amplitude
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The amplitude (in K) of the initial spherical temperature perturbation at the bottom of the model domain. This perturbation will be added to the adiabatic temperature profile, but not to the bottom thermal boundary layer. Instead, the maximum of the perturbation and the bottom boundary layer temperature will be used. 

(parameters:Initial_20temperature_20model:Adiabatic:Position)=
### __Parameter name:__ Position
**Default value:** center 

**Pattern:** [Selection center ] 

**Documentation:** Where the initial temperature perturbation should be placed. If `center' is given, then the perturbation will be centered along a `midpoint' of some sort of the bottom boundary. For example, in the case of a box geometry, this is the center of the bottom face; in the case of a spherical shell geometry, it is along the inner surface halfway between the bounding radial lines. 

(parameters:Initial_20temperature_20model:Adiabatic:Radius)=
### __Parameter name:__ Radius
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The Radius (in m) of the initial spherical temperature perturbation at the bottom of the model domain. 

(parameters:Initial_20temperature_20model:Adiabatic:Subadiabaticity)=
### __Parameter name:__ Subadiabaticity
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** If this value is larger than 0, the initial temperature profile will not be adiabatic, but subadiabatic. This value gives the maximal deviation from adiabaticity. Set to 0 for an adiabatic temperature profile. Units: \si{\kelvin}.

The function object in the Function subsection represents the compositional fields that will be used as a reference profile for calculating the thermal diffusivity. This function is one-dimensional and depends only on depth. The format of this functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}. 

(parameters:Initial_20temperature_20model:Adiabatic:Function)=
## **Parameters in section** Initial temperature model/Adiabatic/Function
(parameters:Initial_20temperature_20model:Adiabatic:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Initial_20temperature_20model:Adiabatic:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Initial_20temperature_20model:Adiabatic:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Initial_20temperature_20model:Adiabatic_20boundary)=
## **Parameters in section** Initial temperature model/Adiabatic boundary
(parameters:Initial_20temperature_20model:Adiabatic_20boundary:Adiabatic_20temperature_20gradient)=
### __Parameter name:__ Adiabatic temperature gradient
**Default value:** 0.0005 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the adiabatic temperature gradient. Units: \si{\kelvin\per\meter}. 

(parameters:Initial_20temperature_20model:Adiabatic_20boundary:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Initial_20temperature_20model:Adiabatic_20boundary:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** adiabatic_boundary.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20temperature_20model:Adiabatic_20boundary:Isotherm_20temperature)=
### __Parameter name:__ Isotherm temperature
**Default value:** 1673.15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the isothermal boundary temperature. Units: \si{\kelvin}. 

(parameters:Initial_20temperature_20model:Adiabatic_20boundary:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20temperature_20model:Adiabatic_20boundary:Surface_20temperature)=
### __Parameter name:__ Surface temperature
**Default value:** 273.15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the surface temperature. Units: \si{\kelvin}. 

(parameters:Initial_20temperature_20model:Ascii_20data_20model)=
## **Parameters in section** Initial temperature model/Ascii data model
(parameters:Initial_20temperature_20model:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.  

(parameters:Initial_20temperature_20model:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** initial_isotherm_500K_box_3d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20temperature_20model:Ascii_20data_20model:Data_20file_20names)=
### __Parameter name:__ Data file names
**Default value:** initial_isotherm_500K_box_3d.txt 

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** The file names of the model data (comma separated).  

(parameters:Initial_20temperature_20model:Ascii_20data_20model:Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** linear 

**Pattern:** [Selection piecewise constant|linear ] 

**Documentation:** Method to interpolate between layer boundaries. Select from piecewise constant or linear. Piecewise constant takes the value from the nearest layer boundary above the data point. The linear option interpolates linearly between layer boundaries. Above and below the domain given by the layer boundaries, the values aregiven by the top and bottom layer boundary. 

(parameters:Initial_20temperature_20model:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20temperature_20model:Ascii_20profile)=
## **Parameters in section** Initial temperature model/Ascii profile
(parameters:Initial_20temperature_20model:Ascii_20profile:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/ascii-profile/tests/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Initial_20temperature_20model:Ascii_20profile:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** simple_test.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20temperature_20model:Ascii_20profile:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20temperature_20model:Continental_20geotherm)=
## **Parameters in section** Initial temperature model/Continental geotherm
(parameters:Initial_20temperature_20model:Continental_20geotherm:Layer_20thicknesses)=
### __Parameter name:__ Layer thicknesses
**Default value:** 30000. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of the 3 thicknesses of the lithospheric layers 'upper\_crust', 'lower\_crust' and 'mantle\_lithosphere'. If only one thickness is given, then the same thickness is used for all layers. Units: \si{meter}. 

(parameters:Initial_20temperature_20model:Continental_20geotherm:Lithosphere_2dAsthenosphere_20boundary_20isotherm)=
### __Parameter name:__ Lithosphere_2dAsthenosphere boundary isotherm
**Default value:** 1673.15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the isotherm that is assumed at the Lithosphere-Asthenosphere boundary. Units: \si{\kelvin}. 

(parameters:Initial_20temperature_20model:Continental_20geotherm:Surface_20temperature)=
### __Parameter name:__ Surface temperature
**Default value:** 273.15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the surface temperature. Units: \si{\kelvin}. 

(parameters:Initial_20temperature_20model:Function)=
## **Parameters in section** Initial temperature model/Function
(parameters:Initial_20temperature_20model:Function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Initial_20temperature_20model:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Initial_20temperature_20model:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Initial_20temperature_20model:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Initial_20temperature_20model:Harmonic_20perturbation)=
## **Parameters in section** Initial temperature model/Harmonic perturbation
(parameters:Initial_20temperature_20model:Harmonic_20perturbation:Lateral_20wave_20number_20one)=
### __Parameter name:__ Lateral wave number one
**Default value:** 3 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Doubled first lateral wave number of the harmonic perturbation. Equals the spherical harmonic degree in 3D spherical shells. In all other cases one equals half of a sine period over the model domain. This allows for single up-/downswings. Negative numbers reverse the sign of the perturbation but are not allowed for the spherical harmonic case. 

(parameters:Initial_20temperature_20model:Harmonic_20perturbation:Lateral_20wave_20number_20two)=
### __Parameter name:__ Lateral wave number two
**Default value:** 2 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Doubled second lateral wave number of the harmonic perturbation. Equals the spherical harmonic order in 3D spherical shells. In all other cases one equals half of a sine period over the model domain. This allows for single up-/downswings. Negative numbers reverse the sign of the perturbation. 

(parameters:Initial_20temperature_20model:Harmonic_20perturbation:Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The magnitude of the Harmonic perturbation. 

(parameters:Initial_20temperature_20model:Harmonic_20perturbation:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 1600.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature that is perturbed by the harmonic function. Only used in incompressible models. 

(parameters:Initial_20temperature_20model:Harmonic_20perturbation:Vertical_20wave_20number)=
### __Parameter name:__ Vertical wave number
**Default value:** 1 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** Doubled radial wave number of the harmonic perturbation.  One equals half of a sine period over the model domain.  This allows for single up-/downswings. Negative numbers  reverse the sign of the perturbation. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation)=
## **Parameters in section** Initial temperature model/Inclusion shape perturbation
(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Ambient_20temperature)=
### __Parameter name:__ Ambient temperature
**Default value:** 1.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The background temperature for the temperature field. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Center_20X)=
### __Parameter name:__ Center X
**Default value:** 0.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The X coordinate for the center of the shape. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Center_20Y)=
### __Parameter name:__ Center Y
**Default value:** 0.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The Y coordinate for the center of the shape. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Center_20Z)=
### __Parameter name:__ Center Z
**Default value:** 0.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The Z coordinate for the center of the shape. This is only necessary for three-dimensional fields. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Inclusion_20gradient)=
### __Parameter name:__ Inclusion gradient
**Default value:** constant 

**Pattern:** [Selection gaussian|linear|constant ] 

**Documentation:** The gradient of the inclusion to be generated. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Inclusion_20shape)=
### __Parameter name:__ Inclusion shape
**Default value:** circle 

**Pattern:** [Selection square|circle ] 

**Documentation:** The shape of the inclusion to be generated. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Inclusion_20temperature)=
### __Parameter name:__ Inclusion temperature
**Default value:** 0.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature of the inclusion shape. This is only the true temperature in the case of the constant gradient. In all other cases, it gives one endpoint of the temperature gradient for the shape. 

(parameters:Initial_20temperature_20model:Inclusion_20shape_20perturbation:Shape_20radius)=
### __Parameter name:__ Shape radius
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The radius of the inclusion to be generated. For shapes with no radius (e.g. square), this will be the width, and for shapes with no width, this gives a general guideline for the size of the shape. 

(parameters:Initial_20temperature_20model:Lithosphere_20Mask)=
## **Parameters in section** Initial temperature model/Lithosphere Mask
(parameters:Initial_20temperature_20model:Lithosphere_20Mask:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/lithosphere-mask/ 

**Pattern:** [DirectoryName] 

**Documentation:** The path to the LAB depth data file 

(parameters:Initial_20temperature_20model:Lithosphere_20Mask:Depth_20specification_20method)=
### __Parameter name:__ Depth specification method
**Default value:** Value 

**Pattern:** [Selection File|Value ] 

**Documentation:** Method that is used to specify the depth of the lithosphere-asthenosphere boundary. 

(parameters:Initial_20temperature_20model:Lithosphere_20Mask:LAB_20depth_20filename)=
### __Parameter name:__ LAB depth filename
**Default value:** LAB_CAM2016.txt 

**Pattern:** [FileName (Type: input)] 

**Documentation:** File from which the lithosphere-asthenosphere boundary depth data is read. 

(parameters:Initial_20temperature_20model:Lithosphere_20Mask:Lithosphere_20temperature)=
### __Parameter name:__ Lithosphere temperature
**Default value:** 1600. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The initial temperature within lithosphere, applied abovethe maximum lithosphere depth. 

(parameters:Initial_20temperature_20model:Lithosphere_20Mask:Maximum_20lithosphere_20depth)=
### __Parameter name:__ Maximum lithosphere depth
**Default value:** 200000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Units: \si{\meter}.The maximum depth of the lithosphere. The model will be NaNs below this depth. 

(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS)=
## **Parameters in section** Initial temperature model/Patch on S40RTS
(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Maximum_20grid_20depth)=
### __Parameter name:__ Maximum grid depth
**Default value:** 700000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum depth of the Vs ascii grid. The model will read in  Vs from S40RTS below this depth. 

(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Remove_20temperature_20heterogeneity_20down_20to_20specified_20depth)=
### __Parameter name:__ Remove temperature heterogeneity down to specified depth
**Default value:** -1.7976931348623157e+308 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** This will set the heterogeneity prescribed by the Vs ascii grid and S40RTS to zero down to the specified depth (in meters). Note that your resolution has to be adequate to capture this cutoff. For example if you specify a depth of 660km, but your closest spherical depth layers are only at 500km and 750km (due to a coarse resolution) it will only zero out heterogeneities down to 500km. Similar caution has to be taken when using adaptive meshing. 

(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Smoothing_20length_20scale)=
### __Parameter name:__ Smoothing length scale
**Default value:** 200000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The depth range (above maximum grid depth) over which to smooth. The boundary is smoothed using a depth weighted combination of Vs values from the ascii grid and S40RTS at each point in the region of smoothing. 

(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Ascii_20data_20model)=
## **Parameters in section** Initial temperature model/Patch on S40RTS/Ascii data model
(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/patch-on-S40RTS/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** upper_shell_3d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20temperature_20model:Patch_20on_20S40RTS:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation)=
## **Parameters in section** Initial temperature model/S40RTS perturbation
(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/ 

**Pattern:** [DirectoryName] 

**Documentation:** The path to the model data.  

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Initial_20condition_20file_20name)=
### __Parameter name:__ Initial condition file name
**Default value:** S40RTS.sph 

**Pattern:** [Anything] 

**Documentation:** The file name of the spherical harmonics coefficients from Ritsema et al. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Maximum_20order)=
### __Parameter name:__ Maximum order
**Default value:** 20 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The maximum order the users specify when reading the data file of spherical harmonic coefficients, which must be smaller than the maximum order the data file stored. This parameter will be used only if 'Specify a lower maximum order' is set to true. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 1600.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature that is perturbed by the spherical harmonic functions. Only used in incompressible models. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Remove_20degree_200_20from_20perturbation)=
### __Parameter name:__ Remove degree 0 from perturbation
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Option to remove the degree zero component from the perturbation, which will ensure that the laterally averaged temperature for a fixed depth is equal to the background temperature. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Remove_20temperature_20heterogeneity_20down_20to_20specified_20depth)=
### __Parameter name:__ Remove temperature heterogeneity down to specified depth
**Default value:** -1.7976931348623157e+308 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** This will set the heterogeneity prescribed by S20RTS or S40RTS to zero down to the specified depth (in meters). Note that your resolution has to be adequate to capture this cutoff. For example if you specify a depth of 660km, but your closest spherical depth layers are only at 500km and 750km (due to a coarse resolution) it will only zero out heterogeneities down to 500km. Similar caution has to be taken when using adaptive meshing. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Specify_20a_20lower_20maximum_20order)=
### __Parameter name:__ Specify a lower maximum order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to use a lower maximum order when reading the data file of spherical harmonic coefficients. This is probably used for the faster tests or when the users only want to see the spherical harmonic pattern up to a certain order. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Spline_20knots_20depth_20file_20name)=
### __Parameter name:__ Spline knots depth file name
**Default value:** Spline_knots.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the spline knot locations from Ritsema et al. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Thermal_20expansion_20coefficient_20in_20initial_20temperature_20scaling)=
### __Parameter name:__ Thermal expansion coefficient in initial temperature scaling
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Use_20thermal_20expansion_20coefficient_20from_20material_20model)=
### __Parameter name:__ Use thermal expansion coefficient from material model
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to take the thermal expansion coefficient from the material model instead of from what is specified in this section. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Vs_20to_20density_20scaling)=
### __Parameter name:__ Vs to density scaling
**Default value:** 0.25 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** This parameter specifies how the perturbation in shear wave velocity as prescribed by S20RTS or S40RTS is scaled into a density perturbation. See the general description of this model for more detailed information. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Vs_20to_20density_20scaling_20method)=
### __Parameter name:__ Vs to density scaling method
**Default value:** constant 

**Pattern:** [Selection file|constant ] 

**Documentation:** Method that is used to specify how the vs-to-density scaling varies with depth. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Ascii_20data_20vs_20to_20density_20model)=
## **Parameters in section** Initial temperature model/S40RTS perturbation/Ascii data vs to density model
(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Ascii_20data_20vs_20to_20density_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Ascii_20data_20vs_20to_20density_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** vs_to_density_Steinberger.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20temperature_20model:S40RTS_20perturbation:Ascii_20data_20vs_20to_20density_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation)=
## **Parameters in section** Initial temperature model/SAVANI perturbation
(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/SAVANI/ 

**Pattern:** [DirectoryName] 

**Documentation:** The path to the model data. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Initial_20condition_20file_20name)=
### __Parameter name:__ Initial condition file name
**Default value:** savani.dlnvs.60.m.ab 

**Pattern:** [Anything] 

**Documentation:** The file name of the spherical harmonics coefficients from Auer et al. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Maximum_20order)=
### __Parameter name:__ Maximum order
**Default value:** 20 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The maximum order the users specify when reading the data file of spherical harmonic coefficients, which must be smaller than the maximum order the data file stored. This parameter will be used only if 'Specify a lower maximum order' is set to true. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 1600.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature that is perturbed by the spherical harmonic functions. Only used in incompressible models. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Remove_20degree_200_20from_20perturbation)=
### __Parameter name:__ Remove degree 0 from perturbation
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Option to remove the degree zero component from the perturbation, which will ensure that the laterally averaged temperature for a fixed depth is equal to the background temperature. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Remove_20temperature_20heterogeneity_20down_20to_20specified_20depth)=
### __Parameter name:__ Remove temperature heterogeneity down to specified depth
**Default value:** -1.7976931348623157e+308 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** This will set the heterogeneity prescribed by SAVANI to zero down to the specified depth (in meters). Note that your resolution has to be adequate to capture this cutoff. For example if you specify a depth of 660km, but your closest spherical depth layers are only at 500km and 750km (due to a coarse resolution) it will only zero out heterogeneities down to 500km. Similar caution has to be taken when using adaptive meshing. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Specify_20a_20lower_20maximum_20order)=
### __Parameter name:__ Specify a lower maximum order
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to use a lower maximum order when reading the data file of spherical harmonic coefficients. This is probably used for the faster tests or when the users only want to see the spherical harmonic pattern up to a certain order. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Spline_20knots_20depth_20file_20name)=
### __Parameter name:__ Spline knots depth file name
**Default value:** Spline_knots.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the spline knots taken from the 28 spherical layers of SAVANI tomography model. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Thermal_20expansion_20coefficient_20in_20initial_20temperature_20scaling)=
### __Parameter name:__ Thermal expansion coefficient in initial temperature scaling
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Use_20thermal_20expansion_20coefficient_20from_20material_20model)=
### __Parameter name:__ Use thermal expansion coefficient from material model
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to take the thermal expansion coefficient from the material model instead of from what is specified in this section. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Vs_20to_20density_20scaling)=
### __Parameter name:__ Vs to density scaling
**Default value:** 0.25 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** This parameter specifies how the perturbation in shear wave velocity as prescribed by SAVANI is scaled into a density perturbation. See the general description of this model for more detailed information. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Vs_20to_20density_20scaling_20method)=
### __Parameter name:__ Vs to density scaling method
**Default value:** constant 

**Pattern:** [Selection file|constant ] 

**Documentation:** Method that is used to specify how the vs-to-density scaling varies with depth. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Ascii_20data_20vs_20to_20density_20model)=
## **Parameters in section** Initial temperature model/SAVANI perturbation/Ascii data vs to density model
(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Ascii_20data_20vs_20to_20density_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Ascii_20data_20vs_20to_20density_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** vs_to_density_Steinberger.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Initial_20temperature_20model:SAVANI_20perturbation:Ascii_20data_20vs_20to_20density_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation)=
## **Parameters in section** Initial temperature model/Spherical gaussian perturbation
(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation:Amplitude)=
### __Parameter name:__ Amplitude
**Default value:** 0.01 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The amplitude of the perturbation. 

(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation:Angle)=
### __Parameter name:__ Angle
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The angle where the center of the perturbation is placed. 

(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation:Filename_20for_20initial_20geotherm_20table)=
### __Parameter name:__ Filename for initial geotherm table
**Default value:** initial-geotherm-table 

**Pattern:** [FileName (Type: input)] 

**Documentation:** The file from which the initial geotherm table is to be read. The format of the file is defined by what is read in source/initial\_temperature/spherical\_shell.cc. 

(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation:Non_2ddimensional_20depth)=
### __Parameter name:__ Non_2ddimensional depth
**Default value:** 0.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The non-dimensional radial distance where the center of the perturbation is placed. 

(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation:Sigma)=
### __Parameter name:__ Sigma
**Default value:** 0.2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The standard deviation of the Gaussian perturbation. 

(parameters:Initial_20temperature_20model:Spherical_20gaussian_20perturbation:Sign)=
### __Parameter name:__ Sign
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The sign of the perturbation. 

(parameters:Initial_20temperature_20model:Spherical_20hexagonal_20perturbation)=
## **Parameters in section** Initial temperature model/Spherical hexagonal perturbation
(parameters:Initial_20temperature_20model:Spherical_20hexagonal_20perturbation:Angular_20mode)=
### __Parameter name:__ Angular mode
**Default value:** 6 

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)] 

**Documentation:** The number of convection cells with which to perturb the system. 

(parameters:Initial_20temperature_20model:Spherical_20hexagonal_20perturbation:Rotation_20offset)=
### __Parameter name:__ Rotation offset
**Default value:** -45. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Amount of clockwise rotation in degrees to apply to the perturbations. Default is set to -45 in order to provide backwards compatibility. 

(parameters:Material_20model)=
## **Parameters in section** Material model
(parameters:Material_20model:Material_20averaging)=
### __Parameter name:__ Material averaging
**Default value:** none 

**Pattern:** [Selection none|arithmetic average|harmonic average|geometric average|pick largest|project to Q1|log average|harmonic average only viscosity|geometric average only viscosity|project to Q1 only viscosity ] 

**Documentation:** Whether or not (and in the first case, how) to do any averaging of material model output data when constructing the linear systems for velocity/pressure, temperature, and compositions in each time step, as well as their corresponding preconditioners.

Possible choices: none|arithmetic average|harmonic average|geometric average|pick largest|project to Q1|log average|harmonic average only viscosity|geometric average only viscosity|project to Q1 only viscosity

The process of averaging, and where it may be used, is discussed in more detail in Section~\ref{sec:sinker-with-averaging}.

More averaging schemes are available in the averaging material model. This material model is a ``compositing material model'' which can be used in combination with other material models. 

(parameters:Material_20model:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** The name of the material model to be used in this simulation. There are many material models you can choose from, as listed below. They generally fall into two category: (i) models that implement a particular case of material behavior, (ii) models that modify other models in some way. We sometimes call the latter ``compositing models''. An example of a compositing model is the ``depth dependent'' model below in that it takes another, freely choosable model as its base and then modifies that model's output in some way.

You can select one of the following models:

`Steinberger': This material model looks up the viscosity from the tables that correspond to the paper of Steinberger and Calderwood 2006 (``Models of large-scale viscous flow in the Earth's mantle with constraints from mineral physics and surface observations'', Geophys. J. Int., 167, 1461-1481, \url{http://dx.doi.org/10.1111/j.1365-246X.2006.03131.x}) and material data from a database generated by the thermodynamics code \texttt{Perplex}, see \url{http://www.perplex.ethz.ch/}. The default example data builds upon the thermodynamic database by Stixrude 2011 and assumes a pyrolitic composition by Ringwood 1988 but is easily replaceable by other data files. 

`ascii reference profile': A material model that reads in a reference state for density, thermal expansivity, compressibility and specific heat from a text file. 
Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of points in the reference state as for example `# POINTS: 3'. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named `density', `thermal\_expansivity', `specific\_heat', and `compressibility'. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

The viscosity $\eta$ is computed as \begin{equation}\eta(z,T) = \eta_r(z) \eta_0 \exp\left(-A \frac{T - T_{\text{adi}}}{T_{\text{adi}}}\right),\end{equation}where $\eta_r(z)$ is the depth-dependence, which is a piecewise constant function computed according to the list of ``Viscosity prefactors'' and ``Transition depths'', $\eta_0$ is the reference viscosity specified by the parameter ``Viscosity'' and $A$ describes the dependence on temperature and corresponds to the parameter ``Thermal viscosity exponent''.

`averaging': The ``averaging'' Material model applies an averaging of the quadrature points within a cell. The values to average are supplied by any of the other available material models. In other words, it is a ``compositing material model''. Parameters related to the average model are read from a subsection ``Material model/Averaging''.

The user must specify a ``Base model'' from which material properties are derived. Furthermore an averaging operation must be selected, where the Choice should be from the list none|arithmetic average|harmonic average|geometric average|pick largest|log average|NWD arithmetic average|NWD harmonic average|NWD geometric average.

NWD stands for Normalized Weighed Distance. The models with this in front of their name work with a weighed average, which means each quadrature point requires an individual weight. The weight is determined by the distance, where the exact relation is determined by a bell shaped curve. A bell shaped curve is a continuous function which is one at its maximum and exactly zero at and beyond its limit. This bell shaped curve is spanned around each quadrature point to determine the weighting map for each quadrature point. The used bell shape comes from Lucy (1977). The distance is normalized so the largest distance becomes one. This means that if variable ''Bell shape limit'' is exactly one, the farthest quadrature point is just on the limit and its weight will be exactly zero. In this plugin it is not implemented as larger and equal than the limit, but larger than, to ensure the quadrature point at distance zero is always included.

`compositing': The ``compositing'' Material model selects material model properties from a given set of other material models, and is intended to make mixing different material models easier. This is useful, for example, when wanting to use the melting parameterization of the ``melt simple'' model (which has a relatively simple viscosity model that only allows for a temperature- but not strain rate-dependent viscosity) with a more realistic viscosity model such as that provided by the ``diffusion dislocation'' model.

Specifically, this material model works by allowing to specify the name of another material model for each coefficient that material models are asked for (such as the viscosity, density, etc.). Whenever the material model is asked for the values of coefficients, it then evaluates all of the ``base models'' that were listed for the various coefficients, and copies the values returned by these base models into the output structure.

The implementation of this material model is somewhat expensive because it has to evaluate all material coefficients of all underlying material models. Consequently, if performance of assembly and postprocessing is important, then implementing a separate material model is a better choice than using this material model.

`composition reaction': A material model that behaves in the same way as the simple material model, but includes two compositional fields and a reaction between them. Above a depth given in the input file, the first fields gets converted to the second field. 

`depth dependent': The ``depth dependent'' Material model applies a depth-dependent scaling to any of the other available material models. In other words, it is a ``compositing material model''.

Parameters related to the depth dependent model are read from a subsection ``Material model/Depth dependent model''. The user must specify a ``Base model'' from which material properties are derived. Currently the depth dependent model only allows depth dependence of viscosity - other material properties are taken from the ``Base model''. Viscosity $\eta$ at depth $z$ is calculated according to:\begin{equation}\eta(z,p,T,X,...) = \eta(z) \eta_b(p,T,X,..)/\eta_{r}\end{equation}where $\eta(z)$ is the depth-dependence specified by the depth dependent model, $\eta_b(p,T,X,...)$ is the viscosity calculated from the base model, and $\eta_{r}$ is the reference viscosity. In addition to the specification of the ``Base model'', the user must specify the method to be used to calculate the depth-dependent viscosity $\eta(z)$ as ``Material model/Depth dependent model/Depth dependence method'', which can be chosen among ``None|Function|File|List''. Each method and the associated parameters are as follows:

``Function'': read a user-specified parsed function from the input file in a subsection ``Material model/Depth dependent model/Viscosity depth function''. By default, this function is uniformly equal to 1.0e21. Specifying a function that returns a value less than or equal to 0.0 anywhere in the model domain will produce an error. 

``File'': read a user-specified file containing viscosity values at specified depths. The file containing depth-dependent viscosities is read from a directory specified by the user as ``Material model/Depth dependent model/Data directory'', from a file with name specified as ``Material model/Depth dependent model/Viscosity depth file''. The format of this file is ascii text and contains two columns with one header line:

example Viscosity depth file:\\Depth (m)    Viscosity (Pa-s)\\0.0000000e+00     1.0000000e+21\\6.7000000e+05     1.0000000e+22\\

Viscosity is interpolated from this file using linear interpolation. ``None'': no depth-dependence. Viscosity is taken directly from ``Base model''

``List:'': read a comma-separated list of depth values corresponding to the maximum depths of layers having constant depth-dependence $\eta(z)$. The layers must be specified in order of increasing depth, and the last layer in the list must have a depth greater than or equal to the maximal depth of the model. The list of layer depths is specified as ``Material model/Depth dependent model/Depth list'' and the corresponding list of layer viscosities is specified as ``Material model/Depth dependent model/Viscosity list''

`diffusion dislocation': An implementation of a viscous rheology including diffusion and dislocation creep. Compositional fields can each be assigned individual activation energies, reference densities, thermal expansivities, and stress exponents. The effective viscosity is defined as 

\[\eta_{\text{eff}} = \left(\frac{1}{\eta_{\text{eff}}^\text{diff}}+ \frac{1}{\eta_{\text{eff}}^\text{dis}}\right)^{-1}\] where \[\eta_{\text{i}} = \frac{1}{2} A^{-\frac{1}{n_i}} d^\frac{m_i}{n_i} \dot{\varepsilon_i}^{\frac{1-n_i}{n_i}} \exp\left(\frac{E_i^\ast + PV_i^\ast}{n_iRT}\right)\] 

where $d$ is grain size, $i$ corresponds to diffusion or dislocation creep, $\dot{\varepsilon}$ is the square root of the second invariant of the strain rate tensor, $R$ is the gas constant, $T$ is temperature, and $P$ is pressure. $A_i$ are prefactors, $n_i$ and $m_i$ are stress and grain size exponents $E_i$ are the activation energies and $V_i$ are the activation volumes. 

This form of the viscosity equation is commonly used in geodynamic simulations See, for example, Billen and Hirth (2007), G3, 8, Q08012. Significantly, other studies may use slightly different forms of the viscosity equation leading to variations in how specific terms are defined or combined. For example, the grain size exponent should always be positive in the diffusion viscosity equation used here, while other studies place the grain size term in the denominator and invert the sign of the grain size exponent. When examining previous work, one should carefully check how the viscous prefactor and grain size terms are defined.  

The ratio of diffusion to dislocation strain rate is found by Newton's method, iterating to find the stress which satisfies the above equations. The value for the components of this formula and additional parameters are read from the parameter file in subsection 'Material model/DiffusionDislocation'.

`drucker prager': A material model that has constant values for all coefficients but the density and viscosity. The defaults for all coefficients are chosen to be similar to what is believed to be correct for Earth's mantle. All of the values that define this model are read from a section ``Material model/Drucker Prager'' in the input file, see Section~\ref{parameters:Material_20model/Drucker_20Prager}. Note that the model does not take into account any dependencies of material properties on compositional fields. 

The viscosity is computed according to the Drucker Prager frictional plasticity criterion (non-associative) based on a user-defined internal friction angle $\phi$ and cohesion $C$. In 3D:  $\sigma_y = \frac{6 C \cos(\phi)}{\sqrt{3} (3+\sin(\phi))} + \frac{6 P \sin(\phi)}{\sqrt{3} (3+\sin(\phi))}$, where $P$ is the pressure. See for example Zienkiewicz, O. C., Humpheson, C. and Lewis, R. W. (1975), G\'{e}otechnique 25, No. 4, 671-689. With this formulation we circumscribe instead of inscribe the Mohr Coulomb yield surface. In 2D the Drucker Prager yield surface is the same as the Mohr Coulomb surface:  $\sigma_y = P \sin(\phi) + C \cos(\phi)$. Note that in 2D for $\phi=0$, these criteria revert to the von Mises criterion (no pressure dependence). See for example Thieulot, C. (2011), PEPI 188, 47-68. 

Note that we enforce the pressure to be positive to prevent negative yield strengths and viscosities. 

We then use the computed yield strength to scale back the viscosity on to the yield surface using the Viscosity Rescaling Method described in Kachanov, L. M. (2004), Fundamentals of the Theory of Plasticity, Dover Publications, Inc. (Not Radial Return.)A similar implementation can be found in GALE (https://geodynamics.org/cig/software/gale/gale-manual.pdf). 

To avoid numerically unfavourably large (or even negative) viscosity ranges, we cut off the viscosity with a user-defined minimum and maximum viscosity: $\eta_eff = \frac{1}{\frac{1}{\eta_min + \eta}+ \frac{1}{\eta_max}}$. 

Note that this model uses the formulation that assumes an incompressible medium despite the fact that the density follows the law $\rho(T)=\rho_0(1-\beta(T-T_{\text{ref}}))$. 

`grain size': A material model that relies on compositional fields that correspond to the average grain sizes of a mineral phase and source terms that determine the grain size evolution in terms of the strain rate, temperature, phase transitions, and the creep regime. This material model only works if a compositional field named 'grain_size' is present. In the diffusion creep regime, the viscosity depends on this grain size field. We use the grain size evolution laws described in Behn et al., 2009. Implications of grain size evolution on the seismic structure of the oceanic upper mantle, Earth Planet. Sci. Letters, 282, 178–189. Other material parameters are either prescribed similar to the 'simple' material model, or read from data files that were generated by the Perplex or Hefesto software. This material model is described in more detail in Dannberg, J., Z. Eilon, U. Faul, R. Gassmoeller, P. Moulik, and R. Myhill (2017), The importance of grain size to mantle dynamics and seismological observations, Geochem. Geophys. Geosyst., 18, 3034–3061, doi:10.1002/2017GC006944.

`latent heat': A material model that includes phase transitions and the possibility that latent heat is released or absorbed when material crosses one of the phase transitions of up to two different materials (compositional fields). This model implements a standard approximation of the latent heat terms following Christensen \& Yuen, 1985. The change of entropy is calculated as $Delta S = \gamma \frac{\Delta\rho}{\rho^2}$ with the Clapeyron slope $\gamma$ and the density change $\Delta\rho$ of the phase transition being input parameters. The model employs an analytic phase function in the form $X=\frac{1}{2} \left( 1 + \tanh \left( \frac{\Delta p}{\Delta p_0} \right) \right)$ with $\Delta p = p - p_{\text{transition}} - \gamma \left( T - T_{\text{transition}} \right)$ and $\Delta p_0$ being the pressure difference over the width of the phase transition (specified as input parameter).

`latent heat melt': A material model that includes the latent heat of melting for two materials: peridotite and pyroxenite. The melting model for peridotite is taken from Katz et al., 2003 (A new parameterization of hydrous mantle melting) and the one for pyroxenite from Sobolev et al., 2011 (Linking mantle plumes, large igneous provinces and environmental catastrophes). The model assumes a constant entropy change for melting 100\% of the material, which can be specified in the input file. The partial derivatives of entropy with respect to temperature and pressure required for calculating the latent heat consumption are then calculated as product of this constant entropy change, and the respective derivative of the function the describes the melt fraction. This is linearly averaged with respect to the fractions of the two materials present. If no compositional fields are specified in the input file, the model assumes that the material is peridotite. If compositional fields are specified, the model assumes that the first compositional field is the fraction of pyroxenite and the rest of the material is peridotite. 

Otherwise, this material model has a temperature- and pressure-dependent density and viscosity and the density and thermal expansivity depend on the melt fraction present. It is possible to extent this model to include a melt fraction dependence of all the material parameters by calling the function melt_fraction in the calculation of the respective parameter. However, melt and solid move with the same velocity and melt extraction is not taken into account (batch melting). 

`melt global': A material model that implements a simple formulation of the material parameters required for the modelling of melt transport, including a source term for the porosity according to a simplified linear melting model similar to \cite{schmeling2006}:
$\phi_{\text{equilibrium}} = \frac{T-T_{\text{sol}}}{T_{\text{liq}}-T_{\text{sol}}}$
with $T_{\text{sol}} = T_{\text{sol,0}} + \Delta T_p \, p + \Delta T_c \, C$ 
$T_{\text{liq}} = T_{\text{sol}}  + \Delta T_{\text{sol-liq}}$.

`melt simple': A material model that implements a simple formulation of the material parameters required for the modelling of melt transport, including a source term for the porosity according to the melting model for dry peridotite of \cite{KSL2003}. This also includes a computation of the latent heat of melting (if the `latent heat' heating model is active).

Most of the material properties are constant, except for the shear, viscosity $\eta$, the compaction viscosity $\xi$, and the permeability $k$, which depend on the porosity; and the solid and melt densities, which depend on temperature and pressure:
 $\eta(\phi,T) = \eta_0 e^{\alpha(\phi-\phi_0)} e^{-\beta(T-T_0)/T_0}$, $\xi(\phi,T) = \xi_0 \frac{\phi_0}{\phi} e^{-\beta(T-T_0)/T_0}$, $k=k_0 \phi^n (1-\phi)^m$, $\rho=\rho_0 (1 - \alpha (T - T_{\text{adi}})) e^{\kappa p}$.

The model is compressible only if this is specified in the input file, and contains compressibility for both solid and melt.

`modified tait': A material model that implements the thermal modified Tait equation of state as written in \cite{HP2011}. Constant values are used for the thermal conductivity and viscosity. The defaults for all coefficients are chosen to be similar to what is believed to be correct for Earth's mantle. All of the values that define this model are read from a section ``Material model/Modified Tait model'' in the input file, see Section~\ref{parameters:Material_20model/Modified_20Tait_20model}.

`multicomponent': This incompressible model is for use with an arbitrary number of compositional fields, where each field represents a rock type which can have completely different properties from the others. However, each rock type itself has constant material properties.  The value of the  compositional field is interpreted as a volume fraction. If the sum of the fields is greater than one, they are renormalized.  If it is less than one, material properties  for ``background mantle'' make up the rest. When more than one field is present, the material properties are averaged arithmetically.  An exception is the viscosity, where the averaging should make more of a difference.  For this, the user selects between arithmetic, harmonic, geometric, or maximum composition averaging.

`multicomponent compressible': This model is for use with an arbitrary number of compositional fields, where each field represents a rock type which can have completely different properties from the others. Each rock type is described by a self-consistent equation of state.  The value of the  compositional field is interpreted as a mass fraction. If the sum of the fields is greater than one, they are renormalized.  If it is less than one, material properties  for ``background mantle'' make up the rest. When more than one field is present, the material properties are averaged arithmetically by mass fraction (for specific heat), or volume fraction (for density, thermal expansivity and compressibility). The thermal conductivity is also arithmetically averaged by volume fraction. Finally, the viscosity is averaged by volume fraction, but the user can choose between arithmetic, harmonic, geometric or maximum composition averaging.

`nondimensional': A material model for nondimensionalized computations for compressible or incompressible computations defined through Rayleigh number 	ext{Ra} and Dissipation number Di. This model is made to be used with the Boussinesq, ALA, or TALA formulation.

The viscosity is defined as \[\eta = \text{Di} / \text{Ra} \cdot \exp(-b T' + c z)\] where $T'$ is the temperature variation from the adiabatic temperature, $z$ is the depth, $b$ is given by ``Viscosity temperature prefactor'', and $c$ by ``Viscosity depth prefactor''. If $\text{Di}$ is zero, it will be replaced by 1.0 in $\eta$.

The density is defined as \[\rho = \exp(\text{Di}/\gamma \cdot z)  (1.0 - \alpha T' + \text{Di} \gamma p'),\] where $\alpha=\text{Di}$ is the thermal expansion coefficient, $\gamma$ is the Grueneisen parameter, and $p'$ is the pressure variation from the adiabatic pressure. The pressure dependent term is not present if ``TALA'' is enabled.

`perplex lookup': A material model that has constant values for viscosity and thermal conductivity, and calculates other properties on-the-fly using PerpleX meemum. Compositional fields correspond to the individual components in the order given in the PerpleX file.

`replace lithosphere viscosity': The ``replace lithosphere viscosity'' Material model sets viscosity to a prescribed constant above the lithosphere-asthenosphere boundary (specified by an ascii file or maximum lithosphere depth). Below the lithosphere-asthenosphereboundary the viscosity is taken from any of the other available material model. In other words, it is a ``compositing material model''.
Parameters related to the replace lithosphere viscosity model are read from a subsection ``Material model/Replace lithosphere viscosity''. The user must specify a ``Base model'' from which other material properties are derived.  
Note the required format of the input data file: The first lines may contain any number of comments if they begin with ‘#’, but one of these lines needs to contain the number of grid points in each dimension as for example‘# POINTS: 3 3’. For a spherical model, the order of the data columns has to be'phi', 'theta','depth (m)', where phi is the  azimuth angle and theta is the polar angle measured positive from the north pole.

`simple': A material model that has constant values for all coefficients but the density and viscosity. The defaults for all coefficients are chosen to be similar to what is believed to be correct for Earth's mantle. All of the values that define this model are read from a section ``Material model/Simple model'' in the input file, see Section~\ref{parameters:Material_20model/Simple_20model}.

This model uses the following set of equations for the two coefficients that are non-constant: \begin{align}  \eta(p,T,\mathfrak c) &= \tau(T) \zeta(\mathfrak c) \eta_0, \\  \rho(p,T,\mathfrak c) &= \left(1-\alpha (T-T_0)\right)\rho_0 + \Delta\rho \; c_0,\end{align}where $c_0$ is the first component of the compositional vector $\mathfrak c$ if the model uses compositional fields, or zero otherwise. 

The temperature pre-factor for the viscosity formula above is defined as \begin{align}  \tau(T) &= H\left(e^{-\beta (T-T_0)/T_0}\right),\intertext{with}   \qquad\qquad H(x) &= \begin{cases}                            \tau_{\text{min}} & \text{if}\; x<\tau_{\text{min}}, \\                            x & \text{if}\; 10^{-2}\le x \le 10^2, \\                            \tau_{\text{max}} & \text{if}\; x>\tau_{\text{max}}, \\                         \end{cases}\end{align} where $x=e^{-\beta (T-T_0)/T_0}$, $\beta$ corresponds to the input parameter ``Thermal viscosity exponent'', and $T_0$ to the parameter ``Reference temperature''. If you set $T_0=0$ in the input file, the thermal pre-factor $\tau(T)=1$. The parameters $\tau_{\text{min}}$ and $\tau_{\text{max}}$ set the minimum and maximum values of the temperature pre-factor and are set using ``Maximum thermal prefactor'' and ``Minimum thermal prefactor''. Specifying a value of 0.0 for the minimum or maximum values will disable pre-factor limiting.

The compositional pre-factor for the viscosity is defined as \begin{align}  \zeta(\mathfrak c) &= \xi^{c_0}\end{align} if the model has compositional fields and equals one otherwise. $\xi$ corresponds to the parameter ``Composition viscosity prefactor'' in the input file.

Finally, in the formula for the density, $\alpha$ corresponds to the ``Thermal expansion coefficient'' and $\Delta\rho$ corresponds to the parameter ``Density differential for compositional field 1''.

Note that this model uses the formulation that assumes an incompressible medium despite the fact that the density follows the law $\rho(T)=\rho_0(1-\alpha(T-T_{\text{ref}}))$. 

\note{Despite its name, this material model is not exactly ``simple'', as indicated by the formulas above. While it was originally intended to be simple, it has over time acquired all sorts of temperature and compositional dependencies that weren't initially intended. Consequently, there is now a ``simpler'' material model that now fills the role the current model was originally intended to fill.}

`simple compressible': A material model that has constant values for all coefficients but the density. The defaults for all coefficients are chosen to be similar to what is believed to be correct for Earth's mantle. All of the values that define this model are read from a section ``Material model/Simple compressible model'' in the input file, see Section~\ref{parameters:Material_20model/Simple_20compressible_20model}.

This model uses the following equations for the density: \begin{align}  \rho(p,T) = \rho_0              \left(1-\alpha (T-T_a)\right)               \exp{\beta (P-P_0))}\end{align}This formulation for the density assumes that the compressibility provided by the user is the adiabatic compressibility ($\beta_S$). The thermal expansivity and isentropic compressibility implied by the pressure and temperature dependence are equal to the user-defined constant values only along the reference isentrope, and there is also an implicit pressure dependence to the heat capacity $C_p$ via Maxwell's relations.

`simpler': A material model that has constant values except for density, which depends linearly on temperature: \begin{align}  \rho(p,T) &= \left(1-\alpha (T-T_0)\right)\rho_0.\end{align}

\note{This material model fills the role the ``simple'' material model was originally intended to fill, before the latter acquired all sorts of complicated temperature and compositional dependencies.}

`visco plastic': An implementation of an incompressible visco(elastic)-plastic rheology with options for selecting dislocation creep, diffusion creep or composite viscous flow laws. Prior to yielding, one may select to modify the viscosity to account for viscoelastic effects. Plasticity limits viscous stresses through a Drucker Prager yield criterion. The implementation of this material model is based heavily on the `DiffusionDislocation' (Bob Myhill), `DruckerPrager' (Anne Glerum), and `Viscoelastic' (John Naliboff) material models. 

 The viscosity for dislocation or diffusion creep is defined as \[v = \frac 12 A^{-\frac{1}{n}} d^{\frac{m}{n}} \dot{\varepsilon}_{ii}^{\frac{1-n}{n}} \exp\left(\frac{E + PV}{nRT}\right)\] where $A$ is the prefactor, $n$ is the stress exponent, $\dot{\varepsilon}_{ii}$ is the square root of the deviatoric strain rate tensor second invariant, $d$ is grain size, $m$ is the grain size exponent, $E$ is activation energy, $V$ is activation volume, $P$ is pressure, $R$ is the gas exponent and $T$ is temperature. This form of the viscosity equation is commonly used in geodynamic simulations. See, for example, Billen and Hirth (2007), G3, 8, Q08012. Significantly, other studies may use slightly different forms of the viscosity equation leading to variations in how specific terms are defined or combined. For example, the grain size exponent should always be positive in the diffusion viscosity equation used here, while other studies place the grain size term in the denominator and invert the sign of the grain size exponent. When examining previous work, one should carefully check how the viscous prefactor and grain size terms are defined. 

 One may select to use the diffusion ($v_{\text{diff}}$; $n=1$, $m
eq 0$), dislocation ($v_{\text{disl}}$, $n>1$, $m=0$) or composite $\frac{v_{\text{diff}} v_{\text{disl}}}{v_{\text{diff}}+v_{\text{disl}}}$ equation form. 

 The diffusion and dislocation prefactors can be weakened with a factor between 0 and 1 according to the total or the viscous strain only. 

 Viscosity is limited through one of two different `yielding' mechanisms. 

The first plasticity mechanism limits viscous stress through a Drucker Prager yield criterion, where the yield stress in 3D is  $\sigma_y = \frac{6C\cos(\phi) + 2P\sin(\phi)} {\sqrt{3}(3+\sin(\phi))}$ and $\sigma_y = C\cos(\phi) + P\sin(\phi)$ in 2D. Above, $C$ is cohesion and $\phi$  is the angle of internal friction.  Note that the 2D form is equivalent to the Mohr Coulomb yield surface.  If $\phi$ is 0, the yield stress is fixed and equal to the cohesion (Von Mises yield criterion). When the viscous stress ($2v{\varepsilon}_{ii}$) exceeds the yield stress, the viscosity is rescaled back to the yield surface: $v_{y}=\sigma_{y}/(2{\varepsilon}_{ii})$. This form of plasticity is commonly used in geodynamic models. See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. 

The user has the option to linearly reduce the cohesion and internal friction angle as a function of the finite strain magnitude. The finite strain invariant or full strain tensor is calculated through compositional fields within the material model. This implementation is identical to the compositional field finite strain plugin and cookbook described in the manual (author: Gassmoeller, Dannberg). If the user selects to track the finite strain invariant ($e_{ii}$), a single compositional field tracks the value derived from $e_{ii}^t = (e_{ii})^{(t-1)} + \dot{e}_{ii}\; dt$, where $t$ and $t-1$ are the current and prior time steps, $\dot{e}_{ii}$ is the second invariant of the strain rate tensor and $dt$ is the time step size. In the case of the full strain tensor $F$, the finite strain magnitude is derived from the second invariant of the symmetric stretching tensor $L$, where $L = F [F]^T$. The user must specify a single compositional field for the finite strain invariant or multiple fields (4 in 2D, 9 in 3D) for the finite strain tensor. These field(s) must be the first listed compositional fields in the parameter file. Note that one or more of the finite strain tensor components must be assigned a non-zero value initially. This value can be be quite small (e.g., 1.e-8), but still non-zero. While the option to track and use the full finite strain tensor exists, tracking the associated compositional fields is computationally expensive in 3D. Similarly, the finite strain magnitudes may in fact decrease if the orientation of the deformation field switches through time. Consequently, the ideal solution is track the finite strain invariant (single compositional) field within the material and track the full finite strain tensor through particles.When only the second invariant of the strain is tracked, one has the option to track the full strain or only the plastic strain. In the latter case, strain is only tracked in case the material is plastically yielding, i.e. the viscous stress > yield stress. 

Viscous stress may also be limited by a non-linear stress limiter that has a form similar to the Peierls creep mechanism. This stress limiter assigns an effective viscosity $\sigma_{\text{eff}} = \frac{\tau_y}{2\varepsilon_y} {\frac{\varepsilon_{ii}}{\varepsilon_y}}^{\frac{1}{n_y}-1}$ Above $\tau_y$ is a yield stress, $\varepsilon_y$ is the reference strain rate, $\varepsilon_{ii}$ is the strain rate and $n_y$ is the stress limiter exponent.  The yield stress, $\tau_y$, is defined through the Drucker Prager yield criterion formulation. This method of limiting viscous stress has been used in various forms within the geodynamic literature \cite{chri92,vavv02,cibi13,cibi15}.When $n_y$ is 1, it essentially becomes a linear viscosity model, and in the limit $n_y\rightarrow \infty$ it converges to the standard viscosity rescaling method (concretely, values $n_y>20$ are large enough).

 The visco-plastic rheology described above may also be modified to include viscoelastic deformation, thus producing a viscoelastic plastic constitutive relationship. 

 The viscoelastic rheology behavior takes into account the elastic shear strength (e.g., shear modulus), while the tensile and volumetric strength (e.g., Young's and bulk modulus) are not considered. The model is incompressible and allows specifying an arbitrary number of compositional fields, where each field represents a different rock type or component of the viscoelastic stress tensor. The stress tensor in 2D and 3D, respectively, contains 3 or 6 components. The compositional fields representing these components must be named and listed in a very specific format, which is designed to minimize mislabeling stress tensor components as distinct 'compositional rock types' (or vice versa). For 2D models, the first three compositional fields must be labeled 'stress\_xx', 'stress\_yy' and 'stress\_xy'. In 3D, the first six compositional fields must be labeled 'stress\_xx', 'stress\_yy', 'stress\_zz', 'stress\_xy', 'stress\_xz', 'stress\_yz'. 

 Combining this viscoelasticity implementation with non-linear viscous flow and plasticity produces a constitutive relationship commonly referred to as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. While extensively discussed and applied within the geodynamics literature, notable references include: Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. Gerya (2010), Introduction to Numerical Geodynamic Modeling. Kaus (2010), Tectonophysics, v. 484, p. 36-47. Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. 

 The overview below directly follows Moresi et al. (2003) eqns. 23-38. However, an important distinction between this material model and the studies above is the use of compositional fields, rather than particles, to track individual components of the viscoelastic stress tensor. The material model will be updated when an option to track and calculate viscoelastic stresses with particles is implemented. 

 Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric rate of deformation ($\hat{D}$) as the sum of elastic ($\hat{D_{e}}$) and viscous ($\hat{D_{v}}$) components: $\hat{D} = \hat{D_{e}} + \hat{D_{v}}$.  These terms further decompose into $\hat{D_{v}} = \frac{\tau}{2\eta}$ and $\hat{D_{e}} = \frac{\overset{\nabla}{\tau}}{2\mu}$, where $\tau$ is the viscous deviatoric stress, $\eta$ is the shear viscosity, $\mu$ is the shear modulus and $\overset{\nabla}{\tau}$ is the Jaumann corotational stress rate. This later term (eqn. 24) contains the time derivative of the deviatoric stress ($\dot{\tau}$) and terms that account for material spin (e.g., rotation) due to advection: $\overset{\nabla}{\tau} = \dot{\tau} + {\tau}W -W\tau$. Above, $W$ is the material spin tensor (eqn. 25): $W_{ij} = \frac{1}{2} \left (\frac{\partial V_{i}}{\partial x_{j}} - \frac{\partial V_{j}}{\partial x_{i}} \right )$. 

 If plasticity is included, the deviatoric rate of deformation may be written as: $\hat{D} = \hat{D_{e}} + \hat{D_{v}} + \hat{D_{p}}$, where $\hat{D_{p}}$ is the plastic component. $\hat{D_{p}}$ decomposes to $\frac{\tau_{y}}{2\eta_{y}}$, where $\tau_{y}$ is the yield stress and $\eta_{y}$ is the viscosity rescaled to the yield surface. The Jaumann stress-rate can also be approximated using terms from the previous time step ($t$) and current time step ($t + \Delta t^{e}$): $\smash[t]{\overset{\nabla}{\tau}}^{t + \Delta t^{e}} \approx \frac{\tau^{t + \Delta t^{e} - \tau^{t}}}{\Delta t^{e}} - W^{t}\tau^{t} + \tau^{t}W^{t}$. In this material model, the size of the time step above ($\Delta t^{e}$) can be specified as the numerical time step size or an independent fixed time step. If the latter case is selected, the user has an option to apply a stress averaging scheme to account for the differences between the numerical and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time step throughout the model run, this can still be achieved by using CFL and maximum time step values that restrict the numerical time step to a specific time.

 The formulation above allows rewriting the total rate of deformation (eqn. 29) as
 $\tau^{t + \Delta t^{e}} = \eta_{eff} \left ( 2\hat{D}^{t + \triangle t^{e}} + \frac{\tau^{t}}{\mu \Delta t^{e}} + \frac{W^{t}\tau^{t} - \tau^{t}W^{t}}{\mu}  \right )$. 

 The effective viscosity (eqn. 28) is a function of the viscosity ($\eta$), elastic time step size ($\Delta t^{e}$) and shear relaxation time ($ \alpha = \frac{\eta}{\mu} $): $\eta_{eff} = \eta \frac{\Delta t^{e}}{\Delta t^{e} + \alpha}$ The magnitude of the shear modulus thus controls how much the effective viscosity is reduced relative to the initial viscosity. 

 Elastic effects are introduced into the governing Stokes equations through an elastic force term (eqn. 30) using stresses from the previous time step: $F^{e,t} = -\frac{\eta_{eff}}{\mu \Delta t^{e}} \tau^{t}$. This force term is added onto the right-hand side force vector in the system of equations. 

 When plastic yielding occurs, the effective viscosity in equation 29 and 30 is the plastic viscosity (equation 36). If the current stress is below the plastic yield stress, the effective viscosity is still as defined in equation 28. During non-linear iterations, we define the current stress prior to yielding (e.g., value compared to yield stress) as $\tau^{t + \Delta t^{e}} = \eta_{eff} \left ( 2\hat{D}^{t + \triangle t^{e}} + \frac{\tau^{t}}{\mu \Delta t^{e}} \right ) $

 Compositional fields can each be assigned individual values of thermal diffusivity, heat capacity, density, thermal expansivity and rheological parameters. 

 If more than one compositional field is present at a given point, viscosities are averaged with an arithmetic, geometric harmonic (default) or maximum composition scheme. 

 The value for the components of this formula and additional parameters are read from the parameter file in subsection  'Material model/Visco Plastic'.

`viscoelastic': An implementation of a simple linear viscoelastic rheology that only includes the deviatoric components of elasticity. Specifically, the viscoelastic rheology only takes into account the elastic shear strength (e.g., shear modulus), while the tensile and volumetric strength (e.g., Young's and bulk modulus) are not considered. The model is incompressible and allows specifying an arbitrary number of compositional fields, where each field represents a different rock type or component of the viscoelastic stress tensor. The stress tensor in 2D and 3D, respectively, contains 3 or 6 components. The compositional fields representing these components must be named and listed in a very specific format, which is designed to minimize mislabeling stress tensor components as distinct 'compositional rock types' (or vice versa). For 2D models, the first three compositional fields must be labeled 'stress\_xx', 'stress\_yy' and 'stress\_xy'. In 3D, the first six compositional fields must be labeled 'stress\_xx', 'stress\_yy', 'stress\_zz', 'stress\_xy', 'stress\_xz', 'stress\_yz'. 

 Expanding the model to include non-linear viscous flow (e.g., diffusion/dislocation creep) and plasticity would produce a constitutive relationship commonly referred to as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. While extensively discussed and applied within the geodynamics literature, notable references include: Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. Gerya (2010), Introduction to Numerical Geodynamic Modeling. Kaus (2010), Tectonophysics, v. 484, p. 36-47. Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. 

 The overview below directly follows Moresi et al. (2003) eqns. 23-32. However, an important distinction between this material model and the studies above is the use of compositional fields, rather than particles, to track individual components of the viscoelastic stress tensor. The material model will be updated when an option to track and calculate viscoelastic stresses with particles is implemented. 

 Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric rate of deformation ($\hat{D}$) as the sum of elastic ($\hat{D_{e}}$) and viscous ($\hat{D_{v}}$) components: $\hat{D} = \hat{D_{e}} + \hat{D_{v}}$.  These terms further decompose into $\hat{D_{v}} = \frac{\tau}{2\eta}$ and $\hat{D_{e}} = \frac{\overset{\nabla}{\tau}}{2\mu}$, where $\tau$ is the viscous deviatoric stress, $\eta$ is the shear viscosity, $\mu$ is the shear modulus and $\overset{\nabla}{\tau}$ is the Jaumann corotational stress rate. This later term (eqn. 24) contains the time derivative of the deviatoric stress ($\dot{\tau}$) and terms that account for material spin (e.g., rotation) due to advection: $\overset{\nabla}{\tau} = \dot{\tau} + {\tau}W -W\tau$. Above, $W$ is the material spin tensor (eqn. 25): $W_{ij} = \frac{1}{2} \left (\frac{\partial V_{i}}{\partial x_{j}} - \frac{\partial V_{j}}{\partial x_{i}} \right )$. 

 The Jaumann stress-rate can also be approximated using terms from the previous time step ($t$) and current time step ($t + \Delta t^{e}$): $\smash[t]{\overset{\nabla}{\tau}}^{t + \Delta t^{e}} \approx \frac{\tau^{t + \Delta t^{e} - \tau^{t}}}{\Delta t^{e}} - W^{t}\tau^{t} + \tau^{t}W^{t}$. In this material model, the size of the time step above ($\Delta t^{e}$) can be specified as the numerical time step size or an independent fixed time step. If the latter case is selected, the user has an option to apply a stress averaging scheme to account for the differences between the numerical and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time step throughout the model run, this can still be achieved by using CFL and maximum time step values that restrict the numerical time step to a specific time.

 The formulation above allows rewriting the total deviatoric stress (eqn. 29) as
 $\tau^{t + \Delta t^{e}} = \eta_\text{eff} \left ( 2\hat{D}^{t + \triangle t^{e}} + \frac{\tau^{t}}{\mu \Delta t^{e}} + \frac{W^{t}\tau^{t} - \tau^{t}W^{t}}{\mu}  \right )$. 

 The effective viscosity (eqn. 28) is a function of the viscosity ($\eta$), elastic time step size ($\Delta t^{e}$) and shear relaxation time ($ \alpha = \frac{\eta}{\mu} $): $\eta_\text{eff} = \eta \frac{\Delta t^{e}}{\Delta t^{e} + \alpha}$ The magnitude of the shear modulus thus controls how much the effective viscosity is reduced relative to the initial viscosity. 

 Elastic effects are introduced into the governing Stokes equations through an elastic force term (eqn. 30) using stresses from the previous time step: $F^{e,t} = -\frac{\eta_\text{eff}}{\mu \Delta t^{e}} \tau^{t}$. This force term is added onto the right-hand side force vector in the system of equations. 

 The value of each compositional field representing distinct rock types at a point is interpreted to be a volume fraction of that rock type. If the sum of the compositional field volume fractions is less than one, then the remainder of the volume is assumed to be 'background material'.

 Several model parameters (densities, elastic shear moduli, thermal expansivities, thermal conductivies, specific heats) can be defined per-compositional field. For each material parameter the user supplies a comma delimited list of length N+1, where N is the number of compositional fields. The additional field corresponds to the value for background material. They should be ordered ''background, composition1, composition2...''. However, the first 3 (2D) or 6 (3D) composition fields correspond to components of the elastic stress tensor and their material values will not contribute to the volume fractions. If a single value is given, then all the compositional fields are given that value. Other lengths of lists are not allowed. For a given compositional field the material parameters are treated as constant, except density, which varies linearly with temperature according to the thermal expansivity. 

 When more than one compositional field is present at a point, they are averaged arithmetically. An exception is viscosity, which may be averaged arithmetically, harmonically, geometrically, or by selecting the viscosity of the composition field with the greatest volume fraction. 

(parameters:Material_20model:Ascii_20reference_20profile)=
## **Parameters in section** Material model/Ascii reference profile
(parameters:Material_20model:Ascii_20reference_20profile:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference conductivity 

(parameters:Material_20model:Ascii_20reference_20profile:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of viscosity. Dimensionless exponent. 

(parameters:Material_20model:Ascii_20reference_20profile:Transition_20depths)=
### __Parameter name:__ Transition depths
**Default value:** 1.5e5, 4.1e5, 6.6e5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of depths where the viscosity changes. Values must monotonically increase. Units: \si{\meter}. 

(parameters:Material_20model:Ascii_20reference_20profile:Use_20TALA)=
### __Parameter name:__ Use TALA
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use the TALA instead of the ALA approximation. 

(parameters:Material_20model:Ascii_20reference_20profile:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 1e21 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Viscosity 

(parameters:Material_20model:Ascii_20reference_20profile:Viscosity_20prefactors)=
### __Parameter name:__ Viscosity prefactors
**Default value:** 10., 0.1, 1., 10. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of prefactors for the viscosity that determine the viscosity profile. Each prefactor is applied in a depth range specified by the list of `Transition depths', i.e. the first prefactor is applied above the first transition depth, the second one between the first and second transition depth, and so on. To compute the viscosity profile, this prefactor is multiplied by the reference viscosity specified through the parameter `Viscosity'. List must have one more entry than Transition depths. Units: non-dimensional. 

(parameters:Material_20model:Ascii_20reference_20profile:Ascii_20data_20model)=
## **Parameters in section** Material model/Ascii reference profile/Ascii data model
(parameters:Material_20model:Ascii_20reference_20profile:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/adiabatic-conditions/ascii-data/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Material_20model:Ascii_20reference_20profile:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Material_20model:Ascii_20reference_20profile:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Material_20model:Averaging)=
## **Parameters in section** Material model/Averaging
(parameters:Material_20model:Averaging:Averaging_20operation)=
### __Parameter name:__ Averaging operation
**Default value:** none 

**Pattern:** [Selection none|arithmetic average|harmonic average|geometric average|pick largest|log average|nwd arithmetic average|nwd harmonic average|nwd geometric average ] 

**Documentation:** Choose the averaging operation to use. 

(parameters:Material_20model:Averaging:Base_20model)=
### __Parameter name:__ Base model
**Default value:** simple 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic ] 

**Documentation:** The name of a material model that will be modified by an averaging operation. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Averaging:Bell_20shape_20limit)=
### __Parameter name:__ Bell shape limit
**Default value:** 1. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The limit normalized distance between 0 and 1 where the bell shape becomes zero. See the manual for a more information. 

(parameters:Material_20model:Compositing)=
## **Parameters in section** Material model/Compositing
(parameters:Material_20model:Compositing:Compressibility)=
### __Parameter name:__ Compressibility
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Compressibility. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Density)=
### __Parameter name:__ Density
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Density. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Entropy_20derivative_20pressure)=
### __Parameter name:__ Entropy derivative pressure
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Entropy derivative pressure. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Entropy_20derivative_20temperature)=
### __Parameter name:__ Entropy derivative temperature
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Entropy derivative temperature. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Reaction_20terms)=
### __Parameter name:__ Reaction terms
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Reaction terms. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Specific_20heat)=
### __Parameter name:__ Specific heat
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Specific heat. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Thermal conductivity. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Thermal expansion coefficient. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Compositing:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** unspecified 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic|unspecified ] 

**Documentation:** Material model to use for Viscosity. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Composition_20reaction_20model)=
## **Parameters in section** Material model/Composition reaction model
(parameters:Material_20model:Composition_20reaction_20model:Composition_20viscosity_20prefactor_201)=
### __Parameter name:__ Composition viscosity prefactor 1
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A linear dependency of viscosity on the first compositional field. Dimensionless prefactor. With a value of 1.0 (the default) the viscosity does not depend on the composition. 

(parameters:Material_20model:Composition_20reaction_20model:Composition_20viscosity_20prefactor_202)=
### __Parameter name:__ Composition viscosity prefactor 2
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A linear dependency of viscosity on the second compositional field. Dimensionless prefactor. With a value of 1.0 (the default) the viscosity does not depend on the composition. 

(parameters:Material_20model:Composition_20reaction_20model:Density_20differential_20for_20compositional_20field_201)=
### __Parameter name:__ Density differential for compositional field 1
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** If compositional fields are used, then one would frequently want to make the density depend on these fields. In this simple material model, we make the following assumptions: if no compositional fields are used in the current simulation, then the density is simply the usual one with its linear dependence on the temperature. If there are compositional fields, then the material model determines how many of them influence the density. The composition-dependence adds a term of the kind $+\Delta \rho \; c_1(\mathbf x)$. This parameter describes the value of $\Delta \rho$. Units: \si{\kilogram\per\meter\cubed}/unit change in composition. 

(parameters:Material_20model:Composition_20reaction_20model:Density_20differential_20for_20compositional_20field_202)=
### __Parameter name:__ Density differential for compositional field 2
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** If compositional fields are used, then one would frequently want to make the density depend on these fields. In this simple material model, we make the following assumptions: if no compositional fields are used in the current simulation, then the density is simply the usual one with its linear dependence on the temperature. If there are compositional fields, then the material model determines how many of them influence the density. The composition-dependence adds a term of the kind $+\Delta \rho \; c_2(\mathbf x)$. This parameter describes the value of $\Delta \rho$. Units: \si{\kilogram\per\meter\cubed}/unit change in composition. 

(parameters:Material_20model:Composition_20reaction_20model:Reaction_20depth)=
### __Parameter name:__ Reaction depth
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Above this depth the compositional fields react: The first field gets converted to the second field. Units: \si{\meter}. 

(parameters:Material_20model:Composition_20reaction_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Composition_20reaction_20model:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Composition_20reaction_20model:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Composition_20reaction_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Composition_20reaction_20model:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Composition_20reaction_20model:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of viscosity. Dimensionless exponent. 

(parameters:Material_20model:Composition_20reaction_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity. Units: \si{\kilogram\per\meter\per\second}. 

(parameters:Material_20model:Depth_20dependent_20model)=
## **Parameters in section** Material model/Depth dependent model
(parameters:Material_20model:Depth_20dependent_20model:Base_20model)=
### __Parameter name:__ Base model
**Default value:** simple 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic ] 

**Documentation:** The name of a material model that will be modified by a depth dependent viscosity. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for that for more information. 

(parameters:Material_20model:Depth_20dependent_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/material-model/rheology/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Material_20model:Depth_20dependent_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** ascii_depth_profile.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Material_20model:Depth_20dependent_20model:Depth_20dependence_20method)=
### __Parameter name:__ Depth dependence method
**Default value:** None 

**Pattern:** [Selection Function|File|List|None ] 

**Documentation:** Method that is used to specify how the viscosity should vary with depth. 

(parameters:Material_20model:Depth_20dependent_20model:Depth_20list)=
### __Parameter name:__ Depth list
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma-separated list of depth values for use with the ``List'' ``Depth dependence method''. The list must be provided in order of increasing depth, and the last value must be greater than or equal to the maximal depth of the model. The depth list is interpreted as a layered viscosity structure and the depth values specify the maximum depths of each layer. 

(parameters:Material_20model:Depth_20dependent_20model:Reference_20viscosity)=
### __Parameter name:__ Reference viscosity
**Default value:** 1.7976931348623157e+308 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant reference viscosity $\eta_r$ that is used to scale the non-dimenional depth-dependent viscosity prefactor. Units: \si{\pascal\second}. 

(parameters:Material_20model:Depth_20dependent_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Material_20model:Depth_20dependent_20model:Viscosity_20depth_20file)=
### __Parameter name__: Viscosity depth file
**Alias:** [Data file name](parameters:Material_20model:Depth_20dependent_20model:Data_20file_20name)

**Deprecation Status:** false

(parameters:Material_20model:Depth_20dependent_20model:Viscosity_20list)=
### __Parameter name:__ Viscosity list
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma-separated list of viscosity values, corresponding to the depth values provided in ``Depth list''. The number of viscosity values specified here must be the same as the number of depths provided in ``Depth list''. 

(parameters:Material_20model:Depth_20dependent_20model:Viscosity_20depth_20function)=
## **Parameters in section** Material model/Depth dependent model/Viscosity depth function
(parameters:Material_20model:Depth_20dependent_20model:Viscosity_20depth_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Material_20model:Depth_20dependent_20model:Viscosity_20depth_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 1.0e21 

**Pattern:** [Anything] 

**Documentation:**  

(parameters:Material_20model:Depth_20dependent_20model:Viscosity_20depth_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Material_20model:Diffusion_20dislocation)=
## **Parameters in section** Material model/Diffusion dislocation
(parameters:Material_20model:Diffusion_20dislocation:Activation_20energies_20for_20diffusion_20creep)=
### __Parameter name:__ Activation energies for diffusion creep
**Default value:** 375e3 

**Pattern:** [Anything] 

**Documentation:** List of activation energies, $E_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Diffusion_20dislocation:Activation_20energies_20for_20dislocation_20creep)=
### __Parameter name:__ Activation energies for dislocation creep
**Default value:** 530e3 

**Pattern:** [Anything] 

**Documentation:** List of activation energies, $E_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Diffusion_20dislocation:Activation_20volumes_20for_20diffusion_20creep)=
### __Parameter name:__ Activation volumes for diffusion creep
**Default value:** 6e-6 

**Pattern:** [Anything] 

**Documentation:** List of activation volumes, $V_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Diffusion_20dislocation:Activation_20volumes_20for_20dislocation_20creep)=
### __Parameter name:__ Activation volumes for dislocation creep
**Default value:** 1.4e-5 

**Pattern:** [Anything] 

**Documentation:** List of activation volumes, $V_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Diffusion_20dislocation:Densities)=
### __Parameter name:__ Densities
**Default value:** 3300. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of densities, $\rho$, for background mantle and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Diffusion_20dislocation:Effective_20viscosity_20coefficient)=
### __Parameter name:__ Effective viscosity coefficient
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Scaling coefficient for effective viscosity. 

(parameters:Material_20model:Diffusion_20dislocation:Grain_20size)=
### __Parameter name:__ Grain size
**Default value:** 1e-3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Units: \si{\meter}. 

(parameters:Material_20model:Diffusion_20dislocation:Grain_20size_20exponents_20for_20diffusion_20creep)=
### __Parameter name:__ Grain size exponents for diffusion creep
**Default value:** 3. 

**Pattern:** [Anything] 

**Documentation:** List of grain size exponents, $m_{\text{diffusion}}$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: None. 

(parameters:Material_20model:Diffusion_20dislocation:Heat_20capacity)=
### __Parameter name:__ Heat capacity
**Default value:** 1.25e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Diffusion_20dislocation:Maximum_20strain_20rate_20ratio_20iterations)=
### __Parameter name:__ Maximum strain rate ratio iterations
**Default value:** 40 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Maximum number of iterations to find the correct diffusion/dislocation strain rate ratio. 

(parameters:Material_20model:Diffusion_20dislocation:Maximum_20viscosity)=
### __Parameter name:__ Maximum viscosity
**Default value:** 1e28 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Upper cutoff for effective viscosity. Units: \si{\pascal\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Minimum_20strain_20rate)=
### __Parameter name:__ Minimum strain rate
**Default value:** 1.4e-20 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Stabilizes strain dependent viscosity. Units: \si{\per\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Minimum_20viscosity)=
### __Parameter name:__ Minimum viscosity
**Default value:** 1e17 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Lower cutoff for effective viscosity. Units: \si{\pascal\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Prefactors_20for_20diffusion_20creep)=
### __Parameter name:__ Prefactors for diffusion creep
**Default value:** 1.5e-15 

**Pattern:** [Anything] 

**Documentation:** List of viscosity prefactors, $A$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\per\pascal\meter}$^{m_{\text{diffusion}}}$\si{\per\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Prefactors_20for_20dislocation_20creep)=
### __Parameter name:__ Prefactors for dislocation creep
**Default value:** 1.1e-16 

**Pattern:** [Anything] 

**Documentation:** List of viscosity prefactors, $A$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\pascal}$^{-n_{\text{dislocation}}}$ \si{\per\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** For calculating density by thermal expansivity. Units: \si{\kelvin}. 

(parameters:Material_20model:Diffusion_20dislocation:Reference_20viscosity)=
### __Parameter name:__ Reference viscosity
**Default value:** 1e22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference viscosity that is used for pressure scaling. To understand how pressure scaling works, take a look at \cite{KHB12}. In particular, the value of this parameter would not affect the solution computed by \aspect{} if we could do arithmetic exactly; however, computers do arithmetic in finite precision, and consequently we need to scale quantities in ways so that their magnitudes are roughly the same. As explained in \cite{KHB12}, we scale the pressure during some computations (never visible by users) by a factor that involves a reference viscosity. This parameter describes this reference viscosity.

For problems with a constant viscosity, you will generally want to choose the reference viscosity equal to the actual viscosity. For problems with a variable viscosity, the reference viscosity should be a value that adequately represents the order of magnitude of the viscosities that appear, such as an average value or the value one would use to compute a Rayleigh number.

Units: \si{\pascal\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Strain_20rate_20residual_20tolerance)=
### __Parameter name:__ Strain rate residual tolerance
**Default value:** 1e-22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Tolerance for correct diffusion/dislocation strain rate ratio. 

(parameters:Material_20model:Diffusion_20dislocation:Stress_20exponents_20for_20diffusion_20creep)=
### __Parameter name:__ Stress exponents for diffusion creep
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of stress exponents, $n_{\text{diffusion}}$, for background mantle and compositional fields, for a total of N+1 values, where N is the number of compositional fields. The stress exponent for diffusion creep is almost always equal to one. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Diffusion_20dislocation:Stress_20exponents_20for_20dislocation_20creep)=
### __Parameter name:__ Stress exponents for dislocation creep
**Default value:** 3.5 

**Pattern:** [Anything] 

**Documentation:** List of stress exponents, $n_{\text{dislocation}}$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Diffusion_20dislocation:Thermal_20diffusivity)=
### __Parameter name:__ Thermal diffusivity
**Default value:** 0.8e-6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Units: \si{\meter\squared\per\second}. 

(parameters:Material_20model:Diffusion_20dislocation:Thermal_20expansivities)=
### __Parameter name:__ Thermal expansivities
**Default value:** 3.5e-5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of thermal expansivities for background mantle and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: \si{\per\kelvin}. 

(parameters:Material_20model:Diffusion_20dislocation:Viscosity_20averaging_20scheme)=
### __Parameter name:__ Viscosity averaging scheme
**Default value:** harmonic 

**Pattern:** [Selection arithmetic|harmonic|geometric|maximum composition ] 

**Documentation:** When more than one compositional field is present at a point with different viscosities, we need to come up with an average viscosity at that point.  Select a weighted harmonic, arithmetic, geometric, or maximum composition. 

(parameters:Material_20model:Drucker_20Prager)=
## **Parameters in section** Material model/Drucker Prager
(parameters:Material_20model:Drucker_20Prager:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Drucker_20Prager:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Drucker_20Prager:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. The reference temperature is used in the density calculation. Units: \si{\kelvin}. 

(parameters:Material_20model:Drucker_20Prager:Reference_20viscosity)=
### __Parameter name:__ Reference viscosity
**Default value:** 1e22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference viscosity that is used for pressure scaling. To understand how pressure scaling works, take a look at \cite{KHB12}. In particular, the value of this parameter would not affect the solution computed by \aspect{} if we could do arithmetic exactly; however, computers do arithmetic in finite precision, and consequently we need to scale quantities in ways so that their magnitudes are roughly the same. As explained in \cite{KHB12}, we scale the pressure during some computations (never visible by users) by a factor that involves a reference viscosity. This parameter describes this reference viscosity.

For problems with a constant viscosity, you will generally want to choose the reference viscosity equal to the actual viscosity. For problems with a variable viscosity, the reference viscosity should be a value that adequately represents the order of magnitude of the viscosities that appear, such as an average value or the value one would use to compute a Rayleigh number.

Units: \si{\pascal\second}. 

(parameters:Material_20model:Drucker_20Prager:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Drucker_20Prager:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Drucker_20Prager:Viscosity)=
## **Parameters in section** Material model/Drucker Prager/Viscosity
(parameters:Material_20model:Drucker_20Prager:Viscosity:Angle_20of_20internal_20friction)=
### __Parameter name:__ Angle of internal friction
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the angle of internal friction $\phi$. For a value of zero, in 2D the von Mises criterion is retrieved. Angles higher than 30 degrees are harder to solve numerically. Units: degrees. 

(parameters:Material_20model:Drucker_20Prager:Viscosity:Cohesion)=
### __Parameter name:__ Cohesion
**Default value:** 2e7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the cohesion $C$. Units: \si{\pascal}. 

(parameters:Material_20model:Drucker_20Prager:Viscosity:Maximum_20viscosity)=
### __Parameter name:__ Maximum viscosity
**Default value:** 1e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the maximum viscosity cutoff $\eta_max$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Drucker_20Prager:Viscosity:Minimum_20viscosity)=
### __Parameter name:__ Minimum viscosity
**Default value:** 1e19 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the minimum viscosity cutoff $\eta_min$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Drucker_20Prager:Viscosity:Reference_20strain_20rate)=
### __Parameter name:__ Reference strain rate
**Default value:** 1e-15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the initial strain rate prescribed during the first nonlinear iteration $\dot{\epsilon}_ref$. Units: \si{\per\second}. 

(parameters:Material_20model:Grain_20size_20model)=
## **Parameters in section** Material model/Grain size model
(parameters:Material_20model:Grain_20size_20model:Advect_20logarithm_20of_20grain_20size)=
### __Parameter name:__ Advect logarithm of grain size
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** This parameter determines whether to advect the logarithm of the grain size or the grain size itself. The equation and the physics are the same, but for problems with high grain size gradients it might be preferable to advect the logarithm.  

(parameters:Material_20model:Grain_20size_20model:Average_20specific_20grain_20boundary_20energy)=
### __Parameter name:__ Average specific grain boundary energy
**Default value:** 1.0 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The average specific grain boundary energy $\gamma$. Units: \si{\joule\per\meter\squared}. 

(parameters:Material_20model:Grain_20size_20model:Bilinear_20interpolation)=
### __Parameter name:__ Bilinear interpolation
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** This parameter determines whether to use bilinear interpolation to compute material properties (slower but more accurate). 

(parameters:Material_20model:Grain_20size_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/material-model/steinberger/ 

**Pattern:** [DirectoryName] 

**Documentation:** The path to the model data. The path may also include the special text '$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the 'data/' subdirectory of ASPECT.  

(parameters:Material_20model:Grain_20size_20model:Derivatives_20file_20names)=
### __Parameter name:__ Derivatives file names
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** The file names of the enthalpy derivatives data. List with as many components as active compositional fields (material data is assumed to be in order with the ordering of the fields).  

(parameters:Material_20model:Grain_20size_20model:Diffusion_20activation_20energy)=
### __Parameter name:__ Diffusion activation energy
**Default value:** 3.35e5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The activation energy for diffusion creep $E_{diff}$. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Grain_20size_20model:Diffusion_20activation_20volume)=
### __Parameter name:__ Diffusion activation volume
**Default value:** 4e-6 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The activation volume for diffusion creep $V_{diff}$. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Grain_20size_20model:Diffusion_20creep_20exponent)=
### __Parameter name:__ Diffusion creep exponent
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The power-law exponent $n_{diff}$ for diffusion creep. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Diffusion_20creep_20grain_20size_20exponent)=
### __Parameter name:__ Diffusion creep grain size exponent
**Default value:** 3. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The diffusion creep grain size exponent $p_{diff}$ that determines the dependence of viscosity on grain size. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Diffusion_20creep_20prefactor)=
### __Parameter name:__ Diffusion creep prefactor
**Default value:** 7.4e-15 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The prefactor for the diffusion creep law $A_{diff}$. Units: \si{\meter}$^{p_{diff}}$\si{\pascal}$^{-n_{diff}}$\si{\per\second}. 

(parameters:Material_20model:Grain_20size_20model:Dislocation_20activation_20energy)=
### __Parameter name:__ Dislocation activation energy
**Default value:** 4.8e5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The activation energy for dislocation creep $E_{dis}$. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Grain_20size_20model:Dislocation_20activation_20volume)=
### __Parameter name:__ Dislocation activation volume
**Default value:** 1.1e-5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The activation volume for dislocation creep $V_{dis}$. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Grain_20size_20model:Dislocation_20creep_20exponent)=
### __Parameter name:__ Dislocation creep exponent
**Default value:** 3.5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The power-law exponent $n_{dis}$ for dislocation creep. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Dislocation_20creep_20prefactor)=
### __Parameter name:__ Dislocation creep prefactor
**Default value:** 4.5e-15 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The prefactor for the dislocation creep law $A_{dis}$. Units: \si{\pascal}$^{-n_{dis}}$\si{\per\second}. 

(parameters:Material_20model:Grain_20size_20model:Dislocation_20viscosity_20iteration_20number)=
### __Parameter name:__ Dislocation viscosity iteration number
**Default value:** 100 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** We need to perform an iteration inside the computation of the dislocation viscosity, because it depends on the dislocation strain rate, which depends on the dislocation viscosity itself. This number determines the maximum number of iterations that are performed.  

(parameters:Material_20model:Grain_20size_20model:Dislocation_20viscosity_20iteration_20threshold)=
### __Parameter name:__ Dislocation viscosity iteration threshold
**Default value:** 1e-3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** We need to perform an iteration inside the computation of the dislocation viscosity, because it depends on the dislocation strain rate, which depends on the dislocation viscosity itself. This number determines the termination accuracy, i.e. if the dislocation viscosity changes by less than this factor we terminate the iteration. 

(parameters:Material_20model:Grain_20size_20model:Geometric_20constant)=
### __Parameter name:__ Geometric constant
**Default value:** 3. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The geometric constant $c$ used in the paleowattmeter grain size reduction law. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Grain_20growth_20activation_20energy)=
### __Parameter name:__ Grain growth activation energy
**Default value:** 3.5e5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The activation energy for grain growth $E_g$. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Grain_20size_20model:Grain_20growth_20activation_20volume)=
### __Parameter name:__ Grain growth activation volume
**Default value:** 8e-6 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The activation volume for grain growth $V_g$. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Grain_20size_20model:Grain_20growth_20exponent)=
### __Parameter name:__ Grain growth exponent
**Default value:** 3. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The exponent of the grain growth law $p_g$. This is an experimentally determined grain growth constant. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Grain_20growth_20rate_20constant)=
### __Parameter name:__ Grain growth rate constant
**Default value:** 1.5e-5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The prefactor for the Ostwald ripening grain growth law $G_0$. This is dependent on water content, which is assumed to be 50 H/$10^6$ Si for the default value. Units: \si{\meter}$^{p_g}$\si{\per\second}. 

(parameters:Material_20model:Grain_20size_20model:Lower_20mantle_20grain_20size_20scaling)=
### __Parameter name:__ Lower mantle grain size scaling
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A scaling factor for the grain size in the lower mantle. In models where the high grain size contrast between the upper and lower mantle causes numerical problems, the grain size in the lower mantle can be scaled to a larger value, simultaneously scaling the viscosity prefactors and grain growth parameters to keep the same physical behavior. Differences to the original formulation only occur when material with a smaller grain size than the recrystallization grain size cross the upper-lower mantle boundary. The real grain size can be obtained by dividing the model grain size by this value. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Material_20file_20format)=
### __Parameter name:__ Material file format
**Default value:** perplex 

**Pattern:** [Selection perplex|hefesto ] 

**Documentation:** The material file format to be read in the property tables. 

(parameters:Material_20model:Grain_20size_20model:Material_20file_20names)=
### __Parameter name:__ Material file names
**Default value:** pyr-ringwood88.txt 

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** The file names of the material data. List with as many components as active compositional fields (material data is assumed to be in order with the ordering of the fields).  

(parameters:Material_20model:Grain_20size_20model:Maximum_20latent_20heat_20substeps)=
### __Parameter name:__ Maximum latent heat substeps
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The maximum number of substeps over the temperature pressure range to calculate the averaged enthalpy gradient over a cell. 

(parameters:Material_20model:Grain_20size_20model:Maximum_20specific_20heat)=
### __Parameter name:__ Maximum specific heat
**Default value:** 6000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum specific heat that is allowed in the whole model domain. Units: J/kg/K. 

(parameters:Material_20model:Grain_20size_20model:Maximum_20temperature_20dependence_20of_20viscosity)=
### __Parameter name:__ Maximum temperature dependence of viscosity
**Default value:** 100. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The factor by which viscosity at adiabatic temperature and ambient temperature are allowed to differ (a value of x means that the viscosity can be x times higher or x times lower compared to the value at adiabatic temperature. This parameter is introduced to limit local viscosity contrasts, but still allow for a widely varying viscosity over the whole mantle range. Units: none. 

(parameters:Material_20model:Grain_20size_20model:Maximum_20thermal_20expansivity)=
### __Parameter name:__ Maximum thermal expansivity
**Default value:** 1e-3 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum thermal expansivity that is allowed in the whole model domain. Units: 1/K. 

(parameters:Material_20model:Grain_20size_20model:Maximum_20viscosity)=
### __Parameter name:__ Maximum viscosity
**Default value:** 1e26 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum viscosity that is allowed in the whole model domain. Units: Pa \, s. 

(parameters:Material_20model:Grain_20size_20model:Minimum_20grain_20size)=
### __Parameter name:__ Minimum grain size
**Default value:** 1e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum grain size that is used for the material model. This parameter is introduced to limit local viscosity contrasts, but still allows for a widely varying viscosity over the whole mantle range. Units: \si{\meter}. 

(parameters:Material_20model:Grain_20size_20model:Minimum_20specific_20heat)=
### __Parameter name:__ Minimum specific heat
**Default value:** 500. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum specific heat that is allowed in the whole model domain. Units: J/kg/K. 

(parameters:Material_20model:Grain_20size_20model:Minimum_20thermal_20expansivity)=
### __Parameter name:__ Minimum thermal expansivity
**Default value:** 1e-5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum thermal expansivity that is allowed in the whole model domain. Units: 1/K. 

(parameters:Material_20model:Grain_20size_20model:Minimum_20viscosity)=
### __Parameter name:__ Minimum viscosity
**Default value:** 1e18 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum viscosity that is allowed in the whole model domain. Units: Pa \, s. 

(parameters:Material_20model:Grain_20size_20model:Phase_20transition_20Clapeyron_20slopes)=
### __Parameter name:__ Phase transition Clapeyron slopes
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of Clapeyron slopes for each phase transition. A positive Clapeyron slope indicates that the phase transition will occur in a greater depth, if the temperature is higher than the one given in Phase transition temperatures and in a smaller depth, if the temperature is smaller than the one given in Phase transition temperatures. For negative slopes the other way round. List must have the same number of entries as Phase transition depths. Units: \si{\pascal\per\kelvin}. 

(parameters:Material_20model:Grain_20size_20model:Phase_20transition_20depths)=
### __Parameter name:__ Phase transition depths
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of depths where phase transitions occur. Values must monotonically increase. Units: \si{\meter}. 

(parameters:Material_20model:Grain_20size_20model:Phase_20transition_20temperatures)=
### __Parameter name:__ Phase transition temperatures
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of temperatures where phase transitions occur. Higher or lower temperatures lead to phase transition occurring in smaller or greater depths than given in Phase transition depths, depending on the Clapeyron slope given in Phase transition Clapeyron slopes. List must have the same number of entries as Phase transition depths. Units: \si{\kelvin}. 

(parameters:Material_20model:Grain_20size_20model:Phase_20transition_20widths)=
### __Parameter name:__ Phase transition widths
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of widths for each phase transition. This is only use to specify the region where the recrystallized grain size is assigned after material has crossed a phase transition and should accordingly be chosen similar to the maximum cell width expected at the phase transition.List must have the same number of entries as Phase transition depths. Units: \si{\meter}. 

(parameters:Material_20model:Grain_20size_20model:Reciprocal_20required_20strain)=
### __Parameter name:__ Reciprocal required strain
**Default value:** 10. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** This parameter ($\lambda$) gives an estimate of the strain necessary to achieve a new grain size.  

(parameters:Material_20model:Grain_20size_20model:Recrystallized_20grain_20size)=
### __Parameter name:__ Recrystallized grain size
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. When set to zero, grain size will not be reduced. Units: \si{\meter}. 

(parameters:Material_20model:Grain_20size_20model:Reference_20compressibility)=
### __Parameter name:__ Reference compressibility
**Default value:** 4e-12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the reference compressibility. Units: \si{\per\pascal}. 

(parameters:Material_20model:Grain_20size_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Grain_20size_20model:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $cp$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Grain_20size_20model:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Grain_20size_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Grain_20size_20model:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Grain_20size_20model:Use_20enthalpy_20for_20material_20properties)=
### __Parameter name:__ Use enthalpy for material properties
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** This parameter determines whether to use the enthalpy to calculate the thermal expansivity and specific heat (if true) or use the thermal expansivity and specific heat values from the material properties table directly (if false). 

(parameters:Material_20model:Grain_20size_20model:Use_20paleowattmeter)=
### __Parameter name:__ Use paleowattmeter
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** A flag indicating whether the computation should be use the paleowattmeter approach of Austin and Evans (2007) for grain size reduction in the dislocation creep regime (if true) or the paleopiezometer approach from Hall and Parmetier (2003) (if false). 

(parameters:Material_20model:Grain_20size_20model:Use_20table_20properties)=
### __Parameter name:__ Use table properties
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** This parameter determines whether to use the table properties also for density, thermal expansivity and specific heat. If false the properties are generated as in the simple compressible plugin. 

(parameters:Material_20model:Grain_20size_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity. Units: \si{\pascal\second}. 

(parameters:Material_20model:Grain_20size_20model:Work_20fraction_20for_20boundary_20area_20change)=
### __Parameter name:__ Work fraction for boundary area change
**Default value:** 0.1 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The fraction $\chi$ of work done by dislocation creep to change the grain boundary area. Units: \si{\joule\per\meter\squared}. 

(parameters:Material_20model:Latent_20heat)=
## **Parameters in section** Material model/Latent heat
(parameters:Material_20model:Latent_20heat:Composition_20viscosity_20prefactor)=
### __Parameter name:__ Composition viscosity prefactor
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A linear dependency of viscosity on composition. Dimensionless prefactor. 

(parameters:Material_20model:Latent_20heat:Compressibility)=
### __Parameter name:__ Compressibility
**Default value:** 5.124e-12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the compressibility $\kappa$. Units: \si{\per\pascal}. 

(parameters:Material_20model:Latent_20heat:Corresponding_20phase_20for_20density_20jump)=
### __Parameter name:__ Corresponding phase for density jump
**Default value:**  

**Pattern:** [List of <[Integer range 0...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of phases, which correspond to the Phase transition density jumps. The density jumps occur only in the phase that is given by this phase value. 0 stands for the 1st compositional fields, 1 for the second compositional field and -1 for none of them. List must have the same number of entries as Phase transition depths. Units: \si{\pascal\per\kelvin}. 

(parameters:Material_20model:Latent_20heat:Define_20transition_20by_20depth_20instead_20of_20pressure)=
### __Parameter name:__ Define transition by depth instead of pressure
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to list phase transitions by depth or pressure. If this parameter is true, then the input file will use Phase transitions depths and Phase transition widths to define the phase transition. If it is false, the parameter file will read in phase transition data from Phase transition pressures and Phase transition pressure widths. 

(parameters:Material_20model:Latent_20heat:Density_20differential_20for_20compositional_20field_201)=
### __Parameter name:__ Density differential for compositional field 1
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** If compositional fields are used, then one would frequently want to make the density depend on these fields. In this simple material model, we make the following assumptions: if no compositional fields are used in the current simulation, then the density is simply the usual one with its linear dependence on the temperature. If there are compositional fields, then the density only depends on the first one in such a way that the density has an additional term of the kind $+\Delta \rho \; c_1(\mathbf x)$. This parameter describes the value of $\Delta \rho$. Units: \si{\kilogram\per\meter\cubed}/unit change in composition. 

(parameters:Material_20model:Latent_20heat:Maximum_20viscosity)=
### __Parameter name:__ Maximum viscosity
**Default value:** 1e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Limit for the maximum viscosity in the model. Units: Pa \, s. 

(parameters:Material_20model:Latent_20heat:Minimum_20viscosity)=
### __Parameter name:__ Minimum viscosity
**Default value:** 1e19 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Limit for the minimum viscosity in the model. Units: Pa \, s. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20Clapeyron_20slopes)=
### __Parameter name:__ Phase transition Clapeyron slopes
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of Clapeyron slopes for each phase transition. A positive Clapeyron slope indicates that the phase transition will occur in a greater depth, if the temperature is higher than the one given in Phase transition temperatures and in a smaller depth, if the temperature is smaller than the one given in Phase transition temperatures. For negative slopes the other way round. List must have the same number of entries as Phase transition depths. Units: \si{\pascal\per\kelvin}. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20density_20jumps)=
### __Parameter name:__ Phase transition density jumps
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of density jumps at each phase transition. A positive value means that the density increases with depth. The corresponding entry in Corresponding phase for density jump determines if the density jump occurs in peridotite, eclogite or none of them.List must have the same number of entries as Phase transition depths. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20depths)=
### __Parameter name:__ Phase transition depths
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of depths where phase transitions occur. Values must monotonically increase. Units: \si{\meter}. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20pressure_20widths)=
### __Parameter name:__ Phase transition pressure widths
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of widths for each phase transition, in terms of pressure. The phase functions are scaled with these values, leading to a jump between phases for a value of zero and a gradual transition for larger values. List must have the same number of entries as Phase transition pressures. Define transition by depth instead of pressure must be set to false to use this parameter. Units: \si{\pascal}. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20pressures)=
### __Parameter name:__ Phase transition pressures
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of pressures where phase transitions occur. Values must monotonically increase. Define transition by depth instead of pressure must be set to false to use this parameter. Units: \si{\pascal}. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20temperatures)=
### __Parameter name:__ Phase transition temperatures
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of temperatures where phase transitions occur. Higher or lower temperatures lead to phase transition occurring in smaller or greater depths than given in Phase transition depths, depending on the Clapeyron slope given in Phase transition Clapeyron slopes. List must have the same number of entries as Phase transition depths. Units: \si{\kelvin}. 

(parameters:Material_20model:Latent_20heat:Phase_20transition_20widths)=
### __Parameter name:__ Phase transition widths
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of widths for each phase transition, in terms of depth. The phase functions are scaled with these values, leading to a jump between phases for a value of zero and a gradual transition for larger values. List must have the same number of entries as Phase transition depths. Units: \si{\meter}. 

(parameters:Material_20model:Latent_20heat:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Latent_20heat:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Latent_20heat:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Latent_20heat:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 2.38 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Latent_20heat:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 4e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Latent_20heat:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of viscosity. Dimensionless exponent. 

(parameters:Material_20model:Latent_20heat:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity. Units: \si{\pascal\second}. 

(parameters:Material_20model:Latent_20heat:Viscosity_20prefactors)=
### __Parameter name:__ Viscosity prefactors
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of prefactors for the viscosity for each phase. The reference viscosity will be multiplied by this factor to get the corresponding viscosity for each phase. List must have one more entry than Phase transition depths. Units: non-dimensional. 

(parameters:Material_20model:Latent_20heat_20melt)=
## **Parameters in section** Material model/Latent heat melt
(parameters:Material_20model:Latent_20heat_20melt:A1)=
### __Parameter name:__ A1
**Default value:** 1085.7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Latent_20heat_20melt:A2)=
### __Parameter name:__ A2
**Default value:** 1.329e-7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Latent_20heat_20melt:A3)=
### __Parameter name:__ A3
**Default value:** -5.1e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Latent_20heat_20melt:B1)=
### __Parameter name:__ B1
**Default value:** 1475.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Latent_20heat_20melt:B2)=
### __Parameter name:__ B2
**Default value:** 8.0e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Latent_20heat_20melt:B3)=
### __Parameter name:__ B3
**Default value:** -3.2e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Latent_20heat_20melt:C1)=
### __Parameter name:__ C1
**Default value:** 1780.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Latent_20heat_20melt:C2)=
### __Parameter name:__ C2
**Default value:** 4.50e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Latent_20heat_20melt:C3)=
### __Parameter name:__ C3
**Default value:** -2.0e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Latent_20heat_20melt:Composition_20viscosity_20prefactor)=
### __Parameter name:__ Composition viscosity prefactor
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A linear dependency of viscosity on composition. Dimensionless prefactor. 

(parameters:Material_20model:Latent_20heat_20melt:Compressibility)=
### __Parameter name:__ Compressibility
**Default value:** 5.124e-12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the compressibility $\kappa$. Units: \si{\per\pascal}. 

(parameters:Material_20model:Latent_20heat_20melt:D1)=
### __Parameter name:__ D1
**Default value:** 976.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of pyroxenite. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Latent_20heat_20melt:D2)=
### __Parameter name:__ D2
**Default value:** 1.329e-7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of pyroxenite. Note that this factor is different from the value given in Sobolev, 2011, because they use the potential temperature whereas we use the absolute temperature. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Latent_20heat_20melt:D3)=
### __Parameter name:__ D3
**Default value:** -5.1e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of pyroxenite. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Latent_20heat_20melt:Density_20differential_20for_20compositional_20field_201)=
### __Parameter name:__ Density differential for compositional field 1
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** If compositional fields are used, then one would frequently want to make the density depend on these fields. In this simple material model, we make the following assumptions: if no compositional fields are used in the current simulation, then the density is simply the usual one with its linear dependence on the temperature. If there are compositional fields, then the density only depends on the first one in such a way that the density has an additional term of the kind $+\Delta \rho \; c_1(\mathbf x)$. This parameter describes the value of $\Delta \rho$. Units: \si{\kilogram\per\meter\cubed}/unit change in composition. 

(parameters:Material_20model:Latent_20heat_20melt:E1)=
### __Parameter name:__ E1
**Default value:** 663.8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear depletion term in the quadratic function that approximates the melt fraction of pyroxenite. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Latent_20heat_20melt:E2)=
### __Parameter name:__ E2
**Default value:** -611.4 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic depletion term in the quadratic function that approximates the melt fraction of pyroxenite. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Latent_20heat_20melt:Mass_20fraction_20cpx)=
### __Parameter name:__ Mass fraction cpx
**Default value:** 0.15 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Mass fraction of clinopyroxene in the peridotite to be molten. Units: non-dimensional. 

(parameters:Material_20model:Latent_20heat_20melt:Maximum_20pyroxenite_20melt_20fraction)=
### __Parameter name:__ Maximum pyroxenite melt fraction
**Default value:** 0.5429 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximum melt fraction of pyroxenite in this parameterization. At higher temperatures peridotite begins to melt. 

(parameters:Material_20model:Latent_20heat_20melt:Peridotite_20melting_20entropy_20change)=
### __Parameter name:__ Peridotite melting entropy change
**Default value:** -300. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The entropy change for the phase transition from solid to melt of peridotite. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Latent_20heat_20melt:Pyroxenite_20melting_20entropy_20change)=
### __Parameter name:__ Pyroxenite melting entropy change
**Default value:** -400. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The entropy change for the phase transition from solid to melt of pyroxenite. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Latent_20heat_20melt:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Latent_20heat_20melt:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Latent_20heat_20melt:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Latent_20heat_20melt:Relative_20density_20of_20melt)=
### __Parameter name:__ Relative density of melt
**Default value:** 0.9 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The relative density of melt compared to the solid material. This means, the density change upon melting is this parameter times the density of solid material.Units: non-dimensional. 

(parameters:Material_20model:Latent_20heat_20melt:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 2.38 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Latent_20heat_20melt:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 4e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha_s$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Latent_20heat_20melt:Thermal_20expansion_20coefficient_20of_20melt)=
### __Parameter name:__ Thermal expansion coefficient of melt
**Default value:** 6.8e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha_f$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Latent_20heat_20melt:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of viscosity. Dimensionless exponent. 

(parameters:Material_20model:Latent_20heat_20melt:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity. Units: \si{\pascal\second}. 

(parameters:Material_20model:Latent_20heat_20melt:beta)=
### __Parameter name:__ beta
**Default value:** 1.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Exponent of the melting temperature in the melt fraction calculation. Units: non-dimensional. 

(parameters:Material_20model:Latent_20heat_20melt:r1)=
### __Parameter name:__ r1
**Default value:** 0.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant in the linear function that approximates the clinopyroxene reaction coefficient. Units: non-dimensional. 

(parameters:Material_20model:Latent_20heat_20melt:r2)=
### __Parameter name:__ r2
**Default value:** 8e-11 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the linear function that approximates the clinopyroxene reaction coefficient. Units: \si{\per\pascal}. 

(parameters:Material_20model:Melt_20global)=
## **Parameters in section** Material model/Melt global
(parameters:Material_20model:Melt_20global:Depletion_20density_20change)=
### __Parameter name:__ Depletion density change
**Default value:** 0.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The density contrast between material with a depletion of 1 and a depletion of zero. Negative values indicate lower densities of depleted material. Depletion is indicated by the compositional field with the name peridotite. Not used if this field does not exist in the model. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Melt_20global:Depletion_20solidus_20change)=
### __Parameter name:__ Depletion solidus change
**Default value:** 200.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The solidus temperature change for a depletion of 100\%. For positive values, the solidus gets increased for a positive peridotite field (depletion) and lowered for a negative peridotite field (enrichment). Scaling with depletion is linear. Only active when fractional melting is used. Units: \si{\kelvin}. 

(parameters:Material_20model:Melt_20global:Exponential_20depletion_20strengthening_20factor)=
### __Parameter name:__ Exponential depletion strengthening factor
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** $\alpha_F$: exponential dependency of viscosity on the depletion field $F$ (called peridotite). Dimensionless factor. With a value of 0.0 (the default) the viscosity does not depend on the depletion. The effective viscosity increasedue to depletion is defined as $exp( \alpha_F * F)$. Rationale: melting dehydrates the source rock by removing most of the volatiles,and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a factor 100 to 1000 viscosity contrast between wet and dry rocks, although some experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013). 

(parameters:Material_20model:Melt_20global:Exponential_20melt_20weakening_20factor)=
### __Parameter name:__ Exponential melt weakening factor
**Default value:** 27. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The porosity dependence of the viscosity. Units: dimensionless. 

(parameters:Material_20model:Melt_20global:Include_20melting_20and_20freezing)=
### __Parameter name:__ Include melting and freezing
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to include melting and freezing (according to a simplified linear melting approximation in the model (if true), or not (if false). 

(parameters:Material_20model:Melt_20global:Maximum_20Depletion_20viscosity_20change)=
### __Parameter name:__ Maximum Depletion viscosity change
**Default value:** 1.0e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** $\Delta \eta_{F,max}$: maximum depletion strengthening of viscosity. Rationale: melting dehydrates the source rock by removing most of the volatiles,and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a factor 100 to 1000 viscosity contrast between wet and dry rocks, although some experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013). 

(parameters:Material_20model:Melt_20global:Melt_20bulk_20modulus_20derivative)=
### __Parameter name:__ Melt bulk modulus derivative
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the pressure derivative of the melt bulk modulus. Units: None. 

(parameters:Material_20model:Melt_20global:Melt_20compressibility)=
### __Parameter name:__ Melt compressibility
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the compressibility of the melt. Units: \si{\per\pascal}. 

(parameters:Material_20model:Melt_20global:Melting_20time_20scale_20for_20operator_20splitting)=
### __Parameter name:__ Melting time scale for operator splitting
**Default value:** 1e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** In case the operator splitting scheme is used, the porosity field can not be set to a new equilibrium melt fraction instantly, but the model has to provide a melting time scale instead. This time scale defines how fast melting happens, or more specifically, the parameter defines the time after which the deviation of the porosity from the equilibrium melt fraction will be reduced to a fraction of $1/e$. So if the melting time scale is small compared to the time step size, the reaction will be so fast that the porosity is very close to the equilibrium melt fraction after reactions are computed. Conversely, if the melting time scale is large compared to the time step size, almost no melting and freezing will occur.

Also note that the melting time scale has to be larger than or equal to the reaction time step used in the operator splitting scheme, otherwise reactions can not be computed. If the model does not use operator splitting, this parameter is not used. Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter. 

(parameters:Material_20model:Melt_20global:Pressure_20solidus_20change)=
### __Parameter name:__ Pressure solidus change
**Default value:** 6e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The linear solidus temperature change with pressure. For positive values, the solidus gets increased for positive pressures. Units: \si{\per\pascal}. 

(parameters:Material_20model:Melt_20global:Reference_20bulk_20viscosity)=
### __Parameter name:__ Reference bulk viscosity
**Default value:** 1e22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant bulk viscosity $\xi_0$ of the solid matrix. This viscosity may be modified by both temperature and porosity dependencies. Units: \si{\pascal\second}. 

(parameters:Material_20model:Melt_20global:Reference_20melt_20density)=
### __Parameter name:__ Reference melt density
**Default value:** 2500. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density of the melt/fluid$\rho_{f,0}$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Melt_20global:Reference_20melt_20viscosity)=
### __Parameter name:__ Reference melt viscosity
**Default value:** 10. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant melt viscosity $\eta_f$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Melt_20global:Reference_20permeability)=
### __Parameter name:__ Reference permeability
**Default value:** 1e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference permeability of the solid host rock.Units: \si{\meter\squared}. 

(parameters:Material_20model:Melt_20global:Reference_20shear_20viscosity)=
### __Parameter name:__ Reference shear viscosity
**Default value:** 5e20 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity $\eta_0$ of the solid matrix. This viscosity may be modified by both temperature and porosity dependencies. Units: \si{\pascal\second}. 

(parameters:Material_20model:Melt_20global:Reference_20solid_20density)=
### __Parameter name:__ Reference solid density
**Default value:** 3000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density of the solid $\rho_{s,0}$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Melt_20global:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Melt_20global:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. The reference temperature is used in both the density and viscosity formulas. Units: \si{\kelvin}. 

(parameters:Material_20model:Melt_20global:Solid_20compressibility)=
### __Parameter name:__ Solid compressibility
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the compressibility of the solid matrix. Units: \si{\per\pascal}. 

(parameters:Material_20model:Melt_20global:Surface_20solidus)=
### __Parameter name:__ Surface solidus
**Default value:** 1300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Solidus for a pressure of zero. Units: \si{\kelvin}. 

(parameters:Material_20model:Melt_20global:Thermal_20bulk_20viscosity_20exponent)=
### __Parameter name:__ Thermal bulk viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of the bulk viscosity. Dimensionless exponent. See the general documentation of this model for a formula that states the dependence of the viscosity on this factor, which is called $\beta$ there. 

(parameters:Material_20model:Melt_20global:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Melt_20global:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Melt_20global:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of the shear viscosity. Dimensionless exponent. See the general documentation of this model for a formula that states the dependence of the viscosity on this factor, which is called $\beta$ there. 

(parameters:Material_20model:Melt_20simple)=
## **Parameters in section** Material model/Melt simple
(parameters:Material_20model:Melt_20simple:A1)=
### __Parameter name:__ A1
**Default value:** 1085.7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Melt_20simple:A2)=
### __Parameter name:__ A2
**Default value:** 1.329e-7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Melt_20simple:A3)=
### __Parameter name:__ A3
**Default value:** -5.1e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Melt_20simple:B1)=
### __Parameter name:__ B1
**Default value:** 1475.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Melt_20simple:B2)=
### __Parameter name:__ B2
**Default value:** 8.0e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Melt_20simple:B3)=
### __Parameter name:__ B3
**Default value:** -3.2e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Melt_20simple:C1)=
### __Parameter name:__ C1
**Default value:** 1780.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius}. 

(parameters:Material_20model:Melt_20simple:C2)=
### __Parameter name:__ C2
**Default value:** 4.50e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius\per\pascal}. 

(parameters:Material_20model:Melt_20simple:C3)=
### __Parameter name:__ C3
**Default value:** -2.0e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Material_20model:Melt_20simple:Depletion_20density_20change)=
### __Parameter name:__ Depletion density change
**Default value:** 0.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The density contrast between material with a depletion of 1 and a depletion of zero. Negative values indicate lower densities of depleted material. Depletion is indicated by the compositional field with the name peridotite. Not used if this field does not exist in the model. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Melt_20simple:Depletion_20solidus_20change)=
### __Parameter name:__ Depletion solidus change
**Default value:** 200.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The solidus temperature change for a depletion of 100\%. For positive values, the solidus gets increased for a positive peridotite field (depletion) and lowered for a negative peridotite field (enrichment). Scaling with depletion is linear. Only active when fractional melting is used. Units: \si{\kelvin}. 

(parameters:Material_20model:Melt_20simple:Exponential_20melt_20weakening_20factor)=
### __Parameter name:__ Exponential melt weakening factor
**Default value:** 27. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The porosity dependence of the viscosity. Units: dimensionless. 

(parameters:Material_20model:Melt_20simple:Freezing_20rate)=
### __Parameter name:__ Freezing rate
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Freezing rate of melt when in subsolidus regions. If this parameter is set to a number larger than 0.0, it specifies the fraction of melt that will freeze per year (or per second, depending on the ``Use years in output instead of seconds'' parameter), as soon as the porosity exceeds the equilibrium melt fraction, and the equilibrium melt fraction falls below the depletion. In this case, melt will freeze according to the given rate until one of those conditions is not fulfilled anymore. The reasoning behind this is that there should not be more melt present than the equilibrium melt fraction, as melt production decreases with increasing depletion, but the freezing process of melt also reduces the depletion by the same amount, and as soon as the depletion falls below the equilibrium melt fraction, we expect that material should melt again (no matter how much melt is present). This is quite a simplification and not a realistic freezing parameterization, but without tracking the melt composition, there is no way to compute freezing rates accurately. If this parameter is set to zero, no freezing will occur. Note that freezing can never be faster than determined by the ``Melting time scale for operator splitting''. The product of the ``Freezing rate'' and the ``Melting time scale for operator splitting'' defines how fast freezing occurs with respect to melting (if the product is 0.5, melting will occur twice as fast as freezing). Units: 1/yr or 1/s, depending on the ``Use years in output instead of seconds'' parameter. 

(parameters:Material_20model:Melt_20simple:Mass_20fraction_20cpx)=
### __Parameter name:__ Mass fraction cpx
**Default value:** 0.15 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Mass fraction of clinopyroxene in the peridotite to be molten. Units: non-dimensional. 

(parameters:Material_20model:Melt_20simple:Melt_20bulk_20modulus_20derivative)=
### __Parameter name:__ Melt bulk modulus derivative
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the pressure derivative of the melt bulk modulus. Units: None. 

(parameters:Material_20model:Melt_20simple:Melt_20compressibility)=
### __Parameter name:__ Melt compressibility
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the compressibility of the melt. Units: \si{\per\pascal}. 

(parameters:Material_20model:Melt_20simple:Melt_20extraction_20depth)=
### __Parameter name:__ Melt extraction depth
**Default value:** 1000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Depth above that melt will be extracted from the model, which is done by a negative reaction term proportional to the porosity field. Units: \si{\meter}. 

(parameters:Material_20model:Melt_20simple:Melting_20time_20scale_20for_20operator_20splitting)=
### __Parameter name:__ Melting time scale for operator splitting
**Default value:** 1e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Because the operator splitting scheme is used, the porosity field can not be set to a new equilibrium melt fraction instantly, but the model has to provide a melting time scale instead. This time scale defines how fast melting happens, or more specifically, the parameter defines the time after which the deviation of the porosity from the equilibrium melt fraction will be reduced to a fraction of $1/e$. So if the melting time scale is small compared to the time step size, the reaction will be so fast that the porosity is very close to the equilibrium melt fraction after reactions are computed. Conversely, if the melting time scale is large compared to the time step size, almost no melting and freezing will occur.

Also note that the melting time scale has to be larger than or equal to the reaction time step used in the operator splitting scheme, otherwise reactions can not be computed. Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter. 

(parameters:Material_20model:Melt_20simple:Peridotite_20melting_20entropy_20change)=
### __Parameter name:__ Peridotite melting entropy change
**Default value:** -300. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The entropy change for the phase transition from solid to melt of peridotite. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Melt_20simple:Reference_20bulk_20viscosity)=
### __Parameter name:__ Reference bulk viscosity
**Default value:** 1e22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant bulk viscosity $\xi_0$ of the solid matrix. This viscosity may be modified by both temperature and porosity dependencies. Units: \si{\pascal\second}. 

(parameters:Material_20model:Melt_20simple:Reference_20melt_20density)=
### __Parameter name:__ Reference melt density
**Default value:** 2500. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density of the melt/fluid$\rho_{f,0}$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Melt_20simple:Reference_20melt_20viscosity)=
### __Parameter name:__ Reference melt viscosity
**Default value:** 10. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant melt viscosity $\eta_f$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Melt_20simple:Reference_20permeability)=
### __Parameter name:__ Reference permeability
**Default value:** 1e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference permeability of the solid host rock.Units: \si{\meter\squared}. 

(parameters:Material_20model:Melt_20simple:Reference_20shear_20viscosity)=
### __Parameter name:__ Reference shear viscosity
**Default value:** 5e20 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity $\eta_0$ of the solid matrix. This viscosity may be modified by both temperature and porosity dependencies. Units: \si{\pascal\second}. 

(parameters:Material_20model:Melt_20simple:Reference_20solid_20density)=
### __Parameter name:__ Reference solid density
**Default value:** 3000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density of the solid $\rho_{s,0}$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Melt_20simple:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Melt_20simple:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. The reference temperature is used in both the density and viscosity formulas. Units: \si{\kelvin}. 

(parameters:Material_20model:Melt_20simple:Solid_20compressibility)=
### __Parameter name:__ Solid compressibility
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the compressibility of the solid matrix. Units: \si{\per\pascal}. 

(parameters:Material_20model:Melt_20simple:Thermal_20bulk_20viscosity_20exponent)=
### __Parameter name:__ Thermal bulk viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of the bulk viscosity. Dimensionless exponent. See the general documentation of this model for a formula that states the dependence of the viscosity on this factor, which is called $\beta$ there. 

(parameters:Material_20model:Melt_20simple:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Melt_20simple:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Melt_20simple:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of the shear viscosity. Dimensionless exponent. See the general documentation of this model for a formula that states the dependence of the viscosity on this factor, which is called $\beta$ there. 

(parameters:Material_20model:Melt_20simple:Use_20fractional_20melting)=
### __Parameter name:__ Use fractional melting
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If fractional melting should be used (if true), including a solidus change based on depletion (in this case, the amount of melt that has migrated away from its origin), and freezing of melt when it has moved to a region with temperatures lower than the solidus; or if batch melting should be used (if false), assuming that the melt fraction only depends on temperature and pressure, and how much melt has already been generated at a given point, but not considering movement of melt in the melting parameterization.

Note that melt does not freeze unless the 'Freezing rate' parameter is set to a value larger than 0. 

(parameters:Material_20model:Melt_20simple:Use_20full_20compressibility)=
### __Parameter name:__ Use full compressibility
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If the compressibility should be used everywhere in the code (if true), changing the volume of material when the density changes, or only in the momentum conservation and advection equations (if false). 

(parameters:Material_20model:Melt_20simple:beta)=
### __Parameter name:__ beta
**Default value:** 1.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Exponent of the melting temperature in the melt fraction calculation. Units: non-dimensional. 

(parameters:Material_20model:Melt_20simple:r1)=
### __Parameter name:__ r1
**Default value:** 0.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant in the linear function that approximates the clinopyroxene reaction coefficient. Units: non-dimensional. 

(parameters:Material_20model:Melt_20simple:r2)=
### __Parameter name:__ r2
**Default value:** 8e-11 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the linear function that approximates the clinopyroxene reaction coefficient. Units: \si{\per\pascal}. 

(parameters:Material_20model:Modified_20Tait_20model)=
## **Parameters in section** Material model/Modified Tait model
(parameters:Material_20model:Modified_20Tait_20model:Einstein_20temperature)=
### __Parameter name:__ Einstein temperature
**Default value:** 600. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The Einstein temperature at the reference pressure and temperature. Units: \si{\kelvin}. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20bulk_20modulus_20derivative)=
### __Parameter name:__ Reference bulk modulus derivative
**Default value:** 4. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the first pressure derivative of the isothermal bulk modulus at the reference pressure and temperature. Units: None. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The density at the reference pressure and temperature. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20isothermal_20bulk_20modulus)=
### __Parameter name:__ Reference isothermal bulk modulus
**Default value:** 125e9 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The isothermal bulk modulus at the reference pressure and temperature. Units: \si{\pascal}. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20pressure)=
### __Parameter name:__ Reference pressure
**Default value:** 1e5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference pressure $P_0$. Units: \si{\pascal}. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 298.15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20thermal_20expansivity)=
### __Parameter name:__ Reference thermal expansivity
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The thermal expansion coefficient at the reference pressure and temperature. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Modified_20Tait_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Modified_20Tait_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 1e21 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity $\eta_0$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20heat_20capacity_20function)=
## **Parameters in section** Material model/Modified Tait model/Reference heat capacity function
(parameters:Material_20model:Modified_20Tait_20model:Reference_20heat_20capacity_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Material_20model:Modified_20Tait_20model:Reference_20heat_20capacity_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 1.25e3 

**Pattern:** [Anything] 

**Documentation:**  

(parameters:Material_20model:Modified_20Tait_20model:Reference_20heat_20capacity_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Material_20model:Multicomponent)=
## **Parameters in section** Material model/Multicomponent
(parameters:Material_20model:Multicomponent:Densities)=
### __Parameter name:__ Densities
**Default value:** 3300. 

**Pattern:** [Anything] 

**Documentation:** List of densities for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Multicomponent:Heat_20capacities)=
### __Parameter name:__ Heat capacities
**Default value:** 1250. 

**Pattern:** [Anything] 

**Documentation:** List of specific heats $C_p$ for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Multicomponent:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Multicomponent:Specific_20heats)=
### __Parameter name__: Specific heats
**Alias:** [Heat capacities](parameters:Material_20model:Multicomponent:Heat_20capacities)

**Deprecation Status:** false

(parameters:Material_20model:Multicomponent:Thermal_20conductivities)=
### __Parameter name:__ Thermal conductivities
**Default value:** 4.7 

**Pattern:** [Anything] 

**Documentation:** List of thermal conductivities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Multicomponent:Thermal_20expansivities)=
### __Parameter name:__ Thermal expansivities
**Default value:** 0.000040 

**Pattern:** [Anything] 

**Documentation:** List of thermal expansivities for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Multicomponent:Viscosities)=
### __Parameter name:__ Viscosities
**Default value:** 1.e21 

**Pattern:** [Anything] 

**Documentation:** List of viscosities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\pascal\second}. 

(parameters:Material_20model:Multicomponent:Viscosity_20averaging_20scheme)=
### __Parameter name:__ Viscosity averaging scheme
**Default value:** harmonic 

**Pattern:** [Selection arithmetic|harmonic|geometric|maximum composition ] 

**Documentation:** When more than one compositional field is present at a point with different viscosities, we need to come up with an average viscosity at that point.  Select a weighted harmonic, arithmetic, geometric, or maximum composition. 

(parameters:Material_20model:Multicomponent_20compressible)=
## **Parameters in section** Material model/Multicomponent compressible
(parameters:Material_20model:Multicomponent_20compressible:Isochoric_20specific_20heats)=
### __Parameter name:__ Isochoric specific heats
**Default value:** 1250. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of isochoric specific heats $C_v$ for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Multicomponent_20compressible:Isothermal_20bulk_20modulus_20pressure_20derivatives)=
### __Parameter name:__ Isothermal bulk modulus pressure derivatives
**Default value:** 4. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of isothermal pressure derivatives of the bulk moduli for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: []. 

(parameters:Material_20model:Multicomponent_20compressible:Reference_20densities)=
### __Parameter name:__ Reference densities
**Default value:** 3300. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of densities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Multicomponent_20compressible:Reference_20isothermal_20compressibilities)=
### __Parameter name:__ Reference isothermal compressibilities
**Default value:** 4e-12 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of isothermal compressibilities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\per\pascal}. 

(parameters:Material_20model:Multicomponent_20compressible:Reference_20temperatures)=
### __Parameter name:__ Reference temperatures
**Default value:** 298.15 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of reference temperatures $T_0$ for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\kelvin}. 

(parameters:Material_20model:Multicomponent_20compressible:Reference_20thermal_20expansivities)=
### __Parameter name:__ Reference thermal expansivities
**Default value:** 4.e-5 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of thermal expansivities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Multicomponent_20compressible:Thermal_20conductivities)=
### __Parameter name:__ Thermal conductivities
**Default value:** 4.7 

**Pattern:** [Anything] 

**Documentation:** List of thermal conductivities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Multicomponent_20compressible:Viscosities)=
### __Parameter name:__ Viscosities
**Default value:** 1.e21 

**Pattern:** [Anything] 

**Documentation:** List of viscosities for background mantle and compositional fields,for a total of N+1 values, where N is the number of compositional fields.If only one value is given, then all use the same value. Units: \si{\pascal\second}. 

(parameters:Material_20model:Multicomponent_20compressible:Viscosity_20averaging_20scheme)=
### __Parameter name:__ Viscosity averaging scheme
**Default value:** harmonic 

**Pattern:** [Selection arithmetic|harmonic|geometric|maximum composition ] 

**Documentation:** When more than one compositional field is present at a point with different viscosities, we need to come up with an average viscosity at that point.  Select a weighted harmonic, arithmetic, geometric, or maximum composition. 

(parameters:Material_20model:Nondimensional_20model)=
## **Parameters in section** Material model/Nondimensional model
(parameters:Material_20model:Nondimensional_20model:Di)=
### __Parameter name:__ Di
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Dissipation number. Pick 0.0 for incompressible computations. 

(parameters:Material_20model:Nondimensional_20model:Ra)=
### __Parameter name:__ Ra
**Default value:** 1e4 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Rayleigh number Ra 

(parameters:Material_20model:Nondimensional_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Nondimensional_20model:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Nondimensional_20model:Use_20TALA)=
### __Parameter name:__ Use TALA
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use the TALA instead of the ALA approximation. 

(parameters:Material_20model:Nondimensional_20model:Viscosity_20depth_20prefactor)=
### __Parameter name:__ Viscosity depth prefactor
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Exponential depth prefactor for viscosity. 

(parameters:Material_20model:Nondimensional_20model:Viscosity_20temperature_20prefactor)=
### __Parameter name:__ Viscosity temperature prefactor
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Exponential temperature prefactor for viscosity. 

(parameters:Material_20model:Nondimensional_20model:gamma)=
### __Parameter name:__ gamma
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Grueneisen parameter 

(parameters:Material_20model:PerpleX_20lookup_20model)=
## **Parameters in section** Material model/PerpleX lookup model
(parameters:Material_20model:PerpleX_20lookup_20model:Maximum_20material_20pressure)=
### __Parameter name:__ Maximum material pressure
**Default value:** 1.e12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the maximum pressure used to query PerpleX. Units: \si{\pascal}. 

(parameters:Material_20model:PerpleX_20lookup_20model:Maximum_20material_20temperature)=
### __Parameter name:__ Maximum material temperature
**Default value:** 6000. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the maximum temperature used to query PerpleX. Units: \si{\kelvin}. 

(parameters:Material_20model:PerpleX_20lookup_20model:Minimum_20material_20pressure)=
### __Parameter name:__ Minimum material pressure
**Default value:** 1.e5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the minimum pressure used to query PerpleX. Units: \si{\pascal}. 

(parameters:Material_20model:PerpleX_20lookup_20model:Minimum_20material_20temperature)=
### __Parameter name:__ Minimum material temperature
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the minimum temperature used to query PerpleX. Units: \si{\kelvin}. 

(parameters:Material_20model:PerpleX_20lookup_20model:PerpleX_20input_20file_20name)=
### __Parameter name:__ PerpleX input file name
**Default value:** rock.dat 

**Pattern:** [Anything] 

**Documentation:** The name of the PerpleX input file (should end with .dat). 

(parameters:Material_20model:PerpleX_20lookup_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:PerpleX_20lookup_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the viscosity $\eta$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Replace_20lithosphere_20viscosity)=
## **Parameters in section** Material model/Replace lithosphere viscosity
(parameters:Material_20model:Replace_20lithosphere_20viscosity:Base_20model)=
### __Parameter name:__ Base model
**Default value:** simple 

**Pattern:** [Selection Steinberger|ascii reference profile|averaging|compositing|composition reaction|depth dependent|diffusion dislocation|drucker prager|grain size|latent heat|latent heat melt|melt global|melt simple|modified tait|multicomponent|multicomponent compressible|nondimensional|perplex lookup|replace lithosphere viscosity|simple|simple compressible|simpler|visco plastic|viscoelastic ] 

**Documentation:** The name of a material model that will be modified by a replacingthe viscosity in the lithosphere by a constant value. Valid values for this parameter are the names of models that are also valid for the ``Material models/Model name'' parameter. See the documentation for more information. 

(parameters:Material_20model:Replace_20lithosphere_20viscosity:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/lithosphere-mask/ 

**Pattern:** [DirectoryName] 

**Documentation:** The path to the LAB depth data file 

(parameters:Material_20model:Replace_20lithosphere_20viscosity:Depth_20specification_20method)=
### __Parameter name:__ Depth specification method
**Default value:** Value 

**Pattern:** [Selection File|Value ] 

**Documentation:** Method that is used to specify the depth of the lithosphere-asthenosphere boundary. 

(parameters:Material_20model:Replace_20lithosphere_20viscosity:LAB_20depth_20filename)=
### __Parameter name:__ LAB depth filename
**Default value:** LAB_CAM2016.txt 

**Pattern:** [FileName (Type: input)] 

**Documentation:** File from which the lithosphere-asthenosphere boundary depth data is read. 

(parameters:Material_20model:Replace_20lithosphere_20viscosity:Lithosphere_20viscosity)=
### __Parameter name:__ Lithosphere viscosity
**Default value:** 1e23 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The viscosity within lithosphere, applied abovethe maximum lithosphere depth. 

(parameters:Material_20model:Replace_20lithosphere_20viscosity:Maximum_20lithosphere_20depth)=
### __Parameter name:__ Maximum lithosphere depth
**Default value:** 200000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Units: \si{\meter}.The maximum depth of the lithosphere. The model will be NaNs below this depth. 

(parameters:Material_20model:Simple_20compressible_20model)=
## **Parameters in section** Material model/Simple compressible model
(parameters:Material_20model:Simple_20compressible_20model:Reference_20compressibility)=
### __Parameter name:__ Reference compressibility
**Default value:** 4e-12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the reference compressibility. Units: \si{\per\pascal}. 

(parameters:Material_20model:Simple_20compressible_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Simple_20compressible_20model:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Simple_20compressible_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Simple_20compressible_20model:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Simple_20compressible_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 1000000000000000000000.000000 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the viscosity $\eta$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Simple_20model)=
## **Parameters in section** Material model/Simple model
(parameters:Material_20model:Simple_20model:Composition_20viscosity_20prefactor)=
### __Parameter name:__ Composition viscosity prefactor
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A linear dependency of viscosity on the first compositional field. Dimensionless prefactor. With a value of 1.0 (the default) the viscosity does not depend on the composition. See the general documentation of this model for a formula that states the dependence of the viscosity on this factor, which is called $\xi$ there. 

(parameters:Material_20model:Simple_20model:Density_20differential_20for_20compositional_20field_201)=
### __Parameter name:__ Density differential for compositional field 1
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** If compositional fields are used, then one would frequently want to make the density depend on these fields. In this simple material model, we make the following assumptions: if no compositional fields are used in the current simulation, then the density is simply the usual one with its linear dependence on the temperature. If there are compositional fields, then the material model determines how many of them influence the density. The composition-dependence adds a term of the kind $+\Delta \rho \; c_1(\mathbf x)$. This parameter describes the value of $\Delta \rho$. Units: \si{\kilogram\per\meter\cubed}/unit change in composition. 

(parameters:Material_20model:Simple_20model:Maximum_20thermal_20prefactor)=
### __Parameter name:__ Maximum thermal prefactor
**Default value:** 1.0e2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum value of the viscosity prefactor associated with temperature dependence. 

(parameters:Material_20model:Simple_20model:Minimum_20thermal_20prefactor)=
### __Parameter name:__ Minimum thermal prefactor
**Default value:** 1.0e-2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum value of the viscosity prefactor associated with temperature dependence. 

(parameters:Material_20model:Simple_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Simple_20model:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Simple_20model:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. The reference temperature is used in both the density and viscosity formulas. Units: \si{\kelvin}. 

(parameters:Material_20model:Simple_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Simple_20model:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Simple_20model:Thermal_20viscosity_20exponent)=
### __Parameter name:__ Thermal viscosity exponent
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The temperature dependence of viscosity. Dimensionless exponent. See the general documentation of this model for a formula that states the dependence of the viscosity on this factor, which is called $\beta$ there. 

(parameters:Material_20model:Simple_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5e24 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the constant viscosity $\eta_0$. This viscosity may be modified by both temperature and compositional dependencies. Units: \si{\pascal\second}. 

(parameters:Material_20model:Simpler_20model)=
## **Parameters in section** Material model/Simpler model
(parameters:Material_20model:Simpler_20model:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference density $\rho_0$. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Simpler_20model:Reference_20specific_20heat)=
### __Parameter name:__ Reference specific heat
**Default value:** 1250. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the specific heat $C_p$. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Simpler_20model:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. The reference temperature is used in both the density and viscosity formulas. Units: \si{\kelvin}. 

(parameters:Material_20model:Simpler_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Simpler_20model:Thermal_20expansion_20coefficient)=
### __Parameter name:__ Thermal expansion coefficient
**Default value:** 2e-5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal expansion coefficient $\alpha$. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Simpler_20model:Viscosity)=
### __Parameter name:__ Viscosity
**Default value:** 5000000000000000452984832.000000 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the viscosity $\eta$. Units: \si{\pascal\second}. 

(parameters:Material_20model:Steinberger_20model)=
## **Parameters in section** Material model/Steinberger model
(parameters:Material_20model:Steinberger_20model:Bilinear_20interpolation)=
### __Parameter name:__ Bilinear interpolation
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to use bilinear interpolation to compute material properties (slower but more accurate).  

(parameters:Material_20model:Steinberger_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/material-model/steinberger/ 

**Pattern:** [DirectoryName] 

**Documentation:** The path to the model data. The path may also include the special text '$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.  

(parameters:Material_20model:Steinberger_20model:Derivatives_20file_20names)=
### __Parameter name:__ Derivatives file names
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** The file names of the enthalpy derivatives data. List with as many components as active compositional fields (material data is assumed to be in order with the ordering of the fields). 

(parameters:Material_20model:Steinberger_20model:Latent_20heat)=
### __Parameter name:__ Latent heat
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to include latent heat effects in the calculation of thermal expansivity and specific heat. If true, ASPECT follows the approach of Nakagawa et al. 2009, using pressure and temperature derivatives of the enthalpy to calculate the thermal expansivity and specific heat. If false, ASPECT uses the thermal expansivity and specific heat values from the material properties table. 

(parameters:Material_20model:Steinberger_20model:Lateral_20viscosity_20file_20name)=
### __Parameter name:__ Lateral viscosity file name
**Default value:** temp-viscosity-prefactor.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the lateral viscosity data.  

(parameters:Material_20model:Steinberger_20model:Material_20file_20format)=
### __Parameter name:__ Material file format
**Default value:** perplex 

**Pattern:** [Selection perplex|hefesto ] 

**Documentation:** The material file format to be read in the property tables. 

(parameters:Material_20model:Steinberger_20model:Material_20file_20names)=
### __Parameter name:__ Material file names
**Default value:** pyr-ringwood88.txt 

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** The file names of the material data (material data is assumed to be in order with the ordering of the compositional fields). Note that there are three options on how many files need to be listed here: 1. If only one file is provided, it is used for the whole model domain, and compositional fields are ignored. 2. If there is one more file name than the number of compositional fields, then the first file is assumed to define a `background composition' that is modified by the compositional fields. If there are exactly as many files as compositional fields, the fields are assumed to represent the fractions of different materials and the average property is computed as a sum of the value of the compositional field times the material property of that field. 

(parameters:Material_20model:Steinberger_20model:Maximum_20latent_20heat_20substeps)=
### __Parameter name:__ Maximum latent heat substeps
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The maximum number of substeps over the temperature pressure range to calculate the averaged enthalpy gradient over a cell. 

(parameters:Material_20model:Steinberger_20model:Maximum_20lateral_20viscosity_20variation)=
### __Parameter name:__ Maximum lateral viscosity variation
**Default value:** 1e2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The relative cutoff value for lateral viscosity variations caused by temperature deviations. The viscosity may vary laterally by this factor squared. 

(parameters:Material_20model:Steinberger_20model:Maximum_20thermal_20conductivity)=
### __Parameter name:__ Maximum thermal conductivity
**Default value:** 1000 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum thermal conductivity that is allowed in the model. Larger values will be cut off. 

(parameters:Material_20model:Steinberger_20model:Maximum_20viscosity)=
### __Parameter name:__ Maximum viscosity
**Default value:** 1e23 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum viscosity that is allowed in the viscosity calculation. Larger values will be cut off. 

(parameters:Material_20model:Steinberger_20model:Minimum_20viscosity)=
### __Parameter name:__ Minimum viscosity
**Default value:** 1e19 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum viscosity that is allowed in the viscosity calculation. Smaller values will be cut off. 

(parameters:Material_20model:Steinberger_20model:Number_20lateral_20average_20bands)=
### __Parameter name:__ Number lateral average bands
**Default value:** 10 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of bands to compute laterally averaged temperature within. 

(parameters:Material_20model:Steinberger_20model:Pressure_20dependencies_20of_20thermal_20conductivity)=
### __Parameter name:__ Pressure dependencies of thermal conductivity
**Default value:** 3.3e-10, 3.4e-10, 3.6e-10, 1.05e-10 

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of values that determine the linear scaling of the thermal conductivity with the pressure in the 'p-T-dependent' Thermal conductivity formulation. Units: \si{\watt\per\meter\per\kelvin\per\pascal}. 

(parameters:Material_20model:Steinberger_20model:Radial_20viscosity_20file_20name)=
### __Parameter name:__ Radial viscosity file name
**Default value:** radial-visc.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the radial viscosity data.  

(parameters:Material_20model:Steinberger_20model:Reference_20temperatures_20for_20thermal_20conductivity)=
### __Parameter name:__ Reference temperatures for thermal conductivity
**Default value:** 300, 300, 300, 1200 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of values of reference temperatures used to determine the temperature-dependence of the thermal conductivity in the 'p-T-dependent' Thermal conductivity formulation. Units: \si{\kelvin}. 

(parameters:Material_20model:Steinberger_20model:Reference_20thermal_20conductivities)=
### __Parameter name:__ Reference thermal conductivities
**Default value:** 2.47, 3.81, 3.52, 4.9 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of base values of the thermal conductivity for each of the horizontal layers in the 'p-T-dependent' Thermal conductivity formulation. Pressure- and temperature-dependence will be appliedon top of this base value, according to the parameters 'Pressure dependencies of thermal conductivity' and 'Reference temperatures for thermal conductivity'. Units: \si{\watt\per\meter\per\kelvin} 

(parameters:Material_20model:Steinberger_20model:Reference_20viscosity)=
### __Parameter name:__ Reference viscosity
**Default value:** 1e23 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference viscosity that is used for pressure scaling. To understand how pressure scaling works, take a look at \cite{KHB12}. In particular, the value of this parameter would not affect the solution computed by \aspect{} if we could do arithmetic exactly; however, computers do arithmetic in finite precision, and consequently we need to scale quantities in ways so that their magnitudes are roughly the same. As explained in \cite{KHB12}, we scale the pressure during some computations (never visible by users) by a factor that involves a reference viscosity. This parameter describes this reference viscosity.

For problems with a constant viscosity, you will generally want to choose the reference viscosity equal to the actual viscosity. For problems with a variable viscosity, the reference viscosity should be a value that adequately represents the order of magnitude of the viscosities that appear, such as an average value or the value one would use to compute a Rayleigh number.

Units: \si{\pascal\second}. 

(parameters:Material_20model:Steinberger_20model:Saturation_20prefactors)=
### __Parameter name:__ Saturation prefactors
**Default value:** 0, 0, 0, 1 

**Pattern:** [List of <[Double 0...1 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of values that indicate how a given layer in the conductivity formulation should take into account the effects of saturation on the temperature-dependence of the thermal conducitivity. This factor is multiplied with a saturation function based on the theory of Roufosse and Klemens, 1974. A value of 1 reproduces the formulation of Stackhouse et al. (2015), a value of 0 reproduces the formulation of Tosi et al., (2013). Units: none. 

(parameters:Material_20model:Steinberger_20model:Thermal_20conductivity)=
### __Parameter name:__ Thermal conductivity
**Default value:** 4.7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The value of the thermal conductivity $k$. Only used in case the 'constant' Thermal conductivity formulation is selected. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Steinberger_20model:Thermal_20conductivity_20exponents)=
### __Parameter name:__ Thermal conductivity exponents
**Default value:** 0.48, 0.56, 0.61, 1.0 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of exponents in the temperature-dependent term of the 'p-T-dependent' Thermal conductivity formulation. Note that this exponent is not used (and should have a value of 1) in the formulation of Stackhouse et al. (2015). Units: none. 

(parameters:Material_20model:Steinberger_20model:Thermal_20conductivity_20formulation)=
### __Parameter name:__ Thermal conductivity formulation
**Default value:** constant 

**Pattern:** [Selection constant|p-T-dependent ] 

**Documentation:** Which law should be used to compute the thermal conductivity. The 'constant' law uses a constant value for the thermal conductivity. The 'p-T-dependent' formulation uses equations from Stackhouse et al. (2015): First-principles calculations of the lattice thermal conductivity of the lower mantle (https://doi.org/10.1016/j.epsl.2015.06.050), and Tosi et al. (2013): Mantle dynamics with pressure- and temperature-dependent thermal expansivity and conductivity (https://doi.org/10.1016/j.pepi.2013.02.004) to compute the thermal conductivity in dependence of temperature and pressure. The thermal conductivity parameter sets can be chosen in such a way that either the Stackhouse or the Tosi relations are used. The conductivity description can consist of several layers with different sets of parameters. Note that the Stackhouse parametrization is only valid for the lower mantle (bridgmanite). 

(parameters:Material_20model:Steinberger_20model:Thermal_20conductivity_20transition_20depths)=
### __Parameter name:__ Thermal conductivity transition depths
**Default value:** 410000, 520000, 660000 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of depth values that indicate where the transitions between the different conductivity parameter sets should occur in the 'p-T-dependent' Thermal conductivity formulation (in most cases, this will be the depths of major mantle phase transitions). Units: \si{\meter}. 

(parameters:Material_20model:Steinberger_20model:Use_20lateral_20average_20temperature_20for_20viscosity)=
### __Parameter name:__ Use lateral average temperature for viscosity
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to use to use the laterally averaged temperature instead of the adiabatic temperature as reference for the viscosity calculation. This ensures that the laterally averaged viscosities remain more or less constant over the model runtime. This behaviour might or might not be desired. 

(parameters:Material_20model:Visco_20Plastic)=
## **Parameters in section** Material model/Visco Plastic
(parameters:Material_20model:Visco_20Plastic:Activation_20energies_20for_20Peierls_20creep)=
### __Parameter name:__ Activation energies for Peierls creep
**Default value:** 320e3 

**Pattern:** [Anything] 

**Documentation:** List of activation energies, $E$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Visco_20Plastic:Activation_20energies_20for_20diffusion_20creep)=
### __Parameter name:__ Activation energies for diffusion creep
**Default value:** 375e3 

**Pattern:** [Anything] 

**Documentation:** List of activation energies, $E_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Visco_20Plastic:Activation_20energies_20for_20dislocation_20creep)=
### __Parameter name:__ Activation energies for dislocation creep
**Default value:** 530e3 

**Pattern:** [Anything] 

**Documentation:** List of activation energies, $E_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\joule\per\mole}. 

(parameters:Material_20model:Visco_20Plastic:Activation_20volumes_20for_20Peierls_20creep)=
### __Parameter name:__ Activation volumes for Peierls creep
**Default value:** 1.4e-5 

**Pattern:** [Anything] 

**Documentation:** List of activation volumes, $V$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Visco_20Plastic:Activation_20volumes_20for_20diffusion_20creep)=
### __Parameter name:__ Activation volumes for diffusion creep
**Default value:** 6e-6 

**Pattern:** [Anything] 

**Documentation:** List of activation volumes, $V_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Visco_20Plastic:Activation_20volumes_20for_20dislocation_20creep)=
### __Parameter name:__ Activation volumes for dislocation creep
**Default value:** 1.4e-5 

**Pattern:** [Anything] 

**Documentation:** List of activation volumes, $V_a$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\meter\cubed\per\mole}. 

(parameters:Material_20model:Visco_20Plastic:Adiabat_20temperature_20gradient_20for_20viscosity)=
### __Parameter name:__ Adiabat temperature gradient for viscosity
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Add an adiabatic temperature gradient to the temperature used in the flow law so that the activation volume is consistent with what one would use in a earth-like (compressible) model. Default is set so this is off. Note that this is a linear approximation of the real adiabatic gradient, which is okay for the upper mantle, but is not really accurate for the lower mantle. Using a pressure gradient of 32436 Pa/m, then a value of 0.3 K/km = 0.0003 K/m = 9.24e-09 K/Pa gives an earth-like adiabat.Units: \si{\kelvin\per\pascal}. 

(parameters:Material_20model:Visco_20Plastic:Allow_20negative_20pressures_20in_20plasticity)=
### __Parameter name:__ Allow negative pressures in plasticity
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to allow negative pressures to be used in the computation of plastic yield stresses and viscosities. Setting this parameter to true may be advantageous in models without gravity where the dynamic stresses are much higher than the lithostatic pressure. If false, the minimum pressure in the plasticity formulation will be set to zero. 

(parameters:Material_20model:Visco_20Plastic:Angles_20of_20internal_20friction)=
### __Parameter name:__ Angles of internal friction
**Default value:** 0. 

**Pattern:** [Anything] 

**Documentation:** List of angles of internal friction, $\phi$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. For a value of zero, in 2D the von Mises criterion is retrieved. Angles higher than 30 degrees are harder to solve numerically. Units: degrees. 

(parameters:Material_20model:Visco_20Plastic:Cohesion_20strain_20weakening_20factors)=
### __Parameter name:__ Cohesion strain weakening factors
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of cohesion strain weakening factors for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Cohesions)=
### __Parameter name:__ Cohesions
**Default value:** 1e20 

**Pattern:** [Anything] 

**Documentation:** List of cohesions, $C$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from exceeding the yield stress. Units: \si{\pascal}. 

(parameters:Material_20model:Visco_20Plastic:Constant_20viscosity_20prefactors)=
### __Parameter name:__ Constant viscosity prefactors
**Default value:** 1.0 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of constant viscosity prefactors (i.e., multiplicative factors) for background material and compositional fields, for a total of N+1 where N is the number of compositional fields. Units: none. 

(parameters:Material_20model:Visco_20Plastic:Cutoff_20stresses_20for_20Peierls_20creep)=
### __Parameter name:__ Cutoff stresses for Peierls creep
**Default value:** 0.0 

**Pattern:** [Anything] 

**Documentation:** List of the Stress thresholds below which the strain rate is solved for as a quadratic function of stress to aid with convergence when stress exponent n=0. Units: \si{\pascal} 

(parameters:Material_20model:Visco_20Plastic:Define_20thermal_20conductivities)=
### __Parameter name:__ Define thermal conductivities
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to directly define thermal conductivities for each compositional field instead of calculating the values through the specified thermal diffusivities, densities, and heat capacities.  

(parameters:Material_20model:Visco_20Plastic:Define_20transition_20by_20depth_20instead_20of_20pressure)=
### __Parameter name:__ Define transition by depth instead of pressure
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to list phase transitions by depth or pressure. If this parameter is true, then the input file will use Phase transitions depths and Phase transition widths to define the phase transition. If it is false, the parameter file will read in phase transition data from Phase transition pressures and Phase transition pressure widths. 

(parameters:Material_20model:Visco_20Plastic:Densities)=
### __Parameter name:__ Densities
**Default value:** 3300. 

**Pattern:** [Anything] 

**Documentation:** List of densities for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Visco_20Plastic:Dynamic_20angles_20of_20internal_20friction)=
### __Parameter name:__ Dynamic angles of internal friction
**Default value:** 2 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of dynamic angles of internal friction, $\phi$, for background material and compositional fields, for a total of N$+$1 values, where N is the number of compositional fields. Dynamic angles of friction are used as the current friction angle when the effective strain rate is well above the 'dynamic characteristic strain rate'. Units: \si{\degree}. 

(parameters:Material_20model:Visco_20Plastic:Dynamic_20characteristic_20strain_20rate)=
### __Parameter name:__ Dynamic characteristic strain rate
**Default value:** 1e-12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The characteristic strain rate value at which the angle of friction is equal to $\mu = (\mu_s+\mu_d)/2$. When the effective strain rate is very high, the dynamic angle of friction is taken, when it is very low, the static angle of internal friction is used. Around the dynamic characteristic strain rate, there is a smooth gradient from the static to the dynamic angle of internal friction. Units: \si{\per\second}. 

(parameters:Material_20model:Visco_20Plastic:Dynamic_20friction_20smoothness_20exponent)=
### __Parameter name:__ Dynamic friction smoothness exponent
**Default value:** 1 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** An exponential factor in the equation for the calculation of the friction angle when a static and a dynamic angle of internal friction are specified. A factor of 1 returns the equation to Equation (13) in \cite{van_dinther_seismic_2013}. A factor between 0 and 1 makes the curve of the friction angle vs. the strain rate smoother, while a factor $>$ 1 makes the change between static and dynamic friction angle more steplike. Units: none. 

(parameters:Material_20model:Visco_20Plastic:Elastic_20damper_20viscosity)=
### __Parameter name:__ Elastic damper viscosity
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Viscosity of a viscous damper that acts in parallel with the elastic element to stabilize behavior. Units: \si{\pascal\second} 

(parameters:Material_20model:Visco_20Plastic:Elastic_20shear_20moduli)=
### __Parameter name:__ Elastic shear moduli
**Default value:** 75.0e9 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of elastic shear moduli, $G$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. The default value of 75 GPa is representative of mantle rocks. Units: Pa. 

(parameters:Material_20model:Visco_20Plastic:End_20plasticity_20strain_20weakening_20intervals)=
### __Parameter name:__ End plasticity strain weakening intervals
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of strain weakening interval final strains for the cohesion and friction angle parameters of the background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:End_20prefactor_20strain_20weakening_20intervals)=
### __Parameter name:__ End prefactor strain weakening intervals
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of strain weakening interval final strains for the diffusion and dislocation prefactor parameters of the background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Fixed_20elastic_20time_20step)=
### __Parameter name:__ Fixed elastic time step
**Default value:** 1.e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The fixed elastic time step $dte$. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Material_20model:Visco_20Plastic:Friction_20mechanism)=
### __Parameter name:__ Friction mechanism
**Default value:** none 

**Pattern:** [Selection none|dynamic friction|function ] 

**Documentation:** Whether to make the friction angle dependent on strain rate or not. This rheology is intended to be used together with the visco-plastic rheology model.

\item ``none'': No dependence of the friction angle is applied. 

\item ``dynamic friction'': The friction angle is rate dependent.When 'dynamic angles of internal friction' are specified, the friction angle will be weakened for high strain rates with: $\mu = \mu_d + \frac{\mu_s-\mu_d}{1+\frac{\dot{\epsilon}_{ii}}{\dot{\epsilon}_C}}^x$  where $\mu_s$ and $\mu_d$ are the friction angles at low and high strain rates, respectively. $\dot{\epsilon}_{ii}$ is the second invariant of the strain rate and $\dot{\epsilon}_C$ is the 'dynamic characteristic strain rate' where $\mu = (\mu_s+\mu_d)/2$. The 'dynamic friction smoothness exponent' x controls how smooth or step-like the change from $\mu_s$ to $\mu_d$ is. The equation is modified after Equation (13) in \cite{van_dinther_seismic_2013}. $\mu_s$ and $\mu_d$ can be specified by setting 'Angles of internal friction' and 'Dynamic angles of internal friction', respectively. This relationship is similar to rate-and-state friction constitutive relationships, which are applicable to the strength of rocks during earthquakes.

\item ``function'': Specify the friction angle as a function of space and time for each compositional field. 

(parameters:Material_20model:Visco_20Plastic:Friction_20strain_20weakening_20factors)=
### __Parameter name:__ Friction strain weakening factors
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of friction strain weakening factors for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Grain_20size)=
### __Parameter name:__ Grain size
**Default value:** 1e-3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Units: \si{\meter}. 

(parameters:Material_20model:Visco_20Plastic:Grain_20size_20exponents_20for_20diffusion_20creep)=
### __Parameter name:__ Grain size exponents for diffusion creep
**Default value:** 3. 

**Pattern:** [Anything] 

**Documentation:** List of grain size exponents, $m_{\text{diffusion}}$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: None. 

(parameters:Material_20model:Visco_20Plastic:Heat_20capacities)=
### __Parameter name:__ Heat capacities
**Default value:** 1250. 

**Pattern:** [Anything] 

**Documentation:** List of specific heats $C_p$ for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Visco_20Plastic:Include_20Peierls_20creep)=
### __Parameter name:__ Include Peierls creep
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to include Peierls creep in the rheological formulation. 

(parameters:Material_20model:Visco_20Plastic:Include_20viscoelasticity)=
### __Parameter name:__ Include viscoelasticity
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to include elastic effects in the rheological formulation. 

(parameters:Material_20model:Visco_20Plastic:Maximum_20Peierls_20strain_20rate_20iterations)=
### __Parameter name:__ Maximum Peierls strain rate iterations
**Default value:** 40 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Maximum number of iterations to find the correct Peierls strain rate. 

(parameters:Material_20model:Visco_20Plastic:Maximum_20viscosity)=
### __Parameter name:__ Maximum viscosity
**Default value:** 1e28 

**Pattern:** [Anything] 

**Documentation:** Upper cutoff for effective viscosity. Units: \si{\pascal\second}. List with as many components as active compositional fields (material data is assumed to be in order with the ordering of the fields).  

(parameters:Material_20model:Visco_20Plastic:Maximum_20yield_20stress)=
### __Parameter name:__ Maximum yield stress
**Default value:** 1e12 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Limits the maximum value of the yield stress determined by the Drucker-Prager plasticity parameters. Default value is chosen so this is not automatically used. Values of 100e6--1000e6 $Pa$ have been used in previous models. Units: \si{\pascal}. 

(parameters:Material_20model:Visco_20Plastic:Minimum_20strain_20rate)=
### __Parameter name:__ Minimum strain rate
**Default value:** 1.0e-20 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Stabilizes strain dependent viscosity. Units: \si{\per\second}. 

(parameters:Material_20model:Visco_20Plastic:Minimum_20viscosity)=
### __Parameter name:__ Minimum viscosity
**Default value:** 1e17 

**Pattern:** [Anything] 

**Documentation:** Lower cutoff for effective viscosity. Units: \si{\pascal\second}. List with as many components as active compositional fields (material data is assumed to be in order with the ordering of the fields).  

(parameters:Material_20model:Visco_20Plastic:Peierls_20creep_20flow_20law)=
### __Parameter name:__ Peierls creep flow law
**Default value:** viscosity approximation 

**Pattern:** [Selection viscosity approximation|exact ] 

**Documentation:** Select what type of Peierls creep flow law to use. Currently, the available options are 'exact', which uses a Newton-Raphson iterative method to find the stress and then compute viscosity, and 'viscosity approximation', in which viscosity is an explicit function of the strain rate invariant, rather than stress.  

(parameters:Material_20model:Visco_20Plastic:Peierls_20fitting_20parameters)=
### __Parameter name:__ Peierls fitting parameters
**Default value:** 0.17 

**Pattern:** [Anything] 

**Documentation:** List of fitting parameters $\gamma$ between stress $\sigma$ and the Peierls stress $\sigma_{\text{peierls}}$ for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: none 

(parameters:Material_20model:Visco_20Plastic:Peierls_20glide_20parameters_20p)=
### __Parameter name:__ Peierls glide parameters p
**Default value:** 0.5 

**Pattern:** [Anything] 

**Documentation:** List of the first Peierls creep glide parameters, $p$, for background and compositional fields for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: none 

(parameters:Material_20model:Visco_20Plastic:Peierls_20glide_20parameters_20q)=
### __Parameter name:__ Peierls glide parameters q
**Default value:** 1.0 

**Pattern:** [Anything] 

**Documentation:** List of the second Peierls creep glide parameters, $q$, for background and compositional fields for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: none 

(parameters:Material_20model:Visco_20Plastic:Peierls_20strain_20rate_20residual_20tolerance)=
### __Parameter name:__ Peierls strain rate residual tolerance
**Default value:** 1e-22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Tolerance for the iterative solve to find the correct Peierls creep strain rate. 

(parameters:Material_20model:Visco_20Plastic:Peierls_20stresses)=
### __Parameter name:__ Peierls stresses
**Default value:** 5.e9 

**Pattern:** [Anything] 

**Documentation:** List of stress limits for Peierls creep $\sigma_{\text{peierls}}$ for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\pascal} 

(parameters:Material_20model:Visco_20Plastic:Phase_20transition_20Clapeyron_20slopes)=
### __Parameter name:__ Phase transition Clapeyron slopes
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of Clapeyron slopes for each phase transition. A positive Clapeyron slope indicates that the phase transition will occur in a greater depth, if the temperature is higher than the one given in Phase transition temperatures and in a smaller depth, if the temperature is smaller than the one given in Phase transition temperatures. For negative slopes the other way round. List must have the same number of entries as Phase transition depths. Units: \si{\pascal\per\kelvin}. 

(parameters:Material_20model:Visco_20Plastic:Phase_20transition_20depths)=
### __Parameter name:__ Phase transition depths
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of depths where phase transitions occur. Values must monotonically increase. Units: \si{\meter}. 

(parameters:Material_20model:Visco_20Plastic:Phase_20transition_20pressure_20widths)=
### __Parameter name:__ Phase transition pressure widths
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of widths for each phase transition, in terms of pressure. The phase functions are scaled with these values, leading to a jump between phases for a value of zero and a gradual transition for larger values. List must have the same number of entries as Phase transition pressures. Define transition by depth instead of pressure must be set to false to use this parameter. Units: \si{\pascal}. 

(parameters:Material_20model:Visco_20Plastic:Phase_20transition_20pressures)=
### __Parameter name:__ Phase transition pressures
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of pressures where phase transitions occur. Values must monotonically increase. Define transition by depth instead of pressure must be set to false to use this parameter. Units: \si{\pascal}. 

(parameters:Material_20model:Visco_20Plastic:Phase_20transition_20temperatures)=
### __Parameter name:__ Phase transition temperatures
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of temperatures where phase transitions occur. Higher or lower temperatures lead to phase transition occurring in smaller or greater depths than given in Phase transition depths, depending on the Clapeyron slope given in Phase transition Clapeyron slopes. List must have the same number of entries as Phase transition depths. Units: \si{\kelvin}. 

(parameters:Material_20model:Visco_20Plastic:Phase_20transition_20widths)=
### __Parameter name:__ Phase transition widths
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of widths for each phase transition, in terms of depth. The phase functions are scaled with these values, leading to a jump between phases for a value of zero and a gradual transition for larger values. List must have the same number of entries as Phase transition depths. Units: \si{\meter}. 

(parameters:Material_20model:Visco_20Plastic:Plastic_20damper_20viscosity)=
### __Parameter name:__ Plastic damper viscosity
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Viscosity of the damper that acts in parallel with the plastic viscosity to produce mesh-independent behavior at sufficient resolutions. Units: \si{\pascal\second} 

(parameters:Material_20model:Visco_20Plastic:Prefactor_20strain_20weakening_20factors)=
### __Parameter name:__ Prefactor strain weakening factors
**Default value:** 1. 

**Pattern:** [List of <[Double 0...1 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of viscous strain weakening factors for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Prefactors_20for_20Frank_20Kamenetskii)=
### __Parameter name:__ Prefactors for Frank Kamenetskii
**Default value:** 1.e21 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A viscosity prefactor for the viscosity approximation, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None 

(parameters:Material_20model:Visco_20Plastic:Prefactors_20for_20Peierls_20creep)=
### __Parameter name:__ Prefactors for Peierls creep
**Default value:** 1.4e-19 

**Pattern:** [Anything] 

**Documentation:** List of viscosity prefactors, $A$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\pascal}$^{-n_{\text{peierls}}}$ \si{\per\second} 

(parameters:Material_20model:Visco_20Plastic:Prefactors_20for_20diffusion_20creep)=
### __Parameter name:__ Prefactors for diffusion creep
**Default value:** 1.5e-15 

**Pattern:** [Anything] 

**Documentation:** List of viscosity prefactors, $A$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\per\pascal\meter}$^{m_{\text{diffusion}}}$\si{\per\second}. 

(parameters:Material_20model:Visco_20Plastic:Prefactors_20for_20dislocation_20creep)=
### __Parameter name:__ Prefactors for dislocation creep
**Default value:** 1.1e-16 

**Pattern:** [Anything] 

**Documentation:** List of viscosity prefactors, $A$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\pascal}$^{-n_{\text{dislocation}}}$ \si{\per\second}. 

(parameters:Material_20model:Visco_20Plastic:Reference_20strain_20rate)=
### __Parameter name:__ Reference strain rate
**Default value:** 1.0e-15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference strain rate for first time step. Units: \si{\per\second}. 

(parameters:Material_20model:Visco_20Plastic:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Visco_20Plastic:Reference_20viscosity)=
### __Parameter name:__ Reference viscosity
**Default value:** 1e22 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Reference viscosity for nondimensionalization. To understand how pressure scaling works, take a look at \cite{KHB12}. In particular, the value of this parameter would not affect the solution computed by \aspect{} if we could do arithmetic exactly; however, computers do arithmetic in finite precision, and consequently we need to scale quantities in ways so that their magnitudes are roughly the same. As explained in \cite{KHB12}, we scale the pressure during some computations (never visible by users) by a factor that involves a reference viscosity. This parameter describes this reference viscosity.

For problems with a constant viscosity, you will generally want to choose the reference viscosity equal to the actual viscosity. For problems with a variable viscosity, the reference viscosity should be a value that adequately represents the order of magnitude of the viscosities that appear, such as an average value or the value one would use to compute a Rayleigh number.

Units: \si{\pascal\second}. 

(parameters:Material_20model:Visco_20Plastic:Specific_20heats)=
### __Parameter name__: Specific heats
**Alias:** [Heat capacities](parameters:Material_20model:Visco_20Plastic:Heat_20capacities)

**Deprecation Status:** false

(parameters:Material_20model:Visco_20Plastic:Stabilization_20time_20scale_20factor)=
### __Parameter name:__ Stabilization time scale factor
**Default value:** 1. 

**Pattern:** [Double 1...MAX_DOUBLE (inclusive)] 

**Documentation:** A stabilization factor for the elastic stresses that influences how fast elastic stresses adjust to deformation. 1.0 is equivalent to no stabilization and may lead to oscillatory motion. Setting the factor to 2 avoids oscillations, but still enables an immediate elastic response. However, in complex models this can lead to problems of convergence, in which case the factor needs to be increased slightly. Setting the factor to infinity is equivalent to not applying elastic stresses at all. The factor is multiplied with the computational time step to create a time scale.  

(parameters:Material_20model:Visco_20Plastic:Start_20plasticity_20strain_20weakening_20intervals)=
### __Parameter name:__ Start plasticity strain weakening intervals
**Default value:** 0. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of strain weakening interval initial strains for the cohesion and friction angle parameters of the background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: None. 

(parameters:Material_20model:Visco_20Plastic:Start_20prefactor_20strain_20weakening_20intervals)=
### __Parameter name:__ Start prefactor strain weakening intervals
**Default value:** 0. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of strain weakening interval initial strains for the diffusion and dislocation prefactor parameters of the background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Strain_20healing_20mechanism)=
### __Parameter name:__ Strain healing mechanism
**Default value:** no healing 

**Pattern:** [Selection no healing|temperature dependent ] 

**Documentation:** Whether to apply strain healing to plastic yielding and viscosity terms, and if yes, which method to use. The following methods are available:

\item ``no healing'': No strain healing is applied. 

\item ``temperature dependent'': Purely temperature dependent strain healing applied to plastic yielding and viscosity terms, similar to the temperature-dependent Frank Kamenetskii formulation, computes strain healing as removing strain as a function of temperature, time, and a user-defined healing rate and prefactor as done in Fuchs and Becker, 2019, for mantle convection 

(parameters:Material_20model:Visco_20Plastic:Strain_20healing_20temperature_20dependent_20prefactor)=
### __Parameter name:__ Strain healing temperature dependent prefactor
**Default value:** 15. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor for temperature dependent strain healing. Units: None 

(parameters:Material_20model:Visco_20Plastic:Strain_20healing_20temperature_20dependent_20recovery_20rate)=
### __Parameter name:__ Strain healing temperature dependent recovery rate
**Default value:** 1.e-15 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Recovery rate prefactor for temperature dependent strain healing. Units: $1/s$ 

(parameters:Material_20model:Visco_20Plastic:Strain_20weakening_20mechanism)=
### __Parameter name:__ Strain weakening mechanism
**Default value:** default 

**Pattern:** [Selection none|finite strain tensor|total strain|plastic weakening with plastic strain only|plastic weakening with total strain only|plastic weakening with plastic strain and viscous weakening with viscous strain|viscous weakening with viscous strain only|default ] 

**Documentation:** Whether to apply strain weakening to viscosity, cohesion and internal angleof friction based on accumulated finite strain, and if yes, which method to use. The following methods are available:

\item ``none'': No strain weakening is applied. 

\item ``finite strain tensor'': The full finite strain tensor is tracked, and its second invariant is used to weaken both the plastic yield stress (specifically, the cohesion and friction angle) and the pre-yield viscosity that arises from diffusion and/or dislocation creep.

\item ``total strain'': The finite strain is approximated as the product of the second invariant of the strain rate in each time step and the time step size, and this quantity is integrated and tracked over time. It is used to weaken both the plastic yield stress (specifically, the cohesion and friction angle) and the pre-yield viscosity.

\item ``plastic weakening with plastic strain only'': The finite strain is approximated as the product of the second invariant of the strain ratein each time step and the time step size in regions where material is plastically yielding. This quantity is integrated and tracked over time, and used to weaken the cohesion and friction angle. The pre-yield viscosity is not weakened.

\item ``plastic weakening with total strain only'': The finite strain is approximated as the product of the second invariant of the strain rate in each time step and the time step size, and this quantity is integrated and tracked over time. It is used to weaken the plastic yield stress (specifically, the cohesion and internal friction angle). The pre-yield viscosity is not weakened.

\item ``plastic weakening with plastic strain and viscous weakening with viscous strain'': Both the finite strain accumulated by plastic deformation and by viscous deformation are computed separately (each approximated as the product of the second invariant of the corresponding strain rate in each time step and the time step size). The plastic strain is used to weaken the plastic yield stress (specifically, the cohesion and yield angle), and the viscous strain is used to weaken the pre-yield viscosity.

\item ``viscous weakening with viscous strain only'': The finite strain is approximated as the product of the second invariant of the strain rate in each time step and the time step size in regions where material is not plastically yielding. This quantity is integrated and tracked over time, and used to weaken the pre-yield viscosity. The cohesion and friction angle are not weakened.

\item ``default'': The default option has the same behavior as ``none'', but is there to make sure that the original parameters for specifying the strain weakening mechanism (``Use plastic/viscous strain weakening'') are still allowed, but to guarantee that one uses either the old parameter names or the new ones, never both.

If a compositional field named 'noninitial\_plastic\_strain' is included in the parameter file, this field will automatically be excluded from from volume fraction calculation and track the cumulative plastic strain with the initial plastic strain values removed. 

(parameters:Material_20model:Visco_20Plastic:Stress_20exponents_20for_20Peierls_20creep)=
### __Parameter name:__ Stress exponents for Peierls creep
**Default value:** 2.0 

**Pattern:** [Anything] 

**Documentation:** List of stress exponents, $n_{\text{peierls}}$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Stress_20exponents_20for_20diffusion_20creep)=
### __Parameter name:__ Stress exponents for diffusion creep
**Default value:** 1. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of stress exponents, $n_{\text{diffusion}}$, for background mantle and compositional fields, for a total of N+1 values, where N is the number of compositional fields. The stress exponent for diffusion creep is almost always equal to one. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Stress_20exponents_20for_20dislocation_20creep)=
### __Parameter name:__ Stress exponents for dislocation creep
**Default value:** 3.5 

**Pattern:** [Anything] 

**Documentation:** List of stress exponents, $n_{\text{dislocation}}$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: None. 

(parameters:Material_20model:Visco_20Plastic:Stress_20limiter_20exponents)=
### __Parameter name:__ Stress limiter exponents
**Default value:** 1.0 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of stress limiter exponents, $n_{\text{lim}}$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. Units: none. 

(parameters:Material_20model:Visco_20Plastic:Thermal_20conductivities)=
### __Parameter name:__ Thermal conductivities
**Default value:** 3.0 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of thermal conductivities, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Visco_20Plastic:Thermal_20diffusivities)=
### __Parameter name:__ Thermal diffusivities
**Default value:** 0.8e-6 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of thermal diffusivities, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value.  Units: \si{\meter\squared\per\second}. 

(parameters:Material_20model:Visco_20Plastic:Thermal_20expansivities)=
### __Parameter name:__ Thermal expansivities
**Default value:** 0.000035 

**Pattern:** [Anything] 

**Documentation:** List of thermal expansivities for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Visco_20Plastic:Use_20adiabatic_20pressure_20in_20creep_20viscosity)=
### __Parameter name:__ Use adiabatic pressure in creep viscosity
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use the adiabatic pressure instead of the full pressure (default) when calculating creep (diffusion, dislocation, and peierls) viscosity. This may be helpful in models where the full pressure has an unusually large negative value arising from large negative dynamic pressure, resulting in solver convergence issue and in some cases a viscosity of zero. 

(parameters:Material_20model:Visco_20Plastic:Use_20fixed_20elastic_20time_20step)=
### __Parameter name:__ Use fixed elastic time step
**Default value:** unspecified 

**Pattern:** [Selection true|false|unspecified ] 

**Documentation:** Select whether the material time scale in the viscoelastic constitutive relationship uses the regular numerical time step or a separate fixed elastic time step throughout the model run. The fixed elastic time step is always used during the initial time step. If a fixed elastic time step is used throughout the model run, a stress averaging scheme can be applied to account for differences with the numerical time step. An alternative approach is to limit the maximum time step size so that it is equal to the elastic time step. The default value of this parameter is 'unspecified', which throws an exception during runtime. In order for the model to run the user must select 'true' or 'false'. 

(parameters:Material_20model:Visco_20Plastic:Use_20plastic_20damper)=
### __Parameter name:__ Use plastic damper
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a plastic damper when computing the Drucker-Prager plastic viscosity. The damper acts to stabilize the plastic shear band width and remove associated mesh-dependent behavior at sufficient resolutions. 

(parameters:Material_20model:Visco_20Plastic:Use_20stress_20averaging)=
### __Parameter name:__ Use stress averaging
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to apply a stress averaging scheme to account for differences between the fixed elastic time step and numerical time step.  

(parameters:Material_20model:Visco_20Plastic:Viscosity_20averaging_20scheme)=
### __Parameter name:__ Viscosity averaging scheme
**Default value:** harmonic 

**Pattern:** [Selection arithmetic|harmonic|geometric|maximum composition ] 

**Documentation:** When more than one compositional field is present at a point with different viscosities, we need to come up with an average viscosity at that point.  Select a weighted harmonic, arithmetic, geometric, or maximum composition. 

(parameters:Material_20model:Visco_20Plastic:Viscosity_20ratios_20for_20Frank_20Kamenetskii)=
### __Parameter name:__ Viscosity ratios for Frank Kamenetskii
**Default value:** 15. 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** An adjusted viscosity ratio, $E$, for the viscosity approximation, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: None 

(parameters:Material_20model:Visco_20Plastic:Viscous_20flow_20law)=
### __Parameter name:__ Viscous flow law
**Default value:** composite 

**Pattern:** [Selection diffusion|dislocation|frank kamenetskii|composite ] 

**Documentation:** Select what type of viscosity law to use between diffusion, dislocation, frank kamenetskii, and composite options. Soon there will be an option to select a specific flow law for each assigned composition  

(parameters:Material_20model:Visco_20Plastic:Yield_20mechanism)=
### __Parameter name:__ Yield mechanism
**Default value:** drucker 

**Pattern:** [Selection drucker|limiter ] 

**Documentation:** Select what type of yield mechanism to use between Drucker Prager and stress limiter options. 

(parameters:Material_20model:Visco_20Plastic:Friction_20function)=
## **Parameters in section** Material model/Visco Plastic/Friction function
(parameters:Material_20model:Visco_20Plastic:Friction_20function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian 

**Pattern:** [Selection cartesian|spherical|depth ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point. 

(parameters:Material_20model:Visco_20Plastic:Friction_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Material_20model:Visco_20Plastic:Friction_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Material_20model:Visco_20Plastic:Friction_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Material_20model:Viscoelastic)=
## **Parameters in section** Material model/Viscoelastic
(parameters:Material_20model:Viscoelastic:Densities)=
### __Parameter name:__ Densities
**Default value:** 3300. 

**Pattern:** [Anything] 

**Documentation:** List of densities for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Material_20model:Viscoelastic:Elastic_20damper_20viscosity)=
### __Parameter name:__ Elastic damper viscosity
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Viscosity of a viscous damper that acts in parallel with the elastic element to stabilize behavior. Units: \si{\pascal\second} 

(parameters:Material_20model:Viscoelastic:Elastic_20shear_20moduli)=
### __Parameter name:__ Elastic shear moduli
**Default value:** 75.0e9 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of elastic shear moduli, $G$, for background material and compositional fields, for a total of N+1 values, where N is the number of compositional fields. The default value of 75 GPa is representative of mantle rocks. Units: Pa. 

(parameters:Material_20model:Viscoelastic:Fixed_20elastic_20time_20step)=
### __Parameter name:__ Fixed elastic time step
**Default value:** 1.e3 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The fixed elastic time step $dte$. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Material_20model:Viscoelastic:Heat_20capacities)=
### __Parameter name:__ Heat capacities
**Default value:** 1250. 

**Pattern:** [Anything] 

**Documentation:** List of specific heats $C_p$ for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\joule\per\kelvin\per\kilogram}. 

(parameters:Material_20model:Viscoelastic:Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 293. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The reference temperature $T_0$. Units: \si{\kelvin}. 

(parameters:Material_20model:Viscoelastic:Specific_20heats)=
### __Parameter name__: Specific heats
**Alias:** [Heat capacities](parameters:Material_20model:Viscoelastic:Heat_20capacities)

**Deprecation Status:** false

(parameters:Material_20model:Viscoelastic:Stabilization_20time_20scale_20factor)=
### __Parameter name:__ Stabilization time scale factor
**Default value:** 1. 

**Pattern:** [Double 1...MAX_DOUBLE (inclusive)] 

**Documentation:** A stabilization factor for the elastic stresses that influences how fast elastic stresses adjust to deformation. 1.0 is equivalent to no stabilization and may lead to oscillatory motion. Setting the factor to 2 avoids oscillations, but still enables an immediate elastic response. However, in complex models this can lead to problems of convergence, in which case the factor needs to be increased slightly. Setting the factor to infinity is equivalent to not applying elastic stresses at all. The factor is multiplied with the computational time step to create a time scale.  

(parameters:Material_20model:Viscoelastic:Thermal_20conductivities)=
### __Parameter name:__ Thermal conductivities
**Default value:** 4.7 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of thermal conductivities for background mantle and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\watt\per\meter\per\kelvin}. 

(parameters:Material_20model:Viscoelastic:Thermal_20expansivities)=
### __Parameter name:__ Thermal expansivities
**Default value:** 0.000035 

**Pattern:** [Anything] 

**Documentation:** List of thermal expansivities for background mantle and compositional fields,for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. If only one value is given, then all use the same value. Units: \si{\per\kelvin}. 

(parameters:Material_20model:Viscoelastic:Use_20fixed_20elastic_20time_20step)=
### __Parameter name:__ Use fixed elastic time step
**Default value:** unspecified 

**Pattern:** [Selection true|false|unspecified ] 

**Documentation:** Select whether the material time scale in the viscoelastic constitutive relationship uses the regular numerical time step or a separate fixed elastic time step throughout the model run. The fixed elastic time step is always used during the initial time step. If a fixed elastic time step is used throughout the model run, a stress averaging scheme can be applied to account for differences with the numerical time step. An alternative approach is to limit the maximum time step size so that it is equal to the elastic time step. The default value of this parameter is 'unspecified', which throws an exception during runtime. In order for the model to run the user must select 'true' or 'false'. 

(parameters:Material_20model:Viscoelastic:Use_20stress_20averaging)=
### __Parameter name:__ Use stress averaging
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to apply a stress averaging scheme to account for differences between the fixed elastic time step and numerical time step.  

(parameters:Material_20model:Viscoelastic:Viscosities)=
### __Parameter name:__ Viscosities
**Default value:** 1.e21 

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of viscosities for background mantle and compositional fields, for a total of N+1 values, where N is the number of compositional fields. If only one value is given, then all use the same value. Units: \si{\pascal\second}. 

(parameters:Material_20model:Viscoelastic:Viscosity_20averaging_20scheme)=
### __Parameter name:__ Viscosity averaging scheme
**Default value:** harmonic 

**Pattern:** [Selection arithmetic|harmonic|geometric|maximum composition  ] 

**Documentation:** When more than one compositional field is present at a point with different viscosities, we need to come up with an average viscosity at that point.  Select a weighted harmonic, arithmetic, geometric, or maximum composition. 

(parameters:Melt_20settings)=
## **Parameters in section** Melt settings
(parameters:Melt_20settings:Average_20melt_20velocity)=
### __Parameter name:__ Average melt velocity
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to cell-wise average the material properties that are used to compute the melt velocity or not. The melt velocity is computed as the sum of the solid velocity and the phase separation flux $ - K_D / \phi (\nabla p_f - \rho_f \mathbf g)$. If this parameter is set to true, $K_D$ and $\phi$ will be averaged cell-wise in the computation of the phase separation flux. This is useful because in some models the melt velocity can have spikes close to the interface between regions of melt and no melt, as both $K_D$ and $\phi$ go to zero for vanishing melt fraction. As the melt velocity is used for computing the time step size, and in models that use heat transport by melt or shear heating of melt, setting this parameter to true can speed up the model and make it mode stable. In computations where accuracy and convergence behavior of the melt velocity is important (like in benchmark cases with an analytical solution), this parameter should probably be set to 'false'. 

(parameters:Melt_20settings:Heat_20advection_20by_20melt)=
### __Parameter name:__ Heat advection by melt
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a porosity weighted average of the melt and solid velocity to advect heat in the temperature equation or not. If this is set to true, additional terms are assembled on the left-hand side of the temperature advection equation. Only used if Include melt transport is true. If this is set to false, only the solid velocity is used (as in models without melt migration). 

(parameters:Melt_20settings:Include_20melt_20transport)=
### __Parameter name:__ Include melt transport
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to include the transport of melt into the model or not. If this is set to true, two additional pressures (the fluid pressure and the compaction pressure) will be added to the finite element. Including melt transport in the simulation also requires that there is one compositional field that has the name `porosity'. This field will be used for computing the additional pressures and the melt velocity, and has a different advection equation than other compositional fields, as it is effectively advected with the melt velocity. 

(parameters:Melt_20settings:Melt_20scaling_20factor_20threshold)=
### __Parameter name:__ Melt scaling factor threshold
**Default value:** 1e-7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** The factor by how much the Darcy coefficient K\_D in a cell can be smaller than the reference Darcy coefficient for this cell still to be considered a melt cell (for which the melt transport equations are solved). For smaller Darcy coefficients, the Stokes equations (without melt) are solved instead. Only used if ``Include melt transport'' is true.  

(parameters:Melt_20settings:Use_20discontinuous_20compaction_20pressure)=
### __Parameter name:__ Use discontinuous compaction pressure
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to use a discontinuous element for the compaction pressure or not. From our preliminary tests, continuous elements seem to work better in models where the porosity is > 0 everywhere in the domain, and discontinuous elements work better in models where in parts of the domain the porosity = 0. 

(parameters:Mesh_20deformation)=
## **Parameters in section** Mesh deformation
(parameters:Mesh_20deformation:Additional_20tangential_20mesh_20velocity_20boundary_20indicators)=
### __Parameter name:__ Additional tangential mesh velocity boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries where there the mesh is allowed to move tangential to the boundary. All tangential mesh movements along those boundaries that have tangential material velocity boundary conditions are allowed by default, this parameters allows to generate mesh movements along other boundaries that are open, or have prescribed material velocities or tractions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model. 

(parameters:Mesh_20deformation:Mesh_20deformation_20boundary_20indicators)=
### __Parameter name:__ Mesh deformation boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries where there the mesh is allowed to move according to the specified mesh deformation objects. 

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model. 

The format is id1: object1 \& object2, id2: object3 \& object2, where objects are one of `ascii data': Implementation of a model in which the initial mesh deformation (initial topography) is derived from a file containing data in ascii format. The following geometry models are currently supported: box, chunk, spherical. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `Topography [m]' in a 2d model and  `x', `y', `Topography [m]' in a 3d model, which means that there has to be a single column containing the topography. Note that the data in the input file needs to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle in radians  and `y' by the polar angle in radians measured positive from the north pole. The grid will be assumed to be a longitude-colatitude grid. Note that the order of spherical coordinates is `phi', `theta' and not `theta', `phi', since this allows for dimension independent expressions.

`boundary function': A plugin, which prescribes the surface mesh to deform according to an analytically prescribed function. Note that the function prescribes a deformation velocity, i.e. the return value of this plugin is later multiplied by the time step length to compute the displacement increment in this time step. Although the function's time variable is interpreted as years when Use years in output instead of seconds is set to true, the boundary deformation velocity should still be given in m/s. The format of the functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`diffusion': A plugin that computes the deformation of surface vertices according to the solution of the hillslope diffusion problem. Specifically, at the end of each timestep, or after a specific number of timesteps, this plugin solves the following equation: \begin{align*}  \frac{\partial h}{\partial t} = \kappa \left( \frac{\partial^{2} h}{\partial x^{2}} + \frac{\partial^{2} h}{\partial y^{2}} \right), \end{align*} where $\kappa$ is the hillslope diffusion coefficient (diffusivity), and $h(x,y)$ the height of a point along the top boundary with respect to the surface of the unperturbed domain. 

Using this definition, the plugin then solves for one time step, i.e., using as initial condition $h(t_{n-1})$ the current surface elevation, and computing $h(t_n)$ from it by solving the equation above over the time interval $t_n-t_{n-1}$. From this, one can then compute a surface velocity $v = \frac{h(t_n)-h(t_{n-1})}{t_n-t_{n-1}}$. 

This surface velocity is used to deform the surface and as a boundary condition for solving the Laplace equation to determine the mesh velocity in the domain interior. Diffusion can be applied every timestep, mimicking surface processes of erosion and deposition, or at a user-defined timestep interval to purely smooth the surface topography to avoid too great a distortion of mesh elements when a free surface is also used.

`free surface': A plugin that computes the deformation of surface vertices according to the solution of the flow problem. In particular this means if the surface of the domain is left open to flow, this flow will carry the mesh with it. The implementation was described in \cite{rose_freesurface}, with the stabilization of the free surface originally described in \cite{KMM2010}. 

(parameters:Mesh_20deformation:Ascii_20data_20model)=
## **Parameters in section** Mesh deformation/Ascii data model
(parameters:Mesh_20deformation:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Mesh_20deformation:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_3d_%s.0.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Mesh_20deformation:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Mesh_20deformation:Boundary_20function)=
## **Parameters in section** Mesh deformation/Boundary function
(parameters:Mesh_20deformation:Boundary_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Mesh_20deformation:Boundary_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Mesh_20deformation:Boundary_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Mesh_20deformation:Diffusion)=
## **Parameters in section** Mesh deformation/Diffusion
(parameters:Mesh_20deformation:Diffusion:Hillslope_20transport_20coefficient)=
### __Parameter name:__ Hillslope transport coefficient
**Default value:** 1e-6 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The hillslope transport coefficient $\kappa$ used to diffuse the free surface, either as a  stabilization step or to mimic erosional and depositional processes. Units: $\si{m^2/s}$.  

(parameters:Mesh_20deformation:Diffusion:Time_20steps_20between_20diffusion)=
### __Parameter name:__ Time steps between diffusion
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of time steps between each application of diffusion. 

(parameters:Mesh_20deformation:Free_20surface)=
## **Parameters in section** Mesh deformation/Free surface
(parameters:Mesh_20deformation:Free_20surface:Free_20surface_20stabilization_20theta)=
### __Parameter name:__ Free surface stabilization theta
**Default value:** 0.5 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** Theta parameter described in \cite{KMM2010}. An unstabilized free surface can overshoot its equilibrium position quite easily and generate unphysical results.  One solution is to use a quasi-implicit correction term to the forces near the free surface.  This parameter describes how much the free surface is stabilized with this term, where zero is no stabilization, and one is fully implicit. 

(parameters:Mesh_20deformation:Free_20surface:Surface_20velocity_20projection)=
### __Parameter name:__ Surface velocity projection
**Default value:** normal 

**Pattern:** [Selection normal|vertical ] 

**Documentation:** After each time step the free surface must be advected in the direction of the velocity field. Mass conservation requires that the mesh velocity is in the normal direction of the surface. However, for steep topography or large curvature, advection in the normal direction can become ill-conditioned, and instabilities in the mesh can form. Projection of the mesh velocity onto the local vertical direction can preserve the mesh quality better, but at the cost of slightly poorer mass conservation of the domain. 

(parameters:Mesh_20refinement)=
## **Parameters in section** Mesh refinement
(parameters:Mesh_20refinement:Adapt_20by_20fraction_20of_20cells)=
### __Parameter name:__ Adapt by fraction of cells
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Use fraction of the total number of cells instead of fraction of the total error as the limit for refinement and coarsening. 

(parameters:Mesh_20refinement:Additional_20refinement_20times)=
### __Parameter name:__ Additional refinement times
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of times so that if the end time of a time step is beyond this time, an additional round of mesh refinement is triggered. This is mostly useful to make sure we can get through the initial transient phase of a simulation on a relatively coarse mesh, and then refine again when we are in a time range that we are interested in and where we would like to use a finer mesh. Units: Each element of the list has units years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Mesh_20refinement:Coarsening_20fraction)=
### __Parameter name:__ Coarsening fraction
**Default value:** 0.05 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** Cells are sorted from largest to smallest by their total error (determined by the Strategy). Then the cells with the smallest error (bottom of this sorted list) that account for the given fraction of the error are coarsened. 

(parameters:Mesh_20refinement:Initial_20adaptive_20refinement)=
### __Parameter name:__ Initial adaptive refinement
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of adaptive refinement steps performed after initial global refinement but while still within the first time step. These refinement steps (n) are added to the value for initial global refinement (m) so that the final mesh has cells that are at most on refinement level $n+m$. 

(parameters:Mesh_20refinement:Initial_20global_20refinement)=
### __Parameter name:__ Initial global refinement
**Default value:** 2 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of global refinement steps performed on the initial coarse mesh, before the problem is first solved there.

Note that it is possible to supply conflicting refinement and coarsening settings, such as an 'Initial global refinement' of 4 and a 'Maximum refinement function' strategy that limits the refinement locally to 2. In this case, the tagging strategies such as the 'Maximum refinement function' will remove refinement flags in each initial global refinement step, such that the resulting mesh is not necessarily uniform or of the level given by the 'Initial global refinement' parameter. 

(parameters:Mesh_20refinement:Minimum_20refinement_20level)=
### __Parameter name:__ Minimum refinement level
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The minimum refinement level each cell should have, and that can not be exceeded by coarsening. Should not be higher than the 'Initial global refinement' parameter. 

(parameters:Mesh_20refinement:Normalize_20individual_20refinement_20criteria)=
### __Parameter name:__ Normalize individual refinement criteria
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** If multiple refinement criteria are specified in the ``Strategy'' parameter, then they need to be combined somehow to form the final refinement indicators. This is done using the method described by the ``Refinement criteria merge operation'' parameter which can either operate on the raw refinement indicators returned by each strategy (i.e., dimensional quantities) or using normalized values where the indicators of each strategy are first normalized to the interval $[0,1]$ (which also makes them non-dimensional). This parameter determines whether this normalization will happen. 

(parameters:Mesh_20refinement:Refinement_20criteria_20merge_20operation)=
### __Parameter name:__ Refinement criteria merge operation
**Default value:** max 

**Pattern:** [Selection plus|max ] 

**Documentation:** If multiple mesh refinement criteria are computed for each cell (by passing a list of more than element to the \texttt{Strategy} parameter in this section of the input file) then one will have to decide which criteria should win when deciding which cells to refine. The operation that determines how to combine these competing criteria is the one that is selected here. The options are:

\begin{itemize}
\item \texttt{plus}: Add the various error indicators together and refine those cells on which the sum of indicators is largest.
\item \texttt{max}: Take the maximum of the various error indicators and refine those cells on which the maximal indicators is largest.
\end{itemize}The refinement indicators computed by each strategy are modified by the ``Normalize individual refinement criteria'' and ``Refinement criteria scale factors'' parameters. 

(parameters:Mesh_20refinement:Refinement_20criteria_20scaling_20factors)=
### __Parameter name:__ Refinement criteria scaling factors
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of scaling factors by which every individual refinement criterion will be multiplied by. If only a single refinement criterion is selected (using the ``Strategy'' parameter, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other. 

If ``Normalize individual refinement criteria'' is set to true, then the criteria will first be normalized to the interval $[0,1]$ and then multiplied by the factors specified here. You will likely want to choose the factors to be not too far from 1 in that case, say between 1 and 10, to avoid essentially disabling those criteria with small weights. On the other hand, if the criteria are not normalized to $[0,1]$ using the parameter mentioned above, then the factors you specify here need to take into account the relative numerical size of refinement indicators (which in that case carry physical units).

You can experimentally play with these scaling factors by choosing to output the refinement indicators into the graphical output of a run.

If the list of indicators given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are indicators chosen in the ``Strategy'' parameter. 

(parameters:Mesh_20refinement:Refinement_20fraction)=
### __Parameter name:__ Refinement fraction
**Default value:** 0.3 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** Cells are sorted from largest to smallest by their total error (determined by the Strategy). Then the cells with the largest error (top of this sorted list) that account for given fraction of the error are refined. 

(parameters:Mesh_20refinement:Run_20postprocessors_20on_20initial_20refinement)=
### __Parameter name:__ Run postprocessors on initial refinement
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether or not the postprocessors should be executed after each of the initial adaptive refinement cycles that are run at the start of the simulation. This is useful for plotting/analyzing how the mesh refinement parameters are working for a particular model. 

(parameters:Mesh_20refinement:Skip_20setup_20initial_20conditions_20on_20initial_20refinement)=
### __Parameter name:__ Skip setup initial conditions on initial refinement
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether or not the initial conditions should be set up during the adaptive refinement cycles that are run at the start of the simulation. 

(parameters:Mesh_20refinement:Skip_20solvers_20on_20initial_20refinement)=
### __Parameter name:__ Skip solvers on initial refinement
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether or not solvers should be executed during the initial adaptive refinement cycles that are run at the start of the simulation. 

(parameters:Mesh_20refinement:Strategy)=
### __Parameter name:__ Strategy
**Default value:** thermal energy density 

**Pattern:** [MultipleSelection artificial viscosity|boundary|compaction length|composition|composition approximate gradient|composition gradient|composition threshold|density|isosurfaces|maximum refinement function|minimum refinement function|nonadiabatic temperature|particle density|slope|strain rate|temperature|thermal energy density|topography|velocity|viscosity|volume of fluid interface ] 

**Documentation:** A comma separated list of mesh refinement criteria that will be run whenever mesh refinement is required. The results of each of these criteria, i.e., the refinement indicators they produce for all the cells of the mesh will then be normalized to a range between zero and one and the results of different criteria will then be merged through the operation selected in this section.

The following criteria are available:

`artificial viscosity': A mesh refinement criterion that computes refinement indicators from the artificial viscosity of the temperature or compositional fields based on user specified weights.

`boundary': A class that implements a mesh refinement criterion which always flags all cells on specified boundaries for refinement. This is useful to provide high accuracy for processes at or close to the edge of the model domain.

To use this refinement criterion, you may want to combine it with other refinement criteria, setting the 'Normalize individual refinement criteria' flag and using the `max' setting for 'Refinement criteria merge operation'.

`compaction length': A mesh refinement criterion for models with melt transport that computes refinement indicators based on the compaction length, defined as $\delta = \sqrt{\frac{(\xi + 4 \eta/3) k}{\eta_f}}$. $\xi$ is the bulk viscosity, $\eta$ is the shear viscosity, $k$ is the permeability and $\eta_f$ is the melt viscosity. If the cell width or height exceeds a multiple (which is specified as an input parameter) of this compaction length, the cell is marked for refinement.

`composition': A mesh refinement criterion that computes refinement indicators from the compositional fields. If there is more than one compositional field, then it simply takes the sum of the indicators computed from each of the compositional field.

The way these indicators are computed is by evaluating the `Kelly error indicator' on each compositional field. This error indicator takes the finite element approximation of the compositional field and uses it to compute an approximation of the second derivatives of the composition for each cell. This approximation is then multiplied by an appropriate power of the cell's diameter to yield an indicator for how large the error is likely going to be on this cell. This construction rests on the observation that for many partial differential equations, the error on each cell is proportional to some power of the cell's diameter times the second derivatives of the solution on that cell.

For complex equations such as those we solve here, this observation may not be strictly true in the mathematical sense, but it often yields meshes that are surprisingly good.

`composition approximate gradient': A mesh refinement criterion that computes refinement indicators from the gradients of compositional fields. If there is more than one compositional field, then it simply takes the sum of the indicators times a user-specified weight for each field.

In contrast to the `composition gradient' refinement criterion, the current criterion does not compute the gradient at quadrature points on each cell, but by a finite difference approximation between the centers of cells. Consequently, it also works if the compositional fields are computed using discontinuous finite elements.

`composition gradient': A mesh refinement criterion that computes refinement indicators from the gradients of compositional fields. If there is more than one compositional field, then it simply takes the sum of the indicators times a user-specified weight for each field.

This refinement criterion computes the gradient of the compositional field at quadrature points on each cell, and then averages them in some way to obtain a refinement indicator for each cell. This will give a reasonable approximation of the true gradient of the compositional field if you are using a continuous finite element.

On the other hand, for discontinuous finite elements (see the `Use discontinuous composition discretization' parameter in the `Discretization' section), the gradient at quadrature points does not include the contribution of jumps in the compositional field between cells, and consequently will not be an accurate approximation of the true gradient. As an extreme example, consider the case of using piecewise constant finite elements for compositional fields; in that case, the gradient of the solution at quadrature points inside each cell will always be exactly zero, even if the finite element solution is different from each cell to the next. Consequently, the current refinement criterion will likely not be useful in this situation. That said, the `composition approximate gradient' refinement criterion exists for exactly this purpose.

`composition threshold': A mesh refinement criterion that computes refinement indicators from the compositional fields. If any field exceeds the threshold given in the input file, the cell is marked for refinement.

`density': A mesh refinement criterion that computes refinement indicators from a field that describes the spatial variability of the density, $\rho$. Because this quantity may not be a continuous function ($\rho$ and $C_p$ may be discontinuous functions along discontinuities in the medium, for example due to phase changes), we approximate the gradient of this quantity to refine the mesh. The error indicator defined here takes the magnitude of the approximate gradient and scales it by $h_K^{1+d/2}$ where $h_K$ is the diameter of each cell and $d$ is the dimension. This scaling ensures that the error indicators converge to zero as $h_K\rightarrow 0$ even if the energy density is discontinuous, since the gradient of a discontinuous function grows like $1/h_K$.

`isosurfaces': A mesh refinement criterion that computes coarsening and refinement indicators between two isosurfaces of specific field entries (e.g. temperature, composition).

The way these indicators are derived between pairs of isosurfaces is by checking whether the solutions of specific fields are within the ranges of the isosurface values given. If these conditions hold, then coarsening and refinement indicators are set such that the mesh refinement levels lies within the range of levels given. Usage of this plugin allows the user to put a conditional minimum and maximum refinement function onto fields that they are interested in.

For now, only temperature and compositional fields are allowed as field entries. The key words could be 'Temperature' or one of the names of the compositional fields which are either specified by user or set up as C\_0, C\_1, etc.

Usage: A list of isosurfaces separated by semi-colons (;). Each isosurface entry consists of multiple entries separated by a comma. The first two entries indicate the minimum and maximum refinement levels respectively. The entries after the first two describe the fields the isosurface applies to, followed by a colon (:), which again is followed by the minimum and maximum field values separated by a bar (|). An example for two isosurface entries is '0, 2, Temperature: 300 | 600; 2, 2, C\_1: 0.5 | 1'. If both isoterm entries are triggered at the same location and the current refinement level is 1, it means that the first isoline will not set any flag and the second isoline will set a refinement flag. This means the cell will be refined. If both the coarsening and refinement flags are set, preference is given to refinement. 

The minimum and maximum refinement levels per isosurface can be provided in absolute values relative to the global minimum and maximum refinement. This is done with the 'min' and 'max' key words. For example: 'set Isosurfaces = max-2,  max,    Temperature: 0 | 600 ; min + 1,min+2, Temperature: 1600 | 3000,   C\_2 : 0.0 | 0.5'.

`maximum refinement function': A mesh refinement criterion that ensures a maximum refinement level described by an explicit formula with the depth or position as argument. Which coordinate representation is used is determined by an input parameter. Whatever the coordinate system chosen, the function you provide in the input file will by default depend on variables `x', `y' and `z' (if in 3d). However, the meaning of these symbols depends on the coordinate system. In the Cartesian coordinate system, they simply refer to their natural meaning. If you have selected `depth' for the coordinate system, then `x' refers to the depth variable and `y' and `z' will simply always be zero. If you have selected a spherical coordinate system, then `x' will refer to the radial distance of the point to the origin, `y' to the azimuth angle and `z' to the polar angle measured positive from the north pole. Note that the order of spherical coordinates is r,phi,theta and not r,theta,phi, since this allows for dimension independent expressions. Each coordinate system also includes a final `t' variable which represents the model time, evaluated in years if the 'Use years in output instead of seconds' parameter is set, otherwise evaluated in seconds. After evaluating the function, its values are rounded to the nearest integer.

The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`minimum refinement function': A mesh refinement criterion that ensures a minimum refinement level described by an explicit formula with the depth or position as argument. Which coordinate representation is used is determined by an input parameter. Whatever the coordinate system chosen, the function you provide in the input file will by default depend on variables `x', `y' and `z' (if in 3d). However, the meaning of these symbols depends on the coordinate system. In the Cartesian coordinate system, they simply refer to their natural meaning. If you have selected `depth' for the coordinate system, then `x' refers to the depth variable and `y' and `z' will simply always be zero. If you have selected a spherical coordinate system, then `x' will refer to the radial distance of the point to the origin, `y' to the azimuth angle and `z' to the polar angle measured positive from the north pole. Note that the order of spherical coordinates is r,phi,theta and not r,theta,phi, since this allows for dimension independent expressions. Each coordinate system also includes a final `t' variable which represents the model time, evaluated in years if the 'Use years in output instead of seconds' parameter is set, otherwise evaluated in seconds. After evaluating the function, its values are rounded to the nearest integer.

The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`nonadiabatic temperature': A mesh refinement criterion that computes refinement indicators from the excess temperature(difference between temperature and adiabatic temperature.

`particle density': A mesh refinement criterion that computes refinement indicators based on the density of particles. In practice this plugin equilibrates the number of particles per cell, leading to fine cells in high particle density regions and coarse cells in low particle density regions. This plugin is mostly useful for models with inhomogeneous particle density, e.g. when tracking an initial interface with a high particle density, or when the spatial particle density denotes the region of interest. Additionally, this plugin tends to balance the computational load between processes in parallel computations, because the particle and mesh density is more aligned.

`slope': A class that implements a mesh refinement criterion intended for use with deforming mesh boundaries, like the free surface. It calculates a local slope based on the angle between the surface normal and the local gravity vector. Cells with larger angles are marked for refinement.

To use this refinement criterion, you may want to combine it with other refinement criteria, setting the 'Normalize individual refinement criteria' flag and using the `max' setting for 'Refinement criteria merge operation'.

`strain rate': A mesh refinement criterion that computes the refinement indicators equal to the strain rate norm computed at the center of the elements.

`temperature': A mesh refinement criterion that computes refinement indicators from the temperature field.

The way these indicators are computed is by evaluating the `Kelly error indicator' on the temperature field. This error indicator takes the finite element approximation of the temperature field and uses it to compute an approximation of the second derivatives of the temperature for each cell. This approximation is then multiplied by an appropriate power of the cell's diameter to yield an indicator for how large the error is likely going to be on this cell. This construction rests on the observation that for many partial differential equations, the error on each cell is proportional to some power of the cell's diameter times the second derivatives of the solution on that cell.

For complex equations such as those we solve here, this observation may not be strictly true in the mathematical sense, but it often yields meshes that are surprisingly good.

`thermal energy density': A mesh refinement criterion that computes refinement indicators from a field that describes the spatial variability of the thermal energy density, $\rho C_p T$. Because this quantity may not be a continuous function ($\rho$ and $C_p$ may be discontinuous functions along discontinuities in the medium, for example due to phase changes), we approximate the gradient of this quantity to refine the mesh. The error indicator defined here takes the magnitude of the approximate gradient and scales it by $h_K^{1.5}$ where $h_K$ is the diameter of each cell. This scaling ensures that the error indicators converge to zero as $h_K\rightarrow 0$ even if the energy density is discontinuous, since the gradient of a discontinuous function grows like $1/h_K$.

`topography': A class that implements a mesh refinement criterion, which always flags all cells in the uppermost layer for refinement. This is useful to provide high accuracy for processes at or close to the surface.

To use this refinement criterion, you may want to combine it with other refinement criteria, setting the 'Normalize individual refinement criteria' flag and using the `max' setting for 'Refinement criteria merge operation'.

`velocity': A mesh refinement criterion that computes refinement indicators from the velocity field.

The way these indicators are computed is by evaluating the `Kelly error indicator' on the velocity field. This error indicator takes the finite element approximation of the velocity field and uses it to compute an approximation of the second derivatives of the velocity for each cell. This approximation is then multiplied by an appropriate power of the cell's diameter to yield an indicator for how large the error is likely going to be on this cell. This construction rests on the observation that for many partial differential equations, the error on each cell is proportional to some power of the cell's diameter times the second derivatives of the solution on that cell.

For complex equations such as those we solve here, this observation may not be strictly true in the mathematical sense, but it often yields meshes that are surprisingly good.

`viscosity': A mesh refinement criterion that computes refinement indicators from a field that describes the spatial variability of the logarithm of the viscosity, $\log\eta$. (We choose the logarithm of the viscosity because it can vary by orders of magnitude.)Because this quantity may not be a continuous function ($\eta$ may be a discontinuous function along discontinuities in the medium, for example due to phase changes), we approximate the gradient of this quantity to refine the mesh. The error indicator defined here takes the magnitude of the approximate gradient and scales it by $h_K^{1+d/2}$ where $h_K$ is the diameter of each cell and $d$ is the dimension. This scaling ensures that the error indicators converge to zero as $h_K\rightarrow 0$ even if the energy density is discontinuous, since the gradient of a discontinuous function grows like $1/h_K$.

`volume of fluid interface': A class that implements a mesh refinement criterion, which ensures a minimum level of refinement near the volume of fluid interface boundary. 

(parameters:Mesh_20refinement:Time_20steps_20between_20mesh_20refinement)=
### __Parameter name:__ Time steps between mesh refinement
**Default value:** 10 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of time steps after which the mesh is to be adapted again based on computed error indicators. If 0 then the mesh will never be changed. 

(parameters:Mesh_20refinement:Artificial_20viscosity)=
## **Parameters in section** Mesh refinement/Artificial viscosity
(parameters:Mesh_20refinement:Artificial_20viscosity:Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of scaling factors by which every individual compositional field will be multiplied. These factors are used to weigh the various indicators relative to each other and to the temperature. 

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to 0. If the list is not empty then it needs to have as many entries as there are compositional fields. 

(parameters:Mesh_20refinement:Artificial_20viscosity:Temperature_20scaling_20factor)=
### __Parameter name:__ Temperature scaling factor
**Default value:** 0.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A scaling factor for the artificial viscosity  of the temperature equation. Use 0.0 to disable. 

(parameters:Mesh_20refinement:Boundary)=
## **Parameters in section** Mesh refinement/Boundary
(parameters:Mesh_20refinement:Boundary:Boundary_20refinement_20indicators)=
### __Parameter name:__ Boundary refinement indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries where there should be mesh refinement.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model. 

(parameters:Mesh_20refinement:Compaction_20length)=
## **Parameters in section** Mesh refinement/Compaction length
(parameters:Mesh_20refinement:Compaction_20length:Mesh_20cells_20per_20compaction_20length)=
### __Parameter name:__ Mesh cells per compaction length
**Default value:** 1.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The desired ratio between compaction length and size of the mesh cells, or, in other words, how many cells the mesh should (at least) have per compaction length. Every cell where this ratio is smaller than the value specified by this parameter (in places with fewer mesh cells per compaction length) is marked for refinement. 

(parameters:Mesh_20refinement:Composition)=
## **Parameters in section** Mesh refinement/Composition
(parameters:Mesh_20refinement:Composition:Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of scaling factors by which every individual compositional field will be multiplied. If only a single compositional field exists, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other. 

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are compositional fields. 

(parameters:Mesh_20refinement:Composition_20approximate_20gradient)=
## **Parameters in section** Mesh refinement/Composition approximate gradient
(parameters:Mesh_20refinement:Composition_20approximate_20gradient:Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of scaling factors by which every individual compositional field gradient will be multiplied. If only a single compositional field exists, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other. 

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are compositional fields. 

(parameters:Mesh_20refinement:Composition_20gradient)=
## **Parameters in section** Mesh refinement/Composition gradient
(parameters:Mesh_20refinement:Composition_20gradient:Compositional_20field_20scaling_20factors)=
### __Parameter name:__ Compositional field scaling factors
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of scaling factors by which every individual compositional field gradient will be multiplied. If only a single compositional field exists, then this parameter has no particular meaning. On the other hand, if multiple criteria are chosen, then these factors are used to weigh the various indicators relative to each other. 

If the list of scaling factors given in this parameter is empty, then this indicates that they should all be chosen equal to one. If the list is not empty then it needs to have as many entries as there are compositional fields. 

(parameters:Mesh_20refinement:Composition_20threshold)=
## **Parameters in section** Mesh refinement/Composition threshold
(parameters:Mesh_20refinement:Composition_20threshold:Compositional_20field_20thresholds)=
### __Parameter name:__ Compositional field thresholds
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** A list of thresholds that every individual compositional field will be evaluated against. 

(parameters:Mesh_20refinement:Isosurfaces)=
## **Parameters in section** Mesh refinement/Isosurfaces
(parameters:Mesh_20refinement:Isosurfaces:Isosurfaces)=
### __Parameter name:__ Isosurfaces
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of isosurfaces separated by semi-colons (;). Each isosurface entry consists of multiple entries separated by a comma. The first two entries indicate the minimum and maximum refinement levels respectively. The entries after the first two describe the fields the isosurface applies to, followed by a colon (:), which is again followed by the minimum and maximum property values separated by bar (|). An example for an isosurface is '0, 2, Temperature: 300 | 600; 2, 2, C\_1: 0.5 | 1'. In this example the mesh refinement is kept between level 0 and level 2 if the temperature is between 300 and 600 and at level 2 when the compositional field C\_1 is between 0.5 and 1. If both happen at the same location and the current refinement level is 1, it means that the first isoline will not set any flag and the second isoline will set a refinement flag. This means the cell will be refined. If both the coarsening and refinement flags are set, preference is given to refinement. 

The first two entries for each isosurface, describing the minimum and maximum grid levels, can be two numbers or contain one of the key values 'min' and 'max'. This indicates the key will be replaced with the global minimum and maximum refinement levels. The 'min' and 'max' keys also accept adding values to be added or substracted from them respectively. This is done by adding a '+' or '-' and a number behind them (e.g. min+2 or max-1). Note that you can't substract a value from a minimum value or add a value to the maximum value. If, for example, `max-4` drops below the minimum or `min+4` goes above the maximum, it will simply use the global minimum and maximum values respectively. The same holds for any mesh refinement level below the global minimum or above the global maximum. 

(parameters:Mesh_20refinement:Maximum_20refinement_20function)=
## **Parameters in section** Mesh refinement/Maximum refinement function
(parameters:Mesh_20refinement:Maximum_20refinement_20function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** depth 

**Pattern:** [Selection depth|cartesian|spherical ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `depth', `cartesian' and `spherical'. `depth' will create a function, in which only the first variable is non-zero, which is interpreted to be the depth of the point. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. 

(parameters:Mesh_20refinement:Maximum_20refinement_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Mesh_20refinement:Maximum_20refinement_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Mesh_20refinement:Maximum_20refinement_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Mesh_20refinement:Minimum_20refinement_20function)=
## **Parameters in section** Mesh refinement/Minimum refinement function
(parameters:Mesh_20refinement:Minimum_20refinement_20function:Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** depth 

**Pattern:** [Selection depth|cartesian|spherical ] 

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `depth', `cartesian' and `spherical'. `depth' will create a function, in which only the first variable is non-zero, which is interpreted to be the depth of the point. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. 

(parameters:Mesh_20refinement:Minimum_20refinement_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Mesh_20refinement:Minimum_20refinement_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Mesh_20refinement:Minimum_20refinement_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Mesh_20refinement:Volume_20of_20fluid_20interface)=
## **Parameters in section** Mesh refinement/Volume of fluid interface
(parameters:Mesh_20refinement:Volume_20of_20fluid_20interface:Strict_20coarsening)=
### __Parameter name:__ Strict coarsening
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If true, then explicitly coarsen any cells not neighboring the VolumeOfFluid interface. 

(parameters:Nullspace_20removal)=
## **Parameters in section** Nullspace removal
(parameters:Nullspace_20removal:Remove_20nullspace)=
### __Parameter name:__ Remove nullspace
**Default value:**  

**Pattern:** [MultipleSelection net rotation|angular momentum|net surface rotation|net translation|linear momentum|net x translation|net y translation|net z translation|linear x momentum|linear y momentum|linear z momentum ] 

**Documentation:** Choose none, one or several from 

\begin{itemize} \item net rotation \item angular momentum \item net translation \item net surface rotation\item linear momentum \item net x translation \item net y translation \item net z translation \item linear x momentum \item linear y momentum \item linear z momentum \end{itemize}

These are a selection of operations to remove certain parts of the nullspace from the velocity after solving. For some geometries and certain boundary conditions the velocity field is not uniquely determined but contains free translations and/or rotations. Depending on what you specify here, these non-determined modes will be removed from the velocity field at the end of the Stokes solve step.


The ``angular momentum'' option removes a rotation such that the net angular momentum is zero. The ``linear * momentum'' options remove translations such that the net momentum in the relevant direction is zero.  The ``net rotation'' option removes the net rotation of the whole domain, the ``net surface rotation'' option removes the net rotation of the top surface, and the ``net * translation'' options remove the net translations in the relevant directions.  For most problems there should not be a significant difference between the momentum and rotation/translation versions of nullspace removal, although the momentum versions are more physically motivated. They are equivalent for constant density simulations, and approximately equivalent when the density variations are small.

Note that while more than one operation can be selected it only makes sense to pick one rotational and one translational operation. 

(parameters:Postprocess)=
## **Parameters in section** Postprocess
(parameters:Postprocess:List_20of_20postprocessors)=
### __Parameter name:__ List of postprocessors
**Default value:**  

**Pattern:** [MultipleSelection Stokes residual|basic statistics|boundary densities|boundary pressures|boundary strain rate residual statistics|boundary velocity residual statistics|command|composition statistics|core statistics|depth average|dynamic topography|entropy viscosity statistics|geoid|global statistics|gravity calculation|heat flux densities|heat flux map|heat flux statistics|heating statistics|load balance statistics|mass flux statistics|material statistics|matrix statistics|maximum depth of field|melt statistics|memory statistics|mobility statistics|particle count statistics|particles|point values|pressure statistics|rotation statistics|spherical velocity statistics|temperature statistics|topography|velocity boundary statistics|velocity statistics|viscous dissipation statistics|visualization|volume of fluid statistics ] 

**Documentation:** A comma separated list of postprocessor objects that should be run at the end of each time step. Some of these postprocessors will declare their own parameters which may, for example, include that they will actually do something only every so many time steps or years. Alternatively, the text `all' indicates that all available postprocessors should be run after each time step.

The following postprocessors are available:

`Stokes residual': A postprocessor that outputs the Stokes residuals during the iterative solver algorithm into a file stokes_residuals.txt in the output directory.

`basic statistics': A postprocessor that outputs some simplified statistics like the Rayleigh number and other quantities that only make sense in certain model setups. The output is written after completing initial adaptive refinement steps. The postprocessor assumes a point at the surface at the adiabatic surface temperature and pressure is a reasonable reference condition for computing these properties. Furthermore, the Rayleigh number is computed using the model depth (i.e. not the radius of the Earth), as we need a definition that is geometry independent. Take care when comparing these values to published studies and make sure they use the same definitions.

`boundary densities': A postprocessor that computes the laterally averaged density at the top and bottom of the domain.

`boundary pressures': A postprocessor that computes the laterally averaged pressure at the top and bottom of the domain.

`boundary strain rate residual statistics': A postprocessor that computes some statistics about the surface strain rate residual along the top boundary. The residual is the difference between the second invariant of the model strain rate and the second strain rate invariant read from the input data file. Currently, the strain residual statistics, i.e., min, max and the rms magnitude, are computed at the top suface.

`boundary velocity residual statistics': A postprocessor that computes some statistics about the velocity residual along the top boundary. The velocity residual is the difference between the model solution velocities and the input velocities (GPlates model or ascii data). Currently, the velocity residual statistics, i.e., min, max and the rms magnitude, is computed at the top suface.

`command': A postprocessor that executes a command line process.

`composition statistics': A postprocessor that computes some statistics about the compositional fields, if present in this simulation. In particular, it computes maximal and minimal values of each field, as well as the total mass contained in this field as defined by the integral $m_i(t) = \int_\Omega c_i(\mathbf x,t) \; \text{d}x$.

`core statistics': A postprocessor that computes some statistics about the core evolution. (Working only with dynamic core boundary temperature plugin)

`depth average': A postprocessor that computes depth averaged quantities and writes them into a file <depth_average.ext> in the output directory, where the extension of the file is determined by the output format you select. In addition to the output format, a number of other parameters also influence this postprocessor, and they can be set in the section \texttt{Postprocess/Depth average} in the input file.

In the output files, the $x$-value of each data point corresponds to the depth, whereas the $y$-value corresponds to the simulation time. The time is provided in seconds or, if the global ``Use years in output instead of seconds'' parameter is set, in years.

`dynamic topography': A postprocessor that computes a measure of dynamic topography based on the stress at the surface and bottom. The data is written into text files named `dynamic\_topography.NNNNN' in the output directory, where NNNNN is the number of the time step.

The exact approach works as follows: At the centers of all cells that sit along the top surface, we evaluate the stress and evaluate the component of it in the direction in which gravity acts. In other words, we compute $\sigma_{rr}={\hat g}^T(2 \eta \varepsilon(\mathbf u)- \frac 13 (\textrm{div}\;\mathbf u)I)\hat g - p_d$ where $\hat g = \mathbf g/\|\mathbf g\|$ is the direction of the gravity vector $\mathbf g$ and $p_d=p-p_a$ is the dynamic pressure computed by subtracting the adiabatic pressure $p_a$ from the total pressure $p$ computed as part of the Stokes solve. From this, the dynamic topography is computed using the formula $h=\frac{\sigma_{rr}}{(\mathbf g \cdot \mathbf n)  \rho}$ where $\rho$ is the density at the cell center. For the bottom surface we chose the convection that positive values are up (out) and negative values are in (down), analogous to the deformation of the upper surface. Note that this implementation takes the direction of gravity into account, which means that reversing the flow in backward advection calculations will not reverse the instantaneous topography because the reverse flow will be divided by the reverse surface gravity.  
The file format then consists of lines with Euclidean coordinates followed by the corresponding topography value.

(As a side note, the postprocessor chooses the cell center instead of the center of the cell face at the surface, where we really are interested in the quantity, since this often gives better accuracy. The results should in essence be the same, though.)

`entropy viscosity statistics': A postprocessor that computes the maximum and volume averagedentropy viscosity stabilization for the temperature field.

`geoid': A postprocessor that computes a representation of the geoid based on the density structure in the mantle, as well as the topography at the surface and core mantle boundary (CMB) if desired. The topography is based on the dynamic topography postprocessor in case of no free surface, and based on the real surface from the geometry model in case of a free surface. The geoid is computed from a spherical harmonic expansion, so the geometry of the domain must be a 3D spherical shell.

`global statistics': A postprocessor that outputs all the global statistics information, e.g. the time of the simulation, the timestep number, number of degrees of freedom and solver iterations for each timestep. The postprocessor can output different formats, the first printing one line in the statistics file per nonlinear solver iteration (if a nonlinear solver scheme is selected). The second prints one line per timestep, summing the information about all nonlinear iterations in this line. Note that this postprocessor is always active independent on whether or not it is selected in the parameter file.

`gravity calculation': A postprocessor that computes gravity, gravity anomalies, gravity potential and gravity gradients for a set of points (e.g. satellites) in or above the model surface for either a user-defined range of latitudes, longitudes and radius or a list of point coordinates.Spherical coordinates in the output file are radius, colatitude and colongitude. Gravity is here based on the density distribution from the material model (and non adiabatic). This means that the density may come directly from an ascii file. This postprocessor also computes theoretical gravity and its derivatives, which corresponds to the analytical solution of gravity in the same geometry but filled with a reference density. The reference density is also used to determine density anomalies for computing gravity anomalies. Thus one must carefully evaluate the meaning of the gravity anomaly output, because the solution may not reflect the actual gravity anomaly (due to differences in the assumed reference density). On way to guarantee correct gravity anomalies is to subtract gravity of a certain point from the average gravity on the map. Another way is to directly use density anomalies for this postprocessor.The average- minimum- and maximum gravity acceleration and potential are written into the statistics file.

`heat flux densities': A postprocessor that computes some statistics about the heat flux density for each boundary id. The heat flux density across each boundary is computed in outward direction, i.e., from the domain to the outside. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in ``Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.''

Note that the ``heat flux statistics'' postprocessor computes the same quantity as the one here, but not divided by the area of the surface. In other words, it computes the \textit{total} heat flux through each boundary.

`heat flux map': A postprocessor that computes the heat flux density across each boundary in outward direction, i.e., from the domain to the outside. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in ``Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.''

`heat flux statistics': A postprocessor that computes some statistics about the heat flux density across each boundary in outward direction, i.e., from the domain to the outside. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in ``Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.''The point-wise heat flux can be obtained from the heat flux map postprocessor, which outputs the heat flux to a file, or the heat flux map visualization postprocessor, which outputs the heat flux for visualization. 

As stated, this postprocessor computes the \textit{outbound} heat flux. If you are interested in the opposite direction, for example from the core into the mantle when the domain describes the mantle, then you need to multiply the result by -1.

\note{In geodynamics, the term ``heat flux'' is often understood to be the quantity $- k \nabla T$, which is really a heat flux \textit{density}, i.e., a vector-valued field. In contrast to this, the current postprocessor only computes the integrated flux over each part of the boundary. Consequently, the units of the quantity computed here are $W=\frac{J}{s}$.}

The ``heat flux densities'' postprocessor computes the same quantity as the one here, but divided by the area of the surface.

`heating statistics': A postprocessor that computes some statistics about heating, averaged by volume. 

`load balance statistics': A postprocessor that computes statistics about the distribution of cells, and if present particles across subdomains. In particular, it computes maximal, average and minimal number of cells across all ranks. If there are particles it also computes the maximal, average, and minimum number of particles across all ranks, and maximal, average, and minimal ratio between local number of particles and local number of cells across all processes. All of these numbers can be useful to assess the load balance between different MPI ranks, as the difference between the mimimal and maximal load should be as small as possible.

`mass flux statistics': A postprocessor that computes some statistics about the mass flux across boundaries. For each boundary indicator (see your geometry description for which boundary indicators are used), the mass flux is computed in outward direction, i.e., from the domain to the outside, using the formula $\int_{\Gamma_i} \rho \mathbf v \cdot \mathbf n$ where $\Gamma_i$ is the part of the boundary with indicator $i$, $\rho$ is the density as reported by the material model, $\mathbf v$ is the velocity, and $\mathbf n$ is the outward normal. 

As stated, this postprocessor computes the \textit{outbound} mass flux. If you are interested in the opposite direction, for example from the core into the mantle when the domain describes the mantle, then you need to multiply the result by -1.

\note{In geodynamics, the term ``mass flux'' is often understood to be the quantity $\rho \mathbf v$, which is really a mass flux \textit{density}, i.e., a vector-valued field. In contrast to this, the current postprocessor only computes the integrated flux over each part of the boundary. Consequently, the units of the quantity computed here are $\frac{kg}{s}$.}

`material statistics': A postprocessor that computes some statistics about the material properties. In particular, it computes the volume-averages of the density and viscosity, and the total mass in the model. Specifically, it implements the following formulas in each time step: $\left<\rho\right> = \frac{1}{|\Omega|} \int_\Omega \rho(\mathbf x) \, \text{d}x$, $\left<\eta\right> = \frac{1}{|\Omega|} \int_\Omega \eta(\mathbf x) \, \text{d}x$, $M = \int_\Omega \rho(\mathbf x) \, \text{d}x$, where $|\Omega|$ is the volume of the domain.

`matrix statistics': A postprocessor that computes some statistics about the matrices. In particular, it outputs total memory consumption, total non-zero elements, and non-zero elements per block, for system matrix and system preconditioner matrix.

`maximum depth of field': A postprocessor that for each compositional field outputs the largest depth at which a quadrature point is found where the field has a value of 0.5 or larger. For fields that do not represent materials, but for example track a certain quantity like strain, this value of 0.5 does not necessarily make sense. 

`melt statistics': A postprocessor that computes some statistics about the melt (volume) fraction. If the material model does not implement a melt fraction function, the output is set to zero.

`memory statistics': A postprocessor that computes some statistics about the memory consumption. In particular, it computes the memory usage of the system matrix, triangulation, p4est, DoFHandler, current constraints, solution vector, and peak virtual memory usage, all in MB. It also outputs the memory usage of the system matrix to the screen.

`mobility statistics': A postprocessor that computes some statistics about mobility following Tackley (2000) and Lourenco et al. (2020).

`particle count statistics': A postprocessor that computes some statistics about the particle distribution, if present in this simulation. In particular, it computes minimal, average and maximal values of particles per cell in the global domain.

`particles': A Postprocessor that creates particles that follow the velocity field of the simulation. The particles can be generated and propagated in various ways and they can carry a number of constant or time-varying properties. The postprocessor can write output positions and properties of all particles at chosen intervals, although this is not mandatory. It also allows other parts of the code to query the particles for information.

`point values': A postprocessor that evaluates the solution (i.e., velocity, pressure, temperature, and compositional fields along with other fields that are treated as primary variables) at the end of every time step or after a user-specified time interval at a given set of points and then writes this data into the file <point\_values.txt> in the output directory. The points at which the solution should be evaluated are specified in the section \texttt{Postprocess/Point values} in the input file.

In the output file, data is organized as (i) time, (ii) the 2 or 3 coordinates of the evaluation points, and (iii) followed by the values of the solution vector at this point. The time is provided in seconds or, if the global ``Use years in output instead of seconds'' parameter is set, in years. In the latter case, the velocity is also converted to meters/year, instead of meters/second.

\note{Evaluating the solution of a finite element field at arbitrarily chosen points is an expensive process. Using this postprocessor will only be efficient if the number of evaluation points or output times is relatively small. If you need a very large number of evaluation points, you should consider extracting this information from the visualization program you use to display the output of the `visualization' postprocessor.}

`pressure statistics': A postprocessor that computes some statistics about the pressure field.

`rotation statistics': A postprocessor that computes some statistics about the rotational velocity of the model (i.e. integrated net rotation and angular momentum). In 2D we assume the model to be a cross-section through an infinite domain in z direction, with a zero z-velocity. Thus, the z-axis is the only possible rotation axis and both moment of inertia and angular momentum are scalar instead of tensor quantities.

`spherical velocity statistics': A postprocessor that computes radial, tangential and total RMS velocity.

`temperature statistics': A postprocessor that computes some statistics about the temperature field.

`topography': A postprocessor intended for use with a deforming top surface.  After every step it loops over all the vertices on the top surface and determines the maximum and minimum topography relative to a reference datum (initial box height for a box geometry model or initial radius for a sphere/spherical shell geometry model). If 'Topography.Output to file' is set to true, also outputs topography into text files named `topography.NNNNN' in the output directory, where NNNNN is the number of the time step.
The file format then consists of lines with Euclidean coordinates followed by the corresponding topography value.Topography is printed/written in meters.

`velocity boundary statistics': A postprocessor that computes some statistics about the velocity along the boundaries. For each boundary indicator (see your geometry description for which boundary indicators are used), the min and max velocity magnitude is computed.

`velocity statistics': A postprocessor that computes some statistics about the velocity field.

`viscous dissipation statistics': A postprocessor that outputs the viscous rate of dissipation of energy for each compositional field (where the field has a value of 0.5 or more) as well as over the whole domain. When all the fields represent lithologies and there is no background field, the sum of the individual field's dissipation should equal that over the whole domain. The viscous dissipation is computed as: $\int_{V}\left(\sigma' \dot{\epsilon}' \right)$, where $\sigma'$  is the deviatoric stress and $\dot{\epsilon}'$ the deviatoric strain rate.Note then when shear heating is included in the temperature equation, it is better to use the 'heating statistics' postprocessor.

`visualization': A postprocessor that takes the solution and writes it into files that can be read by a graphical visualization program. Additional run time parameters are read from the parameter subsection 'Visualization'.

`volume of fluid statistics': A postprocessor that computes some statistics about the volume-of-fluid fields. 

(parameters:Postprocess:Run_20postprocessors_20on_20nonlinear_20iterations)=
### __Parameter name:__ Run postprocessors on nonlinear iterations
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether or not the postprocessors should be executed after each of the nonlinear iterations done within one time step. As this is mainly an option for the purposes of debugging, it is not supported when the 'Time between graphical output' is larger than zero, or when the postprocessor is not intended to be run more than once per timestep. 

(parameters:Postprocess:Boundary_20strain_20rate_20residual_20statistics)=
## **Parameters in section** Postprocess/Boundary strain rate residual statistics
(parameters:Postprocess:Boundary_20strain_20rate_20residual_20statistics:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/postprocess/boundary-strain-rate-residual/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the ascii data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Postprocess:Boundary_20strain_20rate_20residual_20statistics:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_3d_boundary_strain_rate.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the input surface strain rate an ascii data. The file has one column in addition to the coordinate columns corresponding to the second invariant of strain rate.  

(parameters:Postprocess:Boundary_20strain_20rate_20residual_20statistics:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. 

(parameters:Postprocess:Boundary_20velocity_20residual_20statistics)=
## **Parameters in section** Postprocess/Boundary velocity residual statistics
(parameters:Postprocess:Boundary_20velocity_20residual_20statistics:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the GPlates model or the ascii data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Postprocess:Boundary_20velocity_20residual_20statistics:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** current_day.gpml 

**Pattern:** [Anything] 

**Documentation:** The file name of the input velocity as a GPlates model or an ascii data. For the GPlates model, provide file in the same format as described in the 'gplates' boundary velocity plugin. For the ascii data, provide file in the same format as described in  'ascii data' initial composition plugin. 

(parameters:Postprocess:Boundary_20velocity_20residual_20statistics:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/year set this factor to 0.01. 

(parameters:Postprocess:Boundary_20velocity_20residual_20statistics:Use_20ascii_20data)=
### __Parameter name:__ Use ascii data
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Use ascii data files (e.g., GPS) for computing residual velocities instead of GPlates data. 

(parameters:Postprocess:Boundary_20velocity_20residual_20statistics:Use_20spherical_20unit_20vectors)=
### __Parameter name:__ Use spherical unit vectors
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Specify velocity as r, phi, and theta components instead of x, y, and z. Positive velocities point up, east, and north (in 3D) or out and clockwise (in 2D). This setting only makes sense for spherical geometries.GPlates data is always interpreted to be in east and north directions and is not affected by this parameter. 

(parameters:Postprocess:Command)=
## **Parameters in section** Postprocess/Command
(parameters:Postprocess:Command:Command)=
### __Parameter name:__ Command
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Command to execute. 

(parameters:Postprocess:Command:Run_20on_20all_20processes)=
### __Parameter name:__ Run on all processes
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to run command from all processes (true), or only on process 0 (false). 

(parameters:Postprocess:Command:Terminate_20on_20failure)=
### __Parameter name:__ Terminate on failure
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Select whether \aspect{} should terminate if the command returns a non-zero exit status. 

(parameters:Postprocess:Depth_20average)=
## **Parameters in section** Postprocess/Depth average
(parameters:Postprocess:Depth_20average:Depth_20boundaries_20of_20zones)=
### __Parameter name:__ Depth boundaries of zones
**Default value:**  

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The depth boundaries of zones within which we are to compute averages. By default this list is empty and we subdivide the entire domain into equidistant depth zones and compute averages within each of these zones. If this list is not empty it has to contain one more entry than the 'Number of zones' parameter, representing the upper and lower depth boundary of each zone. It is not necessary to cover the whole depth-range (i.e. you can select to only average in a single layer by choosing 2 arbitrary depths as the boundaries of that layer). 

(parameters:Postprocess:Depth_20average:List_20of_20output_20variables)=
### __Parameter name:__ List of output variables
**Default value:** all 

**Pattern:** [MultipleSelection all|temperature|composition|adiabatic temperature|adiabatic pressure|adiabatic density|adiabatic density derivative|velocity magnitude|sinking velocity|rising velocity|Vs|Vp|viscosity|vertical heat flux|vertical mass flux|composition mass ] 

**Documentation:** A comma separated list which specifies which quantities to average in each depth slice. It defaults to averaging all available quantities, but this can be an expensive operation, so you may want to select only a few.

Specifically, the sinking velocity is defined as the scalar product of the velocity and a unit vector in the direction of gravity, if positive (being zero if this product is negative, which would correspond to an upward velocity). The rising velocity is the opposite: the scalar product of the velocity and a unit vector in the direction opposite of gravity, if positive (being zero for downward velocities). 

List of options:
all|temperature|composition|adiabatic temperature|adiabatic pressure|adiabatic density|adiabatic density derivative|velocity magnitude|sinking velocity|rising velocity|Vs|Vp|viscosity|vertical heat flux|vertical mass flux|composition mass 

(parameters:Postprocess:Depth_20average:Number_20of_20zones)=
### __Parameter name:__ Number of zones
**Default value:** 10 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The number of zones in depth direction within which we are to compute averages. By default, we subdivide the entire domain into 10 depth zones and compute temperature and other averages within each of these zones. However, if you have a very coarse mesh, it may not make much sense to subdivide the domain into so many zones and you may wish to choose less than this default. It may also make computations slightly faster. On the other hand, if you have an extremely highly resolved mesh, choosing more zones might also make sense. 

(parameters:Postprocess:Depth_20average:Output_20format)=
### __Parameter name:__ Output format
**Default value:** gnuplot, txt 

**Pattern:** [MultipleSelection none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate|txt ] 

**Documentation:** A list of formats in which the output shall be produced. The format in which the output is generated also determines the extension of the file into which data is written. The list of possible output formats that can be given here is documented in the appendix of the manual where the current parameter is described. By default the output is written as gnuplot file (for plotting), and as a simple text file. 

(parameters:Postprocess:Depth_20average:Time_20between_20graphical_20output)=
### __Parameter name:__ Time between graphical output
**Default value:** 1e8 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The time interval between each generation of graphical output files. A value of zero indicates that output should be generated in each time step. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Postprocess:Dynamic_20core_20statistics)=
## **Parameters in section** Postprocess/Dynamic core statistics
(parameters:Postprocess:Dynamic_20core_20statistics:Excess_20entropy_20only)=
### __Parameter name:__ Excess entropy only
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Output the excess entropy only instead the each entropy terms. 

(parameters:Postprocess:Dynamic_20topography)=
## **Parameters in section** Postprocess/Dynamic topography
(parameters:Postprocess:Dynamic_20topography:Density_20above)=
### __Parameter name:__ Density above
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. This value depends on the density of material that is moved up or down, i.e. crustal rock, and the density of the material that is displaced (generally water or air). While the density of crustal rock is part of the material model, this parameter `Density above' allows the user to specify the density value of material that is displaced above the solid surface. By default this material is assumed to be air, with a density of 0. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Postprocess:Dynamic_20topography:Density_20below)=
### __Parameter name:__ Density below
**Default value:** 9900. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. This value depends on the density of material that is moved up or down, i.e. mantle above CMB, and the density of the material that is displaced (generally outer core material). While the density of mantle rock is part of the material model, this parameter `Density below' allows the user to specify the density value of material that is displaced below the solid surface. By default this material is assumed to be outer core material with a density of 9900. Units: \si{\kilogram\per\meter\cubed}. 

(parameters:Postprocess:Dynamic_20topography:Output_20bottom)=
### __Parameter name:__ Output bottom
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to output a file containing the bottom (i.e., CMB) dynamic topography. 

(parameters:Postprocess:Dynamic_20topography:Output_20surface)=
### __Parameter name:__ Output surface
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to output a file containing the surface dynamic topography. 

(parameters:Postprocess:Geoid)=
## **Parameters in section** Postprocess/Geoid
(parameters:Postprocess:Geoid:Also_20output_20the_20gravity_20anomaly)=
### __Parameter name__: Also output the gravity anomaly
**Alias:** [Output gravity anomaly](parameters:Postprocess:Geoid:Output_20gravity_20anomaly)

**Deprecation Status:** false

(parameters:Postprocess:Geoid:Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20CMB_20dynamic_20topography_20contribution)=
### __Parameter name__: Also output the spherical harmonic coefficients of CMB dynamic topography contribution
**Alias:** [Output CMB topography contribution coefficients](parameters:Postprocess:Geoid:Output_20CMB_20topography_20contribution_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess:Geoid:Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20density_20anomaly_20contribution)=
### __Parameter name__: Also output the spherical harmonic coefficients of density anomaly contribution
**Alias:** [Output density anomaly contribution coefficients](parameters:Postprocess:Geoid:Output_20density_20anomaly_20contribution_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess:Geoid:Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20geoid_20anomaly)=
### __Parameter name__: Also output the spherical harmonic coefficients of geoid anomaly
**Alias:** [Output geoid anomaly coefficients](parameters:Postprocess:Geoid:Output_20geoid_20anomaly_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess:Geoid:Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20surface_20dynamic_20topography_20contribution)=
### __Parameter name__: Also output the spherical harmonic coefficients of surface dynamic topography contribution
**Alias:** [Output surface topography contribution coefficients](parameters:Postprocess:Geoid:Output_20surface_20topography_20contribution_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess:Geoid:Density_20above)=
### __Parameter name:__ Density above
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The density value above the surface boundary. 

(parameters:Postprocess:Geoid:Density_20below)=
### __Parameter name:__ Density below
**Default value:** 9900. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The density value below the CMB boundary. 

(parameters:Postprocess:Geoid:Include_20CMB_20topography_20contribution)=
### __Parameter name:__ Include CMB topography contribution
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Option to include the contribution from CMB topography on geoid. The default is true. 

(parameters:Postprocess:Geoid:Include_20surface_20topography_20contribution)=
### __Parameter name:__ Include surface topography contribution
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Option to include the contribution from surface topography on geoid. The default is true. 

(parameters:Postprocess:Geoid:Include_20the_20contributon_20from_20dynamic_20topography)=
### __Parameter name:__ Include the contributon from dynamic topography
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Option to include the contribution from dynamic topography on geoid. The default is true. 

(parameters:Postprocess:Geoid:Maximum_20degree)=
### __Parameter name:__ Maximum degree
**Default value:** 20 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** This parameter can be a random positive integer. However, the value normally should not exceed the maximum degree of the initial perturbed temperature field. For example, if the initial temperature uses S40RTS, the maximum degree should not be larger than 40. 

(parameters:Postprocess:Geoid:Minimum_20degree)=
### __Parameter name:__ Minimum degree
**Default value:** 2 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** This parameter normally is set to 2 since the perturbed gravitational potential at degree 1 always vanishes in a reference frame with the planetary center of mass same as the center of figure. 

(parameters:Postprocess:Geoid:Output_20CMB_20topography_20contribution_20coefficients)=
### __Parameter name:__ Output CMB topography contribution coefficients
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to output the spherical harmonic coefficients of the CMB topography contribution to the maximum degree. The default is false.  

(parameters:Postprocess:Geoid:Output_20data_20in_20geographical_20coordinates)=
### __Parameter name:__ Output data in geographical coordinates
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to output the geoid anomaly in geographical coordinates (latitude and longitude). The default is false, so postprocess will output the data in geocentric coordinates (x,y,z) as normally. 

(parameters:Postprocess:Geoid:Output_20density_20anomaly_20contribution_20coefficients)=
### __Parameter name:__ Output density anomaly contribution coefficients
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to output the spherical harmonic coefficients of the density anomaly contribution to the maximum degree. The default is false.  

(parameters:Postprocess:Geoid:Output_20geoid_20anomaly_20coefficients)=
### __Parameter name:__ Output geoid anomaly coefficients
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to output the spherical harmonic coefficients of the geoid anomaly up to the maximum degree. The default is false, so postprocess will only output the geoid anomaly in grid format.  

(parameters:Postprocess:Geoid:Output_20gravity_20anomaly)=
### __Parameter name:__ Output gravity anomaly
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to output the free-air gravity anomaly up to the maximum degree. The unit of the output is in SI, hence $m/s^2$ ($1mgal = 10^-5 m/s^2$). The default is false.  

(parameters:Postprocess:Geoid:Output_20surface_20topography_20contribution_20coefficients)=
### __Parameter name:__ Output surface topography contribution coefficients
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Option to output the spherical harmonic coefficients of the surface topography contribution to the maximum degree. The default is false.  

(parameters:Postprocess:Global_20statistics)=
## **Parameters in section** Postprocess/Global statistics
(parameters:Postprocess:Global_20statistics:Write_20statistics_20for_20each_20nonlinear_20iteration)=
### __Parameter name:__ Write statistics for each nonlinear iteration
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to put every nonlinear iteration into a separate line in the statistics file (if true), or to output only one line per time step that contains the total number of iterations of the Stokes and advection linear system solver. 

(parameters:Postprocess:Gravity_20calculation)=
## **Parameters in section** Postprocess/Gravity calculation
(parameters:Postprocess:Gravity_20calculation:List_20of_20latitude)=
### __Parameter name:__ List of latitude
**Default value:**  

**Pattern:** [List of <[Double -90...90 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Parameter for the list of points sampling scheme: List of satellite latitude coordinates. 

(parameters:Postprocess:Gravity_20calculation:List_20of_20longitude)=
### __Parameter name:__ List of longitude
**Default value:**  

**Pattern:** [List of <[Double -180...180 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Parameter for the list of points sampling scheme: List of satellite longitude coordinates. 

(parameters:Postprocess:Gravity_20calculation:List_20of_20radius)=
### __Parameter name:__ List of radius
**Default value:**  

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** Parameter for the list of points sampling scheme: List of satellite radius coordinates. Just specify one radius if all points values have the same radius. If not, make sure there are as many radius as longitude and latitude 

(parameters:Postprocess:Gravity_20calculation:Maximum_20latitude)=
### __Parameter name:__ Maximum latitude
**Default value:** 90 

**Pattern:** [Double -90...90 (inclusive)] 

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the latitude between a minimum and maximum latitude. 

(parameters:Postprocess:Gravity_20calculation:Maximum_20longitude)=
### __Parameter name:__ Maximum longitude
**Default value:** 180. 

**Pattern:** [Double -180...180 (inclusive)] 

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the longitude between a minimum and maximum longitude. 

(parameters:Postprocess:Gravity_20calculation:Maximum_20radius)=
### __Parameter name:__ Maximum radius
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Parameter for the map sampling scheme: Maximum radius can be defined in or outside the model. 

(parameters:Postprocess:Gravity_20calculation:Minimum_20latitude)=
### __Parameter name:__ Minimum latitude
**Default value:** -90. 

**Pattern:** [Double -90...90 (inclusive)] 

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the latitude between a minimum and maximum latitude. 

(parameters:Postprocess:Gravity_20calculation:Minimum_20longitude)=
### __Parameter name:__ Minimum longitude
**Default value:** -180. 

**Pattern:** [Double -180...180 (inclusive)] 

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the longitude between a minimum and maximum longitude. 

(parameters:Postprocess:Gravity_20calculation:Minimum_20radius)=
### __Parameter name:__ Minimum radius
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Parameter for the map sampling scheme: Minimum radius may be defined in or outside the model. Prescribe a minimum radius for a sampling coverage at a specific height. 

(parameters:Postprocess:Gravity_20calculation:Number_20points_20fibonacci_20spiral)=
### __Parameter name:__ Number points fibonacci spiral
**Default value:** 200 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Parameter for the fibonacci spiral sampling scheme: This specifies the desired number of satellites per radius layer. The default value is 200. Note that sampling becomes more uniform with increasing number of satellites 

(parameters:Postprocess:Gravity_20calculation:Number_20points_20latitude)=
### __Parameter name:__ Number points latitude
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Parameter for the map sampling scheme: This specifies the number of points along the latitude (e.g. gravity map) between a minimum and maximum latitude. 

(parameters:Postprocess:Gravity_20calculation:Number_20points_20longitude)=
### __Parameter name:__ Number points longitude
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Parameter for the map sampling scheme: This specifies the number of points along the longitude (e.g. gravity map) between a minimum and maximum longitude. 

(parameters:Postprocess:Gravity_20calculation:Number_20points_20radius)=
### __Parameter name:__ Number points radius
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Parameter for the map sampling scheme: This specifies the number of points along the radius (e.g. depth profile) between a minimum and maximum radius. 

(parameters:Postprocess:Gravity_20calculation:Precision_20in_20gravity_20output)=
### __Parameter name:__ Precision in gravity output
**Default value:** 12 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Set the precision of gravity acceleration, potential and gradients in the gravity output and statistics file. 

(parameters:Postprocess:Gravity_20calculation:Quadrature_20degree_20increase)=
### __Parameter name:__ Quadrature degree increase
**Default value:** 0 

**Pattern:** [Integer range -1...2147483647 (inclusive)] 

**Documentation:** Quadrature degree increase over the velocity element degree may be required when gravity is calculated near the surface or inside the model. An increase in the quadrature element adds accuracy to the gravity solution from noise due to the model grid. 

(parameters:Postprocess:Gravity_20calculation:Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Gravity anomalies may be computed using density anomalies relative to a reference density. 

(parameters:Postprocess:Gravity_20calculation:Sampling_20scheme)=
### __Parameter name:__ Sampling scheme
**Default value:** map 

**Pattern:** [Selection map|list|list of points|fibonacci spiral ] 

**Documentation:** Choose the sampling scheme. By default, the map produces a grid of equally angled points between a minimum and maximum radius, longitude, and latitude. A list of points contains the specific coordinates of the satellites. The fibonacci spiral sampling scheme produces a uniformly distributed map on the surface of sphere defined by a minimum and/or maximum radius. 

(parameters:Postprocess:Gravity_20calculation:Time_20between_20gravity_20output)=
### __Parameter name:__ Time between gravity output
**Default value:** 1e8 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The time interval between each generation of gravity output files. A value of 0 indicates that output should be generated in each time step. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Postprocess:Gravity_20calculation:Time_20steps_20between_20gravity_20output)=
### __Parameter name:__ Time steps between gravity output
**Default value:** 2147483647 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The maximum number of time steps between each generation of gravity output files. 

(parameters:Postprocess:Memory_20statistics)=
## **Parameters in section** Postprocess/Memory statistics
(parameters:Postprocess:Memory_20statistics:Output_20peak_20virtual_20memory_20_28VmPeak_29)=
### __Parameter name:__ Output peak virtual memory _28VmPeak_29
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** If set to 'true', also output the peak virtual memory usage (computed as the maximum over all processors). 

(parameters:Postprocess:Particles)=
## **Parameters in section** Postprocess/Particles
(parameters:Postprocess:Particles:Allow_20cells_20without_20particles)=
### __Parameter name:__ Allow cells without particles
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** By default, every cell needs to contain particles to use this interpolator plugin. If this parameter is set to true, cells are allowed to have no particles, In case both the current cell and its neighbors are empty, the interpolator will return 0 for the current cell's properties. 

(parameters:Postprocess:Particles:Data_20output_20format)=
### __Parameter name:__ Data output format
**Default value:** vtu 

**Pattern:** [MultipleSelection none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate|ascii ] 

**Documentation:** A comma separated list of file formats to be used for graphical output. The list of possible output formats that can be given here is documented in the appendix of the manual where the current parameter is described. 

(parameters:Postprocess:Particles:Exclude_20output_20properties)=
### __Parameter name:__ Exclude output properties
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A comma separated list of particle properties that should \textit{not} be output. If this list contains the entry `all', only the id of particles will be provided in graphical output files. 

(parameters:Postprocess:Particles:Integration_20scheme)=
### __Parameter name:__ Integration scheme
**Default value:** rk2 

**Pattern:** [Selection euler|rk2|rk4 ] 

**Documentation:** This parameter is used to decide which method to use to solve the equation that describes the position of particles, i.e., $\frac{d}{dt}\mathbf x_k(t) = \mathbf u(\mathbf x_k(t),t)$, where $k$ is an index that runs over all particles, and $\mathbf u(\mathbf x,t)$ is the velocity field that results from the Stokes equations.

In practice, the exact velocity $\mathbf u(\mathbf x,t)$ is of course not available, but only a numerical approximation $\mathbf u_h(\mathbf x,t)$. Furthermore, this approximation is only available at discrete time steps, $\mathbf u^n(\mathbf x)=\mathbf u(\mathbf x,t^n)$, and these need to be interpolated between time steps if the integrator for the equation above requires an evaluation at time points between the discrete time steps. If we denote this interpolation in time by $\tilde{\mathbf u}_h(\mathbf x,t)$ where $\tilde{\mathbf u}_h(\mathbf x,t^n)=\mathbf u^n(\mathbf x)$, then the equation the differential equation solver really tries to solve is $\frac{d}{dt}\tilde{\mathbf x}_k(t) =  \tilde{\mathbf u}_h(\mathbf x_k(t),t)$.

As a consequence of these considerations, if you try to assess convergence properties of an ODE integrator -- for example to verify that the RK4 integrator converges with fourth order --, it is important to recall that the integrator may not solve the equation you think it solves. If, for example, we call the numerical solution of the ODE $\tilde{\mathbf x}_{k,h}(t)$, then the error will typically satisfy a relationship like \[  \| \tilde{\mathbf x}_k(T) - \tilde{\mathbf x}_{k,h}(T) \|  \le  C(T) \Delta t^p\] where $\Delta t$ is the time step and $p$ the convergence order of the method, and $C(T)$ is a (generally unknown) constant that depends on the end time $T$ at which one compares the solutions. On the other hand, an analytically computed trajectory would likely use the \textit{exact} velocity, and one may be tempted to compute $\| \mathbf x_k(T) - \tilde{\mathbf x}_{k,h}(T) \|$, but this quantity will, in the best case, only satisfy an estimate of the form \[  \| \mathbf x_k(T) - \tilde{\mathbf x}_{k,h}(T) \|  \le  C_1(T) \Delta t^p  + C_2(T) \| \mathbf u-\mathbf u_h \|  + C_3(T) \| \mathbf u_h-\tilde{\mathbf u}_h \|\] with appropriately chosen norms for the second and third term. These second and third terms typically converge to zero at relatively low rates (compared to the order $p$ of the integrator, which can often be chosen relatively high) in the mesh size $h$ and the time step size $\\Delta t$, limiting the overall accuracy of the ODE integrator.

Select one of the following models:

`euler': Explicit Euler scheme integrator, where $y_{n+1} = y_n + \Delta t \, v(y_n)$. This requires only one integration substep per timestep.

`rk2': Second Order Runge Kutta integrator $y_{n+1} = y_n + \Delta t\, v(t_{n+1/2}, y_{n} + \frac{1}{2} k_1)$ where $k_1 = \Delta t\, v(t_{n}, y_{n})$

`rk4': Runge Kutta fourth order integrator, where $y_{n+1} = y_n + \frac{1}{6} k_1 + \frac{1}{3} k_2 + \frac{1}{3} k_3 + \frac{1}{6} k_4$ and $k_1$, $k_2$, $k_3$, $k_4$ are defined as usual. 

(parameters:Postprocess:Particles:Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** cell average 

**Pattern:** [Selection bilinear least squares|cell average|harmonic average|nearest neighbor|quadratic least squares ] 

**Documentation:** Select one of the following models:

`bilinear least squares': Uses linear least squares to obtain the slopes and center of a 2D or 3D plane from the particle positions and a particular property value on those particles. Interpolate this property onto a vector of points. If the limiter is enabled then it will ensure the interpolated properties do not exceed the range of the minimum and maximum of the values of the property on the particles. Note that deal.II must be configured with BLAS and LAPACK to support this operation.

`cell average': Return the arithmetic average of all particle properties in the given cell, or in the neighboring cells if the given cell is empty. In case the neighboring cells are also empty, and 'Allow cells without particles' is set to true, the interpolator returns 0. Otherwise, an exception is thrown. 

`harmonic average': Return the harmonic average of all particle properties in the given cell. If the cell contains no particles, return the harmonic average of the properties in the neighboring cells. In case the neighboring cells are also empty, and 'Allow cells without particles' is set to true, the interpolator returns 0. Otherwise, an exception is thrown. 

`nearest neighbor': Return the properties of the nearest neighboring particle in the current cell, or nearest particle in nearest neighboring cell if current cell is empty. In case the neighboring cells are also empty, and 'Allow cells without particles' is set to true, the interpolator returns 0. Otherwise, an exception is thrown. 

`quadratic least squares': Interpolates particle properties onto a vector of points using a quadratic least squares method. Note that deal.II must be configured with BLAS/LAPACK. 

(parameters:Postprocess:Particles:List_20of_20particle_20properties)=
### __Parameter name:__ List of particle properties
**Default value:**  

**Pattern:** [MultipleSelection composition|elastic stress|function|initial composition|initial position|integrated strain|integrated strain invariant|melt particle|pT path|position|reference position|velocity|viscoplastic strain invariants ] 

**Documentation:** A comma separated list of particle properties that should be tracked. By default none is selected, which means only position, velocity and id of the particles are output. 

The following properties are available:

`composition': Implementation of a plugin in which the particle property is defined by the compositional fields in the model. This can be used to track solid compositionevolution over time.

`elastic stress': A plugin in which the particle property tensor is defined as the total elastic stress a particle has accumulated. See the viscoelastic material model documentation for more detailed information.

`function': Implementation of a model in which the particle property is set by evaluating an explicit function at the initial position of each particle. The function is defined in the parameters in section ``Particles|Function''. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`initial composition': Implementation of a plugin in which the particle property is given as the initial composition at the particle's initial position. The particle gets as many properties as there are compositional fields.

`initial position': Implementation of a plugin in which the particle property is given as the initial position of the particle. This property is vector-valued with as many components as there are space dimensions. In practice, it is often most useful to only visualize one of the components of this vector, or the magnitude of the vector. For example, in a spherical mantle simulation, the magnitude of this property equals the starting radius of a particle, and is thereby indicative of which part of the mantle a particle comes from.

`integrated strain': A plugin in which the particle property tensor is defined as the deformation gradient tensor $\mathbf F$ this particle has experienced. $\mathbf F$ can be polar-decomposed into the left stretching tensor $\mathbf L$ (the finite strain we are interested in), and the rotation tensor $\mathbf Q$. See the corresponding cookbook in the manual for more detailed information.

`integrated strain invariant': A plugin in which the particle property is defined as the finite strain invariant ($\varepsilon_{ii}$). This property is calculated with the timestep ($dt$) and the second invariant of the deviatoric strain rate tensor ($\dot{\varepsilon}_{ii}$), where the value at time step $n$ is $\varepsilon_{ii}^{n} = \varepsilon_{ii}^{n-1} + dt\dot{\varepsilon}_{ii}$.

`melt particle': Implementation of a plugin in which the particle property is defined as presence of melt above a threshold, which can be set as an input parameter. This property is set to 0 if melt is not present and set to 1 if melt is present.

`pT path': Implementation of a plugin in which the particle property is defined as the current pressure and temperature at this position. This can be used to generate pressure-temperature paths of material points over time.

`position': Implementation of a plugin in which the particle property is defined as the current position.

`reference position': Implementation of a plugin in which the particle property is defined as the current reference position.

`velocity': Implementation of a plugin in which the particle property is defined as the recent velocity at this position.

`viscoplastic strain invariants': A plugin that calculates the finite strain invariant a particle has experienced and assigns it to either the plastic and/or viscous strain field based on whether the material is plastically yielding, or the total strain field used in the visco plastic material model. The implementation of this property is equivalent to the implementation for compositional fields that is located in the plugin in \texttt{benchmarks/buiter\_et\_al\_2008\_jgr/plugin/},and is effectively the same as what the visco plastic material model uses for compositional fields. 

(parameters:Postprocess:Particles:Load_20balancing_20strategy)=
### __Parameter name:__ Load balancing strategy
**Default value:** repartition 

**Pattern:** [MultipleSelection none|remove particles|add particles|remove and add particles|repartition ] 

**Documentation:** Strategy that is used to balance the computational load across processors for adaptive meshes. 

(parameters:Postprocess:Particles:Maximum_20particles_20per_20cell)=
### __Parameter name:__ Maximum particles per cell
**Default value:** 100 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Upper limit for particle number per cell. This limit is useful for adaptive meshes to prevent coarse cells from slowing down the whole model. It will be checked and enforced after mesh refinement, after MPI transfer of particles and after particle movement. If there are \texttt{n\_number\_of\_particles} $>$ \texttt{max\_particles\_per\_cell} particles in one cell then \texttt{n\_number\_of\_particles} - \texttt{max\_particles\_per\_cell} particles in this cell are randomly chosen and destroyed. 

(parameters:Postprocess:Particles:Minimum_20particles_20per_20cell)=
### __Parameter name:__ Minimum particles per cell
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Lower limit for particle number per cell. This limit is useful for adaptive meshes to prevent fine cells from being empty of particles. It will be checked and enforced after mesh refinement and after particle movement. If there are \texttt{n\_number\_of\_particles} $<$ \texttt{min\_particles\_per\_cell} particles in one cell then \texttt{min\_particles\_per\_cell} - \texttt{n\_number\_of\_particles} particles are generated and randomly placed in this cell. If the particles carry properties the individual property plugins control how the properties of the new particles are initialized. 

(parameters:Postprocess:Particles:Number_20of_20grouped_20files)=
### __Parameter name:__ Number of grouped files
**Default value:** 16 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** VTU file output supports grouping files from several CPUs into a given number of files using MPI I/O when writing on a parallel filesystem. Select 0 for no grouping. This will disable parallel file output and instead write one file per processor. A value of 1 will generate one big file containing the whole solution, while a larger value will create that many files (at most as many as there are MPI ranks). 

(parameters:Postprocess:Particles:Number_20of_20particles)=
### __Parameter name:__ Number of particles
**Default value:** 1000 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Total number of particles to create (not per processor or per element). The number is parsed as a floating point number (so that one can specify, for example, '1e4' particles) but it is interpreted as an integer, of course. 

(parameters:Postprocess:Particles:Particle_20generator_20name)=
### __Parameter name:__ Particle generator name
**Default value:** random uniform 

**Pattern:** [Selection ascii file|probability density function|quadrature points|random uniform|reference cell|uniform box|uniform radial ] 

**Documentation:** Select one of the following models:

`ascii file': Generates a distribution of particles from coordinates specified in an Ascii data file. The file format is a simple text file, with as many columns as spatial dimensions and as many lines as particles to be generated. Initial comment lines starting with `#' will be discarded. Note that this plugin always generates as many particles as there are coordinates in the data file, the ``Postprocess/Particles/Number of particles'' parameter has no effect on this plugin. All of the values that define this generator are read from a section ``Postprocess/Particles/Generator/Ascii file'' in the input file, see Section~\ref{parameters:Postprocess/Particles/Generator/Ascii_20file}.

`probability density function': Generate a random distribution of particles over the entire simulation domain. The probability density is prescribed in the form of a user-prescribed function. The format of this function follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}. The return value of the function is always checked to be a non-negative probability density but it can be zero in parts of the domain.

`quadrature points': Generates particles at the quadrature points of each active cell of the triangulation. Here, Gauss quadrature of degree (velocity\_degree + 1), is used similarly to the assembly of Stokes matrix.

`random uniform': Generates a random uniform distribution of particles over the entire simulation domain.

`reference cell': Generates a uniform distribution of particles per cell and spatial direction in the unit cell and transforms each of the particles back to real region in the model domain. Uniform here means the particles will be generated with an equal spacing in each spatial dimension.

`uniform box': Generate a uniform distribution of particles over a rectangular domain in 2D or 3D. Uniform here means the particles will be generated with an equal spacing in each spatial dimension. Note that in order to produce a regular distribution the number of generated particles might not exactly match the one specified in the input file.

`uniform radial': Generate a uniform distribution of particles over a spherical domain in 2D or 3D. Uniform here means the particles will be generated with an equal spacing in each spherical spatial dimension, i.e., the particles are created at positions that increase linearly with equal spacing in radius, colatitude and longitude around a certain center point. Note that in order to produce a regular distribution the number of generated particles might not exactly match the one specified in the input file. 

(parameters:Postprocess:Particles:Particle_20weight)=
### __Parameter name:__ Particle weight
**Default value:** 10 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Weight that is associated with the computational load of a single particle. The sum of particle weights will be added to the sum of cell weights to determine the partitioning of the mesh if the `repartition' particle load balancing strategy is selected. The optimal weight depends on the used integrator and particle properties. In general for a more expensive integrator and more expensive properties a larger particle weight is recommended. Before adding the weights of particles, each cell already carries a weight of 1000 to account for the cost of field-based computations. 

(parameters:Postprocess:Particles:Temporary_20output_20location)=
### __Parameter name:__ Temporary output location
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** On large clusters it can be advantageous to first write the output to a temporary file on a local file system and later move this file to a network file system. If this variable is set to a non-empty string it will be interpreted as a temporary storage location. 

(parameters:Postprocess:Particles:Time_20between_20data_20output)=
### __Parameter name:__ Time between data output
**Default value:** 1e8 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The time interval between each generation of output files. A value of zero indicates that output should be generated every time step.

Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Postprocess:Particles:Update_20ghost_20particles)=
### __Parameter name:__ Update ghost particles
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Some particle interpolation algorithms require knowledge about particles in neighboring cells. To allow this, particles in ghost cells need to be exchanged between the processes neighboring this cell. This parameter determines whether this transport is happening. 

(parameters:Postprocess:Particles:Write_20in_20background_20thread)=
### __Parameter name:__ Write in background thread
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** File operations can potentially take a long time, blocking the progress of the rest of the model run. Setting this variable to `true' moves this process into a background thread, while the rest of the model continues. 

(parameters:Postprocess:Particles:Function)=
## **Parameters in section** Postprocess/Particles/Function
(parameters:Postprocess:Particles:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Postprocess:Particles:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Postprocess:Particles:Function:Number_20of_20components)=
### __Parameter name:__ Number of components
**Default value:** 1 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of function components where each component is described by a function expression delimited by a ';'. 

(parameters:Postprocess:Particles:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Postprocess:Particles:Generator)=
## **Parameters in section** Postprocess/Particles/Generator
(parameters:Postprocess:Particles:Generator:Ascii_20file)=
## **Parameters in section** Postprocess/Particles/Generator/Ascii file
(parameters:Postprocess:Particles:Generator:Ascii_20file:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/particle/generator/ascii/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the particle data. This path may either be absolute (if starting with a '/') or relative to the current directory. The path may also include the special text '$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.  

(parameters:Postprocess:Particles:Generator:Ascii_20file:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** particle.dat 

**Pattern:** [Anything] 

**Documentation:** The name of the particle file. 

(parameters:Postprocess:Particles:Generator:Probability_20density_20function)=
## **Parameters in section** Postprocess/Particles/Generator/Probability density function
(parameters:Postprocess:Particles:Generator:Probability_20density_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Postprocess:Particles:Generator:Probability_20density_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Postprocess:Particles:Generator:Probability_20density_20function:Random_20cell_20selection)=
### __Parameter name:__ Random cell selection
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** If true, particle numbers per cell are calculated randomly according to their respective probability density. This means particle numbers per cell can deviate statistically from the integral of the probability density. If false, first determine how many particles each cell should have based on the integral of the density over each of the cells, and then once we know how many particles we want on each cell, choose their locations randomly within each cell. 

(parameters:Postprocess:Particles:Generator:Probability_20density_20function:Random_20number_20seed)=
### __Parameter name:__ Random number seed
**Default value:** 5432 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The seed for the random number generator that controls the particle generation. Keep constant to generate identical particle distributions in subsequent model runs. Change to get a different distribution. In parallel computations the seed is further modified on each process to ensure different particle patterns on different processes. Note that the number of particles per processor is not affected by the seed. 

(parameters:Postprocess:Particles:Generator:Probability_20density_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Postprocess:Particles:Generator:Reference_20cell)=
## **Parameters in section** Postprocess/Particles/Generator/Reference cell
(parameters:Postprocess:Particles:Generator:Reference_20cell:Number_20of_20particles_20per_20cell_20per_20direction)=
### __Parameter name:__ Number of particles per cell per direction
**Default value:** 2 

**Pattern:** [List of <[Integer range 1...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** List of number of particles to create per cell and spatial dimension. The size of the list is the number of spatial dimensions. If only one value is given, then each spatial dimension is set to the same value. The list of numbers are parsed as a floating point number (so that one can specify, for example, '1e4' particles) but it is interpreted as an integer, of course. 

(parameters:Postprocess:Particles:Generator:Uniform_20box)=
## **Parameters in section** Postprocess/Particles/Generator/Uniform box
(parameters:Postprocess:Particles:Generator:Uniform_20box:Maximum_20x)=
### __Parameter name:__ Maximum x
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximum x coordinate for the region of particles. 

(parameters:Postprocess:Particles:Generator:Uniform_20box:Maximum_20y)=
### __Parameter name:__ Maximum y
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximum y coordinate for the region of particles. 

(parameters:Postprocess:Particles:Generator:Uniform_20box:Maximum_20z)=
### __Parameter name:__ Maximum z
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximum z coordinate for the region of particles. 

(parameters:Postprocess:Particles:Generator:Uniform_20box:Minimum_20x)=
### __Parameter name:__ Minimum x
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimum x coordinate for the region of particles. 

(parameters:Postprocess:Particles:Generator:Uniform_20box:Minimum_20y)=
### __Parameter name:__ Minimum y
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimum y coordinate for the region of particles. 

(parameters:Postprocess:Particles:Generator:Uniform_20box:Minimum_20z)=
### __Parameter name:__ Minimum z
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimum z coordinate for the region of particles. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial)=
## **Parameters in section** Postprocess/Particles/Generator/Uniform radial
(parameters:Postprocess:Particles:Generator:Uniform_20radial:Center_20x)=
### __Parameter name:__ Center x
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** x coordinate for the center of the spherical region, where particles are generated. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Center_20y)=
### __Parameter name:__ Center y
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** y coordinate for the center of the spherical region, where particles are generated. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Center_20z)=
### __Parameter name:__ Center z
**Default value:** 0. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** z coordinate for the center of the spherical region, where particles are generated. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Maximum_20latitude)=
### __Parameter name:__ Maximum latitude
**Default value:** 180. 

**Pattern:** [Double 0...180 (inclusive)] 

**Documentation:** Maximum latitude coordinate for the region of particles in degrees. Measured from the center position, and from the north pole. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Maximum_20longitude)=
### __Parameter name:__ Maximum longitude
**Default value:** 360. 

**Pattern:** [Double -180...360 (inclusive)] 

**Documentation:** Maximum longitude coordinate for the region of particles in degrees. Measured from the center position. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Maximum_20radius)=
### __Parameter name:__ Maximum radius
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Maximum radial coordinate for the region of particles. Measured from the center position. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Minimum_20latitude)=
### __Parameter name:__ Minimum latitude
**Default value:** 0. 

**Pattern:** [Double 0...180 (inclusive)] 

**Documentation:** Minimum latitude coordinate for the region of particles in degrees. Measured from the center position, and from the north pole. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Minimum_20longitude)=
### __Parameter name:__ Minimum longitude
**Default value:** 0. 

**Pattern:** [Double -180...360 (inclusive)] 

**Documentation:** Minimum longitude coordinate for the region of particles in degrees. Measured from the center position. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Minimum_20radius)=
### __Parameter name:__ Minimum radius
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Minimum radial coordinate for the region of particles. Measured from the center position. 

(parameters:Postprocess:Particles:Generator:Uniform_20radial:Radial_20layers)=
### __Parameter name:__ Radial layers
**Default value:** 1 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** The number of radial shells of particles that will be generated around the central point. 

(parameters:Postprocess:Particles:Interpolator)=
## **Parameters in section** Postprocess/Particles/Interpolator
(parameters:Postprocess:Particles:Interpolator:Bilinear_20least_20squares)=
## **Parameters in section** Postprocess/Particles/Interpolator/Bilinear least squares
(parameters:Postprocess:Particles:Interpolator:Bilinear_20least_20squares:Use_20linear_20least_20squares_20limiter)=
### __Parameter name:__ Use linear least squares limiter
**Default value:** false 

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)] 

**Documentation:** Limit the interpolation of particle properties onto the cell so the value of each property is no smaller than its minimum and no larger than its maximum on the particles in each cell. If more than one value is specified, they will be treated as a list. 

(parameters:Postprocess:Particles:Interpolator:Quadratic_20least_20squares)=
## **Parameters in section** Postprocess/Particles/Interpolator/Quadratic least squares
(parameters:Postprocess:Particles:Interpolator:Quadratic_20least_20squares:Global_20particle_20property_20maximum)=
### __Parameter name:__ Global particle property maximum
**Default value:** 1.7976931348623157e+308 

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The maximum global particle property values that will be used as a limiter for the quadratic least squares interpolation. The number of the input 'Global particle property maximum' values separated by ',' has to be the same as the number of particle properties. 

(parameters:Postprocess:Particles:Interpolator:Quadratic_20least_20squares:Global_20particle_20property_20minimum)=
### __Parameter name:__ Global particle property minimum
**Default value:** -1.7976931348623157e+308 

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)] 

**Documentation:** The minimum global particle property that will be used as a limiter for the quadratic least squares interpolation. The number of the input 'Global particle property minimum' values separated by ',' has to be the same as the number of particle properties. 

(parameters:Postprocess:Particles:Interpolator:Quadratic_20least_20squares:Use_20limiter)=
### __Parameter name:__ Use limiter
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to apply a global particle property limiting scheme to the interpolated particle properties. 

(parameters:Postprocess:Particles:Melt_20particle)=
## **Parameters in section** Postprocess/Particles/Melt particle
(parameters:Postprocess:Particles:Melt_20particle:Threshold_20for_20melt_20presence)=
### __Parameter name:__ Threshold for melt presence
**Default value:** 1e-3 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** The minimum porosity that has to be present at the position of a particle for it to be considered a melt particle (in the sense that the melt presence property is set to 1). 

(parameters:Postprocess:Point_20values)=
## **Parameters in section** Postprocess/Point values
(parameters:Postprocess:Point_20values:Evaluation_20points)=
### __Parameter name:__ Evaluation points
**Default value:**  

**Pattern:** [List of <[List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 2...2 (inclusive)]> of length 0...4294967295 (inclusive) separated by <;>] 

**Documentation:** The list of points at which the solution should be evaluated. Points need to be separated by semicolons, and coordinates of each point need to be separated by commas. 

(parameters:Postprocess:Point_20values:Time_20between_20point_20values_20output)=
### __Parameter name:__ Time between point values output
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The time interval between each generation of point values output. A value of zero indicates that output should be generated in each time step. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Postprocess:Point_20values:Use_20natural_20coordinates)=
### __Parameter name:__ Use natural coordinates
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether or not the Evaluation points are specified in the natural coordinates of the geometry model, e.g. radius, lon, lat for the chunk model. Currently, natural coordinates for the spherical shell and sphere geometries are not supported.  

(parameters:Postprocess:Rotation_20statistics)=
## **Parameters in section** Postprocess/Rotation statistics
(parameters:Postprocess:Rotation_20statistics:Output_20full_20moment_20of_20inertia_20tensor)=
### __Parameter name:__ Output full moment of inertia tensor
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to write the full moment of inertia tensor into the statistics output instead of its norm for the current rotation axis. This is a second-order symmetric tensor with 6 components in 3D. In 2D this option has no effect, because the rotation axis is fixed and thus the moment of inertia is always a scalar. 

(parameters:Postprocess:Rotation_20statistics:Use_20constant_20density_20of_20one)=
### __Parameter name:__ Use constant density of one
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use a constant density of one for the computation of the angular momentum and moment of inertia. This is an approximation that assumes that the 'volumetric' rotation is equal to the 'mass' rotation. If this parameter is true this postprocessor computes 'net rotation' instead of 'angular momentum'. 

(parameters:Postprocess:Topography)=
## **Parameters in section** Postprocess/Topography
(parameters:Postprocess:Topography:Output_20to_20file)=
### __Parameter name:__ Output to file
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether or not to write topography to a text file named named 'topography.NNNNN' in the output directory 

(parameters:Postprocess:Topography:Time_20between_20text_20output)=
### __Parameter name:__ Time between text output
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The time interval between each generation of text output files. A value of zero indicates that output should be generated in each time step. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Postprocess:Visualization)=
## **Parameters in section** Postprocess/Visualization
(parameters:Postprocess:Visualization:Filter_20output)=
### __Parameter name:__ Filter output
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** deal.II offers the possibility to filter duplicate vertices for HDF5 output files. This merges the vertices of adjacent cells and therefore saves disk space, but misrepresents discontinuous output properties. Activating this function reduces the disk space by about a factor of $2^{dim}$ for HDF5 output, and currently has no effect on other output formats. \note{\textbf{Warning:} Setting this flag to true will result in visualization output that does not accurately represent discontinuous fields. This may be because you are using a discontinuous finite element for the pressure, temperature, or compositional variables, or because you use a visualization postprocessor that outputs quantities as discontinuous fields (e.g., the strain rate, viscosity, etc.). These will then all be visualized as \textit{continuous} quantities even though, internally, \aspect{} considers them as discontinuous fields.} 

(parameters:Postprocess:Visualization:Interpolate_20output)=
### __Parameter name:__ Interpolate output
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** deal.II offers the possibility to linearly interpolate output fields of higher order elements to a finer resolution. This somewhat compensates the fact that most visualization software only offers linear interpolation between grid points and therefore the output file is a very coarse representation of the actual solution field. Activating this option increases the spatial resolution in each dimension by a factor equal to the polynomial degree used for the velocity finite element (usually 2). In other words, instead of showing one quadrilateral or hexahedron in the visualization per cell on which \aspect{} computes, it shows multiple (for quadratic elements, it will describe each cell of the mesh on which we compute as $2\times 2$ or $2\times 2\times 2$ cells in 2d and 3d, respectively; correspondingly more subdivisions are used if you use cubic, quartic, or even higher order elements for the velocity).

The effect of using this option can be seen in the following picture showing a variation of the output produced with the input files from Section~\ref{sec:shell-simple-2d}:

\begin{center}  \includegraphics[width=0.5\textwidth]{viz/parameters/build-patches}\end{center}Here, the left picture shows one visualization cell per computational cell (i.e., the option is switched off), and the right picture shows the same simulation with the option switched on (which is the default). The images show the same data, demonstrating that interpolating the solution onto bilinear shape functions as is commonly done in visualizing data loses information.

Of course, activating this option also greatly increases the amount of data \aspect{} will write to disk: approximately by a factor of 4 in 2d, and a factor of 8 in 3d, when using quadratic elements for the velocity, and correspondingly more for even higher order elements. 

(parameters:Postprocess:Visualization:List_20of_20output_20variables)=
### __Parameter name:__ List of output variables
**Default value:**  

**Pattern:** [MultipleSelection ISA rotation timescale|Vp anomaly|Vs anomaly|adiabat|artificial viscosity|artificial viscosity composition|boundary indicators|boundary strain rate residual|boundary velocity residual|compositional vector|depth|dynamic topography|error indicator|geoid|grain lag angle|gravity|heat flux map|heating|material properties|maximum horizontal compressive stress|melt fraction|melt material properties|named additional outputs|nonadiabatic pressure|nonadiabatic temperature|particle count|partition|principal stress|shear stress|spd factor|spherical velocity components|strain rate|strain rate tensor|stress|stress second invariant|surface dynamic topography|surface stress|temperature anomaly|vertical heat flux|volume of fluid values|volumetric strain rate|density|specific heat|thermal conductivity|thermal diffusivity|thermal expansivity|viscosity ] 

**Documentation:** A comma separated list of visualization objects that should be run whenever writing graphical output. By default, the graphical output files will always contain the primary variables velocity, pressure, and temperature. However, one frequently wants to also visualize derived quantities, such as the thermodynamic phase that corresponds to a given temperature-pressure value, or the corresponding seismic wave speeds. The visualization objects do exactly this: they compute such derived quantities and place them into the output file. The current parameter is the place where you decide which of these additional output variables you want to have in your output file.

The following postprocessors are available:

`ISA rotation timescale': A visualization output object that generates output showing the timescale for the rotation of grains toward the infinite strain axis. Kaminski and Ribe (see \cite{Kaminski2002}) call this quantity $\tau_\text{ISA}$ and define it as $\tau_\text{ISA} \approx \frac{1}{\dot{\epsilon}}$ where $\dot{\epsilon}$ is the largest eigenvalue of the strain rate tensor. It can be used, along with the grain lag angle $\Theta$, to calculate the grain orientation lag parameter.

`Vp anomaly': A visualization output object that generates output showing the percentage anomaly in the seismic compressional wave speed $V_p$ as a spatially variable function with one value per cell. This anomaly is either shown as a percentage anomaly relative to the reference profile given by adiabatic conditions (with the compositions given by the current composition, such that the reference could potentially change through time), or as a percentage change relative to the laterally averaged velocity at the depth of the cell. This velocity is calculated by linear interpolation between average values calculated within equally thick depth slices. The number of depth slices in the domain is user-defined. Typically, the best results will be obtained if the number of depth slices is balanced between being large enough to capture step changes in velocities, but small enough to maintain a reasonable number of evaluation points per slice. Bear in mind that lateral averaging subsamples the finite element mesh. Note that this plugin requires a material model that provides seismic velocities.

`Vs anomaly': A visualization output object that generates output showing the percentage anomaly in the seismic shear wave speed $V_s$ as a spatially variable function with one value per cell. This anomaly is either shown as a percentage anomaly relative to the reference profile given by adiabatic conditions (with the compositions given by the current composition, such that the reference could potentially change through time), or as a percentage change relative to the laterally averaged velocity at the depth of the cell. This velocity is calculated by linear interpolation between average values calculated within equally thick depth slices. The number of depth slices in the domain is user-defined. Typically, the best results will be obtained if the number of depth slices is balanced between being large enough to capture step changes in velocities, but small enough to maintain a reasonable number of evaluation points per slice. Bear in mind that lateral averaging subsamples the finite element mesh. Note that this plugin requires a material model that provides seismic velocities.

`adiabat': A visualization output object that generates adiabatic temperature, pressure, density, and density derivative as produced by AdiabaticConditions.

`artificial viscosity': A visualization output object that generates output showing the value of the artificial viscosity on each cell.

`artificial viscosity composition': A visualization output object that generates output showing the value of the artificial viscosity for a compositional field on each cell.

`boundary indicators': A visualization output object that generates output about the used boundary indicators. In a loop over the active cells, if a cell lies at a domain boundary, the boundary indicator of the face along the boundary is requested. In case the cell does not lie along any domain boundary, the cell is assigned the value of the largest used boundary indicator plus one. When a cell is situated in one of the corners of the domain, multiple faces will have a boundary indicator. This postprocessor returns the value of the first face along a boundary that is encountered in a loop over all the faces. 

`boundary strain rate residual': A visualization output object that generates output for the strain rate residual at the top surface. The residual is computed at each point at the surface as the difference between the strain rate invariant in the model and the input data, where the invariant is computed like in the 'strain rate' postprocessor. The user chooses the input data as ascii data files with coordinate columns and column corresponding to the surface strain rate norm.

`boundary velocity residual': A visualization output object that generates output for the velocity residual at the top surface. The residual is computed at each point at the surface as the difference between the modeled velocities and the input data velocities for each vector component. The user has an option to choose the input data as ascii data files (e.g. GPS velocities) with columns in the same format as described for the 'ascii data' initial temperature plugin or a velocity field computed from the GPlates program as described in the gplates boundary velocity plugin. 

`compositional vector': A visualization output object that outputs vectors whose components are derived from compositional fields. Input parameters for this postprocessor are defined in section Postprocess/Visualization/Compositional fields as vectors

`depth': A visualization output postprocessor that outputs the depth for all points inside the domain, as determined by the geometry model.

`dynamic topography': A visualization output object that generates output for the dynamic topography at the top and bottom of the model space. The approach to determine the dynamic topography requires us to compute the stress tensor and evaluate the component of it in the direction in which gravity acts. In other words, we compute $\sigma_{rr}={\hat g}^T(2 \eta \varepsilon(\mathbf u)-\frac 13 (\textrm{div}\;\mathbf u)I)\hat g - p_d$ where $\hat g = \mathbf g/\|\mathbf g\|$ is the direction of the gravity vector $\mathbf g$ and $p_d=p-p_a$ is the dynamic pressure computed by subtracting the adiabatic pressure $p_a$ from the total pressure $p$ computed as part of the Stokes solve. From this, the dynamic topography is computed using the formula $h=\frac{\sigma_{rr}}{(\mathbf g \cdot \mathbf n)  \rho}$ where $\rho$ is the density at the cell center. For the bottom surface we chose the convection that positive values are up (out) and negative values are in (down), analogous to the deformation of the upper surface. Note that this implementation takes the direction of gravity into account, which means that reversing the flow in backward advection calculations will not reverse the instantaneous topography because the reverse flow will be divided by the reverse surface gravity.

Strictly speaking, the dynamic topography is of course a quantity that is only of interest at the surface. However, we compute it everywhere to make things fit into the framework within which we produce data for visualization. You probably only want to visualize whatever data this postprocessor generates at the surface of your domain and simply ignore the rest of the data generated.

Alternatively, consider using the "surface dynamic topography" visualization postprocessor to only output the dynamic topography at the boundary of the domain.

`error indicator': A visualization output object that generates output showing the estimated error or other mesh refinement indicator as a spatially variable function with one value per cell.

`geoid': Visualization for the geoid solution. The geoid is given by the equivalent water column height due to a gravity perturbation. Units: \si{\meter}.

`grain lag angle': A visualization output object that generates output showing the angle between the ~infinite strain axis and the flow velocity. Kaminski and Ribe (see \cite{Kaminski2002}) call this quantity $\Theta$ and define it as $\Theta = \cos^{-1}(\hat{u}\cdot\hat{e})$  where $\hat{u}=\vec{u}/|{u}|$, $\vec{u}$ is the local flow velocity, and $\hat{e}$ is the local infinite strain axis, which we calculate as the first eigenvector of the 'left stretch' tensor. $\Theta$ can be used to calculate the grain orientation lag parameter.

`gravity': A visualization output object that outputs the gravity vector.

`heat flux map': A visualization output object that generates output for the heat flux density across the top and bottom boundary in outward direction. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in ``Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.'' If only conductive heat flux through Dirichlet boundaries is of interest, the postprocessor can produce output of higher resolution by evaluating the CBF solution vector point-wise instead of computing cell-wise averaged values.

`heating': A visualization output object that generates output for all the heating terms used in the energy equation.

`material properties': A visualization output object that generates output for the material properties given by the material model. The current postprocessor allows to output a (potentially large) subset of all of the information provided by material models at once, with just a single material model evaluation per output point. Although individual properties can still be listed in the ``List of output variables'', this visualization plugin is called internally to avoid duplicated evaluations of the material model. 

In almost all places inside \aspect{}, the program can use ``averaged'' material properties, for example for the assembly of matrices and right hand side vectors. To accurately reflect the material parameters used internally, this visualization postprocessor averages in the same way as is used to do the assembly, and consequently the graphical output will reflect not pointwise properties, but averaged properties.

`maximum horizontal compressive stress': A plugin that computes the direction and magnitude of the maximum horizontal component of the compressive stress as a vector field. The direction of this vector can often be used to visualize the principal mode of deformation (e.g., at normal faults or extensional margins) and can be correlated with seismic anisotropy. Recall that the \textit{compressive} stress is simply the negative stress, $\sigma_c=-\sigma=-\left[     2\eta (\varepsilon(\mathbf u)             - \frac 13 (\nabla \cdot \mathbf u) I)     + pI\right]$.

Following \cite{LundTownend07}, we define the maximum horizontal stress direction as that \textit{horizontal} direction $\mathbf n$ that maximizes $\mathbf n^T \sigma_c \mathbf n$. We call a vector \textit{horizontal} if it is perpendicular to the gravity vector $\mathbf g$.

In two space dimensions, $\mathbf n$ is simply a vector that is horizontal (we choose one of the two possible choices). This direction is then scaled by the size of the horizontal stress in this direction, i.e., the plugin outputs the vector $\mathbf w = (\mathbf n^T \sigma_c \mathbf n) \; \mathbf n$.

In three space dimensions, given two horizontal, perpendicular, unit length, but otherwise arbitrarily chosen vectors $\mathbf u,\mathbf v$, we can express $\mathbf n = (\cos \alpha)\mathbf u + (\sin\alpha)\mathbf v$ where $\alpha$ maximizes the expression \begin{align*}  f(\alpha) = \mathbf n^T \sigma_c \mathbf n  = (\mathbf u^T \sigma_c \mathbf u)(\cos\alpha)^2    +2(\mathbf u^T \sigma_c \mathbf v)(\cos\alpha)(\sin\alpha)    +(\mathbf v^T \sigma_c \mathbf v)(\sin\alpha)^2.\end{align*}

The maximum of $f(\alpha)$ is attained where $f'(\alpha)=0$. Evaluating the derivative and using trigonometric identities, one finds that $\alpha$ has to satisfy the equation \begin{align*}  \tan(2\alpha) = \frac{2.0\mathbf u^T \sigma_c \mathbf v}                          {\mathbf u^T \sigma_c \mathbf u                            - \mathbf v^T \sigma_c \mathbf v}.\end{align*}Since the transform $\alpha\mapsto\alpha+\pi$ flips the direction of $\mathbf n$, we only need to seek a solution to this equation in the interval $\alpha\in[0,\pi)$. These are given by $\alpha_1=\frac 12 \arctan \frac{\mathbf u^T \sigma_c \mathbf v}{\mathbf u^T \sigma_c \mathbf u - \mathbf v^T \sigma_c \mathbf v}$ and $\alpha_2=\alpha_1+\frac{\pi}{2}$, one of which will correspond to a minimum and the other to a maximum of $f(\alpha)$. One checks the sign of $f''(\alpha)=-2(\mathbf u^T \sigma_c \mathbf u - \mathbf v^T \sigma_c \mathbf v)\cos(2\alpha) - 2 (\mathbf u^T \sigma_c \mathbf v) \sin(2\alpha)$ for each of these to determine the $\alpha$ that maximizes $f(\alpha)$, and from this immediately arrives at the correct form for the maximum horizontal stress $\mathbf n$.

The description above computes a 3d \textit{direction} vector $\mathbf n$. If one were to scale this vector the same way as done in 2d, i.e., with the magnitude of the stress in this direction, one will typically get vectors whose length is principally determined by the hydrostatic pressure at a given location simply because the hydrostatic pressure is the largest component of the overall stress. On the other hand, the hydrostatic pressure does not determine any principal direction because it is an isotropic, anti-compressive force. As a consequence, there are often points in simulations (e.g., at the center of convection rolls) where the stress has no dominant horizontal direction, and the algorithm above will then in essence choose a random direction because the stress is approximately equal in all horizontal directions. If one scaled the output by the magnitude of the stress in this direction (i.e., approximately equal to the hydrostatic pressure at this point), one would get randomly oriented vectors at these locations with significant lengths.

To avoid this problem, we scale the maximal horizontal compressive stress direction $\mathbf n$ by the \textit{difference} between the stress in the maximal and minimal horizontal stress directions. In other words, let $\mathbf n_\perp=(\sin \alpha)\mathbf u - (\cos\alpha)\mathbf v$ be the horizontal direction perpendicular to $\mathbf n$, then this plugin outputs the vector quantity $\mathbf w = (\mathbf n^T \sigma_c \mathbf n                -\mathbf n^T_\perp \sigma_c \mathbf n_\perp)               \; \mathbf n$. In other words, the length of the vector produced indicates \textit{how dominant} the direction of maximal horizontal compressive strength is.

Fig.~\ref{fig:max-horizontal-compressive-stress} shows a simple example for this kind of visualization in 3d.

\begin{figure}  \includegraphics[width=0.3\textwidth]    {viz/plugins/maximum_horizontal_compressive_stress/temperature.png}  \hfill  \includegraphics[width=0.3\textwidth]    {viz/plugins/maximum_horizontal_compressive_stress/velocity.png}  \hfill  \includegraphics[width=0.3\textwidth]    {viz/plugins/maximum_horizontal_compressive_stress/horizontal-stress.png}  \caption{\it Illustration of the `maximum horizontal     compressive stress' visualization plugin. The left     figure shows a ridge-like temperature anomaly. Together     with no-slip boundary along all six boundaries, this     results in two convection rolls (center). The maximal     horizontal compressive strength at the bottom center     of the domain is perpendicular to the ridge because     the flow comes together there from the left and right,     yielding a compressive force in left-right direction.     At the top of the model, the flow separates outward,     leading to a \textit{negative} compressive stress     in left-right direction; because there is no flow     in front-back direction, the compressive strength     in front-back direction is zero, making the along-ridge     direction the dominant one. At the center of the     convection rolls, both horizontal directions yield     the same stress; the plugin therefore chooses an     essentially arbitrary horizontal vector, but then     uses a zero magnitude given that the difference     between the maximal and minimal horizontal stress     is zero at these points.}  \label{fig:max-horizontal-compressive-stress}\end{figure}

`melt fraction': A visualization output object that generates output for the melt fraction at the temperature and pressure of the current point. If the material model computes a melt fraction, this is the quantity that will be visualized. Otherwise, a specific parametrization for batch melting (as described in the following) will be used. It does not take into account latent heat. If there are no compositional fields, or no fields called 'pyroxenite',  this postprocessor will visualize the melt fraction of peridotite (calculated using the anhydrous model of Katz, 2003). If there is a compositional field called 'pyroxenite', the postprocessor assumes that this compositional field is the content of pyroxenite, and will visualize the melt fraction for a mixture of peridotite and pyroxenite (using the melting model of Sobolev, 2011 for pyroxenite). All the parameters that were used in these calculations can be changed in the input file, the most relevant maybe being the mass fraction of Cpx in peridotite in the Katz melting model (Mass fraction cpx), which right now has a default of 15\%. The corresponding p-T-diagrams can be generated by running the tests melt\_postprocessor\_peridotite and melt\_postprocessor\_pyroxenite.

`melt material properties': A visualization output object that generates output for melt related properties of the material model. Note that this postprocessor always outputs the compaction pressure, but can output a large range of additional properties, as selected in the ``List of properties'' parameter.

`named additional outputs': Some material models can compute quantities other than those that typically appear in the equations that \aspect{} solves (such as the viscosity, density, etc). Examples of quantities material models may be able to compute are seismic velocities, or other quantities that can be derived from the state variables and the material coefficients such as the stress or stress anisotropies. These quantities are generically referred to as `named outputs' because they are given an explicit name different from the usual outputs of material models.

This visualization postprocessor outputs whatever quantities the material model can compute. What quantities these are is specific to the material model in use for a simulation, and for many models in fact does not contain any named outputs at all.

`nonadiabatic pressure': A visualization output object that generates output for the non-adiabatic component of the pressure.

The variable that is outputted this way is computed by taking the pressure at each point and subtracting from it the adiabatic pressure computed at the beginning of the simulation. Because the adiabatic pressure is one way of defining a static pressure background field, what this visualization postprocessor therefore produces is \textit{one} way to compute a \textit{dynamic pressure}. There are, however, other ways as well, depending on the choice of the ``background pressure''.

`nonadiabatic temperature': A visualization output object that generates output for the non-adiabatic component of the temperature.

`particle count': A visualization output object that generates output about the number of particles per cell.

`partition': A visualization output object that generates output for the parallel partition that every cell of the mesh is associated with.

`principal stress': A visualization output object that outputs the principal stress values and directions, i.e., the eigenvalues and eigenvectors of the stress tensor. The postprocessor can either operate on the full stress tensor or only on the deviatoric stress tensor.

`shear stress': A visualization output object that generates output for the 3 (in 2d) or 6 (in 3d) components of the shear stress tensor, i.e., for the components of the tensor $-2\eta\varepsilon(\mathbf u)$ in the incompressible case and $-2\eta\left[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I\right]$ in the compressible case. If elasticity is used, the elastic contribution is being accounted for. The shear stress differs from the full stress tensor by the absence of the pressure. Note that the convention of positive compressive stress is followed. 

`spd factor': A visualization output object that generates output for the spd factor. The spd factor is a factor which scales a part of the Jacobian used for the Newton solver to make sure that the Jacobian remains positive definite.

`spherical velocity components': A visualization output object that outputs the polar coordinates components $v_r$ and $v_\phi$ of the velocity field in 2D and the spherical coordinates components $v_r$, $v_{\phi}$ and $v_{\theta}$ of the velocity field in 3D.

`strain rate': A visualization output object that generates output for the norm of the strain rate, i.e., for the quantity $\sqrt{\varepsilon(\mathbf u):\varepsilon(\mathbf u)}$ in the incompressible case and $\sqrt{[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I]:[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I]}$ in the compressible case.

`strain rate tensor': A visualization output object that generates output for the 4 (in 2d) or 9 (in 3d) components of the strain rate tensor, i.e., for the components of the tensor $\varepsilon(\mathbf u)$ in the incompressible case and $\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I$ in the compressible case.

`stress': A visualization output object that generates output for the 3 (in 2d) or 6 (in 3d) components of the stress tensor, i.e., for the components of the tensor $-2\eta\varepsilon(\mathbf u)+pI$ in the incompressible case and $-2\eta\left[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I\right]+pI$ in the compressible case. If elasticity is used, the elastic contribution is being accounted for. Note that the convention of positive compressive stress is followed. 

`stress second invariant': A visualization output object that outputs the second moment invariant of the deviatoric stress tensor.

`surface dynamic topography': A visualization output object that generates output for the dynamic topography at the top and bottom of the model space. The approach to determine the dynamic topography requires us to compute the stress tensor and evaluate the component of it in the direction in which gravity acts. In other words, we compute $\sigma_{rr}={\hat g}^T(2 \eta \varepsilon(\mathbf u)-\frac 13 (\textrm{div}\;\mathbf u)I)\hat g - p_d$ where $\hat g = \mathbf g/\|\mathbf g\|$ is the direction of the gravity vector $\mathbf g$ and $p_d=p-p_a$ is the dynamic pressure computed by subtracting the adiabatic pressure $p_a$ from the total pressure $p$ computed as part of the Stokes solve. From this, the dynamic topography is computed using the formula $h=\frac{\sigma_{rr}}{(\mathbf g \cdot \mathbf n)  \rho}$ where $\rho$ is the density at the cell center. For the bottom surface we chose the convection that positive values are up (out) and negative values are in (down), analogous to the deformation of the upper surface. Note that this implementation takes the direction of gravity into account, which means that reversing the flow in backward advection calculations will not reverse the instantaneous topography because the reverse flow will be divided by the reverse surface gravity.

In contrast to the `dynamic topography' visualization postprocessor, this plugin really only evaluates the dynamic topography at faces of cells that are adjacent to `bottom' and `top' boundaries, and only outputs information on the surface of the domain, rather than padding the information with zeros in the interior of the domain.

`surface stress': A visualization output object that generates output on the surface of the domain for the 3 (in 2d) or 6 (in 3d) components of the stress tensor, i.e., for the components of the tensor $-2\eta\varepsilon(\mathbf u)+pI$ in the incompressible case and $-2\eta\left[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I\right]+pI$ in the compressible case. If elasticity is included, its contribution is accounted for. Note that the convention of positive compressive stress is followed. 

`temperature anomaly': A visualization output postprocessor that outputs the temperature minus the depth-average of the temperature.The average temperature is calculated using the lateral averaging function from the ``depth average'' postprocessor and interpolated linearly between the layers specified through ``Number of depth slices''

`vertical heat flux': A visualization output object that generates output for the heat flux in the vertical direction, which is the sum of the advective and the conductive heat flux, with the sign convention of positive flux upwards.

`volume of fluid values': A visualization output object that outputs the volume fraction and optionally a level set field and the interface normal vectors of volume of fluid fields.

`volumetric strain rate': A visualization output object that generates output for the volumetric strain rate, i.e., for the quantity $\nabla\cdot\mathbf u = \textrm{div}\; \mathbf u = \textrm{trace}\; \varepsilon(\mathbf u)$. This should be zero (in some average sense) in incompressible convection models, but can be non-zero in compressible models and models with melt transport. 

(parameters:Postprocess:Visualization:Number_20of_20grouped_20files)=
### __Parameter name:__ Number of grouped files
**Default value:** 16 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** VTU file output supports grouping files from several CPUs into a given number of files using MPI I/O when writing on a parallel filesystem. Select 0 for no grouping. This will disable parallel file output and instead write one file per processor. A value of 1 will generate one big file containing the whole solution, while a larger value will create that many files (at most as many as there are MPI ranks). 

(parameters:Postprocess:Visualization:Output_20format)=
### __Parameter name:__ Output format
**Default value:** vtu 

**Pattern:** [Selection none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate ] 

**Documentation:** The file format to be used for graphical output. The list of possible output formats that can be given here is documented in the appendix of the manual where the current parameter is described. 

(parameters:Postprocess:Visualization:Output_20mesh_20velocity)=
### __Parameter name:__ Output mesh velocity
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** For computations with deforming meshes, ASPECT uses an Arbitrary-Lagrangian-Eulerian formulation to handle deforming the domain, so the mesh has its own velocity field.  This may be written as an output field by setting this parameter to true. 

(parameters:Postprocess:Visualization:Point_2dwise_20stress_20and_20strain)=
### __Parameter name:__ Point_2dwise stress and strain
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If set to true, quantities related to stress and strain are computed in each vertex. Otherwise, an average per cell is computed. 

(parameters:Postprocess:Visualization:Temporary_20output_20location)=
### __Parameter name:__ Temporary output location
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** On large clusters it can be advantageous to first write the output to a temporary file on a local file system and later move this file to a network file system. If this variable is set to a non-empty string it will be interpreted as a temporary storage location. 

(parameters:Postprocess:Visualization:Time_20between_20graphical_20output)=
### __Parameter name:__ Time between graphical output
**Default value:** 1e8 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The time interval between each generation of graphical output files. A value of zero indicates that output should be generated in each time step. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Postprocess:Visualization:Time_20steps_20between_20graphical_20output)=
### __Parameter name:__ Time steps between graphical output
**Default value:** 2147483647 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The maximum number of time steps between each generation of graphical output files. 

(parameters:Postprocess:Visualization:Write_20higher_20order_20output)=
### __Parameter name:__ Write higher order output
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** deal.II offers the possibility to write vtu files with higher order representations of the output data. This means each cell will correctly show the higher order representation of the output data instead of the linear interpolation between vertices that ParaView and Visit usually show. Note that activating this option is safe and recommended, but requires that (i) ``Output format'' is set to ``vtu'', (ii) ``Interpolate output'' is set to true, (iii) you use a sufficiently new version of Paraview or Visit to read the files (Paraview version 5.5 or newer, and Visit version to be determined), and (iv) you use deal.II version 9.1.0 or newer. 
The effect of using this option can be seen in the following picture:

\begin{center}  \includegraphics[width=0.5\textwidth]{viz/parameters/higher-order-output}\end{center}The top figure shows the plain output without interpolation or higher order output. The middle figure shows output that was interpolated as discussed for the ``Interpolate output'' option. The bottom panel shows higher order output that achieves better accuracy than the interpolated output at a lower memory cost. 

(parameters:Postprocess:Visualization:Write_20in_20background_20thread)=
### __Parameter name:__ Write in background thread
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** File operations can potentially take a long time, blocking the progress of the rest of the model run. Setting this variable to `true' moves this process into a background thread, while the rest of the model continues. 

(parameters:Postprocess:Visualization:Artificial_20viscosity_20composition)=
## **Parameters in section** Postprocess/Visualization/Artificial viscosity composition
(parameters:Postprocess:Visualization:Artificial_20viscosity_20composition:Name_20of_20compositional_20field)=
### __Parameter name:__ Name of compositional field
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** The name of the compositional field whose output should be visualized.  

(parameters:Postprocess:Visualization:Compositional_20fields_20as_20vectors)=
## **Parameters in section** Postprocess/Visualization/Compositional fields as vectors
(parameters:Postprocess:Visualization:Compositional_20fields_20as_20vectors:Names_20of_20fields)=
### __Parameter name:__ Names of fields
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** A list of sets of compositional fields which should be output as vectors. Sets are separated from each other by semicolons and vector components within each set are separated by commas (e.g. $vec1_x$, $vec1_y$ ; $vec2_x$, $vec2_y$) where each name must be a defined named compositional field. If only one name is given in a set, it is interpreted as the first in a sequence of dim consecutive compositional fields. 

(parameters:Postprocess:Visualization:Compositional_20fields_20as_20vectors:Names_20of_20vectors)=
### __Parameter name:__ Names of vectors
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** Names of vectors as they will appear in the output. 

(parameters:Postprocess:Visualization:Heat_20flux_20map)=
## **Parameters in section** Postprocess/Visualization/Heat flux map
(parameters:Postprocess:Visualization:Heat_20flux_20map:Output_20point_20wise_20heat_20flux)=
### __Parameter name:__ Output point wise heat flux
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** A boolean flag that controls whether to output the heat flux map as a point wise value, or as a cell-wise averaged value. The point wise output is more accurate, but it currently omits prescribed heat flux values at boundaries and advective heat flux that is caused by velocities non-tangential to boundaries. If you do not use these two features it is recommended to switch this setting on to benefit from the increased output resolution. 

(parameters:Postprocess:Visualization:Material_20properties)=
## **Parameters in section** Postprocess/Visualization/Material properties
(parameters:Postprocess:Visualization:Material_20properties:List_20of_20material_20properties)=
### __Parameter name:__ List of material properties
**Default value:** density,thermal expansivity,specific heat,viscosity 

**Pattern:** [MultipleSelection viscosity|density|thermal expansivity|specific heat|thermal conductivity|thermal diffusivity|compressibility|entropy derivative temperature|entropy derivative pressure|reaction terms|melt fraction ] 

**Documentation:** A comma separated list of material properties that should be written whenever writing graphical output. By default, the material properties will always contain the density, thermal expansivity, specific heat and viscosity. The following material properties are available:

viscosity|density|thermal expansivity|specific heat|thermal conductivity|thermal diffusivity|compressibility|entropy derivative temperature|entropy derivative pressure|reaction terms|melt fraction 

(parameters:Postprocess:Visualization:Melt_20fraction)=
## **Parameters in section** Postprocess/Visualization/Melt fraction
(parameters:Postprocess:Visualization:Melt_20fraction:A1)=
### __Parameter name:__ A1
**Default value:** 1085.7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius}. 

(parameters:Postprocess:Visualization:Melt_20fraction:A2)=
### __Parameter name:__ A2
**Default value:** 1.329e-7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of peridotite. \si{\degreeCelsius\per\pascal}. 

(parameters:Postprocess:Visualization:Melt_20fraction:A3)=
### __Parameter name:__ A3
**Default value:** -5.1e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of peridotite. \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Postprocess:Visualization:Melt_20fraction:B1)=
### __Parameter name:__ B1
**Default value:** 1475.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius}. 

(parameters:Postprocess:Visualization:Melt_20fraction:B2)=
### __Parameter name:__ B2
**Default value:** 8.0e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. \si{\degreeCelsius\per\pascal}. 

(parameters:Postprocess:Visualization:Melt_20fraction:B3)=
### __Parameter name:__ B3
**Default value:** -3.2e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Postprocess:Visualization:Melt_20fraction:C1)=
### __Parameter name:__ C1
**Default value:** 1780.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius}. 

(parameters:Postprocess:Visualization:Melt_20fraction:C2)=
### __Parameter name:__ C2
**Default value:** 4.50e-8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the liquidus of peridotite. \si{\degreeCelsius\per\pascal}. 

(parameters:Postprocess:Visualization:Melt_20fraction:C3)=
### __Parameter name:__ C3
**Default value:** -2.0e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the liquidus of peridotite. \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Postprocess:Visualization:Melt_20fraction:D1)=
### __Parameter name:__ D1
**Default value:** 976.0 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of pyroxenite. Units: \si{\degreeCelsius}. 

(parameters:Postprocess:Visualization:Melt_20fraction:D2)=
### __Parameter name:__ D2
**Default value:** 1.329e-7 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of pyroxenite. Note that this factor is different from the value given in Sobolev, 2011, because they use the potential temperature whereas we use the absolute temperature. \si{\degreeCelsius\per\pascal}. 

(parameters:Postprocess:Visualization:Melt_20fraction:D3)=
### __Parameter name:__ D3
**Default value:** -5.1e-18 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of pyroxenite. \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Postprocess:Visualization:Melt_20fraction:E1)=
### __Parameter name:__ E1
**Default value:** 663.8 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear depletion term in the quadratic function that approximates the melt fraction of pyroxenite. \si{\degreeCelsius\per\pascal}. 

(parameters:Postprocess:Visualization:Melt_20fraction:E2)=
### __Parameter name:__ E2
**Default value:** -611.4 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the quadratic depletion term in the quadratic function that approximates the melt fraction of pyroxenite. \si{\degreeCelsius\per\pascal\squared}. 

(parameters:Postprocess:Visualization:Melt_20fraction:Mass_20fraction_20cpx)=
### __Parameter name:__ Mass fraction cpx
**Default value:** 0.15 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Mass fraction of clinopyroxene in the peridotite to be molten. Units: non-dimensional. 

(parameters:Postprocess:Visualization:Melt_20fraction:beta)=
### __Parameter name:__ beta
**Default value:** 1.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Exponent of the melting temperature in the melt fraction calculation. Units: non-dimensional. 

(parameters:Postprocess:Visualization:Melt_20fraction:r1)=
### __Parameter name:__ r1
**Default value:** 0.5 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Constant in the linear function that approximates the clinopyroxene reaction coefficient. Units: non-dimensional. 

(parameters:Postprocess:Visualization:Melt_20fraction:r2)=
### __Parameter name:__ r2
**Default value:** 8e-11 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Prefactor of the linear pressure term in the linear function that approximates the clinopyroxene reaction coefficient. Units: \si{\per\pascal}. 

(parameters:Postprocess:Visualization:Melt_20material_20properties)=
## **Parameters in section** Postprocess/Visualization/Melt material properties
(parameters:Postprocess:Visualization:Melt_20material_20properties:List_20of_20properties)=
### __Parameter name:__ List of properties
**Default value:** compaction viscosity,permeability 

**Pattern:** [MultipleSelection compaction viscosity|fluid viscosity|permeability|fluid density|fluid density gradient|is melt cell|darcy coefficient|darcy coefficient no cutoff|compaction length ] 

**Documentation:** A comma separated list of melt properties that should be written whenever writing graphical output. The following material properties are available:

compaction viscosity|fluid viscosity|permeability|fluid density|fluid density gradient|is melt cell|darcy coefficient|darcy coefficient no cutoff|compaction length 

(parameters:Postprocess:Visualization:Principal_20stress)=
## **Parameters in section** Postprocess/Visualization/Principal stress
(parameters:Postprocess:Visualization:Principal_20stress:Use_20deviatoric_20stress)=
### __Parameter name:__ Use deviatoric stress
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to use the deviatoric stress tensor instead of the full stress tensor to compute principal stress directions and values. 

(parameters:Postprocess:Visualization:Temperature_20anomaly)=
## **Parameters in section** Postprocess/Visualization/Temperature anomaly
(parameters:Postprocess:Visualization:Temperature_20anomaly:Number_20of_20depth_20slices)=
### __Parameter name:__ Number of depth slices
**Default value:** 20 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of depth slices used to define average temperature. 

(parameters:Postprocess:Visualization:Temperature_20anomaly:Use_20maximal_20temperature_20for_20bottom)=
### __Parameter name:__ Use maximal temperature for bottom
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** If true, use the specified boundary temperatures as average temperatures at the surface. If false, extrapolate the temperature gradient between the first and second cells to the surface. This option will only work for models with a fixed surface temperature.  

(parameters:Postprocess:Visualization:Temperature_20anomaly:Use_20minimal_20temperature_20for_20surface)=
### __Parameter name:__ Use minimal temperature for surface
**Default value:** true 

**Pattern:** [Bool] 

**Documentation:** Whether to use the minimal specified boundary temperature as the bottom boundary temperature. This option will only work for models with a fixed bottom boundary temperature.  

(parameters:Postprocess:Visualization:Volume_20of_20Fluid)=
## **Parameters in section** Postprocess/Visualization/Volume of Fluid
(parameters:Postprocess:Visualization:Volume_20of_20Fluid:Output_20interface_20normals)=
### __Parameter name:__ Output interface normals
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Include the internal data for the interface normal on the unit cells. 

(parameters:Postprocess:Visualization:Volume_20of_20Fluid:Output_20interface_20reconstruction_20contour)=
### __Parameter name:__ Output interface reconstruction contour
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Include fields defined such that the 0 contour is the fluid interface. 

(parameters:Postprocess:Visualization:Vp_20anomaly)=
## **Parameters in section** Postprocess/Visualization/Vp anomaly
(parameters:Postprocess:Visualization:Vp_20anomaly:Average_20velocity_20scheme)=
### __Parameter name:__ Average velocity scheme
**Default value:** reference profile 

**Pattern:** [Selection reference profile|lateral average ] 

**Documentation:** Scheme to compute the average velocity-depth profile. The reference profile option evaluates the conditions along the reference adiabat according to the material model. The lateral average option instead calculates a lateral average from subdivision of the mesh. The lateral average option may produce spurious results where there are sharp velocity changes. 

(parameters:Postprocess:Visualization:Vp_20anomaly:Number_20of_20depth_20slices)=
### __Parameter name:__ Number of depth slices
**Default value:** 50 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of depth slices used to define average seismic compressional wave velocities from which anomalies are calculated. Units: non-dimensional. 

(parameters:Postprocess:Visualization:Vs_20anomaly)=
## **Parameters in section** Postprocess/Visualization/Vs anomaly
(parameters:Postprocess:Visualization:Vs_20anomaly:Average_20velocity_20scheme)=
### __Parameter name:__ Average velocity scheme
**Default value:** reference profile 

**Pattern:** [Selection reference profile|lateral average ] 

**Documentation:** Scheme to compute the average velocity-depth profile. The reference profile option evaluates the conditions along the reference adiabat according to the material model. The lateral average option instead calculates a lateral average from subdivision of the mesh. The lateral average option may produce spurious results where there are sharp velocity changes. 

(parameters:Postprocess:Visualization:Vs_20anomaly:Number_20of_20depth_20slices)=
### __Parameter name:__ Number of depth slices
**Default value:** 50 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of depth slices used to define average seismic shear wave velocities from which anomalies are calculated. Units: non-dimensional. 

(parameters:Prescribed_20Stokes_20solution)=
## **Parameters in section** Prescribed Stokes solution
(parameters:Prescribed_20Stokes_20solution:Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified 

**Pattern:** [Selection ascii data|circle|function|unspecified ] 

**Documentation:** Select one of the following models:

`ascii data': Implementation of a model in which the velocity is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `v${}_x$' , `v${}_y$' in a 2d model and  `x', `y', `z', `v${}_x$' , `v${}_y$' , `v${}_z$' in a 3d model. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the data will still be handled as Cartesian, however the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`circle': This value describes a vector field that rotates around the z-axis with constant angular velocity (i.e., with a velocity that increases with distance from the axis). The pressure is set to zero.

`function': This plugin allows to prescribe the Stokes solution for the velocity and pressure field in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}. 

(parameters:Prescribed_20Stokes_20solution:Ascii_20data_20model)=
## **Parameters in section** Prescribed Stokes solution/Ascii data model
(parameters:Prescribed_20Stokes_20solution:Ascii_20data_20model:Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/prescribed-stokes-solution/ 

**Pattern:** [DirectoryName] 

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT. 

(parameters:Prescribed_20Stokes_20solution:Ascii_20data_20model:Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d.txt 

**Pattern:** [Anything] 

**Documentation:** The file name of the model data. 

(parameters:Prescribed_20Stokes_20solution:Ascii_20data_20model:Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1. 

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)] 

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01. 

(parameters:Prescribed_20Stokes_20solution:Compaction_20pressure_20function)=
## **Parameters in section** Prescribed Stokes solution/Compaction pressure function
(parameters:Prescribed_20Stokes_20solution:Compaction_20pressure_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Prescribed_20Stokes_20solution:Compaction_20pressure_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Prescribed_20Stokes_20solution:Compaction_20pressure_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Prescribed_20Stokes_20solution:Fluid_20pressure_20function)=
## **Parameters in section** Prescribed Stokes solution/Fluid pressure function
(parameters:Prescribed_20Stokes_20solution:Fluid_20pressure_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Prescribed_20Stokes_20solution:Fluid_20pressure_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Prescribed_20Stokes_20solution:Fluid_20pressure_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Prescribed_20Stokes_20solution:Fluid_20velocity_20function)=
## **Parameters in section** Prescribed Stokes solution/Fluid velocity function
(parameters:Prescribed_20Stokes_20solution:Fluid_20velocity_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Prescribed_20Stokes_20solution:Fluid_20velocity_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Prescribed_20Stokes_20solution:Fluid_20velocity_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Prescribed_20Stokes_20solution:Pressure_20function)=
## **Parameters in section** Prescribed Stokes solution/Pressure function
(parameters:Prescribed_20Stokes_20solution:Pressure_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Prescribed_20Stokes_20solution:Pressure_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Prescribed_20Stokes_20solution:Pressure_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Prescribed_20Stokes_20solution:Velocity_20function)=
## **Parameters in section** Prescribed Stokes solution/Velocity function
(parameters:Prescribed_20Stokes_20solution:Velocity_20function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Prescribed_20Stokes_20solution:Velocity_20function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0 

**Pattern:** [Anything] 

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon. 

(parameters:Prescribed_20Stokes_20solution:Velocity_20function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t 

**Pattern:** [Anything] 

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression. 

(parameters:Solver_20parameters)=
## **Parameters in section** Solver parameters
(parameters:Solver_20parameters:Composition_20solver_20tolerance)=
### __Parameter name:__ Composition solver tolerance
**Default value:** 1e-12 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** The relative tolerance up to which the linear system for the composition system gets solved. See `Stokes solver parameters/Linear solver tolerance' for more details. 

(parameters:Solver_20parameters:Temperature_20solver_20tolerance)=
### __Parameter name:__ Temperature solver tolerance
**Default value:** 1e-12 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** The relative tolerance up to which the linear system for the temperature system gets solved. See `Stokes solver parameters/Linear solver tolerance' for more details. 

(parameters:Solver_20parameters:AMG_20parameters)=
## **Parameters in section** Solver parameters/AMG parameters
(parameters:Solver_20parameters:AMG_20parameters:AMG_20aggregation_20threshold)=
### __Parameter name:__ AMG aggregation threshold
**Default value:** 0.001 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** This threshold tells the AMG setup how the coarsening should be performed. In the AMG used by ML, all points that strongly couple with the tentative coarse-level point form one aggregate. The term strong coupling is controlled by the variable aggregation\_threshold, meaning that all elements that are not smaller than aggregation\_threshold times the diagonal element do couple strongly. The default is strongly recommended. There are indications that for the Newton solver a different value might be better. For extensive benchmarking of various settings of the AMG parameters in this section for the Stokes problem and others, see https://github.com/geodynamics/aspect/pull/234. 

(parameters:Solver_20parameters:AMG_20parameters:AMG_20output_20details)=
### __Parameter name:__ AMG output details
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Turns on extra information on the AMG solver. Note that this will generate much more output. 

(parameters:Solver_20parameters:AMG_20parameters:AMG_20smoother_20sweeps)=
### __Parameter name:__ AMG smoother sweeps
**Default value:** 2 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Determines how many sweeps of the smoother should be performed. When the flag elliptic is set to true, (which is true for ASPECT), the polynomial degree of the Chebyshev smoother is set to this value. The term sweeps refers to the number of matrix-vector products performed in the Chebyshev case. In the non-elliptic case, this parameter sets the number of SSOR relaxation sweeps for post-smoothing to be performed. The default is strongly recommended. There are indications that for the Newton solver a different value might be better. For extensive benchmarking of various settings of the AMG parameters in this section for the Stokes problem and others, see https://github.com/geodynamics/aspect/pull/234. 

(parameters:Solver_20parameters:AMG_20parameters:AMG_20smoother_20type)=
### __Parameter name:__ AMG smoother type
**Default value:** Chebyshev 

**Pattern:** [Selection Chebyshev|symmetric Gauss-Seidel ] 

**Documentation:** This parameter sets the type of smoother for the AMG. The default is strongly recommended for any normal runs with ASPECT. There are some indications that the symmetric Gauss-Seidel might be better and more stable for the Newton solver. For extensive benchmarking of various settings of the AMG parameters in this section for the Stokes problem and others, see https://github.com/geodynamics/aspect/pull/234. 

(parameters:Solver_20parameters:Advection_20solver_20parameters)=
## **Parameters in section** Solver parameters/Advection solver parameters
(parameters:Solver_20parameters:Advection_20solver_20parameters:GMRES_20solver_20restart_20length)=
### __Parameter name:__ GMRES solver restart length
**Default value:** 50 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** This is the number of iterations that define the GMRES solver restart length. Increasing this parameter makes the solver more robust and decreases the number of iterations. Be aware that increasing this number increases the memory usage of the advection solver, and makes individual iterations more expensive. 

(parameters:Solver_20parameters:Diffusion_20solver_20parameters)=
## **Parameters in section** Solver parameters/Diffusion solver parameters
(parameters:Solver_20parameters:Diffusion_20solver_20parameters:Diffusion_20length_20scale)=
### __Parameter name:__ Diffusion length scale
**Default value:** 1.e4 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Set a length scale for the diffusion of advection fields if the ``prescribed field with diffusion'' method is selected for a field. More precisely, this length scale represents the square root of the product of diffusivity and time in the diffusion equation, and controls the distance over which features are diffused. Units: \si{\meter}. 

(parameters:Solver_20parameters:Matrix_20Free)=
## **Parameters in section** Solver parameters/Matrix Free
(parameters:Solver_20parameters:Matrix_20Free:Execute_20solver_20timings)=
### __Parameter name:__ Execute solver timings
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Executes different parts of the Stokes solver repeatedly and print timing information. This is for internal benchmarking purposes: It is useful if you want to see how the solver performs. Otherwise, you don't want to enable this, since it adds additional computational cost to get the timing information. 

(parameters:Solver_20parameters:Matrix_20Free:Output_20details)=
### __Parameter name:__ Output details
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Turns on extra information for the matrix free GMG solver to be printed. 

(parameters:Solver_20parameters:Newton_20solver_20parameters)=
## **Parameters in section** Solver parameters/Newton solver parameters
(parameters:Solver_20parameters:Newton_20solver_20parameters:Max_20Newton_20line_20search_20iterations)=
### __Parameter name:__ Max Newton line search iterations
**Default value:** 5 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The maximum number of line search iterations allowed. If the criterion is not reached after this number of iterations, we apply the scaled increment even though it does not satisfy the necessary criteria and simply continue with the next Newton iteration. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Max_20pre_2dNewton_20nonlinear_20iterations)=
### __Parameter name:__ Max pre_2dNewton nonlinear iterations
**Default value:** 10 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** If the 'Nonlinear Newton solver switch tolerance' is reached before the maximal number of Picard iterations, then the solver switches to Newton solves anyway. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Maximum_20linear_20Stokes_20solver_20tolerance)=
### __Parameter name:__ Maximum linear Stokes solver tolerance
**Default value:** 0.9 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** The linear Stokes solver tolerance is dynamically chosen for the Newton solver, based on the Eisenstat Walker (1994) paper (https://doi.org/10.1137/0917003), equation 2.2. Because this value can become larger than one, we limit this value by this parameter. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Nonlinear_20Newton_20solver_20switch_20tolerance)=
### __Parameter name:__ Nonlinear Newton solver switch tolerance
**Default value:** 1e-5 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** A relative tolerance with respect to the residual of the first iteration, up to which the nonlinear Picard solver will iterate, before changing to the Newton solver. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:SPD_20safety_20factor)=
### __Parameter name:__ SPD safety factor
**Default value:** 0.9 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** When stabilizing the Newton matrix, we can encounter situations where the coefficient inside the elliptic (top-left) block becomes negative or zero. This coefficient has the form $1+x$ where $x$ can sometimes be smaller than $-1$. In this case, the top-left block of the matrix is no longer positive definite, and both preconditioners and iterative solvers may fail. To prevent this, the stabilization computes an $\alpha$ so that $1+\alpha x$ is never negative. This $\alpha$ is chosen as $1$ if $x\ge -1$, and $\alpha=-\frac 1x$ otherwise. (Note that this always leads to $0\le \alpha \le 1$.)  On the other hand, we also want to stay away from $1+\alpha x=0$, and so modify the choice of $\alpha$ to be $1$ if $x\ge -c$, and $\alpha=-\frac cx$ with a $c$ between zero and one. This way, if $c<1$, we are assured that $1-\alpha x>c$, i.e., bounded away from zero. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Stabilization_20preconditioner)=
### __Parameter name:__ Stabilization preconditioner
**Default value:** SPD 

**Pattern:** [Selection SPD|PD|symmetric|none ] 

**Documentation:** This parameters allows for the stabilization of the preconditioner. If one derives the Newton method without any modifications, the matrix created for the preconditioning is not necessarily Symmetric Positive Definite. This is problematic (see \cite{FBTGS19}). When `none' is chosen, the preconditioner is not stabilized. The `symmetric' parameters symmetrizes the matrix, and `PD' makes the matrix Positive Definite. `SPD' is the full stabilization, where the matrix is guaranteed Symmetric Positive Definite. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Stabilization_20velocity_20block)=
### __Parameter name:__ Stabilization velocity block
**Default value:** SPD 

**Pattern:** [Selection SPD|PD|symmetric|none ] 

**Documentation:** This parameters allows for the stabilization of the velocity block. If one derives the Newton method without any modifications, the matrix created for the velocity block is not necessarily Symmetric Positive Definite. This is problematic (see \cite{FBTGS19}). When `none' is chosen, the velocity block is not stabilized. The `symmetric' parameters symmetrizes the matrix, and `PD' makes the matrix Positive Definite. `SPD' is the full stabilization, where the matrix is guaranteed Symmetric Positive Definite. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Use_20Eisenstat_20Walker_20method_20for_20Picard_20iterations)=
### __Parameter name:__ Use Eisenstat Walker method for Picard iterations
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If set to true, the Picard iteration uses the Eisenstat Walker method to determine how accurately linear systems need to be solved. The Picard iteration is used, for example, in the first few iterations of the Newton method before the matrix is built including derivatives of the model, since the Picard iteration generally converges even from points where Newton's method does not. 

Once derivatives are used in a Newton method, \aspect{} always uses the Eisenstat Walker method. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Use_20Newton_20failsafe)=
### __Parameter name:__ Use Newton failsafe
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** When this parameter is true and the linear solver fails, we try again, but now with SPD stabilization for both the preconditioner and the velocity block. The SPD stabilization will remain active until the next timestep, when the default values are restored. 

(parameters:Solver_20parameters:Newton_20solver_20parameters:Use_20Newton_20residual_20scaling_20method)=
### __Parameter name:__ Use Newton residual scaling method
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** This method allows to slowly introduce the derivatives based on the improvement of the residual. If set to false, the scaling factor for the Newton derivatives is set to one immediately when switching on the Newton solver. When this is set to true, the derivatives are slowly introduced by the following equation: $\max(0.0, (1.0-(residual/switch\_initial\_residual)))$, where switch\_initial\_residual is the residual at the time when the Newton solver is switched on. 

(parameters:Solver_20parameters:Operator_20splitting_20parameters)=
## **Parameters in section** Solver parameters/Operator splitting parameters
(parameters:Solver_20parameters:Operator_20splitting_20parameters:Reaction_20time_20step)=
### __Parameter name:__ Reaction time step
**Default value:** 1000.0 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Set a time step size for computing reactions of compositional fields and the temperature field in case operator splitting is used. This is only used when the parameter ``Use operator splitting'' is set to true. The reaction time step must be greater than 0. If you want to prescribe the reaction time step only as a relative value compared to the advection time step as opposed to as an absolute value, you should use the parameter ``Reaction time steps per advection step'' and set this parameter to the same (or larger) value as the ``Maximum time step'' (which is 5.69e+300 by default). Units: Years or seconds, depending on the ``Use years in output instead of seconds'' parameter. 

(parameters:Solver_20parameters:Operator_20splitting_20parameters:Reaction_20time_20steps_20per_20advection_20step)=
### __Parameter name:__ Reaction time steps per advection step
**Default value:** 0 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** The number of reaction time steps done within one advection time step in case operator splitting is used. This is only used if the parameter ``Use operator splitting'' is set to true. If set to zero, this parameter is ignored. Otherwise, the reaction time step size is chosen according to this criterion and the ``Reaction time step'', whichever yields the smaller time step. Units: none. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters)=
## **Parameters in section** Solver parameters/Stokes solver parameters
(parameters:Solver_20parameters:Stokes_20solver_20parameters:GMRES_20solver_20restart_20length)=
### __Parameter name:__ GMRES solver restart length
**Default value:** 50 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** This is the number of iterations that define the GMRES solver restart length. Increasing this parameter helps with convergence issues arising from high localized viscosity jumps in the domain. Be aware that increasing this number increases the memory usage of the Stokes solver, and makes individual Stokes iterations more expensive. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:IDR_28s_29_20parameter)=
### __Parameter name:__ IDR_28s_29 parameter
**Default value:** 2 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** This is the sole parameter for the IDR(s) Krylov solver and will dictate the number of matrix-vector products and preconditioner applications per iteration (s+1) and the total number of temporary vectors required (5+3*s). For s=1, this method is analogous to BiCGStab. As s is increased this method is expected to converge to GMRES in terms of matrix-vector/preconditioner applications to solution. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Krylov_20method_20for_20cheap_20solver_20steps)=
### __Parameter name:__ Krylov method for cheap solver steps
**Default value:** GMRES 

**Pattern:** [Selection GMRES|IDR(s) ] 

**Documentation:** This is the Krylov method used to solve the Stokes system. Both options, GMRES and IDR(s), solve non-symmetric, indefinite systems. GMRES guarantees the residual will be reduced in each iteration while IDR(s) has no such property. On the other hand, the vector storage requirement for GMRES is dependent on the restart length and can be quite restrictive (since, for the matrix-free GMG solver, memory is dominated by these vectors) whereas IDR(s) has a short term recurrence. Note that the IDR(s) Krylov method is not available for the AMG solver since it is not a flexible method, i.e., it cannot handle a preconditioner which may change in each iteration (the AMG-based preconditioner contains a CG solve in the pressure space which may have different number of iterations each step). 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Linear_20solver_20A_20block_20tolerance)=
### __Parameter name:__ Linear solver A block tolerance
**Default value:** 1e-2 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** A relative tolerance up to which the approximate inverse of the $A$ block of the Stokes system is computed. This approximate $A$ is used in the preconditioning used in the GMRES solver. The exact definition of this block preconditioner for the Stokes equation can be found in \cite{KHB12}. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Linear_20solver_20S_20block_20tolerance)=
### __Parameter name:__ Linear solver S block tolerance
**Default value:** 1e-6 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** A relative tolerance up to which the approximate inverse of the $S$ block (i.e., the Schur complement matrix $S = BA^{-1}B^{T}$) of the Stokes system is computed. This approximate inverse of the $S$ block is used in the preconditioning used in the GMRES solver. The exact definition of this block preconditioner for the Stokes equation can be found in \cite{KHB12}. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Linear_20solver_20tolerance)=
### __Parameter name:__ Linear solver tolerance
**Default value:** 1e-7 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** A relative tolerance up to which the linear Stokes systems in each time or nonlinear step should be solved. The absolute tolerance will then be $\| M x_0 - F \| \cdot \text{tol}$, where $x_0 = (0,p_0)$ is the initial guess of the pressure, $M$ is the system matrix, $F$ is the right-hand side, and tol is the parameter specified here. We include the initial guess of the pressure to remove the dependency of the tolerance on the static pressure. A given tolerance value of 1 would mean that a zero solution vector is an acceptable solution since in that case the norm of the residual of the linear system equals the norm of the right hand side. A given tolerance of 0 would mean that the linear system has to be solved exactly, since this is the only way to obtain a zero residual.

In practice, you should choose the value of this parameter to be so that if you make it smaller the results of your simulation do not change any more (qualitatively) whereas if you make it larger, they do. For most cases, the default value should be sufficient. In fact, a tolerance of 1e-4 might be accurate enough. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Maximum_20number_20of_20expensive_20Stokes_20solver_20steps)=
### __Parameter name:__ Maximum number of expensive Stokes solver steps
**Default value:** 1000 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** This sets the maximum number of iterations used in the expensive Stokes solver. If this value is set too low for the size of the problem, the Stokes solver will not converge and return an error message pointing out that the user didn't allow a sufficiently large number of iterations for the iterative solver to converge. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Number_20of_20cheap_20Stokes_20solver_20steps)=
### __Parameter name:__ Number of cheap Stokes solver steps
**Default value:** 200 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** As explained in the paper that describes ASPECT (Kronbichler, Heister, and Bangerth, 2012, see \cite{KHB12}) we first try to solve the Stokes system in every time step using a GMRES iteration with a poor but cheap preconditioner. By default, we try whether we can converge the GMRES solver in 200 such iterations before deciding that we need a better preconditioner. This is sufficient for simple problems with variable viscosity and we never need the second phase with the more expensive preconditioner. On the other hand, for more complex problems, and in particular for problems with strongly nonlinear viscosity, the 200 cheap iterations don't actually do very much good and one might skip this part right away. In that case, this parameter can be set to zero, i.e., we immediately start with the better but more expensive preconditioner. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Stokes_20solver_20type)=
### __Parameter name:__ Stokes solver type
**Default value:** block AMG 

**Pattern:** [Selection block AMG|direct solver|block GMG ] 

**Documentation:** This is the type of solver used on the Stokes system. The block geometric multigrid solver currently has a limited implementation and therefore may trigger Asserts in the code when used. If this is the case, please switch to 'block AMG'. Additionally, the block GMG solver requires using material model averaging. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Use_20direct_20solver_20for_20Stokes_20system)=
### __Parameter name:__ Use direct solver for Stokes system
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** If set to true the linear system for the Stokes equation will be solved using Trilinos klu, otherwise an iterative Schur complement solver is used. The direct solver is only efficient for small problems. 

(parameters:Solver_20parameters:Stokes_20solver_20parameters:Use_20full_20A_20block_20as_20preconditioner)=
### __Parameter name:__ Use full A block as preconditioner
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** This parameter determines whether we use an simplified approximation of the $A$ block as preconditioner for the Stokes solver, or the full $A$ block. The simplified approximation only contains the terms that describe the coupling of identical components (plus boundary conditions) as described in \cite{KHB12}. The full block is closer to the description in \cite{rudi2017weighted}.

There is no clear way to determine which preconditioner performs better. The default value (simplified approximation) requires more outer GMRES iterations, but is faster to apply in each iteration. The full block needs less assembly time (because the block is available anyway), converges in less GMRES iterations, but requires more time per iteration. There are also differences in the amount of memory consumption between the two approaches.

The default value should be good for relatively simple models, but in particular for very strong viscosity contrasts the full $A$ block can be advantageous. 

(parameters:Temperature_20field)=
## **Parameters in section** Temperature field
(parameters:Temperature_20field:Temperature_20method)=
### __Parameter name:__ Temperature method
**Default value:** field 

**Pattern:** [Selection field|prescribed field|prescribed field with diffusion|static ] 

**Documentation:** A comma separated list denoting the solution method of the temperature field. Each entry of the list must be one of the currently implemented field types.

These choices correspond to the following methods by which the temperature field gains its values:\begin{itemize}\item ``field'': If the temperature is marked with this method, then its values are computed in each time step by solving the temperature advection-diffusion equation. In other words, this corresponds to the usual notion of a temperature. 
\item ``prescribed field'': The value of the temperature is determined in each time step from the material model. If a compositional field is marked with this method, then the value of a specific additional material model output, called the `PrescribedTemperatureOutputs' is interpolated onto the temperature. This field does not change otherwise, it is not advected with the flow. 
\item ``prescribed field with diffusion'': If the temperature field is marked this way, the value of a specific additional material model output, called the `PrescribedTemperatureOutputs' is interpolated onto the field, as in the ``prescribed field'' method. Afterwards, the field is diffused based on a solver parameter, the diffusion length scale, smoothing the field. Specifically, the field is updated by solving the equation $(I-l^2 \Delta) T_\text{smoothed} = T_\text{prescribed}$, where $l$ is the diffusion length scale. Note that this means that the amount of diffusion is independent of the time step size, and that the field is not advected with the flow.
\item ``static'': If a temperature field is marked this way, then it does not evolve at all. Its values are simply set to the initial conditions, and will then never change.\end{itemize} 

(parameters:Termination_20criteria)=
## **Parameters in section** Termination criteria
(parameters:Termination_20criteria:Checkpoint_20on_20termination)=
### __Parameter name:__ Checkpoint on termination
**Default value:** false 

**Pattern:** [Bool] 

**Documentation:** Whether to checkpoint the simulation right before termination. 

(parameters:Termination_20criteria:End_20step)=
### __Parameter name:__ End step
**Default value:** 100 

**Pattern:** [Integer range 0...2147483647 (inclusive)] 

**Documentation:** Terminate the simulation once the specified timestep has been reached. 

(parameters:Termination_20criteria:Termination_20criteria)=
### __Parameter name:__ Termination criteria
**Default value:** end time 

**Pattern:** [MultipleSelection end step|end time|steady state heat flux|steady state temperature|steady state velocity|user request|wall time ] 

**Documentation:** A comma separated list of termination criteria that will determine when the simulation should end. Whether explicitly stated or not, the ``end time'' termination criterion will always be used.The following termination criteria are available:

`end step': Terminate the simulation once the specified timestep has been reached. 

`end time': Terminate the simulation once the end time specified in the input file has been reached. Unlike all other termination criteria, this criterion is \textit{always} active, whether it has been explicitly selected or not in the input file (this is done to preserve historical behavior of \aspect{}, but it also likely does not inconvenience anyone since it is what would be selected in most cases anyway).

`steady state heat flux': A criterion that terminates the simulation when the integrated heat flux over a given list of boundaries stays within a certain range for a specified period of time.

The criterion considers the total heat flux over all boundaries listed by their boundary indicators, rather than each boundary separately. As a consequence, if the \textit{sum} of heat fluxes over individual parts of the boundary no longer changes, then this criterion recommends termination, even if the heat flux over individual parts of the boundary continues to change.

`steady state temperature': A criterion that terminates the simulation when the global integral of the temperature field stays within a certain range for a specified period of time.

`steady state velocity': A criterion that terminates the simulation when the RMS of the velocity field stays within a certain range for a specified period of time.

`user request': Terminate the simulation gracefully when a file with a specified name appears in the output directory. This allows the user to gracefully exit the simulation at any time by simply creating such a file using, for example, \texttt{touch output/terminate}. The file's location is chosen to be in the output directory, rather than in a generic location such as the ASPECT directory, so that one can run multiple simulations at the same time (which presumably write to different output directories) and can selectively terminate a particular one.

`wall time': Terminate the simulation once the wall time limit has reached. 

(parameters:Termination_20criteria:Wall_20time)=
### __Parameter name:__ Wall time
**Default value:** 24. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The wall time of the simulation. Unit: hours. 

(parameters:Termination_20criteria:Steady_20state_20heat_20flux)=
## **Parameters in section** Termination criteria/Steady state heat flux
(parameters:Termination_20criteria:Steady_20state_20heat_20flux:Boundary_20indicators)=
### __Parameter name:__ Boundary indicators
**Default value:**  

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)] 

**Documentation:** A comma separated list of names denoting those boundaries that should be taken into account for integrating the heat flux. Note that the plugin will compute the integrated heat flux over these boundaries (instead of taking them into account individually).

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model. 

(parameters:Termination_20criteria:Steady_20state_20heat_20flux:Maximum_20relative_20deviation)=
### __Parameter name:__ Maximum relative deviation
**Default value:** 0.05 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum relative deviation of the heat flux in recent simulation time for the system to be considered in steady state. If the actual deviation is smaller than this number, then the simulation will be terminated. 

(parameters:Termination_20criteria:Steady_20state_20heat_20flux:Time_20in_20steady_20state)=
### __Parameter name:__ Time in steady state
**Default value:** 1e7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum length of simulation time that the system should be in steady state before termination. Note that if the time step size is similar to or larger than this value, the termination criterion will only have very few (in the most extreme case, just two) heat flux values to check. To ensure that a larger number of time steps are included in the check for steady state, this value should be much larger than the time step size. Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Termination_20criteria:Steady_20state_20temperature)=
## **Parameters in section** Termination criteria/Steady state temperature
(parameters:Termination_20criteria:Steady_20state_20temperature:Maximum_20relative_20deviation)=
### __Parameter name:__ Maximum relative deviation
**Default value:** 0.05 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum relative deviation of the temperature in recent simulation time for the system to be considered in steady state. If the actual deviation is smaller than this number, then the simulation will be terminated. 

(parameters:Termination_20criteria:Steady_20state_20temperature:Time_20in_20steady_20state)=
### __Parameter name:__ Time in steady state
**Default value:** 1e7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum length of simulation time that the system should be in steady state before termination.Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Termination_20criteria:Steady_20state_20velocity)=
## **Parameters in section** Termination criteria/Steady state velocity
(parameters:Termination_20criteria:Steady_20state_20velocity:Maximum_20relative_20deviation)=
### __Parameter name:__ Maximum relative deviation
**Default value:** 0.05 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The maximum relative deviation of the RMS in recent simulation time for the system to be considered in steady state. If the actual deviation is smaller than this number, then the simulation will be terminated. 

(parameters:Termination_20criteria:Steady_20state_20velocity:Time_20in_20steady_20state)=
### __Parameter name:__ Time in steady state
**Default value:** 1e7 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** The minimum length of simulation time that the system should be in steady state before termination.Units: years if the 'Use years in output instead of seconds' parameter is set; seconds otherwise. 

(parameters:Termination_20criteria:User_20request)=
## **Parameters in section** Termination criteria/User request
(parameters:Termination_20criteria:User_20request:File_20name)=
### __Parameter name:__ File name
**Default value:** terminate-aspect 

**Pattern:** [FileName (Type: input)] 

**Documentation:** The name of a file that, if it exists in the output directory (whose name is also specified in the input file) will lead to termination of the simulation. The file's location is chosen to be in the output directory, rather than in a generic location such as the ASPECT directory, so that one can run multiple simulations at the same time (which presumably write to different output directories) and can selectively terminate a particular one. 

(parameters:Time_20stepping)=
## **Parameters in section** Time stepping
(parameters:Time_20stepping:List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**  

**Pattern:** [MultipleSelection conduction time step|convection time step|function|repeat on cutback ] 

**Documentation:** A comma separated list of time stepping plugins that will be used to calculate the time step size. The minimum of the  result of each plugin will be used.

The following plugins are available:

`conduction time step': This model computes the conduction time step as the minimum over all cells of $ CFL h^2 \cdot \rho C_p / k$, where k is the thermal conductivity. This plugin will always request advancing to the next time step.

`convection time step': This model computes the convection time step as $ CFL / \max \| u \| / h$ over all cells, where $u$ is the velocity and $h$ is the product of mesh size and temperature polynomial degree.

`function': This model uses a time step specified in the parameter file specified as a function of time. This plugin will always request advancing to the next time step.

`repeat on cutback': This time stepping plugin will detect a situation where the computed time step shrinks by more than a user-controlled factor. In that situation, the previous time step will be repeated with a smaller step size.
A large reduction in time step size typically happens when velocities change abruptly. Repeating the time step ensure properly resolving this event. It is useful to consider setting the "Maximum relative increase in time step" option to avoid repeatedly repeating every other time step. 

(parameters:Time_20stepping:Minimum_20time_20step_20size)=
### __Parameter name:__ Minimum time step size
**Default value:** 0. 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** Specifiy a minimum time step size (or 0 to disable). 

(parameters:Time_20stepping:Function)=
## **Parameters in section** Time stepping/Function
(parameters:Time_20stepping:Function:Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**  

**Pattern:** [Anything] 

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.) 

(parameters:Time_20stepping:Function:Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 1.0 

**Pattern:** [Anything] 

**Documentation:** Expression for the time step size as a function of 'time'. 

(parameters:Time_20stepping:Function:Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** time 

**Pattern:** [Anything] 

**Documentation:** Name for the variable representing the current time. 

(parameters:Time_20stepping:Repeat_20on_20cutback)=
## **Parameters in section** Time stepping/Repeat on cutback
(parameters:Time_20stepping:Repeat_20on_20cutback:Cut_20back_20amount)=
### __Parameter name:__ Cut back amount
**Default value:** 0.5 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A factor that controls the size of the time step when repeating. The default of 0.5 corresponds to 50\% of the original step taken. 

(parameters:Time_20stepping:Repeat_20on_20cutback:Relative_20repeat_20threshold)=
### __Parameter name:__ Relative repeat threshold
**Default value:** 0.2 

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)] 

**Documentation:** A factor that controls when a step is going to be repeated. If the newly computed step size is smaller than the last step size multiplied by this factor, the step is repeated. 

(parameters:Volume_20of_20Fluid)=
## **Parameters in section** Volume of Fluid
(parameters:Volume_20of_20Fluid:Number_20initialization_20samples)=
### __Parameter name:__ Number initialization samples
**Default value:** 3 

**Pattern:** [Integer range 1...2147483647 (inclusive)] 

**Documentation:** Number of divisions per dimension when computing the initial volume fractions.If set to the default of 3 for a 2D model, then initialization will be based on the initialization criterion at $3^2=9$ points within each cell. If the initialization based on a composition style initial condition, a larger value may be desired for better approximation of the initial fluid fractions. Smaller values will suffice in the case of level set initializations due to the presence of more information to better approximate the initial fluid fractions. 

(parameters:Volume_20of_20Fluid:Volume_20fraction_20threshold)=
### __Parameter name:__ Volume fraction threshold
**Default value:** 1e-6 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** Minimum significant volume. Fluid fractions below this value are considered to be zero. 

(parameters:Volume_20of_20Fluid:Volume_20of_20Fluid_20solver_20tolerance)=
### __Parameter name:__ Volume of Fluid solver tolerance
**Default value:** 1e-12 

**Pattern:** [Double 0...1 (inclusive)] 

**Documentation:** The relative tolerance up to which the linear system for the Volume of Fluid system gets solved. See 'Solver parameters/Composition solver tolerance' for more details. 


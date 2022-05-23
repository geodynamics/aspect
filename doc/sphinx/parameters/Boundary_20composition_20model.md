(parameters:Boundary_20composition_20model)=
# Boundary composition model


## **Subsection:** Boundary composition model


(parameters:Boundary_20composition_20model/Allow_20fixed_20composition_20on_20outflow_20boundaries)=
### __Parameter name:__ Allow fixed composition on outflow boundaries
**Default value:** false for models without melt

**Pattern:** [Selection true|false|false for models without melt ]

**Documentation:** When the composition is fixed on a given boundary as determined by the list of 'Fixed composition boundary indicators', there might be parts of the boundary where material flows out and one may want to prescribe the composition only on those parts of the boundary where there is inflow. This parameter determines if compositions are only prescribed at these inflow parts of the boundary (if false) or everywhere on a given boundary, independent of the flow direction (if true). By default, this parameter is set to false, except in models with melt transport (see below). Note that in this context, `fixed' refers to the fact that these are the boundary indicators where Dirichlet boundary conditions are applied, and does not imply that the boundary composition is time-independent.

Mathematically speaking, the compositional fields satisfy an advection equation that has no diffusion. For this equation, one can only impose Dirichlet boundary conditions (i.e., prescribe a fixed compositional field value at the boundary) at those boundaries where material flows in. This would correspond to the ``false'' setting of this parameter, which is correspondingly the default. On the other hand, on a finite dimensional discretization such as the one one obtains from the finite element method, it is possible to also prescribe values on outflow boundaries, even though this may make no physical sense. This would then correspond to the ``true'' setting of this parameter.

A warning for models with melt transport: In models with fluid flow, some compositional fields (in particular the porosity) might be transported with the fluid velocity, and would need to set the constraints based on the fluid velocity. However, this is currently not possible, because we reuse the same matrix for all compositional fields, and therefore can not use different constraints for different fields. Consequently, we set this parameter to true by default in models where melt transport is enabled. Be aware that if you change this default setting, you will not use the melt velocity, but the solid velocity to determine on which parts of the boundaries there is outflow.

(parameters:Boundary_20composition_20model/Fixed_20composition_20boundary_20indicators)=
### __Parameter name:__ Fixed composition boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries on which the composition is fixed and described by the boundary composition object selected in its own section of this input file. All boundary indicators used by the geometry but not explicitly listed here will end up with no-flux (insulating) boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed composition, but not what composition should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryComposition group, unless an existing implementation in this group already provides what you want.

(parameters:Boundary_20composition_20model/List_20of_20model_20names)=
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

(parameters:Boundary_20composition_20model/List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ]

**Documentation:** A comma-separated list of operators that will be used to append the listed composition models onto the previous models. If only one operator is given, the same operator is applied to all models.

(parameters:Boundary_20composition_20model/Model_20name)=
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

(parameters:Boundary_20composition_20model/Ascii_20data_20model)=
## **Subsection:** Boundary composition model / Ascii data model
(parameters:Boundary_20composition_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-composition/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.

(parameters:Boundary_20composition_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number.

(parameters:Boundary_20composition_20model/Ascii_20data_20model/Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Time step between following data files. Depending on the setting of the global `Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years.

(parameters:Boundary_20composition_20model/Ascii_20data_20model/Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false

**Pattern:** [Bool]

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. `Ma BP'). If this flag is set to `True' the plugin will first load the file with the number `First data file number' and decrease the file number during the model run.

(parameters:Boundary_20composition_20model/Ascii_20data_20model/First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The `First data file model time' parameter has been deactivated and will be removed in a future release. Do not use this paramter and instead provide data files starting from the model start time.

(parameters:Boundary_20composition_20model/Ascii_20data_20model/First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than `First velocity file model time'.

(parameters:Boundary_20composition_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Boundary_20composition_20model/Box)=
## **Subsection:** Boundary composition model / Box
(parameters:Boundary_20composition_20model/Box/Bottom_20composition)=
### __Parameter name:__ Bottom composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the bottom boundary (at minimal $y$-value in 2d, or minimal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box/Left_20composition)=
### __Parameter name:__ Left composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the left boundary (at minimal $x$-value). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box/Right_20composition)=
### __Parameter name:__ Right composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the right boundary (at maximal $x$-value). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box/Top_20composition)=
### __Parameter name:__ Top composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the top boundary (at maximal $y$-value in 2d, or maximal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators)=
## **Subsection:** Boundary composition model / Box with lithosphere boundary indicators
(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators/Bottom_20composition)=
### __Parameter name:__ Bottom composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the bottom boundary (at minimal $y$-value in 2d, or minimal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators/Left_20composition)=
### __Parameter name:__ Left composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the left boundary (at minimal $x$-value). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators/Left_20composition_20lithosphere)=
### __Parameter name:__ Left composition lithosphere
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the left boundary (at minimal $x$-value). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators/Right_20composition)=
### __Parameter name:__ Right composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the right boundary (at maximal $x$-value). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators/Right_20composition_20lithosphere)=
### __Parameter name:__ Right composition lithosphere
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the right boundary (at maximal $x$-value). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Box_20with_20lithosphere_20boundary_20indicators/Top_20composition)=
### __Parameter name:__ Top composition
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of composition boundary values at the top boundary (at maximal $y$-value in 2d, or maximal $z$-value in 3d). This list must have as many entries as there are compositional fields. Units: none.

(parameters:Boundary_20composition_20model/Function)=
## **Subsection:** Boundary composition model / Function
(parameters:Boundary_20composition_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are 'cartesian', 'spherical', and 'depth'. 'spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. 'depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Boundary_20composition_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.)

(parameters:Boundary_20composition_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Boundary_20composition_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression.

(parameters:Boundary_20composition_20model/Initial_20composition)=
## **Subsection:** Boundary composition model / Initial composition
(parameters:Boundary_20composition_20model/Initial_20composition/Maximal_20composition)=
### __Parameter name:__ Maximal composition
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximal composition. Units: none.

(parameters:Boundary_20composition_20model/Initial_20composition/Minimal_20composition)=
### __Parameter name:__ Minimal composition
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimal composition. Units: none.

(parameters:Boundary_20composition_20model/Spherical_20constant)=
## **Subsection:** Boundary composition model / Spherical constant
(parameters:Boundary_20composition_20model/Spherical_20constant/Inner_20composition)=
### __Parameter name:__ Inner composition
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Composition at the inner boundary (core mantle boundary). Units: none.

(parameters:Boundary_20composition_20model/Spherical_20constant/Outer_20composition)=
### __Parameter name:__ Outer composition
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Composition at the outer boundary (lithosphere water/air). For a spherical geometry model, this is the only boundary. Units: none.

(parameters:Boundary_20traction_20model)=
# Boundary traction model


## **Parameters in section** Boundary traction model


(parameters:Boundary_20traction_20model/Prescribed_20traction_20boundary_20indicators)=
### __Parameter name:__ Prescribed traction boundary indicators
**Default value:**

**Pattern:** [Map of <[Anything]>:<[Selection ascii data|function|initial lithostatic pressure|zero traction ]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting those boundaries on which a traction force is prescribed, i.e., where known external forces act, resulting in an unknown velocity. This is often used to model ``open'' boundaries where we only know the pressure. This pressure then produces a force that is normal to the boundary and proportional to the pressure.

The format of valid entries for this parameter is that of a map given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where each key must be a valid boundary indicator (which is either an integer or the symbolic name the geometry model in use may have provided for this part of the boundary) and each value must be one of the currently implemented boundary traction models. ``selector'' is an optional string given as a subset of the letters `xyz' that allows you to apply the boundary conditions only to the components listed. As an example, '1 y: function' applies the type `function' to the y component on boundary 1. Without a selector it will affect all components of the traction.

(parameters:Boundary_20traction_20model/Ascii_20data_20model)=
## **Parameters in section** Boundary traction model/Ascii data model
(parameters:Boundary_20traction_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-traction/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Time step between following data files. Depending on the setting of the global `Use years in output instead of seconds' flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false

**Pattern:** [Bool]

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. `Ma BP'). If this flag is set to `True' the plugin will first load the file with the number `First data file number' and decrease the file number during the model run.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The `First data file model time' parameter has been deactivated and will be removed in a future release. Do not use this paramter and instead provide data files starting from the model start time.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than `First velocity file model time'.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Boundary_20traction_20model/Function)=
## **Parameters in section** Boundary traction model/Function
(parameters:Boundary_20traction_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Boundary_20traction_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.)

(parameters:Boundary_20traction_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Boundary_20traction_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression.

(parameters:Boundary_20traction_20model/Initial_20lithostatic_20pressure)=
## **Parameters in section** Boundary traction model/Initial lithostatic pressure
(parameters:Boundary_20traction_20model/Initial_20lithostatic_20pressure/Number_20of_20integration_20points)=
### __Parameter name:__ Number of integration points
**Default value:** 1000

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of integration points over which we integrate the lithostatic pressure downwards.

(parameters:Boundary_20traction_20model/Initial_20lithostatic_20pressure/Representative_20point)=
### __Parameter name:__ Representative point
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The point where the pressure profile will be calculated. Cartesian coordinates $(x,y,z)$ when geometry is a box, otherwise enter radius, longitude, and in 3D latitude. Note that the coordinate related to the depth ($y$ in 2D cartesian, $z$ in 3D cartesian and radius in spherical coordinates) is not used. Units: \si{\meter} or degrees.

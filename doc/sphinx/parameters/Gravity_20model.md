(parameters:Gravity_20model)=
# **Gravity model**


## **Parameters in section** Gravity model


(parameters:Gravity_20model/Model_20name)=
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

(parameters:Gravity_20model/Ascii_20data_20model)=
## **Parameters in section** Gravity model/Ascii data model
(parameters:Gravity_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/gravity-model/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.

(parameters:Gravity_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** prem.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Gravity_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Gravity_20model/Function)=
## **Parameters in section** Gravity model/Function
(parameters:Gravity_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.)

(parameters:Gravity_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Gravity_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression.

(parameters:Gravity_20model/Radial_20constant)=
## **Parameters in section** Gravity model/Radial constant
(parameters:Gravity_20model/Radial_20constant/Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 9.81

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Magnitude of the gravity vector in $m/s^2$. For positive values the direction is radially inward towards the center of the earth.

(parameters:Gravity_20model/Radial_20linear)=
## **Parameters in section** Gravity model/Radial linear
(parameters:Gravity_20model/Radial_20linear/Magnitude_20at_20bottom)=
### __Parameter name:__ Magnitude at bottom
**Default value:** 10.7

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Magnitude of the radial gravity vector at the bottom of the domain. `Bottom' means themaximum depth in the chosen geometry, and for example represents the core-mantle boundary in the case of the `spherical shell' geometry model, and the center in the case of the `sphere' geometry model. Units: \si{\meter\per\second\squared}.

(parameters:Gravity_20model/Radial_20linear/Magnitude_20at_20surface)=
### __Parameter name:__ Magnitude at surface
**Default value:** 9.8

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Magnitude of the radial gravity vector at the surface of the domain. Units: \si{\meter\per\second\squared}.

(parameters:Gravity_20model/Vertical)=
## **Parameters in section** Gravity model/Vertical
(parameters:Gravity_20model/Vertical/Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Value of the gravity vector in $m/s^2$ directed along negative y (2D) or z (3D) axis (if the magnitude is positive.

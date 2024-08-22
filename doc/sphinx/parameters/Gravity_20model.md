(parameters:Gravity_20model)=
# Gravity model


## **Subsection:** Gravity model


(parameters:Gravity_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified

**Pattern:** [Selection ascii data|function|radial constant|radial earth-like|radial linear|vertical|unspecified ]

**Documentation:** Select one of the following models:

&lsquo;ascii data&rsquo;: Gravity is read from a file that describes the reference state. The default profile follows the preliminary reference Earth model (PREM, Dziewonski and Anderson, 1981). Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of points in the reference state as for example &lsquo;# POINTS: 3&rsquo;. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide a column named &lsquo;gravity&rsquo;. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

&lsquo;function&rsquo;: Gravity is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Gravity model|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;radial constant&rsquo;: A gravity model in which the gravity has a constant magnitude and the direction is radial (pointing inward if the value is positive). The magnitude is read from the parameter file in subsection &rsquo;Radial constant&rsquo;.

&lsquo;radial earth-like&rsquo;: This plugin has been removed due to its misleading name. The included profile was hard-coded and was less earth-like than the &lsquo;ascii data&rsquo; plugin, which uses the profile of the Preliminary Reference Earth Model (PREM). Use &lsquo;ascii data&rsquo; instead of &lsquo;radial earth-like&rsquo;.

&lsquo;radial linear&rsquo;: A gravity model which is radial (pointing inward if the gravity is positive) and the magnitude changes linearly with depth. The magnitude of gravity at the surface and bottom is read from the input file in a section &ldquo;Gravity model/Radial linear&rdquo;.

&lsquo;vertical&rsquo;: A gravity model in which the gravity direction is vertical (pointing downward for positive values) and at a constant magnitude by default equal to one.

(parameters:Gravity_20model/Ascii_20data_20model)=
## **Subsection:** Gravity model / Ascii data model
(parameters:Gravity_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/gravity-model/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

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
## **Subsection:** Gravity model / Function
(parameters:Gravity_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Gravity_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Gravity_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Gravity_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Gravity_20model/Radial_20constant)=
## **Subsection:** Gravity model / Radial constant
(parameters:Gravity_20model/Radial_20constant/Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 9.81

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Magnitude of the gravity vector in $m/s^2$. For positive values the direction is radially inward towards the center of the earth.

(parameters:Gravity_20model/Radial_20linear)=
## **Subsection:** Gravity model / Radial linear
(parameters:Gravity_20model/Radial_20linear/Magnitude_20at_20bottom)=
### __Parameter name:__ Magnitude at bottom
**Default value:** 10.7

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Magnitude of the radial gravity vector at the bottom of the domain. &lsquo;Bottom&rsquo; means themaximum depth in the chosen geometry, and for example represents the core-mantle boundary in the case of the &lsquo;spherical shell&rsquo; geometry model, and the center in the case of the &lsquo;sphere&rsquo; geometry model. Units: \si{\meter\per\second\squared}.

(parameters:Gravity_20model/Radial_20linear/Magnitude_20at_20surface)=
### __Parameter name:__ Magnitude at surface
**Default value:** 9.8

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Magnitude of the radial gravity vector at the surface of the domain. Units: \si{\meter\per\second\squared}.

(parameters:Gravity_20model/Vertical)=
## **Subsection:** Gravity model / Vertical
(parameters:Gravity_20model/Vertical/Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Value of the gravity vector in $m/s^2$ directed along negative y (2d) or z (3d) axis (if the magnitude is positive.

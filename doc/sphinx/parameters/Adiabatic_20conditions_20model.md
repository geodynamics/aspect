(parameters:Adiabatic_20conditions_20model)=
# Adiabatic conditions model


## **Subsection:** Adiabatic conditions model


(parameters:Adiabatic_20conditions_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** compute profile

**Pattern:** [Selection ascii data|compute entropy profile|compute profile|function ]

**Documentation:** Select one of the following models:

&lsquo;ascii data&rsquo;: A model in which the adiabatic profile is read from a file that describes the reference state. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of points in the reference state as for example &lsquo;# POINTS: 3&rsquo;. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide columns named &lsquo;temperature&rsquo;, &lsquo;pressure&rsquo;, and &lsquo;density&rsquo;. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

&lsquo;compute entropy profile&rsquo;: A model in which the adiabatic profile is calculated by solving the hydrostatic equations for pressure and entropy in depth. Of course the entropy along an adiabat is constant. This plugin requires the material model to provide an additional output object of type PrescribedTemperatureOutputs. It also requires that there is a compositional field of type &rsquo;entropy&rsquo; that represents the entropy of the material.

&lsquo;compute profile&rsquo;: A model in which the adiabatic profile is calculated by solving the hydrostatic equations for pressure and temperature in depth. The gravity is assumed to be in depth direction and the composition is either given by the initial composition at reference points or computed as a reference depth-function. All material parameters are computed by the material model plugin. The surface conditions are either constant or changing over time as prescribed by a user-provided function.

&lsquo;function&rsquo;: A model in which the adiabatic profile is specified by a user defined function. The supplied function has to contain temperature, pressure, and density as a function of depth in this order.

(parameters:Adiabatic_20conditions_20model/Ascii_20data_20model)=
## **Subsection:** Adiabatic conditions model / Ascii data model
(parameters:Adiabatic_20conditions_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/tests/adiabatic-conditions/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Adiabatic_20conditions_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:**

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Adiabatic_20conditions_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Adiabatic_20conditions_20model/Compute_20entropy_20profile)=
## **Subsection:** Adiabatic conditions model / Compute entropy profile
(parameters:Adiabatic_20conditions_20model/Compute_20entropy_20profile/Number_20of_20points)=
### __Parameter name:__ Number of points
**Default value:** 1000

**Pattern:** [Integer range 5...2147483647 (inclusive)]

**Documentation:** The number of points we use to compute the adiabatic profile. The higher the number of points, the more accurate the downward integration from the adiabatic surface conditions will be.

(parameters:Adiabatic_20conditions_20model/Compute_20entropy_20profile/Surface_20entropy)=
### __Parameter name:__ Surface entropy
**Default value:** 0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The surface entropy for the profile.

(parameters:Adiabatic_20conditions_20model/Compute_20profile)=
## **Subsection:** Adiabatic conditions model / Compute profile
(parameters:Adiabatic_20conditions_20model/Compute_20profile/Composition_20reference_20profile)=
### __Parameter name:__ Composition reference profile
**Default value:** initial composition

**Pattern:** [Selection initial composition|function ]

**Documentation:** Select how the reference profile for composition is computed. This profile is used to evaluate the material model, when computing the pressure and temperature profile.

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Number_20of_20points)=
### __Parameter name:__ Number of points
**Default value:** 1000

**Pattern:** [Integer range 5...2147483647 (inclusive)]

**Documentation:** The number of points we use to compute the adiabatic profile. The higher the number of points, the more accurate the downward integration from the adiabatic surface temperature will be.

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Use_20surface_20condition_20function)=
### __Parameter name:__ Use surface condition function
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use the &rsquo;Surface condition function&rsquo; to determine surface conditions, or the &rsquo;Adiabatic surface temperature&rsquo; and &rsquo;Surface pressure&rsquo; parameters. If this is set to true the reference profile is updated every timestep. The function expression of the function should be independent of space, but can depend on time &rsquo;t&rsquo;. The function must return two components, the first one being reference surface pressure, the second one being reference surface temperature.

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Surface_20condition_20function)=
## **Subsection:** Adiabatic conditions model / Compute profile / Surface condition function
(parameters:Adiabatic_20conditions_20model/Compute_20profile/Surface_20condition_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Surface_20condition_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Adiabatic_20conditions_20model/Compute_20profile/Surface_20condition_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Adiabatic_20conditions_20model/Function)=
## **Subsection:** Adiabatic conditions model / Function
(parameters:Adiabatic_20conditions_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Adiabatic_20conditions_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0.0; 0.0; 1.0

**Pattern:** [Anything]

**Documentation:** Expression for the adiabatic temperature, pressure, and density separated by semicolons as a function of &lsquo;depth&rsquo;.

(parameters:Adiabatic_20conditions_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** depth

**Pattern:** [Anything]

**Documentation:**

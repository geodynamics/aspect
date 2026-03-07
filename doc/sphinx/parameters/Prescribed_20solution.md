(parameters:Prescribed_20solution)=
# Prescribed solution


## **Subsection:** Prescribed solution


::::{dropdown} __Parameter:__ {ref}`List of model names<parameters:Prescribed_20solution/List_20of_20model_20names>`
:name: parameters:Prescribed_20solution/List_20of_20model_20names
**Default value:**

**Pattern:** [MultipleSelection initial temperature|temperature function|velocity function ]

**Documentation:** A comma-separated list of prescribed solution models that will be used to compute the solution in certain regions. These plugins are loaded in the order given, and are combined via the operators listed in &rsquo;List of model operators&rsquo;.

The following prescribed solution models are available:

&lsquo;initial temperature&rsquo;: Prescribe the temperature in a selected region using the active initial temperature model. The selected region is defined through an indicator function. At locations where the indicator value is greater than 0.5, the temperature is constrained to the initial temperature evaluated at that position.

&lsquo;temperature function&rsquo;: Prescribe the temperature in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library.

&lsquo;velocity function&rsquo;: Prescribe the velocity in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.
::::

(parameters:Prescribed_20solution/Initial_20temperature)=
## **Subsection:** Prescribed solution / Initial temperature
::::{dropdown} __Parameter:__ {ref}`Coordinate system<parameters:Prescribed_20solution/Initial_20temperature/Coordinate_20system>`
:name: parameters:Prescribed_20solution/Initial_20temperature/Coordinate_20system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the indicator function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;.
::::

(parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function)=
## **Subsection:** Prescribed solution / Initial temperature / Indicator function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function/Function_20constants>`
:name: parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function/Function_20expression>`
:name: parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function/Function_20expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function/Variable_20names>`
:name: parameters:Prescribed_20solution/Initial_20temperature/Indicator_20function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Prescribed_20solution/Temperature_20function)=
## **Subsection:** Prescribed solution / Temperature function
::::{dropdown} __Parameter:__ {ref}`Coordinate system<parameters:Prescribed_20solution/Temperature_20function/Coordinate_20system>`
:name: parameters:Prescribed_20solution/Temperature_20function/Coordinate_20system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;.
::::

(parameters:Prescribed_20solution/Temperature_20function/Function)=
## **Subsection:** Prescribed solution / Temperature function / Function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Prescribed_20solution/Temperature_20function/Function/Function_20constants>`
:name: parameters:Prescribed_20solution/Temperature_20function/Function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Prescribed_20solution/Temperature_20function/Function/Function_20expression>`
:name: parameters:Prescribed_20solution/Temperature_20function/Function/Function_20expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Prescribed_20solution/Temperature_20function/Function/Variable_20names>`
:name: parameters:Prescribed_20solution/Temperature_20function/Function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Prescribed_20solution/Temperature_20function/Indicator_20function)=
## **Subsection:** Prescribed solution / Temperature function / Indicator function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Prescribed_20solution/Temperature_20function/Indicator_20function/Function_20constants>`
:name: parameters:Prescribed_20solution/Temperature_20function/Indicator_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Prescribed_20solution/Temperature_20function/Indicator_20function/Function_20expression>`
:name: parameters:Prescribed_20solution/Temperature_20function/Indicator_20function/Function_20expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Prescribed_20solution/Temperature_20function/Indicator_20function/Variable_20names>`
:name: parameters:Prescribed_20solution/Temperature_20function/Indicator_20function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Prescribed_20solution/Velocity_20function)=
## **Subsection:** Prescribed solution / Velocity function
::::{dropdown} __Parameter:__ {ref}`Coordinate system<parameters:Prescribed_20solution/Velocity_20function/Coordinate_20system>`
:name: parameters:Prescribed_20solution/Velocity_20function/Coordinate_20system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.
::::

::::{dropdown} __Parameter:__ {ref}`Use spherical unit vectors<parameters:Prescribed_20solution/Velocity_20function/Use_20spherical_20unit_20vectors>`
:name: parameters:Prescribed_20solution/Velocity_20function/Use_20spherical_20unit_20vectors
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Specify velocity as $r$, $\phi$, and $\theta$ components instead of $x$, $y$, and $z$. Positive velocities point up, east, and north (in 3d) or out and clockwise (in 2d). This setting only makes sense for spherical geometries.
::::

(parameters:Prescribed_20solution/Velocity_20function/Function)=
## **Subsection:** Prescribed solution / Velocity function / Function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Prescribed_20solution/Velocity_20function/Function/Function_20constants>`
:name: parameters:Prescribed_20solution/Velocity_20function/Function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Prescribed_20solution/Velocity_20function/Function/Function_20expression>`
:name: parameters:Prescribed_20solution/Velocity_20function/Function/Function_20expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Prescribed_20solution/Velocity_20function/Function/Variable_20names>`
:name: parameters:Prescribed_20solution/Velocity_20function/Function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Prescribed_20solution/Velocity_20function/Indicator_20function)=
## **Subsection:** Prescribed solution / Velocity function / Indicator function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Prescribed_20solution/Velocity_20function/Indicator_20function/Function_20constants>`
:name: parameters:Prescribed_20solution/Velocity_20function/Indicator_20function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Prescribed_20solution/Velocity_20function/Indicator_20function/Function_20expression>`
:name: parameters:Prescribed_20solution/Velocity_20function/Indicator_20function/Function_20expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Prescribed_20solution/Velocity_20function/Indicator_20function/Variable_20names>`
:name: parameters:Prescribed_20solution/Velocity_20function/Indicator_20function/Variable_20names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.
::::

(parameters:Boundary_20heat_20flux_20model)=
# Boundary heat flux model


## **Subsection:** Boundary heat flux model


(parameters:Boundary_20heat_20flux_20model/Fixed_20heat_20flux_20boundary_20indicators)=
### __Parameter name:__ Fixed heat flux boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries on which the heat flux is fixed and described by the boundary heat flux object selected in the &rsquo;Model name&rsquo; parameter. All boundary indicators used by the geometry but not explicitly listed here or in the list of &rsquo;Fixed temperature boundary indicators&rsquo; in the &rsquo;Boundary temperature model&rsquo; will end up with no-flux (insulating) boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed heat flux, but not what heat flux should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryHeatFlux group, unless an existing implementation in this group already provides what you want.

(parameters:Boundary_20heat_20flux_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** function

**Pattern:** [Selection function ]

**Documentation:** Select one of the following plugins:

&lsquo;function&rsquo;: Implementation of a model in which the boundary heat flux is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Boundary heat flux model|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

The formula you describe in the mentioned section is a scalar value for the heat flux that is assumed to be the flux normal to the boundary, and that has the unit W/(m$^2$) (in 3d) or W/m (in 2d). Negative fluxes are interpreted as the flow of heat into the domain, and positive fluxes are interpreted as heat flowing out of the domain.

The symbol $t$ indicating time that may appear in the formulas for the prescribed heat flux is interpreted as having units seconds unless the global parameter &ldquo;Use years in output instead of seconds&rdquo; has been set.

(parameters:Boundary_20heat_20flux_20model/Function)=
## **Subsection:** Boundary heat flux model / Function
(parameters:Boundary_20heat_20flux_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Boundary_20heat_20flux_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Boundary_20heat_20flux_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Boundary_20heat_20flux_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Boundary_20convective_20heating_20model)=
# Boundary convective heating model


## **Subsection:** Boundary convective heating model


(parameters:Boundary_20convective_20heating_20model/Convective_20heating_20boundary_20indicators)=
### __Parameter name:__ Convective heating boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries on which convective heating should be applied, i.e., where we want to use Robin boundary conditions. On these boundaries, we can prescribe both a heat flux and a boundary temperature, and the heat transfer coefficient input parameter determines the weighting of the two conditions. The boundary temperature is described by the boundary temperature object selected in the &rsquo;List of boundary temperature model names&rsquo; parameter, and the boundary heat flux is prescribed by the boundary heat flux object selected in the &rsquo;List of boundary heat flux model names&rsquo; parameter. All boundary indicators used by the geometry but not explicitly listed here will end up with no-flux (insulating) boundary conditions, or, if they are listed in the &rsquo;Fixed heat flux boundary indicators&rsquo;, with Neumann boundary conditions, or, if they are listed in the &rsquo;Fixed temperature boundary indicators&rsquo;, with Dirichlet boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed convective flux (i.e. Robin boundarry condictions), but not what conditions should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryConvectiveHeating group, unless an existing implementation in this group already provides what you want.

(parameters:Boundary_20convective_20heating_20model/List_20of_20boundary_20heat_20flux_20model_20names)=
### __Parameter name:__ List of boundary heat flux model names
**Default value:**

**Pattern:** [MultipleSelection function ]

**Documentation:** A comma-separated list of boundary heat flux models that will be used to determine the temperature boundary conditions. At the moment, this list can only have one entry.

The following boundary heat flux models are available:

&lsquo;function&rsquo;: Implementation of a model in which the boundary heat flux is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Boundary heat flux model|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

The formula you describe in the mentioned section is a scalar value for the heat flux that is assumed to be the flux normal to the boundary, and that has the unit W/(m$^2$) (in 3d) or W/m (in 2d). Negative fluxes are interpreted as the flow of heat into the domain, and positive fluxes are interpreted as heat flowing out of the domain.

The symbol $t$ indicating time that may appear in the formulas for the prescribed heat flux is interpreted as having units seconds unless the global parameter &ldquo;Use years in output instead of seconds&rdquo; has been set.

(parameters:Boundary_20convective_20heating_20model/List_20of_20boundary_20temperature_20model_20names)=
### __Parameter name:__ List of boundary temperature model names
**Default value:**

**Pattern:** [MultipleSelection ascii data|box|box with lithosphere boundary indicators|constant|dynamic core|function|initial temperature|spherical constant ]

**Documentation:** A comma-separated list of boundary temperature models that will be used to determine the temperature boundary conditions. At the moment, this list can only have one entry.

The following boundary temperature models are available:

&lsquo;ascii data&rsquo;: Implementation of a model in which the boundary data is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;Temperature [K]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Temperature [K]&rsquo; in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;box&rsquo;: A model in which the temperature is chosen constant on the sides of a box which are selected by the parameters Left/Right/Top/Bottom/Front/Back temperature

&lsquo;box with lithosphere boundary indicators&rsquo;: A model in which the temperature is chosen constant on all the sides of a box. Additional boundary indicators are added to the lithospheric parts of the vertical boundaries. This model is to be used with the &rsquo;Two Merged Boxes&rsquo; Geometry Model.

&lsquo;constant&rsquo;: A model in which the temperature is chosen constant on a given boundary indicator.  Parameters are read from the subsection &rsquo;Constant&rsquo;.

&lsquo;dynamic core&rsquo;: This is a boundary temperature model working only with spherical shell geometry and core statistics postprocessor. The temperature at the top is constant, and the core mantle boundary temperature is dynamically evolving through time by calculating the heat flux into the core and solving the core energy balance. The formulation is mainly following {cite}`NPB+04`, and the plugin is used in Zhang et al. [2016]. The energy of core cooling and freeing of the inner core is included in the plugin. However, current plugin can not deal with the energy balance if the core is in the &lsquo;snowing core&rsquo; regime (i.e., the core solidifies from the top instead of bottom).

&lsquo;function&rsquo;: Implementation of a model in which the boundary temperature is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Boundary temperature model|Function&rdquo;.

Since the symbol $t$ indicating time may appear in the formulas for the prescribed temperatures, it is interpreted as having units seconds unless the global input parameter &ldquo;Use years in output instead of seconds&rdquo; is set, in which case we interpret the formula expressions as having units year.

Because this class simply takes what the function calculates, this class can not know certain pieces of information such as the minimal and maximal temperature on the boundary. For operations that require this, for example in post-processing, this boundary temperature model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section &ldquo;Boundary temperature model/Initial temperature&rdquo;.

The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;initial temperature&rsquo;: A model in which the temperature at the boundary is chosen to be the same as given in the initial conditions.

Because this class simply takes what the initial temperature had described, this class can not know certain pieces of information such as the minimal and maximal temperature on the boundary. For operations that require this, for example in post-processing, this boundary temperature model must therefore be told what the minimal and maximal values on the boundary are. This is done using parameters set in section &ldquo;Boundary temperature model/Initial temperature&rdquo;.

&lsquo;spherical constant&rsquo;: A model in which the temperature is chosen constant on the inner and outer boundaries of a spherical shell, ellipsoidal chunk or chunk. Parameters are read from subsection &rsquo;Spherical constant&rsquo;.

(parameters:Boundary_20convective_20heating_20model/List_20of_20heat_20transfer_20coefficient_20model_20names)=
### __Parameter name:__ List of heat transfer coefficient model names
**Default value:**

**Pattern:** [MultipleSelection function ]

**Documentation:** A comma-separated list of boundary convective heating models that will be used to determine the heat transfer coefficient across the boundary. The heat transfer coefficient characterises the heat exchange between the solid model interior and an adjacent fluid. In the context of a Robin boundary condition, the heat transfer coefficient governs the strength of the convective coupling: For heat transfer coefficient --> zero, the boundary approaches insulating (Neumann) behaviour; For heat transfer coefficient --> infinity, the boundary approaches a prescribed-temperature (Dirichlet) condition.The unit of the heat transfer coefficient is \si{\watt\per\meter\squared\per\kelvin}.At the moment, this list can only have one entry.

The following heat transfer coefficient models are available:

&lsquo;function&rsquo;: Implementation of a model in which the boundary heat transfer coefficient is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Boundary convective heating model|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

The formula you describe in the mentioned section is a scalar value for the heat transfer coefficient across the boundary that has the unit W/(m$^2$)/K (in 3d) or W/m/K (in 2d). The heat flux across the boundary is then computed as the sum of a term that is proportional to the product of the heat transfer coefficient and the difference between the temperature given by the boundary temperature model and the current temperature at the boundary and a term that prescribes a fixed heat flux across the boundary.

The symbol $t$ indicating time that may appear in the formulas for the prescribed heat flux is interpreted as having units seconds unless the global parameter &ldquo;Use years in output instead of seconds&rdquo; has been set.

(parameters:Boundary_20convective_20heating_20model/Function)=
## **Subsection:** Boundary convective heating model / Function
(parameters:Boundary_20convective_20heating_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Boundary_20convective_20heating_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Boundary_20convective_20heating_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Boundary_20convective_20heating_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

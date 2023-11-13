(parameters:Boundary_20temperature_20model)=
# Boundary temperature model


## **Subsection:** Boundary temperature model


(parameters:Boundary_20temperature_20model/Allow_20fixed_20temperature_20on_20outflow_20boundaries)=
### __Parameter name:__ Allow fixed temperature on outflow boundaries
**Default value:** true

**Pattern:** [Bool]

**Documentation:** When the temperature is fixed on a given boundary as determined by the list of &rsquo;Fixed temperature boundary indicators&rsquo;, there might be parts of the boundary where material flows out and one may want to prescribe the temperature only on the parts of the boundary where there is inflow. This parameter determines if temperatures are only prescribed at these inflow parts of the boundary (if false) or everywhere on a given boundary, independent of the flow direction (if true).Note that in this context, &lsquo;fixed&rsquo; refers to the fact that these are the boundary indicators where Dirichlet boundary conditions are applied, and does not imply that the boundary temperature is time-independent.

Mathematically speaking, the temperature satisfies an advection-diffusion equation. For this type of equation, one can prescribe the temperature even on outflow boundaries as long as the diffusion coefficient is nonzero. This would correspond to the &ldquo;true&rdquo; setting of this parameter, which is correspondingly the default. In practice, however, this would only make physical sense if the diffusion coefficient is actually quite large to prevent the creation of a boundary layer. In addition, if there is no diffusion, one can only impose Dirichlet boundary conditions (i.e., prescribe a fixed temperature value at the boundary) at those boundaries where material flows in. This would correspond to the &ldquo;false&rdquo; setting of this parameter.

(parameters:Boundary_20temperature_20model/Fixed_20temperature_20boundary_20indicators)=
### __Parameter name:__ Fixed temperature boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries on which the temperature is fixed and described by the boundary temperature object selected in the &rsquo;List of model names&rsquo; parameter. All boundary indicators used by the geometry but not explicitly listed here will end up with no-flux (insulating) boundary conditions, or, if they are listed in the &rsquo;Fixed heat flux boundary indicators&rsquo;, with Neumann boundary conditions.

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

This parameter only describes which boundaries have a fixed temperature, but not what temperature should hold on these boundaries. The latter piece of information needs to be implemented in a plugin in the BoundaryTemperature group, unless an existing implementation in this group already provides what you want.

(parameters:Boundary_20temperature_20model/List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**

**Pattern:** [MultipleSelection ascii data|box|box with lithosphere boundary indicators|constant|dynamic core|function|initial temperature|spherical constant ]

**Documentation:** A comma-separated list of boundary temperature models that will be used to initialize the temperature. These plugins are loaded in the order given, and modify the existing temperature field via the operators listed in &rsquo;List of model operators&rsquo;.

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

(parameters:Boundary_20temperature_20model/List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ]

**Documentation:** A comma-separated list of operators that will be used to append the listed temperature models onto the previous models. If only one operator is given, the same operator is applied to all models.

(parameters:Boundary_20temperature_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified

**Pattern:** [Selection ascii data|box|box with lithosphere boundary indicators|constant|dynamic core|function|initial temperature|spherical constant|unspecified ]

**Documentation:** Select one of the following models:

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

**Warning**: This parameter provides an old and deprecated way of specifying boundary temperature models and shouldn&rsquo;t be used. Please use &rsquo;List of model names&rsquo; instead.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model)=
## **Subsection:** Boundary temperature model / Ascii data model
(parameters:Boundary_20temperature_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-temperature/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model/Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Time step between following data files. Depending on the setting of the global &lsquo;Use years in output instead of seconds&rsquo; flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model/Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false

**Pattern:** [Bool]

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. &lsquo;Ma BP&rsquo;). If this flag is set to &lsquo;True&rsquo; the plugin will first load the file with the number &lsquo;First data file number&rsquo; and decrease the file number during the model run.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model/First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The &lsquo;First data file model time&rsquo; parameter has been deactivated and will be removed in a future release. Do not use this parameter and instead provide data files starting from the model start time.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model/First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than &lsquo;First velocity file model time&rsquo;.

(parameters:Boundary_20temperature_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Boundary_20temperature_20model/Box)=
## **Subsection:** Boundary temperature model / Box
(parameters:Boundary_20temperature_20model/Box/Bottom_20temperature)=
### __Parameter name:__ Bottom temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the bottom boundary (at minimal $z$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box/Left_20temperature)=
### __Parameter name:__ Left temperature
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the left boundary (at minimal $x$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box/Right_20temperature)=
### __Parameter name:__ Right temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the right boundary (at maximal $x$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box/Top_20temperature)=
### __Parameter name:__ Top temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the top boundary (at maximal $x$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators)=
## **Subsection:** Boundary temperature model / Box with lithosphere boundary indicators
(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators/Bottom_20temperature)=
### __Parameter name:__ Bottom temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the bottom boundary (at minimal $z$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators/Left_20temperature)=
### __Parameter name:__ Left temperature
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the left boundary (at minimal $x$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators/Left_20temperature_20lithosphere)=
### __Parameter name:__ Left temperature lithosphere
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the additional left lithosphere boundary (specified by user in Geometry Model). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators/Right_20temperature)=
### __Parameter name:__ Right temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the right boundary (at maximal $x$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators/Right_20temperature_20lithosphere)=
### __Parameter name:__ Right temperature lithosphere
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the additional right lithosphere boundary (specified by user in Geometry Model). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Box_20with_20lithosphere_20boundary_20indicators/Top_20temperature)=
### __Parameter name:__ Top temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the top boundary (at maximal $x$-value). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Constant)=
## **Subsection:** Boundary temperature model / Constant
(parameters:Boundary_20temperature_20model/Constant/Boundary_20indicator_20to_20temperature_20mappings)=
### __Parameter name:__ Boundary indicator to temperature mappings
**Default value:**

**Pattern:** [Map of <[Anything]>:<[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of mappings between boundary indicators and the temperature associated with the boundary indicators. The format for this list is &ldquo;indicator1 : value1, indicator2 : value2, ...&rdquo;, where each indicator is a valid boundary indicator (either a number or the symbolic name of a boundary as provided by the geometry model) and each value is the temperature of that boundary.

(parameters:Boundary_20temperature_20model/Dynamic_20core)=
## **Subsection:** Boundary temperature model / Dynamic core
(parameters:Boundary_20temperature_20model/Dynamic_20core/Alpha)=
### __Parameter name:__ Alpha
**Default value:** 1.35e-5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Core thermal expansivity. Units: \si{\per\kelvin}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Beta_20composition)=
### __Parameter name:__ Beta composition
**Default value:** 1.1

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Compositional expansion coefficient $Beta_c$. See {cite}`NPB+04` for more details.

(parameters:Boundary_20temperature_20model/Dynamic_20core/CMB_20pressure)=
### __Parameter name:__ CMB pressure
**Default value:** 0.14e12

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Pressure at CMB. Units: \si{\pascal}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Core_20conductivity)=
### __Parameter name:__ Core conductivity
**Default value:** 60.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Core heat conductivity $k_c$. Units: \si{\watt\per\meter\per\kelvin}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Core_20density)=
### __Parameter name:__ Core density
**Default value:** 12.5e3

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Density of the core. Units: \si{\kilogram\per\meter\cubed}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Core_20heat_20capacity)=
### __Parameter name:__ Core heat capacity
**Default value:** 840.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Heat capacity of the core. Units: \si{\joule\per\kelvin\per\kilogram}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Delta)=
### __Parameter name:__ Delta
**Default value:** 0.5

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** Partition coefficient of the light element.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Gravity_20acceleration)=
### __Parameter name:__ Gravity acceleration
**Default value:** 9.8

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Gravitation acceleration at CMB. Units: \si{\meter\per\second\squared}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Initial_20light_20composition)=
### __Parameter name:__ Initial light composition
**Default value:** 0.01

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Initial light composition (eg. S,O) concentration in weight fraction.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Inner_20temperature)=
### __Parameter name:__ Inner temperature
**Default value:** 6000.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the inner boundary (core mantle boundary) at the beginning. Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/K0)=
### __Parameter name:__ K0
**Default value:** 4.111e11

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Core compressibility at zero pressure. See {cite}`NPB+04` for more details.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Lh)=
### __Parameter name:__ Lh
**Default value:** 750e3

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The latent heat of core freeze. Units: \si{\joule\per\kilogram}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Max_20iteration)=
### __Parameter name:__ Max iteration
**Default value:** 30000

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The max iterations for nonlinear core energy solver.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Outer_20temperature)=
### __Parameter name:__ Outer temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the outer boundary (lithosphere water/air). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Rh)=
### __Parameter name:__ Rh
**Default value:** -27.7e6

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The heat of reaction. Units: \si{\joule\per\kilogram}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Rho0)=
### __Parameter name:__ Rho0
**Default value:** 7.019e3

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Core density at zero pressure. Units: \si{\kilogram\per\meter\cubed}. See {cite}`NPB+04` for more details.

(parameters:Boundary_20temperature_20model/Dynamic_20core/dR_20over_20dt)=
### __Parameter name:__ dR over dt
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Initial inner core radius changing rate. Units: \si{\kilo\meter}/year.

(parameters:Boundary_20temperature_20model/Dynamic_20core/dT_20over_20dt)=
### __Parameter name:__ dT over dt
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Initial CMB temperature changing rate. Units: \si{\kelvin}/year.

(parameters:Boundary_20temperature_20model/Dynamic_20core/dX_20over_20dt)=
### __Parameter name:__ dX over dt
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Initial light composition changing rate. Units: 1/year.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters)=
## **Subsection:** Boundary temperature model / Dynamic core / Geotherm parameters
(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters/Composition_20dependency)=
### __Parameter name:__ Composition dependency
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If melting curve dependent on composition.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters/Theta)=
### __Parameter name:__ Theta
**Default value:** 0.11

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Melting curve ({cite}`NPB+04` eq. (40)) parameter Theta.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters/Tm0)=
### __Parameter name:__ Tm0
**Default value:** 1695.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Melting curve ({cite}`NPB+04` eq. (40)) parameter Tm0. Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters/Tm1)=
### __Parameter name:__ Tm1
**Default value:** 10.9

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Melting curve ({cite}`NPB+04` eq. (40)) parameter Tm1. Units: \si{\per\tera\pascal}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters/Tm2)=
### __Parameter name:__ Tm2
**Default value:** -8.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Melting curve ({cite}`NPB+04` eq. (40)) parameter Tm2. Units: \si{\per\tera\pascal\squared}.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Geotherm_20parameters/Use_20BW11)=
### __Parameter name:__ Use BW11
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If using the Fe-FeS system solidus from Buono \& Walker (2011) instead.

(parameters:Boundary_20temperature_20model/Dynamic_20core/Other_20energy_20source)=
## **Subsection:** Boundary temperature model / Dynamic core / Other energy source
(parameters:Boundary_20temperature_20model/Dynamic_20core/Other_20energy_20source/File_20name)=
### __Parameter name:__ File name
**Default value:**

**Pattern:** [Anything]

**Documentation:** Data file name for other energy source into the core. The &rsquo;other energy source&rsquo; is used for external core energy source.For example if someone want to test the early lunar core powered by precession (Dwyer, C. A., et al. (2011). A long-lived lunar dynamo driven by continuous mechanical stirring. Nature 479(7372): 212-214.)Format [Time(Gyr)   Energy rate(W)]

(parameters:Boundary_20temperature_20model/Dynamic_20core/Radioactive_20heat_20source)=
## **Subsection:** Boundary temperature model / Dynamic core / Radioactive heat source
(parameters:Boundary_20temperature_20model/Dynamic_20core/Radioactive_20heat_20source/Half_20life_20times)=
### __Parameter name:__ Half life times
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Half decay times of different elements (Ga)

(parameters:Boundary_20temperature_20model/Dynamic_20core/Radioactive_20heat_20source/Heating_20rates)=
### __Parameter name:__ Heating rates
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Heating rates of different elements (W/kg)

(parameters:Boundary_20temperature_20model/Dynamic_20core/Radioactive_20heat_20source/Initial_20concentrations)=
### __Parameter name:__ Initial concentrations
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Initial concentrations of different elements (ppm)

(parameters:Boundary_20temperature_20model/Dynamic_20core/Radioactive_20heat_20source/Number_20of_20radioactive_20heating_20elements)=
### __Parameter name:__ Number of radioactive heating elements
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Number of different radioactive heating elements in core

(parameters:Boundary_20temperature_20model/Function)=
## **Subsection:** Boundary temperature model / Function
(parameters:Boundary_20temperature_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Boundary_20temperature_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Boundary_20temperature_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Boundary_20temperature_20model/Function/Maximal_20temperature)=
### __Parameter name:__ Maximal temperature
**Default value:** 3773.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximal temperature. Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Function/Minimal_20temperature)=
### __Parameter name:__ Minimal temperature
**Default value:** 273.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimal temperature. Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Boundary_20temperature_20model/Initial_20temperature)=
## **Subsection:** Boundary temperature model / Initial temperature
(parameters:Boundary_20temperature_20model/Initial_20temperature/Maximal_20temperature)=
### __Parameter name:__ Maximal temperature
**Default value:** 3773.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximal temperature. Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Initial_20temperature/Minimal_20temperature)=
### __Parameter name:__ Minimal temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimal temperature. Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Spherical_20constant)=
## **Subsection:** Boundary temperature model / Spherical constant
(parameters:Boundary_20temperature_20model/Spherical_20constant/Inner_20temperature)=
### __Parameter name:__ Inner temperature
**Default value:** 6000.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the inner boundary (core mantle boundary). Units: \si{\kelvin}.

(parameters:Boundary_20temperature_20model/Spherical_20constant/Outer_20temperature)=
### __Parameter name:__ Outer temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Temperature at the outer boundary (lithosphere water/air). Units: \si{\kelvin}.

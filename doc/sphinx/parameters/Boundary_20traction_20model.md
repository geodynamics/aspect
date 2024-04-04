(parameters:Boundary_20traction_20model)=
# Boundary traction model


## **Subsection:** Boundary traction model


(parameters:Boundary_20traction_20model/Prescribed_20traction_20boundary_20indicators)=
### __Parameter name:__ Prescribed traction boundary indicators
**Default value:**

**Pattern:** [Map of <[Anything]>:<[Selection ascii data|function|initial lithostatic pressure|zero traction ]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting those boundaries on which the traction is prescribed, i.e., where unknown external forces act to prescribe a particular traction. This is often used to prescribe a traction that equals that of overlying plates.

The format of valid entries for this parameter is that of a map given as &ldquo;key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...&rdquo; where each key must be a valid boundary indicator (which is either an integer or the symbolic name the geometry model in use may have provided for this part of the boundary) and each value must be one of the currently implemented boundary traction models. &ldquo;selector&rdquo; is an optional string given as a subset of the letters &lsquo;xyz&rsquo; that allows you to apply the boundary conditions only to the components listed. As an example, &rsquo;1 y: function&rsquo; applies the type &lsquo;function&rsquo; to the y component on boundary 1. Without a selector it will affect all components of the traction.

Note that traction should be given in N/m^2. The following boundary traction models are available:

&lsquo;ascii data&rsquo;: Implementation of a model in which the boundary traction is derived from files containing pressure data in ascii format. The pressure given in the data file is applied as traction normal to the surface of a given boundary, pointing inwards. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;Pressure [Pa]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Pressure [Pa]&rsquo; in a 3d model, which means that there has to be a single column containing the pressure. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the data will still be handled as Cartesian, however the assumed grid changes. &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;function&rsquo;: Implementation of a model in which the boundary traction is given in terms of an explicit formula that is elaborated in the parameters in section &ldquo;Boundary traction model|Function&rdquo;.

The formula you describe in the mentioned section is a semicolon separated list of traction components for each of the $d$ components of the traction vector. These $d$ formulas are interpreted as having units Pa.

&lsquo;initial lithostatic pressure&rsquo;: Implementation of a model in which the boundary traction is given in terms of a normal traction component set to the lithostatic pressure calculated according to the parameters in section &ldquo;Boundary traction model|Lithostatic pressure&rdquo;.

The lithostatic pressure is calculated by integrating the pressure downward based on the initial composition and temperature along the user-specified depth profile. The user-specified profile is given in terms of a point in Cartesian coordinates for box geometries and in spherical coordinates for all other geometries (radius, longitude, latitude), and the number of integration points. The lateral coordinates of the point are used to calculate the lithostatic pressure profile with depth. This means that the depth coordinate is not used. Note that when initial topography is included, the initial topography at the user-provided representative point is used to compute the profile. If at other points the (initial) topography is higher, the behavior of this plugin at later timesteps depends on the domain geometry. The depth returned by the geometry model does (box geometries) or does not (spherical geometries) include the initial topography. This depth is used to interpolate between the points of the reference pressure profile. Depths outside the reference profile get returned the pressure value of the closest profile depth. When only the bottom boundary is prescribed initial lithostatic pressure, the pressure value of the deepest depth of the profile is returned.

Gravity is expected to point along the depth direction.

&lsquo;zero traction&rsquo;: Implementation of a model in which the boundary traction is zero. This is commonly referred to as an &ldquo;open boundary condition&rdquo;, indicating that the material experiences no forces in response to what might exist on the other side of the boundary. However, this is only true in the case where hydrostatic pressure is not relevant. If hydrostatic pressure is not negligible, for example at the sides of a regional model, the material at the other side of the boundary does exceed a force, namely the force normal to the boundary induced by the hydrostatic pressure.

(parameters:Boundary_20traction_20model/Ascii_20data_20model)=
## **Subsection:** Boundary traction model / Ascii data model
(parameters:Boundary_20traction_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-traction/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.%d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s\%d, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model), and \%d is any sprintf integer qualifier specifying the format of the current file number.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Data_20file_20time_20step)=
### __Parameter name:__ Data file time step
**Default value:** 1e6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Time step between following data files. Depending on the setting of the global &lsquo;Use years in output instead of seconds&rsquo; flag in the input file, this number is either interpreted as seconds or as years. The default is one million, i.e., either one million seconds or one million years.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Decreasing_20file_20order)=
### __Parameter name:__ Decreasing file order
**Default value:** false

**Pattern:** [Bool]

**Documentation:** In some cases the boundary files are not numbered in increasing but in decreasing order (e.g. &lsquo;Ma BP&rsquo;). If this flag is set to &lsquo;True&rsquo; the plugin will first load the file with the number &lsquo;First data file number&rsquo; and decrease the file number during the model run.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/First_20data_20file_20model_20time)=
### __Parameter name:__ First data file model time
**Default value:** 0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The &lsquo;First data file model time&rsquo; parameter has been deactivated and will be removed in a future release. Do not use this parameter and instead provide data files starting from the model start time.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/First_20data_20file_20number)=
### __Parameter name:__ First data file number
**Default value:** 0

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Number of the first velocity file to be loaded when the model time is larger than &lsquo;First velocity file model time&rsquo;.

(parameters:Boundary_20traction_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Boundary_20traction_20model/Function)=
## **Subsection:** Boundary traction model / Function
(parameters:Boundary_20traction_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Boundary_20traction_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Boundary_20traction_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0; 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Boundary_20traction_20model/Function/Use_20spherical_20unit_20vectors)=
### __Parameter name:__ Use spherical unit vectors
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Specify traction as $r$, $\phi$, and $\theta$ components instead of $x$, $y$, and $z$. Positive tractions point up, east, and north (in 3d) or out and clockwise (in 2d). This setting only makes sense for spherical geometries.

(parameters:Boundary_20traction_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Boundary_20traction_20model/Initial_20lithostatic_20pressure)=
## **Subsection:** Boundary traction model / Initial lithostatic pressure
(parameters:Boundary_20traction_20model/Initial_20lithostatic_20pressure/Number_20of_20integration_20points)=
### __Parameter name:__ Number of integration points
**Default value:** 1000

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of integration points over which we integrate the lithostatic pressure downwards.

(parameters:Boundary_20traction_20model/Initial_20lithostatic_20pressure/Representative_20point)=
### __Parameter name:__ Representative point
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The point where the pressure profile will be calculated. Cartesian coordinates $(x,y,z)$ when geometry is a box, otherwise enter radius, longitude, and in 3d latitude. Note that the coordinate related to the depth ($y$ in 2d Cartesian, $z$ in 3d Cartesian and radius in spherical coordinates) is not used. Units: \si{\meter} or degrees.

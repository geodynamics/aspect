(parameters:Initial_20composition_20model)=
# Initial composition model


## **Subsection:** Initial composition model


(parameters:Initial_20composition_20model/List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**

**Pattern:** [MultipleSelection adiabatic density|ascii data|ascii data layered|entropy table lookup|function|porosity|slab model|world builder ]

**Documentation:** A comma-separated list of initial composition models that together describe the initial composition field. These plugins are loaded in the order given, and modify the existing composition field via the operators listed in &rsquo;List of model operators&rsquo;.

The following composition models are available:

&lsquo;adiabatic density&rsquo;: Specify the initial composition as the adiabatic reference density at each position. Note that only the field of the type &rsquo;density&rsquo; will be filled. For all other fields this plugin returns 0.0.

&lsquo;ascii data&rsquo;: Implementation of a model in which the initial composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo;, etc. in a 2d model and &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo;, etc. in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model.Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;ascii data layered&rsquo;: Implementation of a model in which the initial composition is derived from files containing data in ascii format. Each file defines a surface on which compositional fields are defined. Between the surfaces, the fields can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo; etc. in a 2d model and &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo; etc. in a 3d model; i.e. the columns before the compositional field always contains the position of the surface along the vertical direction. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle and &lsquo;y&rsquo; (if 3d) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3d is &lsquo;phi&rsquo;, &lsquo;theta&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo;and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo; as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

&lsquo;entropy table lookup&rsquo;: A class that implements initial conditions for the entropy field by converting the initial temperature field through a look up table. Note that this plugin only works if there is a compositional field called &lsquo;entropy&rsquo;, and an additional look up table that can convert pressure and temperature to entropy. For all compositional fields except entropy this plugin returns 0.0, and they are therefore not changed as long as the default &lsquo;add&rsquo; operator is selected for this plugin.

&lsquo;function&rsquo;: Specify the composition in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;porosity&rsquo;: A class that implements initial conditions for the porosity field by computing the equilibrium melt fraction for the given initial condition and reference pressure profile. Note that this plugin only works if there is a compositional field called &lsquo;porosity&rsquo;, and the used material model implements the &rsquo;MeltFractionModel&rsquo; interface. For all compositional fields except porosity this plugin returns 0.0, and they are therefore not changed as long as the default &lsquo;add&rsquo; operator is selected for this plugin.

&lsquo;slab model&rsquo;: An initial composition model that implements subducted slab geometries as a compositional field determined from an input file. The file defines the depth to the top of the slab and the slab thickness. The computed compositional value is 1 within the slabs and zero elsewhere. An example model that is included is Slab2 described in Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., Flamme, H., Furtney, M., \& Smoczyk, G. M. (2018). Slab2, a comprehensive subduction zone geometry model. Science, 362(6410), 58-61. The script to convert the Slab2 model into an aspect input data file is available in the directory data/initial-composition/slab-model/. Please note that Slab2 and the example data file assume spherical geometry (latitude, longitude coordinates), however, that is not necessary for this plugin, data files in Cartesian coordinates will work with box geometries.

&lsquo;world builder&rsquo;: Specify the initial composition through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter &rsquo;World builder file&rsquo;. It is possible to use the World Builder only for selected compositional fields by specifying the parameter &rsquo;List of relevant compositions&rsquo;.

(parameters:Initial_20composition_20model/List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ]

**Documentation:** A comma-separated list of operators that will be used to append the listed composition models onto the previous models. If only one operator is given, the same operator is applied to all models.

(parameters:Initial_20composition_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified

**Pattern:** [Selection adiabatic density|ascii data|ascii data layered|entropy table lookup|function|porosity|slab model|world builder|unspecified ]

**Documentation:** Select one of the following models:

&lsquo;adiabatic density&rsquo;: Specify the initial composition as the adiabatic reference density at each position. Note that only the field of the type &rsquo;density&rsquo; will be filled. For all other fields this plugin returns 0.0.

&lsquo;ascii data&rsquo;: Implementation of a model in which the initial composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo;, etc. in a 2d model and &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo;, etc. in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model.Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;ascii data layered&rsquo;: Implementation of a model in which the initial composition is derived from files containing data in ascii format. Each file defines a surface on which compositional fields are defined. Between the surfaces, the fields can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo; etc. in a 2d model and &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;composition1&rsquo;, &lsquo;composition2&rsquo; etc. in a 3d model; i.e. the columns before the compositional field always contains the position of the surface along the vertical direction. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle and &lsquo;y&rsquo; (if 3d) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3d is &lsquo;phi&rsquo;, &lsquo;theta&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo;and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo; as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

&lsquo;entropy table lookup&rsquo;: A class that implements initial conditions for the entropy field by converting the initial temperature field through a look up table. Note that this plugin only works if there is a compositional field called &lsquo;entropy&rsquo;, and an additional look up table that can convert pressure and temperature to entropy. For all compositional fields except entropy this plugin returns 0.0, and they are therefore not changed as long as the default &lsquo;add&rsquo; operator is selected for this plugin.

&lsquo;function&rsquo;: Specify the composition in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;porosity&rsquo;: A class that implements initial conditions for the porosity field by computing the equilibrium melt fraction for the given initial condition and reference pressure profile. Note that this plugin only works if there is a compositional field called &lsquo;porosity&rsquo;, and the used material model implements the &rsquo;MeltFractionModel&rsquo; interface. For all compositional fields except porosity this plugin returns 0.0, and they are therefore not changed as long as the default &lsquo;add&rsquo; operator is selected for this plugin.

&lsquo;slab model&rsquo;: An initial composition model that implements subducted slab geometries as a compositional field determined from an input file. The file defines the depth to the top of the slab and the slab thickness. The computed compositional value is 1 within the slabs and zero elsewhere. An example model that is included is Slab2 described in Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., Flamme, H., Furtney, M., \& Smoczyk, G. M. (2018). Slab2, a comprehensive subduction zone geometry model. Science, 362(6410), 58-61. The script to convert the Slab2 model into an aspect input data file is available in the directory data/initial-composition/slab-model/. Please note that Slab2 and the example data file assume spherical geometry (latitude, longitude coordinates), however, that is not necessary for this plugin, data files in Cartesian coordinates will work with box geometries.

&lsquo;world builder&rsquo;: Specify the initial composition through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter &rsquo;World builder file&rsquo;. It is possible to use the World Builder only for selected compositional fields by specifying the parameter &rsquo;List of relevant compositions&rsquo;.

**Warning**: This parameter provides an old and deprecated way of specifying initial composition models and shouldn&rsquo;t be used. Please use &rsquo;List of model names&rsquo; instead.

(parameters:Initial_20composition_20model/Volume_20of_20fluid_20initialization_20type)=
### __Parameter name:__ Volume of fluid initialization type
**Default value:**

**Pattern:** [Map of <[Anything]>:<[Selection composition|level set ]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting the method to be used to initialize a composition field specified to be advected using the volume of fluid method.

The format of valid entries for this parameter is that of a map given as &ldquo;key1:value1, key2:value2&ldquo; where each key must be the name of a compositional field using the volume of fluid advection method, and the value is one of &ldquo;composition&ldquo; or &ldquo;level set&ldquo;. &ldquo;composition&ldquo; is the default

When &ldquo;composition is specified, the initial model is treated as a standard composition field with bounds between 0 and 1 assumed, The initial fluid fractions are then based on an iterated midpoint quadrature. Resultant volume fractions outside of the bounds will be coerced to the nearest valid value (ie 0 or 1). If &ldquo;level set&ldquo; is specified, the initial data will be assumed to be in the form of a signed distance level set function (i.e. a function which is positive when in the fluid, negative outside, and zero on the interface and the magnitude is always the distance to the interface so the gradient is one everywhere).

(parameters:Initial_20composition_20model/Ascii_20data_20model)=
## **Subsection:** Initial composition model / Ascii data model
(parameters:Initial_20composition_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** initial_composition_top_mantle_box_3d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Data_20file_20names)=
### __Parameter name:__ Data file names
**Default value:** initial_composition_top_mantle_box_3d.txt

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** The file names of the model data (comma separated).

(parameters:Initial_20composition_20model/Ascii_20data_20model/First_20point_20on_20slice)=
### __Parameter name:__ First point on slice
**Default value:** 0.0,1.0,0.0

**Pattern:** [Anything]

**Documentation:** Point that determines the plane in which the 2d slice lies in. This variable is only used if &rsquo;Slice dataset in 2d plane&rsquo; is true. The slice will go through this point, the point defined by the parameter &rsquo;Second point on slice&rsquo;, and the center of the model domain. After the rotation, this first point will lie along the (0,1,0) axis of the coordinate system. The coordinates of the point have to be given in Cartesian coordinates.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** linear

**Pattern:** [Selection piecewise constant|linear ]

**Documentation:** Method to interpolate between layer boundaries. Select from piecewise constant or linear. Piecewise constant takes the value from the nearest layer boundary above the data point. The linear option interpolates linearly between layer boundaries. Above and below the domain given by the layer boundaries, the values aregiven by the top and bottom layer boundary.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Second_20point_20on_20slice)=
### __Parameter name:__ Second point on slice
**Default value:** 1.0,0.0,0.0

**Pattern:** [Anything]

**Documentation:** Second point that determines the plane in which the 2d slice lies in. This variable is only used if &rsquo;Slice dataset in 2d plane&rsquo; is true. The slice will go through this point, the point defined by the parameter &rsquo;First point on slice&rsquo;, and the center of the model domain. The coordinates of the point have to be given in Cartesian coordinates.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Slice_20dataset_20in_202D_20plane)=
### __Parameter name:__ Slice dataset in 2D plane
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a 2d data slice of a 3d data file or the entire data file. Slicing a 3d dataset is only supported for 2d models.

(parameters:Initial_20composition_20model/Entropy_20table_20lookup)=
## **Subsection:** Initial composition model / Entropy table lookup
(parameters:Initial_20composition_20model/Entropy_20table_20lookup/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/material-model/entropy-table/pyrtable/

**Pattern:** [DirectoryName]

**Documentation:** The path to the model data. The path may also include the special text &rsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20composition_20model/Entropy_20table_20lookup/Material_20file_20name)=
### __Parameter name:__ Material file name
**Default value:** material_table_temperature_pressure.txt

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** The file name of the material data.

(parameters:Initial_20composition_20model/Function)=
## **Subsection:** Initial composition model / Function
(parameters:Initial_20composition_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Initial_20composition_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Initial_20composition_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Initial_20composition_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Initial_20composition_20model/Slab_20model)=
## **Subsection:** Initial composition model / Slab model
(parameters:Initial_20composition_20model/Slab_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-composition/slab-model/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20composition_20model/Slab_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** shell_3d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data. Provide file in format: (File name).\%s, where \%s is a string specifying the boundary of the model according to the names of the boundary indicators (of the chosen geometry model).

(parameters:Initial_20composition_20model/Slab_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20composition_20model/World_20builder)=
## **Subsection:** Initial composition model / World builder
(parameters:Initial_20composition_20model/World_20builder/List_20of_20relevant_20compositions)=
### __Parameter name:__ List of relevant compositions
**Default value:**

**Pattern:** [Anything]

**Documentation:** A list of names of compositional fields for which to determine the initial composition using the World Builder. As World Builder evaluations can be expensive, this parameter allows to only evaluate the fields that are relevant. This plugin returns 0.0 for all compositions that are not selected in the list. By default the list is empty and the world builder is evaluated for all compositional fields.

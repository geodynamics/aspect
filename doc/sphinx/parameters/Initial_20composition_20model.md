(parameters:Initial_20composition_20model)=
# Initial composition model


## **Subsection:** Initial composition model


(parameters:Initial_20composition_20model/List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**

**Pattern:** [MultipleSelection adiabatic density|ascii data|ascii data layered|function|porosity|world builder ]

**Documentation:** A comma-separated list of initial composition models that together describe the initial composition field. These plugins are loaded in the order given, and modify the existing composition field via the operators listed in 'List of model operators'.

The following composition models are available:

`adiabatic density': Specify the initial composition as the adiabatic reference density at each position. Note that only the field of the type 'density' will be filled. For all other fields this plugin returns 0.0.

`ascii data': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2', etc. in a 2d model and `x', `y', `z', `composition1', `composition2', etc. in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model.Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`ascii data layered': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Each file defines a surface on which compositional fields are defined. Between the surfaces, the fields can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2' etc. in a 2d model and `x', `y', `z', `composition1', `composition2' etc. in a 3d model; i.e. the columns before the compositional field always contains the position of the surface along the vertical direction. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle and `y' (if 3D) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3D is `phi', `theta', `r', `T'and not `theta', `phi', `r', `T' as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

`function': Specify the composition in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`porosity': A class that implements initial conditions for the porosity field by computing the equilibrium melt fraction for the given initial condition and reference pressure profile. Note that this plugin only works if there is a compositional field called `porosity', and the used material model implements the 'MeltFractionModel' interface. For all compositional fields except porosity this plugin returns 0.0, and they are therefore not changed as long as the default `add' operator is selected for this plugin.

`world builder': Specify the initial composition through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter 'World builder file'. It is possible to use the World Builder only for selected compositional fields by specifying the parameter 'List of relevant compositions'.

(parameters:Initial_20composition_20model/List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ]

**Documentation:** A comma-separated list of operators that will be used to append the listed composition models onto the previous models. If only one operator is given, the same operator is applied to all models.

(parameters:Initial_20composition_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified

**Pattern:** [Selection adiabatic density|ascii data|ascii data layered|function|porosity|world builder|unspecified ]

**Documentation:** Select one of the following models:

`adiabatic density': Specify the initial composition as the adiabatic reference density at each position. Note that only the field of the type 'density' will be filled. For all other fields this plugin returns 0.0.

`ascii data': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2', etc. in a 2d model and `x', `y', `z', `composition1', `composition2', etc. in a 3d model, according to the number of compositional fields, which means that there has to be a single column for every composition in the model.Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the radial distance of the point to the bottom of the model, `y' by the azimuth angle and `z' by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is `r', `phi', `theta' and not `r', `theta', `phi', since this allows for dimension independent expressions.

`ascii data layered': Implementation of a model in which the initial composition is derived from files containing data in ascii format. Each file defines a surface on which compositional fields are defined. Between the surfaces, the fields can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with `#', but one of these lines needs to contain the number of grid points in each dimension as for example `# POINTS: 3 3'. The order of the data columns has to be `x', `y', `composition1', `composition2' etc. in a 2d model and `x', `y', `z', `composition1', `composition2' etc. in a 3d model; i.e. the columns before the compositional field always contains the position of the surface along the vertical direction. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. `x' will be replaced by the azimuth angle and `y' (if 3D) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3D is `phi', `theta', `r', `T'and not `theta', `phi', `r', `T' as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

`function': Specify the composition in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

`porosity': A class that implements initial conditions for the porosity field by computing the equilibrium melt fraction for the given initial condition and reference pressure profile. Note that this plugin only works if there is a compositional field called `porosity', and the used material model implements the 'MeltFractionModel' interface. For all compositional fields except porosity this plugin returns 0.0, and they are therefore not changed as long as the default `add' operator is selected for this plugin.

`world builder': Specify the initial composition through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter 'World builder file'. It is possible to use the World Builder only for selected compositional fields by specifying the parameter 'List of relevant compositions'.

\textbf{Warning}: This parameter provides an old and deprecated way of specifying initial composition models and shouldn't be used. Please use 'List of model names' instead.

(parameters:Initial_20composition_20model/Volume_20of_20fluid_20initialization_20type)=
### __Parameter name:__ Volume of fluid initialization type
**Default value:**

**Pattern:** [Map of <[Anything]>:<[Selection composition|level set ]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list denoting the method to be used to initialize a composition field specified to be advected using the volume of fluid method.

The format of valid entries for this parameter is that of a map given as ``key1:value1, key2:value2`` where each key must be the name of a compositional field using the volume of fluid advection method, and the value is one of ``composition`` or ``level set``. ``composition`` is the default

When ``composition is specified, the initial model is treated as a standard composition field with bounds between 0 and 1 assumed, The initial fluid fractions are then based on an iterated midpoint quadrature. Resultant volume fractions outside of the bounds will be coerced to the nearest valid value (ie 0 or 1). If ``level set`` is specified, the intial data will be assumed to be in the form of a signed distance level set function (i.e. a function which is positive when in the fluid, negative outside, and zero on the interface and the magnitude is always the distance to the interface so the gradient is one everywhere).

(parameters:Initial_20composition_20model/Ascii_20data_20model)=
## **Subsection:** Initial composition model / Ascii data model
(parameters:Initial_20composition_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a `/') or relative to the current directory. The path may also include the special text `$ASPECT_SOURCE_DIR' which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the `data/' subdirectory of ASPECT.

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

**Documentation:** Point that determines the plane in which the 2D slice lies in. This variable is only used if 'Slice dataset in 2D plane' is true. The slice will go through this point, the point defined by the parameter 'Second point on slice', and the center of the model domain. After the rotation, this first point will lie along the (0,1,0) axis of the coordinate system. The coordinates of the point have to be given in Cartesian coordinates.

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

**Documentation:** Second point that determines the plane in which the 2D slice lies in. This variable is only used if 'Slice dataset in 2D plane' is true. The slice will go through this point, the point defined by the parameter 'First point on slice', and the center of the model domain. The coordinates of the point have to be given in Cartesian coordinates.

(parameters:Initial_20composition_20model/Ascii_20data_20model/Slice_20dataset_20in_202D_20plane)=
### __Parameter name:__ Slice dataset in 2D plane
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a 2D data slice of a 3D data file or the entire data file. Slicing a 3D dataset is only supported for 2D models.

(parameters:Initial_20composition_20model/Function)=
## **Subsection:** Initial composition model / Function
(parameters:Initial_20composition_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are `cartesian', `spherical', and `depth'. `spherical' coordinates are interpreted as r,phi or r,phi,theta in 2D/3D respectively with theta being the polar angle. `depth' will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Initial_20composition_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form `var1=value1, var2=value2, ...'.

A typical example would be to set this runtime parameter to `pi=3.1415926536' and then use `pi' in the expression of the actual formula. (That said, for convenience this class actually defines both `pi' and `Pi' by default, but you get the idea.)

(parameters:Initial_20composition_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as `sin' or `cos'. In addition, it may contain expressions like `if(x>0, 1, -1)' where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Initial_20composition_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in 3d) for spatial coordinates and `t' for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to `r,phi,theta,t' and then use these variable names in your function expression.

(parameters:Initial_20composition_20model/World_20builder)=
## **Subsection:** Initial composition model / World builder
(parameters:Initial_20composition_20model/World_20builder/List_20of_20relevant_20compositions)=
### __Parameter name:__ List of relevant compositions
**Default value:**

**Pattern:** [Anything]

**Documentation:** A list of names of compositional fields for which to determine the initial composition using the World Builder. As World Builder evaluations can be expensive, this parameter allows to only evaluate the fields that are relevant. This plugin returns 0.0 for all compositions that are not selected in the list. By default the list is empty and the world builder is evaluated for all compositional fields.

:::::::::::::::::::::::::{dropdown} /
:open:
:name: open_

- **type**:object
- **description**:Root object
- **additionalProperties**:false
- **required**:[version, features]

::::::::::::::::::::::::{dropdown} /version
:open:
:name: open_version

- **default value**:
- **type**:string
- **description**:The major and minor version number for which the input file was written.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /$schema
:open:
:name: open_$schema

- **default value**:
- **type**:string
- **description**:The optional filename or https address to a JSON schema file
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /cross section
:open:
:name: open_cross-section

- **type**:array
- **minItems**:2
- **maxItems**:2
- **uniqueItems**:false
- **description**:This is an array of two points along where the cross section is taken
:::::::::::::::::::::::{dropdown} /cross section/items
:open:
:name: open_cross-section_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::::::{dropdown} /cross section/items/items
:open:
:name: open_cross-section_items_items

- **type**:number
::::::::::::::::::::::

:::::::::::::::::::::::

::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /potential mantle temperature
:open:
:name: open_potential-mantle-temperature

- **default value**:1600.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /surface temperature
:open:
:name: open_surface-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the surface in Kelvin.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /force surface temperature
:open:
:name: open_force-surface-temperature

- **default value**:false
- **type**:boolean
- **description**:Force the provided surface temperature to be set at the surface
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /thermal expansion coefficient
:open:
:name: open_thermal-expansion-coefficient

- **default value**:0.000035
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /specific heat
:open:
:name: open_specific-heat

- **default value**:1250.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}.$
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /thermal diffusivity
:open:
:name: open_thermal-diffusivity

- **default value**:8.04e-7
- **type**:number
- **description**:The thermal diffusivity in $m^{2} s^{-1}$.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /maximum distance between coordinates
:open:
:name: open_maximum-distance-between-coordinates

- **default value**:0.0
- **type**:number
- **description**:This enforces a maximum distance (in degree for spherical coordinates or meter in cartesian coordinates) between coordinates in the model. If the distance is larger, extra points are added by interpolation. Requires interpolation to be not 'none'.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /interpolation
:open:
:name: open_interpolation

- **default value**:continuous monotone spline
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are none, linear, monotone spline and continuous monotone spline interpolation.
::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /coordinate system
:open:
:name: open_coordinate-system

- **description**:A coordinate system. Cartesian or spherical.
- **default value**:cartesian
- **type**:object
:::::::::::::::::::::::{dropdown} /coordinate system/oneOf
:open:
:name: open_coordinate-system_oneOf

::::::::::::::::::::::{dropdown} /coordinate system/oneOf/1
:open:
:name: open_coordinate-system_oneOf_1

- **type**:object
- **description**:A Cartesian coordinate system. Coordinates are (x,y,z) and extend infinitely in all directions.
- **additionalProperties**:false
- **required**:[model]

:::::::::::::::::::::{dropdown} /coordinate system/oneOf/1/model
:open:
:name: open_coordinate-system_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the coordinate system to use.
- **enum**:[cartesian]
:::::::::::::::::::::



::::::::::::::::::::::

::::::::::::::::::::::{dropdown} /coordinate system/oneOf/2
:open:
:name: open_coordinate-system_oneOf_2

- **type**:object
- **description**:A spherical coordinate system. The coordinates are (radius, longitude, latitude). The radius is set in this plugin, the longitude extends at least from -360 to 360 degrees, and the latitude extends from -90 to 90. It is required to choose a depth method. Please see the manual for more information.
- **additionalProperties**:false
- **required**:[model, depth method]

:::::::::::::::::::::{dropdown} /coordinate system/oneOf/2/model
:open:
:name: open_coordinate-system_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the coordinate system to use.
- **enum**:[spherical]
:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /coordinate system/oneOf/2/depth method
:open:
:name: open_coordinate-system_oneOf_2_depth-method

- **default value**:
- **type**:string
- **description**:Which depth method to use in the spherical case. The available options are 'starting point', 'begin segment' and 'begin at end segment'. See the manual section on coordinate systems for more info.
- **enum**:[starting point, begin segment, begin at end segment, continuous]
:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /coordinate system/oneOf/2/radius
:open:
:name: open_coordinate-system_oneOf_2_radius

- **default value**:6371000.0
- **type**:number
- **description**:The radius of the sphere.
:::::::::::::::::::::



::::::::::::::::::::::


::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /gravity model
:open:
:name: open_gravity-model

- **description**:A gravity model for the world.
- **default value**:uniform
- **type**:object
:::::::::::::::::::::::{dropdown} /gravity model/oneOf
:open:
:name: open_gravity-model_oneOf

::::::::::::::::::::::{dropdown} /gravity model/oneOf/1
:open:
:name: open_gravity-model_oneOf_1

- **type**:object
- **description**:Uniform gravity model. It returns the gravity vector in a Cartesian coordinate system at a given position, which has a constant magitude for the whole domain. The vector points down in cartesian coordinates and to the center of the sphere in spherical coordinates.
- **additionalProperties**:false
- **required**:[model]

:::::::::::::::::::::{dropdown} /gravity model/oneOf/1/model
:open:
:name: open_gravity-model_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the model for the gravity to use.
- **enum**:[uniform]
:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /gravity model/oneOf/1/magnitude
:open:
:name: open_gravity-model_oneOf_1_magnitude

- **default value**:9.81
- **type**:number
- **description**:The magnitude of the gravity.
:::::::::::::::::::::



::::::::::::::::::::::


::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /features
:open:
:name: open_features

- **description**:A list of features.
- **default value**:
- **type**:array
:::::::::::::::::::::::{dropdown} /features/items
:open:
:name: open_features_items

::::::::::::::::::::::{dropdown} /features/items/oneOf
:open:
:name: open_features_items_oneOf

:::::::::::::::::::::{dropdown} /features/items/oneOf/1
:open:
:name: open_features_items_oneOf_1

- **type**:object
- **description**:Continental plate object. Requires properties `model` and `coordinates`.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::::::{dropdown} /features/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The model name of the feature determining its type.
- **enum**:[continental plate]
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/name
:open:
:name: open_features_items_oneOf_1_name

- **default value**:
- **type**:string
- **description**:The name which the user has given to the feature. This is mostly used for documentation purposes, and should in most cases be unique, although this is not enforced.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/tag
:open:
:name: open_features_items_oneOf_1_tag

- **default value**:
- **type**:string
- **description**:A tag which can be given to a feature. This is meant to categorize different features. If the tag is not provided or empty, it is set to the model name.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/coordinates
:open:
:name: open_features_items_oneOf_1_coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An array of 2d Points representing an array of coordinates where the feature is located.
:::::::::::::::::::{dropdown} /features/items/oneOf/1/coordinates/items
:open:
:name: open_features_items_oneOf_1_coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::{dropdown} /features/items/oneOf/1/coordinates/items/items
:open:
:name: open_features_items_oneOf_1_coordinates_items_items

- **type**:number
::::::::::::::::::

:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/interpolation
:open:
:name: open_features_items_oneOf_1_interpolation

- **default value**:global
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are 'global' and 'continuous monotone spline' interpolation. If this value is set to global, the global value for interpolation is used. This option is deprecated and will be removed in a future release.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_1_min-depth

- **description**:The depth from which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf

::::::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::::::{dropdown} /features/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

:::::::::::::::::

::::::::::::::::::


::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_1_max-depth

- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf

::::::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::::::{dropdown} /features/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

:::::::::::::::::

::::::::::::::::::


::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models
:open:
:name: open_features_items_oneOf_1_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max depth]

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/temperature models/items/oneOf/3/temperature
:open:
:name: open_features_items_oneOf_1_temperature-models_items_oneOf_3_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/composition models
:open:
:name: open_features_items_oneOf_1_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items
:open:
:name: open_features_items_oneOf_1_composition-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[random]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min value
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-value

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:Minimum value of the range within which we want to generate a random compositional value corresponding to the compositional field.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/min value/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_min-value_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max value
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-value

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:Maximum value of the range within which we want to generate a random compositional value corresponding to the compositional field.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/max value/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_max-value_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/fractions
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/fractions/items
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_1_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/1/grains models
:open:
:name: open_features_items_oneOf_1_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items
:open:
:name: open_features_items_oneOf_1_grains-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace, multiply]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/1/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_1_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::



:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /features/items/oneOf/2
:open:
:name: open_features_items_oneOf_2

- **type**:object
- **description**:Fault object. Requires properties `model` and `coordinates`.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::::::{dropdown} /features/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The model name of the feature determining its type.
- **enum**:[fault]
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/name
:open:
:name: open_features_items_oneOf_2_name

- **default value**:
- **type**:string
- **description**:The name which the user has given to the feature. This is mostly used for documentation purposes, and should in most cases be unique, although this is not enforced.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/tag
:open:
:name: open_features_items_oneOf_2_tag

- **default value**:
- **type**:string
- **description**:A tag which can be given to a feature. This is meant to categorize different features. If the tag is not provided or empty, it is set to the model name.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/coordinates
:open:
:name: open_features_items_oneOf_2_coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An array of 2d Points representing an array of coordinates where the feature is located.
:::::::::::::::::::{dropdown} /features/items/oneOf/2/coordinates/items
:open:
:name: open_features_items_oneOf_2_coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::{dropdown} /features/items/oneOf/2/coordinates/items/items
:open:
:name: open_features_items_oneOf_2_coordinates_items_items

- **type**:number
::::::::::::::::::

:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/interpolation
:open:
:name: open_features_items_oneOf_2_interpolation

- **default value**:global
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are 'global' and 'continuous monotone spline' interpolation. If this value is set to global, the global value for interpolation is used. This option is deprecated and will be removed in a future release.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_2_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_2_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/dip point
:open:
:name: open_features_items_oneOf_2_dip-point

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/2/dip point/items
:open:
:name: open_features_items_oneOf_2_dip-point_items

- **type**:number
:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/segments
:open:
:name: open_features_items_oneOf_2_segments

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items
:open:
:name: open_features_items_oneOf_2_segments_items

- **type**:object
- **additionalProperties**:false
- **description**:
- **required**:[length, thickness, angle]

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/length
:open:
:name: open_features_items_oneOf_2_segments_items_length

- **type**:number
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/thickness
:open:
:name: open_features_items_oneOf_2_segments_items_thickness

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/thickness/items
:open:
:name: open_features_items_oneOf_2_segments_items_thickness_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/top truncation
:open:
:name: open_features_items_oneOf_2_segments_items_top-truncation

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/top truncation/items
:open:
:name: open_features_items_oneOf_2_segments_items_top-truncation_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/angle
:open:
:name: open_features_items_oneOf_2_segments_items_angle

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/angle/items
:open:
:name: open_features_items_oneOf_2_segments_items_angle_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items

::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance fault center]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature starts.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature end.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2/center temperature
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2_center-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the center of this feature in degree Kelvin.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/2/side temperature
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_2_side-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the sides of this feature in degree Kelvin. If the value is below zero, an adiabatic temperature is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_3

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/temperature models/items/oneOf/3/temperature
:open:
:name: open_features_items_oneOf_2_segments_items_temperature-models_items_oneOf_3_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items

::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/side distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_side-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance over which the composition is reduced from 1 to 0.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/center fractions
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_center-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the center of the fault.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/center fractions/items
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_center-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/side fractions
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_side-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the sides of this feature.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/side fractions/items
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_side-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/fractions
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/fractions/items
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_segments_items_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items

::::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/2/segments/items/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_2_segments_items_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::



:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models
:open:
:name: open_features_items_oneOf_2_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items
:open:
:name: open_features_items_oneOf_2_temperature-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance fault center]

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature starts.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature end.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2/center temperature
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2_center-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the center of this feature in degree Kelvin.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/2/side temperature
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_2_side-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the sides of this feature in degree Kelvin. If the value is below zero, an adiabatic temperature is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_3

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/temperature models/items/oneOf/3/temperature
:open:
:name: open_features_items_oneOf_2_temperature-models_items_oneOf_3_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/composition models
:open:
:name: open_features_items_oneOf_2_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items
:open:
:name: open_features_items_oneOf_2_composition-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/side distance fault center
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_side-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance over which the composition is reduced from 1 to 0.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/center fractions
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_center-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the center of the fault.
:::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/center fractions/items
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_center-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/side fractions
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_side-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the sides of this feature.
:::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/side fractions/items
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_side-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/fractions
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/fractions/items
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/grains models
:open:
:name: open_features_items_oneOf_2_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items
:open:
:name: open_features_items_oneOf_2_grains-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/2/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_2_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/2/sections
:open:
:name: open_features_items_oneOf_2_sections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of feature properties for a coordinate.
:::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items
:open:
:name: open_features_items_oneOf_2_sections_items

- **description**:
- **default value**:
- **type**:object

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/min depth
:open:
:name: open_features_items_oneOf_2_sections_items_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/max depth
:open:
:name: open_features_items_oneOf_2_sections_items_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/dip point
:open:
:name: open_features_items_oneOf_2_sections_items_dip-point

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:The depth to which this feature is present
:::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/dip point/items
:open:
:name: open_features_items_oneOf_2_sections_items_dip-point_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments
:open:
:name: open_features_items_oneOf_2_sections_items_segments

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The depth to which this feature is present
:::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items

- **type**:object
- **additionalProperties**:false
- **description**:
- **required**:[length, thickness, angle]

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/length
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_length

- **type**:number
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/thickness
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_thickness

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/thickness/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_thickness_items

- **type**:number
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/top truncation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_top-truncation

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/top truncation/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_top-truncation_items

- **type**:number
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/angle
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_angle

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/angle/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_angle_items

- **type**:number
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance fault center]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature starts.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature end.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2/center temperature
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2_center-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the center of this feature in degree Kelvin.If the value is below zero, the an adiabatic temperature is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/2/side temperature
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_2_side-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the sides of this feature in degree Kelvin. If the value is below zero, an adiabatic temperature is used.
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_3

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/temperature models/items/oneOf/3/temperature
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_temperature-models_items_oneOf_3_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::



:::::::::::::


:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/side distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_side-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance over which the composition is reduced from 1 to 0.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/center fractions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_center-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the center of the fault.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/center fractions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_center-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/side fractions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_side-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the sides of this feature.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/side fractions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_side-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/fractions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/fractions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::



:::::::::::::


:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::

::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::

::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::

:::::::::::

::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::

::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/segments/items/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_segments_items_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::



:::::::::::::


:::::::::::::::

::::::::::::::::



:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance fault center]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature starts.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The minimum distance to the center of the fault. This determines where the linear temperature end.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2/center temperature
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2_center-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the center of this feature in degree Kelvin.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/2/side temperature
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_2_side-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the sides of this feature in degree Kelvin. If the value is below zero, an adiabatic temperature is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_3

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/temperature models/items/oneOf/3/temperature
:open:
:name: open_features_items_oneOf_2_sections_items_temperature-models_items_oneOf_3_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/side distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_side-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance over which the composition is reduced from 1 to 0.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/center fractions
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_center-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the center of the fault.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/center fractions/items
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_center-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/side fractions
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_side-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the sides of this feature.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/side fractions/items
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_side-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/fractions
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/fractions/items
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_2_sections_items_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items

::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/min distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_min-distance-fault-center

- **default value**:0.0
- **type**:number
- **description**:The distance from the fault center in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/max distance fault center
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_max-distance-fault-center

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the fault in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/2/sections/items/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_2_sections_items_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/2/sections/items/coordinate
:open:
:name: open_features_items_oneOf_2_sections_items_coordinate

- **default value**:0
- **type**:integer
- **description**:The coordinate which should be overwritten
::::::::::::::::::



:::::::::::::::::::

::::::::::::::::::::



:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /features/items/oneOf/3
:open:
:name: open_features_items_oneOf_3

- **type**:object
- **description**:Mantle layer object. Requires properties `model` and `coordinates`.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::::::{dropdown} /features/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The model name of the feature determining its type.
- **enum**:[mantle layer]
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/name
:open:
:name: open_features_items_oneOf_3_name

- **default value**:
- **type**:string
- **description**:The name which the user has given to the feature. This is mostly used for documentation purposes, and should in most cases be unique, although this is not enforced.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/tag
:open:
:name: open_features_items_oneOf_3_tag

- **default value**:
- **type**:string
- **description**:A tag which can be given to a feature. This is meant to categorize different features. If the tag is not provided or empty, it is set to the model name.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/coordinates
:open:
:name: open_features_items_oneOf_3_coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An array of 2d Points representing an array of coordinates where the feature is located.
:::::::::::::::::::{dropdown} /features/items/oneOf/3/coordinates/items
:open:
:name: open_features_items_oneOf_3_coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::{dropdown} /features/items/oneOf/3/coordinates/items/items
:open:
:name: open_features_items_oneOf_3_coordinates_items_items

- **type**:number
::::::::::::::::::

:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/interpolation
:open:
:name: open_features_items_oneOf_3_interpolation

- **default value**:global
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are 'global' and 'continuous monotone spline' interpolation. If this value is set to global, the global value for interpolation is used. This option is deprecated and will be removed in a future release.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_3_min-depth

- **description**:The depth from which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf

::::::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::::::{dropdown} /features/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

:::::::::::::::::

::::::::::::::::::


::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_3_max-depth

- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf

::::::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::::::{dropdown} /features/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

:::::::::::::::::

::::::::::::::::::


::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models
:open:
:name: open_features_items_oneOf_3_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max depth]

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/temperature models/items/oneOf/3/temperature
:open:
:name: open_features_items_oneOf_3_temperature-models_items_oneOf_3_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/composition models
:open:
:name: open_features_items_oneOf_3_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items
:open:
:name: open_features_items_oneOf_3_composition-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/fractions
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/fractions/items
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_3_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/3/grains models
:open:
:name: open_features_items_oneOf_3_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items
:open:
:name: open_features_items_oneOf_3_grains-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/3/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_3_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::



:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /features/items/oneOf/4
:open:
:name: open_features_items_oneOf_4

- **type**:object
- **description**:Oceanic plate object. Requires properties `model` and `coordinates`.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::::::{dropdown} /features/items/oneOf/4/model
:open:
:name: open_features_items_oneOf_4_model

- **default value**:
- **type**:string
- **description**:The model name of the feature determining its type.
- **enum**:[oceanic plate]
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/name
:open:
:name: open_features_items_oneOf_4_name

- **default value**:
- **type**:string
- **description**:The name which the user has given to the feature. This is mostly used for documentation purposes, and should in most cases be unique, although this is not enforced.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/tag
:open:
:name: open_features_items_oneOf_4_tag

- **default value**:
- **type**:string
- **description**:A tag which can be given to a feature. This is meant to categorize different features. If the tag is not provided or empty, it is set to the model name.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/coordinates
:open:
:name: open_features_items_oneOf_4_coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An array of 2d Points representing an array of coordinates where the feature is located.
:::::::::::::::::::{dropdown} /features/items/oneOf/4/coordinates/items
:open:
:name: open_features_items_oneOf_4_coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::{dropdown} /features/items/oneOf/4/coordinates/items/items
:open:
:name: open_features_items_oneOf_4_coordinates_items_items

- **type**:number
::::::::::::::::::

:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/interpolation
:open:
:name: open_features_items_oneOf_4_interpolation

- **default value**:global
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are 'global' and 'continuous monotone spline' interpolation. If this value is set to global, the global value for interpolation is used. This option is deprecated and will be removed in a future release.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/min depth
:open:
:name: open_features_items_oneOf_4_min-depth

- **description**:The depth from which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf

::::::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items_items

:::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf

::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::::::{dropdown} /features/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

:::::::::::::::::

::::::::::::::::::


::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/max depth
:open:
:name: open_features_items_oneOf_4_max-depth

- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf

::::::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items_items

:::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf

::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::::::{dropdown} /features/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

:::::::::::::::::

::::::::::::::::::


::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models
:open:
:name: open_features_items_oneOf_4_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2

- **type**:object
- **description**:Half space cooling mode
- **additionalProperties**:false
- **required**:[model, ridge coordinates, spreading velocity, max depth]

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[half space model]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.Because half-space reaches background temperature asymptotically, this value should be ~2 times the nominal plate thickness of 100 km
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The actual surface temperature in degree Kelvin for this feature.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The mantle temperature for the half-space cooling modelin degree Kelvin for this feature. If the model has an adiabatic gradientthis should be the mantle potential temperature, and T = Tad + Thalf. 
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity

- **description**:The spreading velocity of the plate in meter per year. This is the velocity with which one side moves away from the ridge.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:18446744073709551615
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.05
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:18446744073709551615
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/spreading velocity/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/ridge coordinates
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_ridge-coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An list of ridges. Each ridge is a lists of at least 2 2d points which define the location of the ridge. You need to define at least one ridge.So the an example with two ridges is [[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/ridge coordinates/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_ridge-coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/ridge coordinates/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_ridge-coordinates_items_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/2/ridge coordinates/items/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_2_ridge-coordinates_items_items_items

- **type**:number
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max depth]

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/top temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/3/bottom temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_3_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4

- **type**:object
- **description**:Plate model.
- **additionalProperties**:false
- **required**:[model, max depth]

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/model
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[plate model]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/operation
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/top temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/bottom temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity

- **description**:The spreading velocity of the plate in meter per year. This is the velocity with which one side moves away from the ridge.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:18446744073709551615
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.05
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:18446744073709551615
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/spreading velocity/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_spreading-velocity_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/ridge coordinates
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_ridge-coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An list of ridges. Each ridge is a lists of at least 2 2d points which define the location of the ridge. You need to define at least one ridge.So the an example with two ridges is [[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/ridge coordinates/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_ridge-coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/ridge coordinates/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_ridge-coordinates_items_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/4/ridge coordinates/items/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_4_ridge-coordinates_items_items_items

- **type**:number
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5

- **type**:object
- **description**:Plate model, but with a fixed age.
- **additionalProperties**:false
- **required**:[model, max depth]

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/model
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[plate model constant age]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/operation
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/top temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/bottom temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/5/plate age
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_5_plate-age

- **default value**:80000.0
- **type**:number
- **description**:The age of the plate in year. This age is assigned to the whole plate. 
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/model
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/operation
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth

- **description**:The depth in meters from which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth

- **description**:The depth in meters to which the temperature of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/temperature models/items/oneOf/6/temperature
:open:
:name: open_features_items_oneOf_4_temperature-models_items_oneOf_6_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/composition models
:open:
:name: open_features_items_oneOf_4_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items
:open:
:name: open_features_items_oneOf_4_composition-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1

- **type**:object
- **description**:TianWaterContent compositional model. Sets bound water content as a compositional field. The returned water content is based on the the temperature and pressure at a point within the world. Currently, the bound water content can be determined for four different lithologies: 'sediment', mid-ocean ridge basalt ('MORB'), 'gabbro', and 'peridotite', using parameterized phase diagrams from Tian et al., 2019 (https://doi.org/10.1029/2019GC008488). The pressure is lithostatic, calculated with a constant user defined density, and is limited by a user defined cutoff pressure (in GPa) for each lithology. This is required because the parameterization breaks down at large pressures. Recommended cutoff pressures are 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[tian water content]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/density
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_density

- **default value**:3000.0
- **type**:number
- **description**:The reference density used for determining the lithostatic pressure for calculating the bound water content.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/lithology
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_lithology

- **default value**:peridotite
- **type**:string
- **description**:The lithology used to determine which polynomials to use for calculating the water content. Valid options are: 'sediment', 'MORB', 'gabbro', and 'peridotite'.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/initial water content
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_initial-water-content

- **default value**:5.0
- **type**:number
- **description**:The value of the initial water content (in wt%) for the lithology at the trench. This represents the max value applied to this lithology.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/cutoff pressure
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_cutoff-pressure

- **default value**:10.0
- **type**:number
- **description**:The upper bound for the pressure, in GPa, for the specified lithology in the Tian parameterization. This is necessary because the parameterization breaks down for high pressures. It is recommended that 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/fractions
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/fractions/items
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_4_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/4/grains models
:open:
:name: open_features_items_oneOf_4_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items
:open:
:name: open_features_items_oneOf_4_grains-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth

- **description**:The depth in meters from which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_1

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.0
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/min depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_min-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth

- **description**:The depth in meters to which the composition of this feature is present.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_1

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:2
- **description**:
::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:1.7976931348623157e308
::::::::::

::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:2
::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/max depth/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_max-depth_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/4/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_4_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::



:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /features/items/oneOf/5
:open:
:name: open_features_items_oneOf_5

- **type**:object
- **description**:Plume object. Requires properties `model` and `coordinates`.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::::::{dropdown} /features/items/oneOf/5/model
:open:
:name: open_features_items_oneOf_5_model

- **default value**:
- **type**:string
- **description**:The model name of the feature determining its type.
- **enum**:[plume]
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/name
:open:
:name: open_features_items_oneOf_5_name

- **default value**:
- **type**:string
- **description**:The name which the user has given to the feature. This is mostly used for documentation purposes, and should in most cases be unique, although this is not enforced.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/tag
:open:
:name: open_features_items_oneOf_5_tag

- **default value**:
- **type**:string
- **description**:A tag which can be given to a feature. This is meant to categorize different features. If the tag is not provided or empty, it is set to the model name.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/coordinates
:open:
:name: open_features_items_oneOf_5_coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An array of 2d Points representing an array of coordinates where the feature is located.
:::::::::::::::::::{dropdown} /features/items/oneOf/5/coordinates/items
:open:
:name: open_features_items_oneOf_5_coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::{dropdown} /features/items/oneOf/5/coordinates/items/items
:open:
:name: open_features_items_oneOf_5_coordinates_items_items

- **type**:number
::::::::::::::::::

:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/interpolation
:open:
:name: open_features_items_oneOf_5_interpolation

- **default value**:global
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are 'global' and 'continuous monotone spline' interpolation. If this value is set to global, the global value for interpolation is used. This option is deprecated and will be removed in a future release.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/min depth
:open:
:name: open_features_items_oneOf_5_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth from which this feature is present, in other words, the depth of the tip of the plume. If the first entry in the cross section depths has a greater depth, an ellipsoidal plume head will be added in between. Units: m.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/max depth
:open:
:name: open_features_items_oneOf_5_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth to which this feature is present. Units: m.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/cross section depths
:open:
:name: open_features_items_oneOf_5_cross-section-depths

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The depths of the elliptic cross section of the plume. Units: m.
:::::::::::::::::::{dropdown} /features/items/oneOf/5/cross section depths/items
:open:
:name: open_features_items_oneOf_5_cross-section-depths_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/semi-major axis
:open:
:name: open_features_items_oneOf_5_semi-major-axis

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The lengths of the semi-major axes of the elliptic cross sections of the plume. In spherical coordinates, this is in degrees, otherwise in meters.
:::::::::::::::::::{dropdown} /features/items/oneOf/5/semi-major axis/items
:open:
:name: open_features_items_oneOf_5_semi-major-axis_items

- **default value**:100000.0
- **type**:number
- **description**:
:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/eccentricity
:open:
:name: open_features_items_oneOf_5_eccentricity

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The eccentricities of the cross sections.
:::::::::::::::::::{dropdown} /features/items/oneOf/5/eccentricity/items
:open:
:name: open_features_items_oneOf_5_eccentricity_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/rotation angles
:open:
:name: open_features_items_oneOf_5_rotation-angles

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The directions that the semi-major axis of the elliptic cross-sections are pointing to, in degrees. This direction is expressed as the angle from geographic North in spherical coordinates, or as the angle from the Y axis (clockwise) in Cartesian coordinates. The angle should be between 0 and 360 degrees.
:::::::::::::::::::{dropdown} /features/items/oneOf/5/rotation angles/items
:open:
:name: open_features_items_oneOf_5_rotation-angles_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models
:open:
:name: open_features_items_oneOf_5_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items
:open:
:name: open_features_items_oneOf_5_temperature-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1

- **type**:object
- **description**:Gaussian temperature model. The temperature is interpolated between the plume center and margin (as defined by the plume feature) using a Gaussian function: T(r) = T_center(z) exp(-r^2/(2 sigma^2). The temperature at the plume centerline T_center can be changed with depth by defining an array of depths and centerline temperatures, and temperature is interpolated linearly with depth. Similarly, the sigma of the Gaussian function (relative to the width of the plume as given by the plume feature) can be changed with depth. Temperature is always interpolated in a horizonzal/radial plane, except for the plume head: If the first depth of the plume centerline and the minimum depth of the plume feature are different, an ellipsoidal plume head is created in this depth range. Within this plume head, temperature is interpolated radially, i.e., depending on the distance from the center of the ellipsoid.
- **additionalProperties**:false
- **required**:[model, centerline temperatures]

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[gaussian]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/depths
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_depths

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The list of depths where both the temperature in the center of the plume and the width of the temperature anomaly in terms of the sigma of a Gaussian function can be provided. Temperature is interpolated linearly in vertical direction between these depths. Units: m.
:::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/depths/items
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_depths_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/centerline temperatures
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_centerline-temperatures

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The temperature at the center of this feature in degree Kelvin.If the value is below zero, then an adiabatic temperature is used.
:::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/centerline temperatures/items
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_centerline-temperatures_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/gaussian sigmas
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_gaussian-sigmas

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The sigma (standard deviation) of the Gaussian function used to compute the temperature distribution within the plume. This sigma is non-dimensional, i.e. it is defined relative to the distance between the plume center and margin as defined by the plume feature. Choosing a sigma of 1 therefore means that the temperature at the plume margin is set to a fraction of 1/sqrt(e) (approx. 0.61) of the centerline temperature. To achieve a smoother transition between the plume temperature and the outside temperature a smaller values has to be chosen for the gaussian sigmas.
:::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/1/gaussian sigmas/items
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_1_gaussian-sigmas_items

- **default value**:0.3
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_2

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_2_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth in meters from which the temperature of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_2_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth in meters to which the temperature of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/temperature models/items/oneOf/2/temperature
:open:
:name: open_features_items_oneOf_5_temperature-models_items_oneOf_2_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/composition models
:open:
:name: open_features_items_oneOf_5_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items
:open:
:name: open_features_items_oneOf_5_composition-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/fractions
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/fractions/items
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_5_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/5/grains models
:open:
:name: open_features_items_oneOf_5_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items
:open:
:name: open_features_items_oneOf_5_grains-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/min depth
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth in meters from which the grains of this feature are present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/max depth
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth in meters to which the grains of this feature are present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace, multiply]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/deflections
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/deflections/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis rotation matrices
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/1/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_1_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/min depth
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth in meters from which the grains of this feature are present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/max depth
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth in meters to which the grains of this feature are present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/rotation matrices
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/rotation matrices/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace, multiply]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/5/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_5_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::



:::::::::::::::::::::

:::::::::::::::::::::{dropdown} /features/items/oneOf/6
:open:
:name: open_features_items_oneOf_6

- **type**:object
- **description**:Subducting slab object. Requires properties `model` and `coordinates`.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::::::{dropdown} /features/items/oneOf/6/model
:open:
:name: open_features_items_oneOf_6_model

- **default value**:
- **type**:string
- **description**:The model name of the feature determining its type.
- **enum**:[subducting plate]
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/name
:open:
:name: open_features_items_oneOf_6_name

- **default value**:
- **type**:string
- **description**:The name which the user has given to the feature. This is mostly used for documentation purposes, and should in most cases be unique, although this is not enforced.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/tag
:open:
:name: open_features_items_oneOf_6_tag

- **default value**:
- **type**:string
- **description**:A tag which can be given to a feature. This is meant to categorize different features. If the tag is not provided or empty, it is set to the model name.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/coordinates
:open:
:name: open_features_items_oneOf_6_coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An array of 2d Points representing an array of coordinates where the feature is located.
:::::::::::::::::::{dropdown} /features/items/oneOf/6/coordinates/items
:open:
:name: open_features_items_oneOf_6_coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
::::::::::::::::::{dropdown} /features/items/oneOf/6/coordinates/items/items
:open:
:name: open_features_items_oneOf_6_coordinates_items_items

- **type**:number
::::::::::::::::::

:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/interpolation
:open:
:name: open_features_items_oneOf_6_interpolation

- **default value**:global
- **type**:string
- **description**:What type of interpolation should be used to enforce the minimum points per distance parameter. Options are 'global' and 'continuous monotone spline' interpolation. If this value is set to global, the global value for interpolation is used. This option is deprecated and will be removed in a future release.
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/min depth
:open:
:name: open_features_items_oneOf_6_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/max depth
:open:
:name: open_features_items_oneOf_6_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/dip point
:open:
:name: open_features_items_oneOf_6_dip-point

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/6/dip point/items
:open:
:name: open_features_items_oneOf_6_dip-point_items

- **type**:number
:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/segments
:open:
:name: open_features_items_oneOf_6_segments

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The depth to which this feature is present
:::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items
:open:
:name: open_features_items_oneOf_6_segments_items

- **type**:object
- **additionalProperties**:false
- **description**:
- **required**:[length, thickness, angle]

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/length
:open:
:name: open_features_items_oneOf_6_segments_items_length

- **type**:number
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/thickness
:open:
:name: open_features_items_oneOf_6_segments_items_thickness

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/thickness/items
:open:
:name: open_features_items_oneOf_6_segments_items_thickness_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/top truncation
:open:
:name: open_features_items_oneOf_6_segments_items_top-truncation

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/top truncation/items
:open:
:name: open_features_items_oneOf_6_segments_items_top-truncation_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/angle
:open:
:name: open_features_items_oneOf_6_segments_items_angle

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/angle/items
:open:
:name: open_features_items_oneOf_6_segments_items_angle_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items

::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance slab top]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the bottom in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3

- **type**:object
- **description**:Mass conserving temperature model. The temperature model uses the heat content (proportional to to thermal mass anomaly) to define a smooth temperature profile that conserves mass along the slab length. An empirical model, using error functions for smooth transitions, is used to  define how the minimum temperature increases with depth and how the location of the minimum temperature shifts into the slab interior. The slab is divided into top and bottom parts, which meet at the location where the minimum temperature occurs in the slab. For the bottom slab, the temperature is defined by a half-space cooling model. For the top of the slab the temperature is defined by one side of a 1D infinite space cooling model: this function was chosen to have a smoother temperature function across the minimum temperature position. The age of the overriding plate is used so the slab temperature at shallow depth smoothly transitions to the temperature of the overriding plate: this is not perfect, and is affected by the value of "top truncation" parameter subducting plate. Notes:1) the parameter "thickness" for the subducting plate segments needs to be defined but is not used. 2) because we use a negative truncation for distance above the slab, it is recommended to usedepth method:begin at end segment, in the main part of the world-builder file.Other methods may lead to gpas in temperatures at the segment boundaries.3)the empirical model used to define how Tmin increases with depth and how the position of Tmin shift with depth is expected to change somewhat after better calibrating with further tests.
- **additionalProperties**:false
- **required**:[model, spreading velocity, subducting velocity]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[mass conserving]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be negative and should be 1.5-2 times larger than the nominal slab thickness to allow the diffusion of cold temperatures from in the slab into the mantle above the slab surface. Also note that the top truncation value for the slab segment needs to have a value of -1, otherwise the temperature above the slab will be cut off at a distance less than the value set here.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be positive and approximately 2.5-3.0 times larger than the nominal slab thickness to allow the diffusion of coldtemperatures from in the slab into the mantle below the slab surface.For example if the slab starts with cold temperatures over a 100 km wide region, thisparameters should be about 250 km.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/density
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity

- **description**:The velocity with which the ridge spreads and create the plate in meters per year. Default is 5 cm/yr
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf

::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:18446744073709551615
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items

:::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf

::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.05
::::::::

::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:18446744073709551615
::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::

:::::::

::::::::


::::::::::

:::::::::::

::::::::::::


::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/subducting velocity
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_subducting-velocity

- **description**:The velocity with which the slab is subducting through time. Default is 5 cm/yr
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf

::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items_items

- **default value**:0.05
- **type**:number
- **description**:
::::::::::

:::::::::::

::::::::::::


::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/coupling depth
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_coupling-depth

- **default value**:100000.0
- **type**:number
- **description**:The depth at which the slab surface first comes in contact with the hot mantle wedge in meters. Default is 100 km.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/forearc cooling factor
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_forearc-cooling-factor

- **default value**:1.0
- **type**:number
- **description**:Increase the value to create thin (~2 km) cold thermal boundary layer above the slab.Any value greater than 1 does NOT meet the instantaneous conservation of mass, but does allow one to account for the history of insulating the forearc from heating up to this point in time. Note younger subducting lithosphere provides less insulation, while thicker, older slabs provide more insulation. Values up to 10 to 30 have been tested and don't cause any other extraneous effects. The larger th value the more you are not meeting the mass conserving criteria, so you don't want to see this affecting the temperature beyond the coupling depth as it will increase the mass of the slab and affect how it sinks.  If you use higher values, you will start to see that this creates a very thick cool layer above the entire slab - if you see this extending beyond the coupling zone reduce the value. You should use a value of 1 first and then only increase as little as possible to cool just the forearc region. Please examine the output temperature carefully. 
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/thermal conductivity
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_thermal-conductivity

- **default value**:3.3
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/specific heat
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/thermal diffusivity
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_thermal-diffusivity

- **default value**:-1.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/adiabatic heating
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/taper distance
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_taper-distance

- **default value**:100000.0
- **type**:number
- **description**:Distance over which to taper the slab tip.tapers the initial heat content to zero and the minimum temperature to the background temperature.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/ridge coordinates
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_ridge-coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An list of ridges. Each ridge is a lists of at least 2 2d points which define the location of the ridge. You need to define at least one ridge.So the an example with two ridges is [[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/ridge coordinates/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_ridge-coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/ridge coordinates/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_ridge-coordinates_items_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/ridge coordinates/items/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_ridge-coordinates_items_items_items

- **type**:number
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/reference model name
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_reference-model-name

- **default value**:half space model
- **type**:string
- **description**:The type of thermal model to use in the mass conserving model of slab temperature. Options are half space model and plate model
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/apply spline
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_apply-spline

- **default value**:false
- **type**:boolean
- **description**:Whether a spline should be applied on the mass conserving model.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/3/number of points in spline
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_3_number-of-points-in-spline

- **default value**:5
- **type**:integer
- **description**:The number of points in the spline
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4

- **type**:object
- **description**:Plate model (based on McKenzie, 1970).
- **additionalProperties**:false
- **required**:[model, plate velocity]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/model
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[plate model]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/operation
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/density
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/plate velocity
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_plate-velocity

- **default value**:NaN
- **type**:number
- **description**:The velocity in meters per year with which the plate subducts in meters per year.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/thermal conductivity
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_thermal-conductivity

- **default value**:2.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/specific heat
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/adiabatic heating
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab. Setting the parameter to false leads to equation 26 from McKenzie (1970),which is the result obtained from McKenzie 1969.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/4/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_4_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/5
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_5

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/5/model
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_5_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/5/operation
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_5_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/5/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_5_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/5/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_5_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/temperature models/items/oneOf/5/temperature
:open:
:name: open_features_items_oneOf_6_segments_items_temperature-models_items_oneOf_5_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items

::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_max-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/top fractions
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_top-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the top of the slab (layer).
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/top fractions/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_top-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/bottom fractions
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_bottom-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the bottom of the slab (layer).
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/bottom fractions/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_bottom-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2

- **type**:object
- **description**:TianWaterContent compositional model. Sets bound water content as a compositional field. The returned water content is based on the the temperature and pressure at a point within the world. Currently, the bound water content can be determined for four different lithologies: 'sediment', mid-ocean ridge basalt ('MORB'), 'gabbro', and 'peridotite', using parameterized phase diagrams from Tian et al., 2019 (https://doi.org/10.1029/2019GC008488). The pressure is lithostatic, calculated with a constant user defined density, and is limited by a user defined cutoff pressure (in GPa) for each lithology. This is required because the parameterization breaks down at large pressures. Recommended cutoff pressures are 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[tian water content]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/density
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_density

- **default value**:3000.0
- **type**:number
- **description**:The reference density used for determining the lithostatic pressure for calculating the bound water content.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/lithology
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_lithology

- **default value**:peridotite
- **type**:string
- **description**:The lithology used to determine which polynomials to use for calculating the water content. Valid options are: 'sediment', 'MORB', 'gabbro', and 'peridotite'.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/initial water content
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_initial-water-content

- **default value**:5.0
- **type**:number
- **description**:The value of the initial water content (in wt%) for the lithology at the trench. This represents the max value applied to this lithology.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/cutoff pressure
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_cutoff-pressure

- **default value**:10.0
- **type**:number
- **description**:The upper bound for the pressure, in GPa, for the specified lithology in the Tian parameterization. This is necessary because the parameterization breaks down for high pressures. It is recommended that 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/fractions
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/fractions/items
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/composition models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_segments_items_composition-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items

::::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/6/segments/items/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_6_segments_items_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::



:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models
:open:
:name: open_features_items_oneOf_6_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance slab top]

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the bottom in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3

- **type**:object
- **description**:Mass conserving temperature model. The temperature model uses the heat content (proportional to to thermal mass anomaly) to define a smooth temperature profile that conserves mass along the slab length. An empirical model, using error functions for smooth transitions, is used to  define how the minimum temperature increases with depth and how the location of the minimum temperature shifts into the slab interior. The slab is divided into top and bottom parts, which meet at the location where the minimum temperature occurs in the slab. For the bottom slab, the temperature is defined by a half-space cooling model. For the top of the slab the temperature is defined by one side of a 1D infinite space cooling model: this function was chosen to have a smoother temperature function across the minimum temperature position. The age of the overriding plate is used so the slab temperature at shallow depth smoothly transitions to the temperature of the overriding plate: this is not perfect, and is affected by the value of "top truncation" parameter subducting plate. Notes:1) the parameter "thickness" for the subducting plate segments needs to be defined but is not used. 2) because we use a negative truncation for distance above the slab, it is recommended to usedepth method:begin at end segment, in the main part of the world-builder file.Other methods may lead to gpas in temperatures at the segment boundaries.3)the empirical model used to define how Tmin increases with depth and how the position of Tmin shift with depth is expected to change somewhat after better calibrating with further tests.
- **additionalProperties**:false
- **required**:[model, spreading velocity, subducting velocity]

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[mass conserving]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be negative and should be 1.5-2 times larger than the nominal slab thickness to allow the diffusion of cold temperatures from in the slab into the mantle above the slab surface. Also note that the top truncation value for the slab segment needs to have a value of -1, otherwise the temperature above the slab will be cut off at a distance less than the value set here.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be positive and approximately 2.5-3.0 times larger than the nominal slab thickness to allow the diffusion of coldtemperatures from in the slab into the mantle below the slab surface.For example if the slab starts with cold temperatures over a 100 km wide region, thisparameters should be about 250 km.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/density
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity

- **description**:The velocity with which the ridge spreads and create the plate in meters per year. Default is 5 cm/yr
:::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf

::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:18446744073709551615
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items

:::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf

::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.05
::::::::::

::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:18446744073709551615
::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::::

:::::::::

::::::::::


::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/subducting velocity
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_subducting-velocity

- **description**:The velocity with which the slab is subducting through time. Default is 5 cm/yr
:::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/subducting velocity/oneOf
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_subducting-velocity_oneOf

::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/subducting velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_subducting-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/subducting velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items_items

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::


::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/coupling depth
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_coupling-depth

- **default value**:100000.0
- **type**:number
- **description**:The depth at which the slab surface first comes in contact with the hot mantle wedge in meters. Default is 100 km.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/forearc cooling factor
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_forearc-cooling-factor

- **default value**:1.0
- **type**:number
- **description**:Increase the value to create thin (~2 km) cold thermal boundary layer above the slab.Any value greater than 1 does NOT meet the instantaneous conservation of mass, but does allow one to account for the history of insulating the forearc from heating up to this point in time. Note younger subducting lithosphere provides less insulation, while thicker, older slabs provide more insulation. Values up to 10 to 30 have been tested and don't cause any other extraneous effects. The larger th value the more you are not meeting the mass conserving criteria, so you don't want to see this affecting the temperature beyond the coupling depth as it will increase the mass of the slab and affect how it sinks.  If you use higher values, you will start to see that this creates a very thick cool layer above the entire slab - if you see this extending beyond the coupling zone reduce the value. You should use a value of 1 first and then only increase as little as possible to cool just the forearc region. Please examine the output temperature carefully. 
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/thermal conductivity
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_thermal-conductivity

- **default value**:3.3
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/specific heat
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/thermal diffusivity
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_thermal-diffusivity

- **default value**:-1.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/adiabatic heating
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/taper distance
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_taper-distance

- **default value**:100000.0
- **type**:number
- **description**:Distance over which to taper the slab tip.tapers the initial heat content to zero and the minimum temperature to the background temperature.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/ridge coordinates
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_ridge-coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An list of ridges. Each ridge is a lists of at least 2 2d points which define the location of the ridge. You need to define at least one ridge.So the an example with two ridges is [[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].
:::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/ridge coordinates/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_ridge-coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/ridge coordinates/items/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_ridge-coordinates_items_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/ridge coordinates/items/items/items
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_ridge-coordinates_items_items_items

- **type**:number
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/reference model name
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_reference-model-name

- **default value**:half space model
- **type**:string
- **description**:The type of thermal model to use in the mass conserving model of slab temperature. Options are half space model and plate model
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/apply spline
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_apply-spline

- **default value**:false
- **type**:boolean
- **description**:Whether a spline should be applied on the mass conserving model.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/3/number of points in spline
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_3_number-of-points-in-spline

- **default value**:5
- **type**:integer
- **description**:The number of points in the spline
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4

- **type**:object
- **description**:Plate model (based on McKenzie, 1970).
- **additionalProperties**:false
- **required**:[model, plate velocity]

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/model
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[plate model]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/operation
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/min distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/max distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/density
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/plate velocity
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_plate-velocity

- **default value**:NaN
- **type**:number
- **description**:The velocity in meters per year with which the plate subducts in meters per year.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/thermal conductivity
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_thermal-conductivity

- **default value**:2.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/specific heat
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/adiabatic heating
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab. Setting the parameter to false leads to equation 26 from McKenzie (1970),which is the result obtained from McKenzie 1969.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/4/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_4_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/5
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_5

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/5/model
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_5_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/5/operation
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_5_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/5/min distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_5_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/5/max distance slab top
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_5_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/temperature models/items/oneOf/5/temperature
:open:
:name: open_features_items_oneOf_6_temperature-models_items_oneOf_5_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/composition models
:open:
:name: open_features_items_oneOf_6_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items
:open:
:name: open_features_items_oneOf_6_composition-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_max-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/top fractions
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_top-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the top of the slab (layer).
:::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/top fractions/items
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_top-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/bottom fractions
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_bottom-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the bottom of the slab (layer).
:::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/bottom fractions/items
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_bottom-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2

- **type**:object
- **description**:TianWaterContent compositional model. Sets bound water content as a compositional field. The returned water content is based on the the temperature and pressure at a point within the world. Currently, the bound water content can be determined for four different lithologies: 'sediment', mid-ocean ridge basalt ('MORB'), 'gabbro', and 'peridotite', using parameterized phase diagrams from Tian et al., 2019 (https://doi.org/10.1029/2019GC008488). The pressure is lithostatic, calculated with a constant user defined density, and is limited by a user defined cutoff pressure (in GPa) for each lithology. This is required because the parameterization breaks down at large pressures. Recommended cutoff pressures are 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[tian water content]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/density
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_density

- **default value**:3000.0
- **type**:number
- **description**:The reference density used for determining the lithostatic pressure for calculating the bound water content.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/lithology
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_lithology

- **default value**:peridotite
- **type**:string
- **description**:The lithology used to determine which polynomials to use for calculating the water content. Valid options are: 'sediment', 'MORB', 'gabbro', and 'peridotite'.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/initial water content
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_initial-water-content

- **default value**:5.0
- **type**:number
- **description**:The value of the initial water content (in wt%) for the lithology at the trench. This represents the max value applied to this lithology.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/cutoff pressure
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_cutoff-pressure

- **default value**:10.0
- **type**:number
- **description**:The upper bound for the pressure, in GPa, for the specified lithology in the Tian parameterization. This is necessary because the parameterization breaks down for high pressures. It is recommended that 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/fractions
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/fractions/items
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/composition models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_composition-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/grains models
:open:
:name: open_features_items_oneOf_6_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items
:open:
:name: open_features_items_oneOf_6_grains-models_items

::::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf

:::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::



:::::::::::::::::

:::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::::

:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::::{dropdown} /features/items/oneOf/6/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_6_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::::

::::::::::::::::



:::::::::::::::::


:::::::::::::::::::

::::::::::::::::::::

::::::::::::::::::::{dropdown} /features/items/oneOf/6/sections
:open:
:name: open_features_items_oneOf_6_sections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of feature properties for a coordinate.
:::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items
:open:
:name: open_features_items_oneOf_6_sections_items

- **description**:
- **default value**:
- **type**:object

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/min depth
:open:
:name: open_features_items_oneOf_6_sections_items_min-depth

- **default value**:0.0
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/max depth
:open:
:name: open_features_items_oneOf_6_sections_items_max-depth

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The depth to which this feature is present
::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/dip point
:open:
:name: open_features_items_oneOf_6_sections_items_dip-point

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:The depth to which this feature is present
:::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/dip point/items
:open:
:name: open_features_items_oneOf_6_sections_items_dip-point_items

- **type**:number
:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments
:open:
:name: open_features_items_oneOf_6_sections_items_segments

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The depth to which this feature is present
:::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items

- **type**:object
- **additionalProperties**:false
- **description**:
- **required**:[length, thickness, angle]

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/length
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_length

- **type**:number
::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/thickness
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_thickness

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/thickness/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_thickness_items

- **type**:number
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/top truncation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_top-truncation

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/top truncation/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_top-truncation_items

- **type**:number
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/angle
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_angle

- **type**:array
- **minItems**:1
- **maxItems**:2
:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/angle/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_angle_items

- **type**:number
:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance slab top]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the bottom in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3

- **type**:object
- **description**:Mass conserving temperature model. The temperature model uses the heat content (proportional to to thermal mass anomaly) to define a smooth temperature profile that conserves mass along the slab length. An empirical model, using error functions for smooth transitions, is used to  define how the minimum temperature increases with depth and how the location of the minimum temperature shifts into the slab interior. The slab is divided into top and bottom parts, which meet at the location where the minimum temperature occurs in the slab. For the bottom slab, the temperature is defined by a half-space cooling model. For the top of the slab the temperature is defined by one side of a 1D infinite space cooling model: this function was chosen to have a smoother temperature function across the minimum temperature position. The age of the overriding plate is used so the slab temperature at shallow depth smoothly transitions to the temperature of the overriding plate: this is not perfect, and is affected by the value of "top truncation" parameter subducting plate. Notes:1) the parameter "thickness" for the subducting plate segments needs to be defined but is not used. 2) because we use a negative truncation for distance above the slab, it is recommended to usedepth method:begin at end segment, in the main part of the world-builder file.Other methods may lead to gpas in temperatures at the segment boundaries.3)the empirical model used to define how Tmin increases with depth and how the position of Tmin shift with depth is expected to change somewhat after better calibrating with further tests.
- **additionalProperties**:false
- **required**:[model, spreading velocity, subducting velocity]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[mass conserving]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be negative and should be 1.5-2 times larger than the nominal slab thickness to allow the diffusion of cold temperatures from in the slab into the mantle above the slab surface. Also note that the top truncation value for the slab segment needs to have a value of -1, otherwise the temperature above the slab will be cut off at a distance less than the value set here.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be positive and approximately 2.5-3.0 times larger than the nominal slab thickness to allow the diffusion of coldtemperatures from in the slab into the mantle below the slab surface.For example if the slab starts with cold temperatures over a 100 km wide region, thisparameters should be about 250 km.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/density
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity

- **description**:The velocity with which the ridge spreads and create the plate in meters per year. Default is 5 cm/yr
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf

::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::

::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:18446744073709551615
- **description**:
::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items

:::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf

::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.05
::::::

::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:18446744073709551615
::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::

:::::

::::::


::::::::

:::::::::

::::::::::


::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/subducting velocity
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_subducting-velocity

- **description**:The velocity with which the slab is subducting through time. Default is 5 cm/yr
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf

::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::

::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items_items

- **default value**:0.05
- **type**:number
- **description**:
::::::::

:::::::::

::::::::::


::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/coupling depth
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_coupling-depth

- **default value**:100000.0
- **type**:number
- **description**:The depth at which the slab surface first comes in contact with the hot mantle wedge in meters. Default is 100 km.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/forearc cooling factor
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_forearc-cooling-factor

- **default value**:1.0
- **type**:number
- **description**:Increase the value to create thin (~2 km) cold thermal boundary layer above the slab.Any value greater than 1 does NOT meet the instantaneous conservation of mass, but does allow one to account for the history of insulating the forearc from heating up to this point in time. Note younger subducting lithosphere provides less insulation, while thicker, older slabs provide more insulation. Values up to 10 to 30 have been tested and don't cause any other extraneous effects. The larger th value the more you are not meeting the mass conserving criteria, so you don't want to see this affecting the temperature beyond the coupling depth as it will increase the mass of the slab and affect how it sinks.  If you use higher values, you will start to see that this creates a very thick cool layer above the entire slab - if you see this extending beyond the coupling zone reduce the value. You should use a value of 1 first and then only increase as little as possible to cool just the forearc region. Please examine the output temperature carefully. 
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/thermal conductivity
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_thermal-conductivity

- **default value**:3.3
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/specific heat
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/thermal diffusivity
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_thermal-diffusivity

- **default value**:-1.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/adiabatic heating
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/taper distance
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_taper-distance

- **default value**:100000.0
- **type**:number
- **description**:Distance over which to taper the slab tip.tapers the initial heat content to zero and the minimum temperature to the background temperature.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/ridge coordinates
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_ridge-coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An list of ridges. Each ridge is a lists of at least 2 2d points which define the location of the ridge. You need to define at least one ridge.So the an example with two ridges is [[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/ridge coordinates/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_ridge-coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/ridge coordinates/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_ridge-coordinates_items_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
:::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/ridge coordinates/items/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_ridge-coordinates_items_items_items

- **type**:number
:::::::::

::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/reference model name
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_reference-model-name

- **default value**:half space model
- **type**:string
- **description**:The type of thermal model to use in the mass conserving model of slab temperature. Options are half space model and plate model
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/apply spline
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_apply-spline

- **default value**:false
- **type**:boolean
- **description**:Whether a spline should be applied on the mass conserving model.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/3/number of points in spline
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_3_number-of-points-in-spline

- **default value**:5
- **type**:integer
- **description**:The number of points in the spline
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4

- **type**:object
- **description**:Plate model (based on McKenzie, 1970).
- **additionalProperties**:false
- **required**:[model, plate velocity]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[plate model]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/density
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/plate velocity
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_plate-velocity

- **default value**:NaN
- **type**:number
- **description**:The velocity in meters per year with which the plate subducts in meters per year.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/thermal conductivity
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_thermal-conductivity

- **default value**:2.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/specific heat
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/adiabatic heating
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab. Setting the parameter to false leads to equation 26 from McKenzie (1970),which is the result obtained from McKenzie 1969.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/4/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_4_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/5
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_5

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/5/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_5_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/5/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_5_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/5/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_5_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/5/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_5_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/temperature models/items/oneOf/5/temperature
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_temperature-models_items_oneOf_5_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::



:::::::::::::


:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_max-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/top fractions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_top-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the top of the slab (layer).
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/top fractions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_top-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/bottom fractions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_bottom-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the bottom of the slab (layer).
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/bottom fractions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_bottom-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2

- **type**:object
- **description**:TianWaterContent compositional model. Sets bound water content as a compositional field. The returned water content is based on the the temperature and pressure at a point within the world. Currently, the bound water content can be determined for four different lithologies: 'sediment', mid-ocean ridge basalt ('MORB'), 'gabbro', and 'peridotite', using parameterized phase diagrams from Tian et al., 2019 (https://doi.org/10.1029/2019GC008488). The pressure is lithostatic, calculated with a constant user defined density, and is limited by a user defined cutoff pressure (in GPa) for each lithology. This is required because the parameterization breaks down at large pressures. Recommended cutoff pressures are 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[tian water content]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/density
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_density

- **default value**:3000.0
- **type**:number
- **description**:The reference density used for determining the lithostatic pressure for calculating the bound water content.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/lithology
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_lithology

- **default value**:peridotite
- **type**:string
- **description**:The lithology used to determine which polynomials to use for calculating the water content. Valid options are: 'sediment', 'MORB', 'gabbro', and 'peridotite'.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/initial water content
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_initial-water-content

- **default value**:5.0
- **type**:number
- **description**:The value of the initial water content (in wt%) for the lithology at the trench. This represents the max value applied to this lithology.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/cutoff pressure
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_cutoff-pressure

- **default value**:10.0
- **type**:number
- **description**:The upper bound for the pressure, in GPa, for the specified lithology in the Tian parameterization. This is necessary because the parameterization breaks down for high pressures. It is recommended that 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/fractions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/fractions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/composition models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_composition-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::



:::::::::::::


:::::::::::::::

::::::::::::::::

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models

- **description**:
- **default value**:
- **type**:array
:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::

::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::

::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::

:::::::::::

::::::::::::



:::::::::::::

:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::

::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::

:::::::::::

::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/segments/items/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_segments_items_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::

::::::::::::



:::::::::::::


:::::::::::::::

::::::::::::::::



:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models

- **description**:A list of temperature models.
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1

- **type**:object
- **description**:Adiabatic temperature model. Uses global values by default.
- **additionalProperties**:false
- **required**:[model]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[adiabatic]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansion coefficient in $K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/1/specific heat
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_1_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat in $J kg^{-1} K^{-1}$. If the value is lower then zero, the global value is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2

- **type**:object
- **description**:Linear temperature model. Can be set to use an adiabatic temperature at the boundaries.
- **additionalProperties**:false
- **required**:[model, max distance slab top]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[linear]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2/top temperature
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2_top-temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature at the top in degree Kelvin of this feature.If the value is below zero, the an adiabatic temperature is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/2/bottom temperature
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_2_bottom-temperature

- **default value**:-1.0
- **type**:number
- **description**:The temperature at the bottom in degree Kelvin of this feature. If the value is below zero, an adiabatic temperature is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3

- **type**:object
- **description**:Mass conserving temperature model. The temperature model uses the heat content (proportional to to thermal mass anomaly) to define a smooth temperature profile that conserves mass along the slab length. An empirical model, using error functions for smooth transitions, is used to  define how the minimum temperature increases with depth and how the location of the minimum temperature shifts into the slab interior. The slab is divided into top and bottom parts, which meet at the location where the minimum temperature occurs in the slab. For the bottom slab, the temperature is defined by a half-space cooling model. For the top of the slab the temperature is defined by one side of a 1D infinite space cooling model: this function was chosen to have a smoother temperature function across the minimum temperature position. The age of the overriding plate is used so the slab temperature at shallow depth smoothly transitions to the temperature of the overriding plate: this is not perfect, and is affected by the value of "top truncation" parameter subducting plate. Notes:1) the parameter "thickness" for the subducting plate segments needs to be defined but is not used. 2) because we use a negative truncation for distance above the slab, it is recommended to usedepth method:begin at end segment, in the main part of the world-builder file.Other methods may lead to gpas in temperatures at the segment boundaries.3)the empirical model used to define how Tmin increases with depth and how the position of Tmin shift with depth is expected to change somewhat after better calibrating with further tests.
- **additionalProperties**:false
- **required**:[model, spreading velocity, subducting velocity]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[mass conserving]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be negative and should be 1.5-2 times larger than the nominal slab thickness to allow the diffusion of cold temperatures from in the slab into the mantle above the slab surface. Also note that the top truncation value for the slab segment needs to have a value of -1, otherwise the temperature above the slab will be cut off at a distance less than the value set here.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance in meters from the top surface of the slab over which the temperature is determined by this feature. This parameter should be positive and approximately 2.5-3.0 times larger than the nominal slab thickness to allow the diffusion of coldtemperatures from in the slab into the mantle below the slab surface.For example if the slab starts with cold temperatures over a 100 km wide region, thisparameters should be about 250 km.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/density
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity

- **description**:The velocity with which the ridge spreads and create the plate in meters per year. Default is 5 cm/yr
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items

- **type**:array
- **additionalProperties**:false
- **minItems**:1
- **maxItems**:18446744073709551615
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items

:::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf

::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_1

- **type**:number
- **default value**:0.05
::::::::

::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
:::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:18446744073709551615
::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/spreading velocity/oneOf/2/items/items/anyOf/2/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_spreading-velocity_oneOf_2_items_items_anyOf_2_items_items

- **type**:number
::::::

:::::::

::::::::


::::::::::

:::::::::::

::::::::::::


::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/subducting velocity
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_subducting-velocity

- **description**:The velocity with which the slab is subducting through time. Default is 5 cm/yr
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/subducting velocity/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/subducting velocity/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_1

- **default value**:0.05
- **type**:number
- **description**:
::::::::::::

::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/subducting velocity/oneOf/2/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_subducting-velocity_oneOf_2_items_items

- **default value**:0.05
- **type**:number
- **description**:
::::::::::

:::::::::::

::::::::::::


::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/coupling depth
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_coupling-depth

- **default value**:100000.0
- **type**:number
- **description**:The depth at which the slab surface first comes in contact with the hot mantle wedge in meters. Default is 100 km.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/forearc cooling factor
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_forearc-cooling-factor

- **default value**:1.0
- **type**:number
- **description**:Increase the value to create thin (~2 km) cold thermal boundary layer above the slab.Any value greater than 1 does NOT meet the instantaneous conservation of mass, but does allow one to account for the history of insulating the forearc from heating up to this point in time. Note younger subducting lithosphere provides less insulation, while thicker, older slabs provide more insulation. Values up to 10 to 30 have been tested and don't cause any other extraneous effects. The larger th value the more you are not meeting the mass conserving criteria, so you don't want to see this affecting the temperature beyond the coupling depth as it will increase the mass of the slab and affect how it sinks.  If you use higher values, you will start to see that this creates a very thick cool layer above the entire slab - if you see this extending beyond the coupling zone reduce the value. You should use a value of 1 first and then only increase as little as possible to cool just the forearc region. Please examine the output temperature carefully. 
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/thermal conductivity
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_thermal-conductivity

- **default value**:3.3
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/specific heat
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/thermal diffusivity
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_thermal-diffusivity

- **default value**:-1.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/adiabatic heating
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/taper distance
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_taper-distance

- **default value**:100000.0
- **type**:number
- **description**:Distance over which to taper the slab tip.tapers the initial heat content to zero and the minimum temperature to the background temperature.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/ridge coordinates
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_ridge-coordinates

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:An list of ridges. Each ridge is a lists of at least 2 2d points which define the location of the ridge. You need to define at least one ridge.So the an example with two ridges is [[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/ridge coordinates/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_ridge-coordinates_items

- **type**:array
- **minItems**:2
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/ridge coordinates/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_ridge-coordinates_items_items

- **type**:array
- **minItems**:2
- **maxItems**:2
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/ridge coordinates/items/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_ridge-coordinates_items_items_items

- **type**:number
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/reference model name
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_reference-model-name

- **default value**:half space model
- **type**:string
- **description**:The type of thermal model to use in the mass conserving model of slab temperature. Options are half space model and plate model
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/apply spline
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_apply-spline

- **default value**:false
- **type**:boolean
- **description**:Whether a spline should be applied on the mass conserving model.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/3/number of points in spline
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_3_number-of-points-in-spline

- **default value**:5
- **type**:integer
- **description**:The number of points in the spline
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4

- **type**:object
- **description**:Plate model (based on McKenzie, 1970).
- **additionalProperties**:false
- **required**:[model, plate velocity]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/model
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[plate model]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/operation
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/density
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_density

- **default value**:3300.0
- **type**:number
- **description**:The reference density of the subducting plate in $kg/m^3$
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/plate velocity
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_plate-velocity

- **default value**:NaN
- **type**:number
- **description**:The velocity in meters per year with which the plate subducts in meters per year.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/thermal conductivity
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_thermal-conductivity

- **default value**:2.0
- **type**:number
- **description**:The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/thermal expansion coefficient
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_thermal-expansion-coefficient

- **default value**:-1.0
- **type**:number
- **description**:The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/specific heat
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_specific-heat

- **default value**:-1.0
- **type**:number
- **description**:The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/adiabatic heating
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_adiabatic-heating

- **default value**:true
- **type**:boolean
- **description**:Whether adiabatic heating should be used for the slab. Setting the parameter to false leads to equation 26 from McKenzie (1970),which is the result obtained from McKenzie 1969.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/4/potential mantle temperature
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_4_potential-mantle-temperature

- **default value**:-1.0
- **type**:number
- **description**:The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/5
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_5

- **type**:object
- **description**:Uniform temperature model. Set the temperature to a constant value.
- **additionalProperties**:false
- **required**:[model, temperature]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/5/model
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_5_model

- **default value**:
- **type**:string
- **description**:The name of the temperature model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/5/operation
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_5_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace), add the value to the previously define value (add) or subtract the value to the previously define value (subtract).
- **enum**:[replace, add, subtract]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/5/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_5_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/5/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_5_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/temperature models/items/oneOf/5/temperature
:open:
:name: open_features_items_oneOf_6_sections_items_temperature-models_items_oneOf_5_temperature

- **default value**:293.15
- **type**:number
- **description**:The temperature in degree Kelvin which this feature should have
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models

- **description**:A list of composition models.
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1

- **type**:object
- **description**:Compositional model object
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[smooth]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_max-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance in meters from which the composition of this layer is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/top fractions
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_top-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the top of the slab (layer).
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/top fractions/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_top-fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/bottom fractions
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_bottom-fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:The composition fraction at the bottom of the slab (layer).
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/bottom fractions/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_bottom-fractions_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/1/operation
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_1_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2

- **type**:object
- **description**:TianWaterContent compositional model. Sets bound water content as a compositional field. The returned water content is based on the the temperature and pressure at a point within the world. Currently, the bound water content can be determined for four different lithologies: 'sediment', mid-ocean ridge basalt ('MORB'), 'gabbro', and 'peridotite', using parameterized phase diagrams from Tian et al., 2019 (https://doi.org/10.1029/2019GC008488). The pressure is lithostatic, calculated with a constant user defined density, and is limited by a user defined cutoff pressure (in GPa) for each lithology. This is required because the parameterization breaks down at large pressures. Recommended cutoff pressures are 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[tian water content]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/density
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_density

- **default value**:3000.0
- **type**:number
- **description**:The reference density used for determining the lithostatic pressure for calculating the bound water content.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/lithology
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_lithology

- **default value**:peridotite
- **type**:string
- **description**:The lithology used to determine which polynomials to use for calculating the water content. Valid options are: 'sediment', 'MORB', 'gabbro', and 'peridotite'.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/initial water content
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_initial-water-content

- **default value**:5.0
- **type**:number
- **description**:The value of the initial water content (in wt%) for the lithology at the trench. This represents the max value applied to this lithology.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/cutoff pressure
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_cutoff-pressure

- **default value**:10.0
- **type**:number
- **description**:The upper bound for the pressure, in GPa, for the specified lithology in the Tian parameterization. This is necessary because the parameterization breaks down for high pressures. It is recommended that 10 GPa is used for 'peridotite', 26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/2/operation
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_2_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3

- **type**:object
- **description**:Uniform compositional model. Sets constant compositional field.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the composition model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:todo The depth in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:todo The depth in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/fractions
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_fractions

- **type**:array
- **minItems**:1
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:TA list of compositional fractions corresponding to the compositions list.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/fractions/items
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_fractions_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/composition models/items/oneOf/3/operation
:open:
:name: open_features_items_oneOf_6_sections_items_composition-models_items_oneOf_3_operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value. Replacing implies that all compositions not explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.
- **enum**:[replace, replace defined only, add, subtract]
::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models

- **description**:A list of grains models.
- **default value**:
- **type**:array
:::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items

::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/model
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/orientation operation
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/1/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_1_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2

- **type**:object
- **description**:Random uniform distribution grains model. The size of the grains can be independently set to a single value or to a random distribution.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/model
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[random uniform distribution deflected]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/orientation operation
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be randomized between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_grain-sizes_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/normalize grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_normalize-grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of whether the sizes of the grains should be normalized or not. If normalized, the total of the grains of a composition will be equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/normalize grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_normalize-grain-sizes_items

- **default value**:true
- **type**:boolean
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/deflections
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_deflections

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the deflections of all of the grains in each composition between 0 and 1.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/deflections/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_deflections_items

- **default value**:1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis rotation matrices
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the rotation matrices of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis rotation matrices/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/2/basis Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_2_basis-Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::



:::::::::::::::

:::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3

- **type**:object
- **description**:Uniform grains model. All grains start exactly the same.
- **additionalProperties**:false
- **required**:[model, compositions]

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/model
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_model

- **default value**:
- **type**:string
- **description**:The name of the grains model.
- **enum**:[uniform]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/min distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_min-distance-slab-top

- **default value**:0.0
- **type**:number
- **description**:The distance from the slab top in meters from which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/max distance slab top
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_max-distance-slab-top

- **default value**:1.7976931348623157e308
- **type**:number
- **description**:The distance from the slab top in meters to which the composition of this feature is present.
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/compositions
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_compositions

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the integer labels of the composition which are present there.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/compositions/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_compositions_items

- **default value**:0
- **type**:integer
- **description**:
:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/rotation matrices
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_rotation-matrices

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the labels of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/rotation matrices/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_rotation-matrices_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/rotation matrices/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_rotation-matrices_items_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
:::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/rotation matrices/items/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_rotation-matrices_items_items_items

- **default value**:0.0
- **type**:number
- **description**:
:::::::::::

::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/Euler angles z-x-z
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_Euler-angles-z-x-z

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list with the z-x-z Euler angles of the grains which are present there for each compositions.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/Euler angles z-x-z/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items

- **type**:array
- **minItems**:3
- **maxItems**:3
- **uniqueItems**:false
- **description**:
::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/Euler angles z-x-z/items/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_Euler-angles-z-x-z_items_items

- **default value**:0.0
- **type**:number
- **description**:
::::::::::::

:::::::::::::

::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/orientation operation
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_orientation-operation

- **default value**:replace
- **type**:string
- **description**:Whether the value should replace any value previously defined at this location (replace) or add the value to the previously define value (add, not implemented). Replacing implies that all values not explicitly defined are set to zero.
- **enum**:[replace]
::::::::::::::

::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/grain sizes
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_grain-sizes

- **type**:array
- **minItems**:0
- **maxItems**:4294967295
- **uniqueItems**:false
- **description**:A list of the size of all of the grains in each composition. If set to <0, the size will be set so that the total is equal to 1.
:::::::::::::{dropdown} /features/items/oneOf/6/sections/items/grains models/items/oneOf/3/grain sizes/items
:open:
:name: open_features_items_oneOf_6_sections_items_grains-models_items_oneOf_3_grain-sizes_items

- **default value**:-1.0
- **type**:number
- **description**:
:::::::::::::

::::::::::::::



:::::::::::::::


:::::::::::::::::

::::::::::::::::::

::::::::::::::::::{dropdown} /features/items/oneOf/6/sections/items/coordinate
:open:
:name: open_features_items_oneOf_6_sections_items_coordinate

- **default value**:0
- **type**:integer
- **description**:The coordinate which should be overwritten
::::::::::::::::::



:::::::::::::::::::

::::::::::::::::::::



:::::::::::::::::::::


:::::::::::::::::::::::

::::::::::::::::::::::::

::::::::::::::::::::::::{dropdown} /random number seed
:open:
:name: open_random-number-seed

- **default value**:-1
- **type**:integer
- **description**:This allows the input of a preferred random number seed to generate random numbers. If no input is given, this value is -1 and triggers the use of default seed = 1.
::::::::::::::::::::::::






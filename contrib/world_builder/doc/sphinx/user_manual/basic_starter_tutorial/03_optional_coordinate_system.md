(part:user_manual:chap:basic_starter_tutorial:sec:03_coordinate_system)=
Optional parameter: coordinate system
=====================================


When the minimal example on the previous page is passed to the World Builder, the World Builder fills in anything else it needs with default values. The default values can be viewed in the Parameter listings part of this manual in {ref}`part:GWB_parameter_listings:chap:world_builder_file:sec:index`.

One of the optional parameters which can most fundamentally change what the resulting world looks like is the {ref}`open_coordinate-system`. Its default value is `cartesian`. Reading this document for the first (or even 10th time) is overwhelming. If you want to learn more about `/coordinate system`, open the drop box below for a detailed explanation about how to read this section of the manual.

:::::{dropdown} Coordinate system Parameter listings explanation
:name: part:user_manual:chap:concepts:sec:03_coordinate_system:parameter_explanation

Parameter listings is a very big collapsible tree in which each heading contains its path. `/` is the root of the tree, and, for example, `/coordinate system` is an object directly on the root. You can collapse objects below `/coordinate system` by clicking the upward arrow on the right, and click it again (now pointing downward) to expand it back. Each object contains its description, type, default value and additional specifications. Here you see the default value for `/coordinate system` is set to `cartesian`. Hence, if your model is in Cartesian coordinates, you don't need to do anything!

Since coordinate system is an object, we need to write something like `"coordinate system": {}` in our input file. 

The next dropbox has the path `/coordinate system/oneOf`. This is the only entry beneath the `/coordinate system` dropbox. The oneOf subpath indicates you need to choose one object out of a few options. These are listed in the oneOf box. The paths of these dropboxes in the oneOf dropbox are numbered and unfortunately not named. This means we will have to inspect each to see what they actually represent.

 In this case there are two dropboxes. Let's look at path number 1 (`/coordinate system/oneOf/1`). We can see that this is an object which has one required key (parameter): `model`. The documentation already states that this is a Cartesian coordinate system. When we look into the model dropbox (`/coordinate system/oneOf/1/model`), we see that there is an enum with only one possible value: `cartesian`. 
 
 This is a very long winded way to say that if we want a Cartesian coordinate system, we need to provide a coordinate system object which has a key called `model` and has an enum value `cartesian`: `"coordinate system": {"model":cartesian}`.

```{code-block} json
---
lineno-start: 3
---
    "coordinate system": {"model":"cartesian"}, 
```

 Now see if you can use the documentation how to create a spherical model.

::::{dropdown} Solution
:name: part:user_manual:chap:concepts:sec:03_coordinate_system:parameter_explanation:spherical

If you managed to do this in one go, great, well done! If you just tried changing `cartesian` to `spherical`, you will get an error. Take a good look at the error message. The error message indicates that a required parameter is missing. `/coordinate system/oneOf/2` requires a key `depth method`. Please see {ref}`part:user_manual:chap:concepts:sec:const_angle_spherical` for more information on why this is needed. The last option in this section is radius. This key has a default value, so it is optional. The default is set to 6371000.0, which should be fine for most models.

So a spherical model with a user defined radius of 1.0 would look like this:

```{code-block} json
---
lineno-start: 3
---
    "coordinate system":{"model":"spherical", "depth method":"begin segment", "radius":1.0}, 
```

::::


:::::

Our previous minimal example looks like this:
```{code-block} json
---
lineno-start: 1
---
{
    "version": "1.0",
    "features":[ ]
}
```


We can be more explicit and add one line setting it to the default value. However, there is no difference between this one and the previous code block. 
```{code-block} json
---
lineno-start: 1
---
{
    "version": "1.0",
    "coordinate system":{"model":"cartesian"},
    "features":[ ]
}
```

````{note}
If you want to have a spherical model, please see {ref}`part:user_manual:chap:concepts:sec:const_angle_spherical` first. An input file for a spherical model would like something like this:

```{code-block} json
---
lineno-start: 1
---
{
    "version": "1.0",
    "coordinate system":{"model":"spherical", "depth method":"begin segment"}, 
    "features":[ ]
}
```

This should be a good default spherical coordinate system. For more information on how to derive this from the parameter listing and what the options are, please expand and read the dropbox above.

In the following sections we will continue using the Cartesian coordinate system. We will also show this in  {ref}`part:user_manual:chap:basic_starter_tutorial:sec:19_spherical_models`, where we show how easy it is to switch between Cartesian and spherical coordinate systems in our finished subduction example.
````

:::{important}
The coordinate convention used in the World Builder for spherical geometries is (Radius, Longitude, Latitude)
:::

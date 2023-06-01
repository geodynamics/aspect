(part:user_manual:chap:basic_starter_tutorial:sec:03_coordinate_system)=
Optional parameter: coordinate system
=====================================


When the minimal example on the previous page is passed to the world builder, the world builder fills in anything else it needs with default values. The default values can be viewed in the Parameter listing part of this manual in {ref}`part:GWB_parameter_listings:chap:world_builder_file:sec:index`.

One of the optional parameters which can most fundamentally change what the resulting world looks like is the {ref}`open_coordinate-system`. Reading this document for the first (or even 10th time) is going to probably make you feel a bit overwhelmed and it is not needed to continue. If you want to know more, you can open the drop box below of a detailed explanation of how to read this part of the parameter section.

:::::{dropdown} Coordinate system parameter listing explanation
:name: part:user_manual:chap:concepts:sec:03_coordinate_system:parameter_explanation
 For now I will just explain what you need to know to read this part of that document.

The parameter listing is a very big collapsible tree, which in the heading has the path. `/` is the root of the tree, and `/coordinate system` is a object directly on the root. You can collapse everything below `/coordinate system` by clicking the upward arrow on the right, and click it again (now pointing downward) to expand it back. These entries contain information about them, such as the type, default value and documentation. The default is set to Cartesian, so of your model is in Cartesian coordinates, you don't need to do anything.

In our case the coordinate system is an object, which means we need to write something like `"coordinate system": {}` in our input file. 

The next dropbox has the path `/coordinate system/oneOf`. You can collapse it to see that it is the only one dropbox in the `/coordinate system` dropbox. The oneOf subpath indicates you need to choose one object out of a few options. They are listed in the oneOf box. The paths of these dropboxes in the oneOf dropbox are numbered and not named, unfortunately. This means we will have to look into them to see what they actually represent.

 In this case there are two dropboxes, so let's start with number 1 (path `/coordinate system/oneOf/1`). We can see that this is an object which has one required key/parameter: `model`. The documentation already states that this is a Cartesian coordinate system. When we look into the model dropbox (`/coordinate system/oneOf/1/model`) we see that that is an enum with only one possible value: `cartesian`. 
 
 This is a very long winded way to say that if we want a Cartesian coordinate system, we need to provide a coordinate system object which has a key called model and as a value an enum `cartesian`: `"coordinate system": {"model":cartesian}`.

```{code-block} json
---
lineno-start: 3
---
    "coordinate system": {"model":"cartesian"}, 
```

 Now see if you can find out from the documentation how to create a spherical model.

::::{dropdown} Solution
:name: part:user_manual:chap:concepts:sec:03_coordinate_system:parameter_explanation:spherical

If you managed to do it in one go, great, well done! If you just tried adding changing `cartesian` to `spherical` you will get an error. It is a nice change to learn how to read the error messages, so take a good look at it. The problem here is that when we look at the required parameters in `/coordinate system/oneOf/2`, we see that it also requires a key `depth method`. Please see {ref}`part:user_manual:chap:concepts:sec:const_angle_spherical` for more info on why this is needed. The last option in this section is radius. This key has a default value, so it is optional. The default is set to 6371000.0, which should be fine for most models.

So a spherical model with a user defined radius of 1.0 would look like this:

```{code-block} json
---
lineno-start: 3
---
    "coordinate system":{"model":"spherical", "depth method":"begin segment", "radius":1.0}, 
```

::::


:::::

This is what we had in our minimal example.
```{code-block} json
---
lineno-start: 1
---
{
    "version":"0.5",
    "features":[ ]
}
```


And now we add one line and set it to the default value. There is no difference between this one and the previous code block. 
```{code-block} json
---
lineno-start: 1
---
{
    "version":"0.5",
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
    "version":"0.5",
    "coordinate system":{"model":"spherical", "depth method":"begin segment"}, 
    "features":[ ]
}
```

This should be a good default spherical coordinate system. For more information on how to derive this from the parameter listing and what the options are, please expand and read the dropbox above.

In the following sections we will continue using the Cartesian coordinate system. We will also show this in  {ref}`part:user_manual:chap:basic_starter_tutorial:sec:17_spherical_models`, where we show how easy it is to switch between cartesian and spherical coordinate systems in our finished subduction example.
````
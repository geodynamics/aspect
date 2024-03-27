(part:user_manual:chap:concepts:sec:boundless_world)=
The world is boundless
======================

A concept that often confuses new and potential users is that the GWB domain doesn't have any bounds and is meshless. To explain how this works, we will first talk about the Cartesian coordinate system and then talk about the spherical coordinate system.

## Boundless Cartesian coordinate system

Users of the World Builder provide the locations where they want to have some tectonic features. They do this through providing map view coordinates (x and y) through the input file. We will explain later how we get the z coordinate. Coordinates can be set anywhere in an infinite 2D plane. For example the origin (0,0) is not restricted to the "center" of the domain but all features share a common coordinate system. This is not an issue for the World Builder since when World Builder computes a property, e.g., temperature, at a location, it just checks whether a feature is present at that location. So there is no need to define or limit the size of the domain! 

The reason the GWB doesn't need a mesh is similar. When the World Builder computes the temperature at a location, it determines whether it is within a feature and then determines its location within the feature. That allows the temperature model to provide different temperatures based on where it is within the feature. 

## Boundless spherical coordinate system

When we are using a spherical coordinate system, the physical domain is limited to a sphere. The input in the World Builder is, however, not restricted to 360&deg; since you may want to cross the 180&deg; or 360&deg; barrier one or many times. To be able to do this, 0&deg; and 360&deg;need to be different coordinates even though they may represent the same physical location. An example of what you can make with this can be seen below.

![feature wrapping around sphere](https://user-images.githubusercontent.com/7631629/136886253-73ac05f5-4df0-420b-9bc3-dce7c53dd14e.png)

```{note}
Although it should be possible exceed the -380&deg; to 380&deg; range, this has not been tested. If you are using a bigger range, please check the results carefully and report any issues.
```




(part:user_manual:chap:concepts:sec:boundless_world)=
The world is boundless
======================

One of the things which often confuses new and potential users is that the GWB domain doesn't have any bounds and is meshless. To explain how this works, we will first talk about the Cartesian coordinate system and then talk about the spherical coordinate system.

## Boundless Cartesian coordinate system

Users of the world builder provide the location where they want to have some tectonic features. They do this through providing map view coordinates (x and y) through the input file. We will explain later how we get the z coordinate. This means that they can set their coordinates wherever they want in an infinite 2D plane. The reason this is not an issue for the world builder is that when the world builder computes what for example the temperature at a location is, it just checks if that feature is present at that location or not. So no need to limit the domain! 

The reason the GWB doesn't need a mesh is similar. When the world builder computes the temperature at a location, and it found that it is within the feature, then it can look where it is within that feature. That allows the temperature model to provide different temperatures based on where it is within the feature. 

## Boundless spherical coordinate system

When we are using a spherical coordinate system, the physical domain is limited to a sphere. The input in the world builder is not restricted by 360 degrees though. The reason for this is that you may want to cross the 180 or 360 barrier one or many times. To be able to do this 0 and 360 degrees need to be different coordinates, even though they may represent the same physical location. An example of what you can make with this can be seen below.

![feature wrapping around sphere](https://user-images.githubusercontent.com/7631629/136886253-73ac05f5-4df0-420b-9bc3-dce7c53dd14e.png)

```{note}
Although it should be possible to do more than the -380 to 380 range, only that range has been tested. If you are using a bigger range, please check the results carefully and report any issues.
```




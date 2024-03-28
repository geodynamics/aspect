(part:user_manual:chap:parameter_documentation:sec:features:subsec:plume)=
Plume
======================

The plume feature allows adding a mantle plume to a GWB world. Specifically, the plume is described in terms of its centerline and a number of cross sections.
The location of the centerline of the plume is prescribed through the `coordinates` and the `cross section depths` parameter. The coordinates express the location of the plume center projected onto the surface of the model. In other words, the coordinates are 2D points, defining the x and y location of the plume center (in a Cartesian geometry) or the latitude and longitude of the plume center (in a spherical geometry). The cross section depths describe the depth of each of these 2D points, and they have to be listed from top to bottom (in order of ascending depth). The plume centerline is then computed through a linear interpolation between these points.

At each centerpoint defined in this way, the shape of the cross section of the plume needs to be described. We assume that the cross section has an elliptical shape, so it can be defined in terms of the length of the semi-major axis, the direction of the semi-major axis, and the eccentricity of the ellipse. This combination of parameters allows it to prescribe a pipe-like structure with changing cross-sections.

In addition, it is also possible to add a plume head at the top of the plume: If the first (shallowest) entry in the `cross section depths` list is located at a greater depth than the `min depth` of the plume feature, an ellipsoidal plume head will be added between these two depths. The ellipsoid is defined by three axes: A vertical axis pointing upwards from the location of the first coordinates on the plume centerline towards the `min depth`, and the semi-major and semi-minor axes of the ellipse definined as the first cross section of the plume. This guarantees a smooth plume head.

```{figure} ./plume_description.png
:name: plume_description
:alt: Plume geometry.
:align: center

Geometry of the plume that can be defined in the plume feature. Black ellipses are the cross sections.
```

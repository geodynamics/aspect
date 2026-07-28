(part:dev_manual:chap:parameter_documentation:sec:features:subsec:subducting_plate)=
Subducting plate
==========

The subducting plate feature provides the user with the ability to add a subducting slab to the GWB world. To define a slab, the user must specify the `coordinates` where the slab intersects with the surface of the world, which for a subducting plate represents the location of the trench. In the Cartesian coordinate system, `coordinates` are represented by `[x, y]` points (in meters), and in the spherical coordinate system `coordinates` are represented by `[longitude, latitude]` points (in degrees). Additionally, because the trench is simply a line on the surface, the user must also specify the direction that the subducting plate dips from the trench via the `dip point`. The dip point is also represented as a surface point (`[x, y]` in Cartesian, `[longitude, latitude]` in spherical).

```{figure} ./map_view.png
:name: view_of_coordinates
:alt: view of coordinates.
:align: center

Map-view of a world where a trench is defined with three `coordinates` (black points), and contains a slab dipping to the left towards the user-specified `dip point` (white point). The slab is coloured by the depth.
```

To create a slab that does not change along-strike, the user can simply specify a number of `segments`. Segments are defined through a `length`, a `thickness`, and an `angle`. When the `angle` and the `thickness` are input as single values and not as arrays, the resulting segment is defined by a straight line segment, and dips away from the trench towards the `dip point` at the specified `angle`, with a uniform `thickness`. However, the segment can vary in `thickness` if the user specifies `thickness` as an array of two values, in which case the segment will still be represented as a straight line segment, but will linearly vary from the first thickness to the second thickness across the `length` of the segment. If the user also inputs `angle` as an array of two values, the segment will be represented as a piece of a circle, where the dip of the segment varies linearly from the first angle to the second angle across the `length` of the segment. Any given segment can vary in `thickness` and in `angle`, and when multiple segments are defined the termination point of a given segment acts as the starting point of the following segment. Combining segments in this way enables a user to define extremely complicated slab geometries down-dip of the trench `coordinates`.

```{figure} ./2D_cross_section.png
:name: view_of_segments
:alt: view of segments.
:align: center

Cross section of a subducting plate feature showing how `segments` can create complicated slab geometries down-dip of the trench. Black region represents where the slab is located. The starting point of "Segment 0" is the trench, defined by the user-specified `coordinates`. The starting point of "Segment 1" is the end point of "Segment 0", and so on.
```

Additional complexity can be introduced by specifying `sections`, which allow the user to vary the shape of the slab along-strike of the trench. If the trench is defined by _N_ surface `coordinates`, then _N_ `sections` can be specified, where each `section` is mapped to the corresponding `coordinates` of the trench. The `sections` are made up of _M_ `segments` (_M_ does not need to be equal to _N_, but each section MUST be composed of _M_ segments), thereby creating _N_ unique curves along-strike of the trench, and these `sections` are interpolated along strike using a Bezier interpolation to form a complex 3-dimensional slab that varies along-strike and down-dip.

```{figure} ./downdip_sections.png
:name: view_of_sections
:alt: view of sections.
:align: center

Side view of the subducting plate feature showing how `sections` can be used to add along-strike complexity to the slab. Each of the `sections` are made up of unique `segments`, which are then interpolated along-strike using a Bezier interpolation. The slab is coloured by the depth.
```
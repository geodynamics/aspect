# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
with the addition of author(s), date of change and optionally the relevant issue. 

Add new entries a the bottom of the current list in the subheading. Item format: 
- Description. [Name; date; relevant github issue tag(s) and or pull requests]

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0]
### Added
- The World Builder Visualizer can now compress vtu files with zlib and output in binary through the vtu11 library. ASCII output is still available.  \[Menno Fraters; 2021-06-26; [#282](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/282)\]
- The option to use a half space cooling model for the temperature of an oceanic plate has been added to the oceanic plate feature. \[Magali Billen; 2021-07-15; [#317](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/317)\] 
- Information about the depth of the reference surface in the PointDistanceFromCurvedPlanes structure. \[Menno Fraters; 2021-07-20; [#320](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/320)\]
- The option for a smooth slab temperature structure that conserves mass with distance along the slab surface. \[Magali Billen; 2021-07-28; [#323](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/323)\]
- The option for a temperature model for an oceanic plate with a constant age using the plate model. \[Haoyuan Li; 2021-10-29; [#344](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/344)\]
- The cmake targets to easily switch between debug and release mode \[Menno Fraters; 2021-10-31; [#361](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/361)\]
- A new depth method for the spherical coordinate system called begin at end segment. This adds the spherical correction to the end of the segment instead of the beginning, resulting in a smoother transition between segments. \[Menno Fraters; 2021-11-06; [#364](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/364), [#365](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/365)\]
- A new input parameter and accociated functions which define the maximum depth of a model. This allows the world builder to create a complete picture of the world. \[Menno Fraters; 2021-11-08; [#367](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/367) and [#331](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/331)\]
- A one of type which can create a JSON schema entry where one of several types can be chosen. \[Menno Fraters; 2022-03-26; [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- A value at points type which can read in an array containing a value and a list of points from the input. \[Menno Fraters; 2022-03-26; [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- A surface object which can take an array of value at points and create an interpolation through delaunay triangulation (with the delaunator library which was added) and a KD-Tree and barycentric linear interplation. \[Menno Fraters, KD-Tree with help of Oliver Kreylos; 2022-03-26; [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- The option to make min and max depth into variable surfaces for all area features (continental plate, oceanic plate and mantle layer) and their temperature, composition and grain plugins. \[Menno Fraters; 2022-03-26; [#366](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/366) and [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- The option to the WordBuilderApp to input 3D spherical coordinates \[Menno Fraters; 2022-03-26; [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- An equal operator (operator==) for the point class, which compares two points with an espilon. \[Menno Fraters; 2022-03-26; [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- New gravity plugin system with a uniform gravity plugin where the gravity magnitude can be set from the input file. This how replaces the gravity input provided through the interface. The interface itself will be updated in a future pull request, where the gravity norm parameter will be removed. \[Menno Fraters; 2022-03-27; [#370](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/370)\]
- ApprovalTest have been added to the unit tests. \[Menno Fraters; 2022-04-09; [#401](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/401)\]
- Added a cmake target to update the reference test results called `update_test_references`. \[Menno Fraters; 2022-04-12; [#404](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/404)\]
- Added a new multi property query interface called properties to the world. This allows to ask for different properties at the same time, which significantly improve performance. Internally all other interface now use this properties function to reduce complexity. \[Menno Fraters; 2022-04-18; [#409](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/409) and [#410](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/410)\]
- Added a new compositional model for fault models such that ensures a smooth transition of compositional value from the fault trace until a particular user-determined distance. This feature can be helpful for model that uses composition of fault to compute other material properties, e.g., viscosity. \[Arushi Saxena; 2022-05-19; [#356](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/356)
- Added an interface to compute the distance of a query point to a feature's plane. For example, the distance of a point to a subducting slab could be duely computed. This interface simply calls the previously implemented interfaces of the feature objects and wrap them up. Thus it only takes variables like coordinates and depth in the model and could be called from ASPECT directly.
\[Haoyuan Li; 2022-12-23; [#453](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/453)
- Added tester options to allow running the tester with gdb and/or valgrind. Also setup a github action which automatically runs gdb and valgrind when running the tester. \[Menno Fraters; 2023-01-26; [#466](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/466)\]
- Added some tips and tricks in the doc/sphinx/developer_manual/contributing_to_the_code/tips_and_tricks.md file. \[Haoyuan Li; 2023-02-09; [#472](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/472)]
- Added operation options `add`, `subtract` and `replace defined only` to the the composition plugins \[Menno Fraters; 2023-02-17; [#474](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/474)\]
- Added a new compositional model for subducting slab models such that ensures a smooth transition of compositional value from one part of a compositional layer to the other side of the layer. This is based on Arushi Saxena's fault composition plugin with the same name ([#356](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/356)) \[Menno Fraters; 2023-02-18; [#477](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/356)\]
- If World Builder is configured with MPI it now reads input files on a single process and distributes them via MPI to other processes to reduce I/O load. This can be extended in the future to other input files. \[Rene Gassmoeller; 2023-04-13; [#480](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/480)\]
- Added a new html manual based on Sphinx, including a new tutorial which was designed from scratch to get new users up to speed quickly. It also contains a much structure to add cookbooks, developer and other documentation in the future. \[Menno Fraters; 2023-05-31; [#379](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/489), [#379](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/379) the related pull request mentioned in that issue.\]

### Changed
- The World Builder Visualizer will now use zlib compression for vtu files by default. If zlib is not available binary output will be used. \[Menno Fraters; 2021-06-26; [#282](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/282)\]
- The return argument type of the distance_point_from_curved_planes function has been converted from a map to a struct, requiring a change in the plugin interfaces for 'fault_models' and 'subducting_plate_models', but significantly increasing the speed of the function and all functions that access its returned values. \[Rene Gassmoeller; 2021-06-27; [#289](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/289)\]
- The plugin systems (temperature, composition and grains) and the distance_point_from_curved_planes function now all pass a precomputed NaturalCoordinate, besides just the cartesian position. It turns out that this can make a significant performance differce. \[Issue found and solution suggested by Wolfgang Bangerth, implemented and tested by Menno Fraters; 2021-07-03; [#300](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/300) and [#219](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/219)\]
- Introduces a bounding box for searching the fault and the subducting plate model properties (temperature, composition, grains) in the lateral direction. This reduces the computation time by constraining the number of points for which the expensive distance_point_from_curved_planes function is called. \[Arushi Saxena; 2021-07-15; [#219](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/219)\]
- Changing the internal interface of the get function to add a new structure (AdditionalParameters) to hold additional parameters (total_local_segment_length and local_thickness). The local_thickness has been moved away from the PointDistanceFromCurvedPlanes structure.  \[Menno Fraters; 2021-09-08; [#330](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/282)\]
- Changed the name of the main branch from master to main \[Menno Fraters; 2021-10-28; [#350](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/350)\]
- Ridge coordinates are now an array of ridges, allowing multiple ridges within a single oceanic plate with transform faults in between. \[Menno Fraters; 2021-11-03; [#362](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/362)\]
- parse entries for all temperature, composition and grain plugins now require the coordinates of the feature to be passed as a parameter. \[Menno Fraters; 2022-03-26; [#396](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/396)\]
- NaturalCoordinate has moved out of utilities to the objects namespace and folder. \[Menno Fraters; 2022-03-26; [#399](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/399)\]
- Breaking Change: Non-continuous interpolation has been removed, a lot of corner cases for the continuous interpolation have been fixed and accuracy has been improved with a new algorithm to compute the closest point on a spline. \[Menno Fraters; 2022-04-09; [#401](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/401)\]
- Breaking change: The temperature interfaces with gravity as a parameter are now deprecated and the App and Visualization programs no longer use or accept gravity. \[Menno Fraters; 2022-04-12; [#404](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/404)\]
- Breaking change: The WorldBuilderApp application has been renamed to gwb-dat. \[Menno Fraters; 2022-06-27; [#379](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/379), [#440](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/440)\]
- Breaking change: The WorldBuilderVisualization application has been renamed to gwb-grid. \[Menno Fraters; 2022-06-27; [#379](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/379), [#440](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/440)\]
- Breaking change: The Cubic monotone spline interpolation has been replaced with a Cubic Bezier spline. \[Menno Fraters; 2023-01-27; [#462](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/462), [#452](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/452), [#449](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/339), [#429](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/429)\]
- The slab/fault is now always on one side of the trench/fault line. The reference point is now used to determine whether it on the right or left side of the line between the first and last point, and the rest of the slab will remain on that side. \[Menno Fraters; 2023-01-27; [#462](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/462)\]
- Modified the mass conserving slab temperature model to use error functions to parameterize the minimum slab temperature and to allow the thickness of the top thermal boundary layer to be increased independent from the mass conserving constraints. \[Magali Billen; 2023-02-16; [#470](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/470)\]
- The interface of the properties functions of features now takes a reference to the properties vector instead of creating a copy. All internal features have been changed, but all user created features will have to be adjusted as well.  \[Rene Gassmoeller; 2023-05-11; [#485](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/485)\]
- The World Builder can now be easily build included in another CMake project by using ADD_SUBDIRECTORY. It now uses `PROJECT_*_DIR` instead of `CMAKE_*_DIR` and if it is build by an external project, the targets will have a `wb_` prefix. \[Menno Fraters; 2023-01-27; [#492](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/492)\]

### Fixed 
- Using slabs and faults beyond the -180 to 180 range gave issues. These are now fixed and it now works and is tested for the -380 to 380 range. \[Menno Fraters; 2021-10-22; [#338](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/338), [#340](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/340) and [#342](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/342)\]
- The ridge coordinates for the mass conservative slab temperature model were not converted to radians internally. \[Menno Fraters; 2021-10-27; [#352](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/352)\]
- Fixed the taper of temperature at the slab tip for the bottom part of the slab, and fixed issue with negative temperatures above slab when there is an existing overririding plate temperature. \[Magali Billen; 2021-11-02; [#353](https:https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/353)
- The input dip point defined for subduction plate and fault models is now in degrees (as opposed radians) consistent with the system for coordinates. \[Arushi Saxena; 2022-10-07; [#448](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/448)\]

## [0.4.0] - 2021-06-22
### Added
- Added basic unity build. \[Menno Fraters; 2020-08-15; [#206](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/206)\]
- Allow to chose interpolation type per feature. \[Menno Fraters; 2021-04-23; [#224](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/224)\]
- Added a benchmarks and let them automatically run by the github actions. A summary of the results are automatically posted on the github pages.  \[Menno Fraters; 2021-05-22; [#238](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/238)\]
- Added an option to not build the unit tests.  \[Menno Fraters; 2021-05-22; [#238](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/238)\]
- Added an option for Continuous Monotone Spline interpolation. This is useful for faults and slabs which need to be continuous and smooth. \[Menno Fraters; 2021-05-24; [#130](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/130),[#237](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/237)\]

### Changed
- Change minimum chame version fom 2.8.12 to 2.8.13. \[Menno Fraters; 2020-11-16; [#215](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/215)\]
- Change minimum xcode version to 11. \[Menno Fraters; 2020-12-10; [#217](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/217)\]
- Change changelog style to markdown based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). \[Menno Fraters; 2021-05-01; [#230](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/230),[#231](https://github.com/GeodynamicWorldBuilder/WorldBuilder/pull/231)\]
- Changed the underlying function which determines when you are inside a slab or a fault. This may slightly change the result, but the new result is more realistic I think. \[Menno Fraters; 2021-05-01; [#229](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/229)\]
- The Vizualizer now uses compressed output by default. This decreases the file size and increases performance. \[Menno Fraters; 2021-05-22; [#239](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/239)\]
- The Vizualizer buffers output before it writes it to a file. This increases performance. \[Menno Fraters; 2021-05-22; [#239](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/239)\]
- In the fault temperature model linear, the entiry top temperature is changed to center temperature and the entry bottom temperature is changed to side temperature, since this represents the actual sides more accurately. \[Menno Fraters; 2021-07-09; [#260](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/260)\]
- A large overhoal of the distance_point_from_curved_planes function improving the accuracy in spherical coordinates. This may slightly change the results, but it should be an improvement both in accuracy and performance. Also makes available a coordinate system independent distrance function and some fast trigonometric functions. \[Menno Fraters; 2021-06-11; [#255](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/255)\]

### Fixed
- Fixed namespaces, adding `WorldBuilder::` where needed. \[Timo Heister; 2020-08-10; [#205](github.com/GeodynamicWorldBuilder/WorldBuilder/issues/205)\] 
- Fixed bug in linear temperature models in several features when using a top/center temperature which computes the background (-1). \[Menno Fraters; 2021-06-22; [#273](github.com/GeodynamicWorldBuilder/WorldBuilder/pull/273)\] 

## [0.3.0] - 2020-06-19
### Added
- An experimental grains interface containing grain sizes and orientations has been added. [Menno Fraters;;]
### Changed
- Fortran compiler is now optional [Menno Fraters;;]
- The world builder no longer requires or uses OpenMP. It was not used in the core  library, but used to speed up the Visualization program. That uses its own threads now. [Menno Fraters;;] 
- The usage of cmake has been moderized allowing for more flexibility [Menno Fraters;;]
### Fixed
- The code has been refactored and cleaned exstinively and new compiler warnings have been  enabled and fixed. [Menno Fraters;;]

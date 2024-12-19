(parameters:Geometry_20model)=
# Geometry model


## **Subsection:** Geometry model


(parameters:Geometry_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified

**Pattern:** [Selection box|box with lithosphere boundary indicators|chunk|chunk with lithosphere boundary indicators|ellipsoidal chunk|sphere|spherical shell|unspecified ]

**Documentation:** Select one of the following models:

&lsquo;box&rsquo;: A box geometry parallel to the coordinate directions. The extent of the box in each coordinate direction is set in the parameter file. The box geometry labels its 2*dim sides as follows: in 2d, boundary indicators 0 through 3 denote the left, right, bottom and top boundaries; in 3d, boundary indicators 0 through 5 indicate left, right, front, back, bottom and top boundaries (see also the documentation of the deal.II class &ldquo;ReferenceCell&rdquo;). You can also use symbolic names &ldquo;left&rdquo;, &ldquo;right&rdquo;, etc., to refer to these boundaries in input files. It is also possible to add initial topography to the box model. Note however that this is done after the last initial adaptive refinement cycle. Also, initial topography is supposed to be small, as it is not taken into account when depth or a representative point is computed.

&lsquo;box with lithosphere boundary indicators&rsquo;: A box geometry parallel to the coordinate directions. The extent of the box in each coordinate direction is set in the parameter file. This geometry model labels its sides with 2*dim+2*(dim-1) boundary indicators: in 2d, boundary indicators 0 through 3 denote the left, right, bottom and top boundaries, while indicators4 and 5 denote the upper part of the left and right vertical boundary, respectively. In 3d, boundary indicators 0 through 5 indicate left, right, front, back, bottom and top boundaries (see also the documentation of the deal.II class &ldquo;ReferenceCell&rdquo;), while indicators 6, 7, 8 and 9 denote the left, right, front and back upper parts of the vertical boundaries, respectively. You can also use symbolic names &ldquo;left&rdquo;, &ldquo;right&rdquo;, &ldquo;left lithosphere&rdquo;, etc., to refer to these boundaries in input files.

Note that for a given &ldquo;Global refinement level&rdquo; and no user-specified &ldquo;Repetitions&rdquo;, the lithosphere part of the mesh will be more refined.

The additional boundary indicators for the lithosphere allow for selecting boundary conditions for the lithosphere different from those for the underlying mantle. An example application of this geometry is to prescribe a velocity on the lithospheric plates, but use open boundary conditions underneath.

&lsquo;chunk&rsquo;: A geometry which can be described as a chunk of a spherical shell, bounded by lines of longitude, latitude and radius. The minimum and maximum longitude, latitude (if in 3d) and depth of the chunk is set in the parameter file. The chunk geometry labels its 2*dim sides as follows: &ldquo;west&rdquo; and &ldquo;east&rdquo;: minimum and maximum longitude, &ldquo;south&rdquo; and &ldquo;north&rdquo;: minimum and maximum latitude, &ldquo;inner&rdquo; and &ldquo;outer&rdquo;: minimum and maximum radii.

The dimensions of the model are specified by parameters of the following form: Chunk (minimum || maximum) (longitude || latitude): edges of geographical quadrangle (in degrees)Chunk (inner || outer) radius: Radii at bottom and top of chunk(Longitude || Latitude || Radius) repetitions: number of cells in each coordinate direction.

When used in 2d, this geometry does not imply the use of a spherical coordinate system. Indeed, in 2d the geometry is simply a sector of an annulus in a Cartesian coordinate system and consequently would correspond to a sector of a cross section of the fluid filled space between two infinite cylinders where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in {ref}`sec:methods:2d-models`. It is also possible to add initial topography to the chunk geometry, based on an ascii data file.

&lsquo;chunk with lithosphere boundary indicators&rsquo;: A geometry which can be described as a chunk of a spherical shell, bounded by lines of longitude, latitude and radius. The side boundaries have two boundary indicators, so the user can prescribe different boundary conditions on these boundaries. The minimum and maximum longitude, (latitude) and depth of the chunk are set in the parameter file. The chunk geometry labels its 2*dim+2*(dim-1) sides as follows: &ldquo;lowerwest&rdquo; and &ldquo;lowereast&rdquo;: minimum and maximum longitude of the lower part of the east and west side boundaries, &ldquo;upperwest&rdquo; and &ldquo;uppereast&rdquo;: minimum and maximum longitude of the upper part of the east and west side boundaries, &ldquo;lowersouth&rdquo; and &ldquo;lowernorth&rdquo;: minimum and maximum latitude of the lower part of the south and north side boundaries, &ldquo;uppersouth&rdquo; and &ldquo;uppernorth&rdquo;: minimum and maximum latitude of the upper part of the south and north side boundaries,

The dimensions of the model are specified by parameters of the following form: Chunk (minimum | maximum) (longitude | latitude): edges of geographical quadrangle (in degrees). Chunk (inner | outer | middle boundary) radius: Radii at bottom and top of chunk and the radius at which the lower boundary indicator along a side boundary transitions into the upper boundary indicator. (Longitude | Latitude) repetitions: number of cells in each coordinate direction.(Inner | Outer) chunk radius repetitions: number of cells in the radial coordinate direction for the lower part of the domain (up to the Middle boundary radius) and for the upper part of the domain.

When used in 2d, this geometry does not imply the use of a spherical coordinate system. Indeed, in 2d the geometry is simply a sector of an annulus in a Cartesian coordinate system and consequently would correspond to a sector of a cross section of the fluid filled space between two infinite cylinders where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in {ref}`sec:methods:2d-models`. It is also possible to add initial topography to the chunk geometry, based on an ascii data file.

&lsquo;ellipsoidal chunk&rsquo;: A 3d chunk geometry that accounts for Earth&rsquo;s ellipticity (default assuming the WGS84 ellipsoid definition) which can be defined in non-coordinate directions. In the description of the ellipsoidal chunk, two of the ellipsoidal axes have the same length so that there is only a semi-major axis and a semi-minor axis. The user has two options for creating an ellipsoidal chunk geometry: 1) by defining two opposing points (SW and NE or NW and SE) a coordinate parallel ellipsoidal chunk geometry will be created. 2) by defining three points a non-coordinate parallel ellipsoidal chunk will be created. The points are defined in the input file by longitude:latitude. It is also possible to define additional subdivisions of the mesh in each direction. The boundary of the domain is formed by linear interpolation in longitude-latitude space between adjacent points (i.e. ${lon, lat}(f) = [lon1 \cdot f + lon2 \cdot(1-f), lat1 \cdot f + lat2 \cdot (1-f)]$, where f is a value between 0 and 1). Faces of the model are defined as 0, west; 1,east; 2, south; 3, north; 4, inner; 5, outer.

This geometry model supports initial topography for deforming the initial mesh.

&lsquo;sphere&rsquo;: A geometry model for a sphere with a user specified radius. This geometry has only a single boundary, so the only valid boundary indicator to specify in input files is &ldquo;0&rdquo;. It can also be referenced by the symbolic name &ldquo;surface&rdquo; in input files.

Despite the name, this geometry does not imply the use of a spherical coordinate system when used in 2d. Indeed, in 2d the geometry is simply a circle in a Cartesian coordinate system and consequently would correspond to a cross section of the fluid filled interior of an infinite cylinder where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in {ref}`sec:methods:2d-models`.

&lsquo;spherical shell&rsquo;: A geometry representing a spherical shell or a piece of it. Inner and outer radii are read from the parameter file in subsection &rsquo;Spherical shell&rsquo;.

The spherical shell may be generated as per the original code (with respect to the inner and outer radius, and an initial number of cells along circumference) or following a custom mesh scheme: list of radial values or number of slices. A surface mesh is first generated and refined as desired, before it is extruded radially. A list of radial values subdivides the spherical shell at specified radii. The number of slices subdivides the spherical shell into N slices of equal thickness. The custom spherical shell only works with an opening angle of 360 degrees.

Despite the name, this geometry does not imply the use of a spherical coordinate system when used in 2d. Indeed, in 2d the geometry is simply an annulus in a Cartesian coordinate system and consequently would correspond to a cross section of the fluid filled space between two infinite cylinders where one has made the assumption that the velocity in direction of the cylinder axes is zero. This is consistent with the definition of what we consider the two-dimension case given in {ref}`sec:methods:2d-models`.

The model assigns boundary indicators as follows: In 2d, inner and outer boundaries get boundary indicators zero and one, and if the opening angle set in the input file is less than 360, then left and right boundaries are assigned indicators two and three. These boundaries can also be referenced using the symbolic names &lsquo;inner&rsquo;, &lsquo;outer&rsquo; and (if applicable) &lsquo;left&rsquo;, &lsquo;right&rsquo;.

In 3d, inner and outer indicators are treated as in 2d. If the opening angle is chosen as 90 degrees, i.e., the domain is the intersection of a spherical shell and the first octant, then indicator 2 is at the face $x=0$, 3 at $y=0$, and 4 at $z=0$. These last three boundaries can then also be referred to as &lsquo;east&rsquo;, &lsquo;west&rsquo; and &lsquo;south&rsquo; symbolically in input files.

(parameters:Geometry_20model/Box)=
## **Subsection:** Geometry model / Box
(parameters:Geometry_20model/Box/Box_20origin_20X_20coordinate)=
### __Parameter name:__ Box origin X coordinate
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** X coordinate of box origin. Units: \si{\meter}.

(parameters:Geometry_20model/Box/Box_20origin_20Y_20coordinate)=
### __Parameter name:__ Box origin Y coordinate
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Y coordinate of box origin. Units: \si{\meter}.

(parameters:Geometry_20model/Box/Box_20origin_20Z_20coordinate)=
### __Parameter name:__ Box origin Z coordinate
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Z coordinate of box origin. This value is ignored if the simulation is in 2d. Units: \si{\meter}.

(parameters:Geometry_20model/Box/X_20extent)=
### __Parameter name:__ X extent
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Extent of the box in x-direction. Units: \si{\meter}.

(parameters:Geometry_20model/Box/X_20periodic)=
### __Parameter name:__ X periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in X direction

(parameters:Geometry_20model/Box/X_20repetitions)=
### __Parameter name:__ X repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in X direction.

(parameters:Geometry_20model/Box/Y_20extent)=
### __Parameter name:__ Y extent
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Extent of the box in y-direction. Units: \si{\meter}.

(parameters:Geometry_20model/Box/Y_20periodic)=
### __Parameter name:__ Y periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in Y direction

(parameters:Geometry_20model/Box/Y_20repetitions)=
### __Parameter name:__ Y repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in Y direction.

(parameters:Geometry_20model/Box/Z_20extent)=
### __Parameter name:__ Z extent
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Extent of the box in z-direction. This value is ignored if the simulation is in 2d. Units: \si{\meter}.

(parameters:Geometry_20model/Box/Z_20periodic)=
### __Parameter name:__ Z periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in Z direction

(parameters:Geometry_20model/Box/Z_20repetitions)=
### __Parameter name:__ Z repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in Z direction.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators)=
## **Subsection:** Geometry model / Box with lithosphere boundary indicators
(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Box_20origin_20X_20coordinate)=
### __Parameter name:__ Box origin X coordinate
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** X coordinate of box origin. Units: \si{\meter}.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Box_20origin_20Y_20coordinate)=
### __Parameter name:__ Box origin Y coordinate
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Y coordinate of box origin. Units: \si{\meter}.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Box_20origin_20Z_20coordinate)=
### __Parameter name:__ Box origin Z coordinate
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Z coordinate of box origin. This value is ignored if the simulation is in 2d. Units: \si{\meter}.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Lithospheric_20thickness)=
### __Parameter name:__ Lithospheric thickness
**Default value:** 0.2

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The thickness of the lithosphere used to create additional boundary indicators to set specific boundary conditions for the lithosphere.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Use_20merged_20grids)=
### __Parameter name:__ Use merged grids
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to make the grid by gluing together two boxes, or just use one chunk to make the grid. Using two grids glued together is a safer option, since it forces the boundary conditions to be always applied to the same depth, but using one grid allows for a more flexible usage of the adaptive refinement. Note that if there is no cell boundary exactly on the boundary between the lithosphere and the mantle, the velocity boundary will not be exactly at that depth. Therefore, using a merged grid is generally recommended over using one grid.When using one grid, the parameter for lower repetitions is used and the upper repetitions are ignored.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/X_20extent)=
### __Parameter name:__ X extent
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Extent of the box in x-direction. Units: \si{\meter}.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/X_20periodic)=
### __Parameter name:__ X periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in X direction.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/X_20periodic_20lithosphere)=
### __Parameter name:__ X periodic lithosphere
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in X direction in the lithosphere.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/X_20repetitions)=
### __Parameter name:__ X repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in X direction of the lower box. The same number of repetitions will be used in the upper box.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Y_20extent)=
### __Parameter name:__ Y extent
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Extent of the box in y-direction. Units: \si{\meter}.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Y_20periodic)=
### __Parameter name:__ Y periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in Y direction.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Y_20periodic_20lithosphere)=
### __Parameter name:__ Y periodic lithosphere
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in Y direction in the lithosphere. This value is ignored if the simulation is in 2d.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Y_20repetitions)=
### __Parameter name:__ Y repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in Y direction of the lower box. If the simulation is in 3d, the same number of repetitions will be used in the upper box.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Y_20repetitions_20lithosphere)=
### __Parameter name:__ Y repetitions lithosphere
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in Y direction in the lithosphere. This value is ignored if the simulation is in 3d.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Z_20extent)=
### __Parameter name:__ Z extent
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Extent of the box in z-direction. This value is ignored if the simulation is in 2d. Units: \si{\meter}.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Z_20periodic)=
### __Parameter name:__ Z periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the box should be periodic in Z direction. This value is ignored if the simulation is in 2d.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Z_20repetitions)=
### __Parameter name:__ Z repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in Z direction of the lower box. This value is ignored if the simulation is in 2d.

(parameters:Geometry_20model/Box_20with_20lithosphere_20boundary_20indicators/Z_20repetitions_20lithosphere)=
### __Parameter name:__ Z repetitions lithosphere
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in Z direction in the lithosphere. This value is ignored if the simulation is in 2d.

(parameters:Geometry_20model/Chunk)=
## **Subsection:** Geometry model / Chunk
(parameters:Geometry_20model/Chunk/Chunk_20inner_20radius)=
### __Parameter name:__ Chunk inner radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Radius at the bottom surface of the chunk. Units: \si{\meter}.

(parameters:Geometry_20model/Chunk/Chunk_20maximum_20latitude)=
### __Parameter name:__ Chunk maximum latitude
**Default value:** 1.

**Pattern:** [Double -90...90 (inclusive)]

**Documentation:** Maximum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees.

(parameters:Geometry_20model/Chunk/Chunk_20maximum_20longitude)=
### __Parameter name:__ Chunk maximum longitude
**Default value:** 1.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Maximum longitude of the chunk. Units: degrees.

(parameters:Geometry_20model/Chunk/Chunk_20minimum_20latitude)=
### __Parameter name:__ Chunk minimum latitude
**Default value:** 0.

**Pattern:** [Double -90...90 (inclusive)]

**Documentation:** Minimum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees.

(parameters:Geometry_20model/Chunk/Chunk_20minimum_20longitude)=
### __Parameter name:__ Chunk minimum longitude
**Default value:** 0.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Minimum longitude of the chunk. Units: degrees.

(parameters:Geometry_20model/Chunk/Chunk_20outer_20radius)=
### __Parameter name:__ Chunk outer radius
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Radius at the top surface of the chunk. Units: \si{\meter}.

(parameters:Geometry_20model/Chunk/Latitude_20repetitions)=
### __Parameter name:__ Latitude repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in latitude. This value is ignored if the simulation is in 2d

(parameters:Geometry_20model/Chunk/Longitude_20repetitions)=
### __Parameter name:__ Longitude repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in longitude.

(parameters:Geometry_20model/Chunk/Radius_20repetitions)=
### __Parameter name:__ Radius repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in radius.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators)=
## **Subsection:** Geometry model / Chunk with lithosphere boundary indicators
(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20inner_20radius)=
### __Parameter name:__ Chunk inner radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Radius at the bottom surface of the chunk. Units: \si{\meter}.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20maximum_20latitude)=
### __Parameter name:__ Chunk maximum latitude
**Default value:** 1.

**Pattern:** [Double -90...90 (inclusive)]

**Documentation:** Maximum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20maximum_20longitude)=
### __Parameter name:__ Chunk maximum longitude
**Default value:** 1.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Maximum longitude of the chunk. Units: degrees.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20middle_20boundary_20radius)=
### __Parameter name:__ Chunk middle boundary radius
**Default value:** 1

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Radius at the top surface of the lower chunk, where it merges with the upper chunk. Units: \si{\meter}.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20minimum_20latitude)=
### __Parameter name:__ Chunk minimum latitude
**Default value:** 0.

**Pattern:** [Double -90...90 (inclusive)]

**Documentation:** Minimum latitude of the chunk. This value is ignored if the simulation is in 2d. Units: degrees.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20minimum_20longitude)=
### __Parameter name:__ Chunk minimum longitude
**Default value:** 0.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Minimum longitude of the chunk. Units: degrees.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Chunk_20outer_20radius)=
### __Parameter name:__ Chunk outer radius
**Default value:** 1.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Radius at the top surface of the chunk. Units: \si{\meter}.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Inner_20chunk_20radius_20repetitions)=
### __Parameter name:__ Inner chunk radius repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in radial direction for the lower chunk.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Latitude_20repetitions)=
### __Parameter name:__ Latitude repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in latitude. This value is ignored if the simulation is in 2d

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Longitude_20repetitions)=
### __Parameter name:__ Longitude repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in longitude.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Outer_20chunk_20radius_20repetitions)=
### __Parameter name:__ Outer chunk radius repetitions
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of cells in radial direction for the upper chunk.

(parameters:Geometry_20model/Chunk_20with_20lithosphere_20boundary_20indicators/Use_20merged_20grids)=
### __Parameter name:__ Use merged grids
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to make the grid by gluing together two boxes, or just use one chunk to make the grid. Using two grids glued together is a safer option, since it forces the boundary conditions to be always applied to the same depth, but using one grid allows for a more flexible usage of the adaptive refinement. Note that if there is no cell boundary exactly on the boundary between the lithosphere and the mantle, the velocity boundary will not be exactly at that depth. Therefore, using a merged grid is generally recommended over using one grid. When using one grid, the parameter for lower repetitions is used and the upper repetitions are ignored.

(parameters:Geometry_20model/Ellipsoidal_20chunk)=
## **Subsection:** Geometry model / Ellipsoidal chunk
(parameters:Geometry_20model/Ellipsoidal_20chunk/Depth)=
### __Parameter name:__ Depth
**Default value:** 500000.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Bottom depth of model region.

(parameters:Geometry_20model/Ellipsoidal_20chunk/Depth_20subdivisions)=
### __Parameter name:__ Depth subdivisions
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of subdivisions of the coarse (initial) mesh in depth.

(parameters:Geometry_20model/Ellipsoidal_20chunk/East_2dWest_20subdivisions)=
### __Parameter name:__ East-West subdivisions
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of subdivisions of the coarse (initial) mesh in the East-West direction.

(parameters:Geometry_20model/Ellipsoidal_20chunk/Eccentricity)=
### __Parameter name:__ Eccentricity
**Default value:** 8.1819190842622e-2

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Eccentricity of the ellipsoid. Zero is a perfect sphere, default (8.1819190842622e-2) is WGS84.

(parameters:Geometry_20model/Ellipsoidal_20chunk/NE_20corner)=
### __Parameter name:__ NE corner
**Default value:**

**Pattern:** [Anything]

**Documentation:** Longitude:latitude in degrees of the North-East corner point of model region.The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel.

(parameters:Geometry_20model/Ellipsoidal_20chunk/NW_20corner)=
### __Parameter name:__ NW corner
**Default value:**

**Pattern:** [Anything]

**Documentation:** Longitude:latitude in degrees of the North-West corner point of model region. The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel.

(parameters:Geometry_20model/Ellipsoidal_20chunk/North_2dSouth_20subdivisions)=
### __Parameter name:__ North-South subdivisions
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of subdivisions of the coarse (initial) mesh in the North-South direction.

(parameters:Geometry_20model/Ellipsoidal_20chunk/SE_20corner)=
### __Parameter name:__ SE corner
**Default value:**

**Pattern:** [Anything]

**Documentation:** Longitude:latitude in degrees of the South-East corner point of model region. The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel.

(parameters:Geometry_20model/Ellipsoidal_20chunk/SW_20corner)=
### __Parameter name:__ SW corner
**Default value:**

**Pattern:** [Anything]

**Documentation:** Longitude:latitude in degrees of the South-West corner point of model region. The North-East direction is positive. If one of the three corners is not provided the missing corner value will be calculated so all faces are parallel.

(parameters:Geometry_20model/Ellipsoidal_20chunk/Semi_2dmajor_20axis)=
### __Parameter name:__ Semi-major axis
**Default value:** 6378137.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The semi-major axis (a) of an ellipsoid. This is the radius for a sphere (eccentricity=0). Default WGS84 semi-major axis.

(parameters:Geometry_20model/Initial_20topography_20model)=
## **Subsection:** Geometry model / Initial topography model
(parameters:Geometry_20model/Initial_20topography_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** zero topography

**Pattern:** [Selection ascii data|function|prm polygon|zero topography ]

**Documentation:** Select one of the following models:

&lsquo;ascii data&rsquo;: Implementation of a model in which the surface topography is derived from a file containing data in ascii format. The following geometry models are currently supported: box, chunk. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;Topography [m]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Topography [m]&rsquo; in a 3d model, which means that there has to be a single column containing the topography. Note that the data in the input file needs to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle (longitude) in radians (between zero and $2\pi$, not between $-\pi$ corresponding to 180 degrees west, and $+\pi$ corresponding to 180 degrees east), and &lsquo;y&rsquo; by the polar angle in radians (between $0$ and $\pi$) measured positive from the north pole. The grid will be assumed to be a longitude-colatitude grid. Note that the order of spherical coordinates is &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;function&rsquo;: Implementation of a model in which the initial topography is described by a function in Cartesian or spherical coordinates.

&lsquo;prm polygon&rsquo;: An initial topography model that defines the initial topography as constant inside each of a set of polygonal parts of the surface. The polygons, and their associated surface elevation, are defined in the &lsquo;Geometry model/Initial topography/Prm polygon&rsquo; section.

&lsquo;zero topography&rsquo;: Implementation of a model in which the initial topography is zero.

(parameters:Geometry_20model/Initial_20topography_20model/Ascii_20data_20model)=
## **Subsection:** Geometry model / Initial topography model / Ascii data model
(parameters:Geometry_20model/Initial_20topography_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Geometry_20model/Initial_20topography_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_2d_%s.0.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Geometry_20model/Initial_20topography_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Geometry_20model/Initial_20topography_20model/Function)=
## **Subsection:** Geometry model / Initial topography model / Function
(parameters:Geometry_20model/Initial_20topography_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo; and &lsquo;spherical&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle.

(parameters:Geometry_20model/Initial_20topography_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Geometry_20model/Initial_20topography_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Geometry_20model/Initial_20topography_20model/Function/Maximum_20topography_20value)=
### __Parameter name:__ Maximum topography value
**Default value:** 2000.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum value the topography given by the function can take.

(parameters:Geometry_20model/Initial_20topography_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Geometry_20model/Initial_20topography_20model/Prm_20polygon)=
## **Subsection:** Geometry model / Initial topography model / Prm polygon
(parameters:Geometry_20model/Initial_20topography_20model/Prm_20polygon/Topography_20parameters)=
### __Parameter name:__ Topography parameters
**Default value:**

**Pattern:** [Anything]

**Documentation:** Set the topography height and the polygon which should be set to that height. The format is : "The topography height   extgreater The point list describing a polygon \& The next topography height   extgreater the next point list describing a polygon." The format for the point list describing the polygon is "x1,y1;x2,y2". For example for two triangular areas of 100 and -100 meters high set: &rsquo;100   extgreater 0,0;5,5;0,10 \& -100   extgreater 10,10;10,15;20,15&rsquo;. Units of the height are always in meters. The units of the coordinates are dependent on the geometry model. In the box model they are in meters, in the chunks they are in degrees, etc. Please refer to the manual of the individual geometry model to so see how the topography is implemented.

(parameters:Geometry_20model/Sphere)=
## **Subsection:** Geometry model / Sphere
(parameters:Geometry_20model/Sphere/Radius)=
### __Parameter name:__ Radius
**Default value:** 6371000.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Radius of the sphere. Units: \si{\meter}.

(parameters:Geometry_20model/Spherical_20shell)=
## **Subsection:** Geometry model / Spherical shell
(parameters:Geometry_20model/Spherical_20shell/Cells_20along_20circumference)=
### __Parameter name:__ Cells along circumference
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of cells in circumferential direction that are created in the coarse mesh in 2d. If zero, this number is chosen automatically in a way that produces meshes in which cells have a reasonable aspect ratio for models in which the depth of the mantle is roughly that of the Earth. For planets with much shallower mantles and larger cores, you may want to chose a larger number to avoid cells that are elongated in tangential and compressed in radial direction.

In 3d, the number of cells is computed differently and does not have an easy interpretation. Valid values for this parameter in 3d are 0 (let this class choose), 6, 12 and 96. Other possible values may be discussed in the documentation of the deal.II function GridGenerator::hyper_shell. The parameter is best left at its default in 3d.

In either case, this parameter is ignored unless the opening angle of the domain is 360 degrees. This parameter is also ignored when using a custom mesh subdivision scheme.

(parameters:Geometry_20model/Spherical_20shell/Custom_20mesh_20subdivision)=
### __Parameter name:__ Custom mesh subdivision
**Default value:** none

**Pattern:** [Selection none|list of radial values|number of slices ]

**Documentation:** Choose how the spherical shell mesh is generated. By default, a coarse mesh is generated with respect to the inner and outer radius, and an initial number of cells along circumference. In the other cases, a surface mesh is first generated and refined as desired, before it is extruded radially following the specified subdivision scheme.

(parameters:Geometry_20model/Spherical_20shell/Initial_20lateral_20refinement)=
### __Parameter name:__ Initial lateral refinement
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Initial lateral refinement for the custom mesh subdivision schemes.The number of refinement steps performed on the initial coarse surface mesh, before the surface is extruded radially. This parameter allows the user more control over the ratio between radial and lateral refinement of the mesh.

(parameters:Geometry_20model/Spherical_20shell/Inner_20radius)=
### __Parameter name:__ Inner radius
**Default value:** 3481000.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Inner radius of the spherical shell. Units: \si{\meter}.

:::{note}
The default value of 3,481,000 m equals the radius of a sphere with equal volume as Earth (i.e., 6371 km) minus the average depth of the core-mantle boundary (i.e., 2890 km).
:::

(parameters:Geometry_20model/Spherical_20shell/List_20of_20radial_20values)=
### __Parameter name:__ List of radial values
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** List of radial values for the custom mesh scheme. Units: $\si{m}$. A list of radial values subdivides the spherical shell at specified radii. The list must be strictly ascending, and the first value must be greater than the inner radius while the last must be less than the outer radius.

(parameters:Geometry_20model/Spherical_20shell/Number_20of_20slices)=
### __Parameter name:__ Number of slices
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of slices for the custom mesh subdivision scheme. The number of slices subdivides the spherical shell into N slices of equal thickness. Must be greater than 0.

(parameters:Geometry_20model/Spherical_20shell/Opening_20angle)=
### __Parameter name:__ Opening angle
**Default value:** 360.

**Pattern:** [Double 0...360 (inclusive)]

**Documentation:** Opening angle in degrees of the section of the shell that we want to build. The only opening angles that are allowed for this geometry are 90, 180, and 360 in 2d; and 90 and 360 in 3d. Units: degrees.

(parameters:Geometry_20model/Spherical_20shell/Outer_20radius)=
### __Parameter name:__ Outer radius
**Default value:** 6336000.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Outer radius of the spherical shell. Units: \si{\meter}.

:::{note}
The default value of 6,336,000 m equals the radius of a sphere with equal volume as Earth (i.e., 6371 km) minus the average depth of the mantle-crust interface (i.e., 35 km).
:::

(parameters:Geometry_20model/Spherical_20shell/Phi_20periodic)=
### __Parameter name:__ Phi periodic
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether the shell should be periodic in the phi direction.

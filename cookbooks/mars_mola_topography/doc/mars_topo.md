# Mars topography

*This section was contributed by Cedric Thieulot and Bart Root.*

From the original MOLA data (<https://attic.gsfc.nasa.gov/mola/>), we have created a topography file
(present in the folder of this cookbook) with a 1 degree resolution in both longitude and latitude
that is suitable for the Initial topography model.

Note that for the time being the Initial topography model is not compatible with the sphere geometry
model so the model geometry for this cookbook is a chunk with inner radius 2000 km and outer 
radius 3390 km, the latter being the actual radius of the planet.

The topography is then projected onto the surface of the planet, and as a first step both viscosity and
density are kept constant in the domain. No-slip boundary conditions are prescribed on the bottom, and
free slip is prescribed on all four sides while the top boundary is free.
The model is isothermal, and there are no other buoyancy forces.

```{figure-md} fig:mars-topo
<img src="topo3.*" width="90%" />

Mars topography.
For visualization purposes the topography has been scaled up by a factor 10.
```

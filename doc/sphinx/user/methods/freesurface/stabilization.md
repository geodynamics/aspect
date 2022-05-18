
# Free surface stabilization

#### Free surface stabilization

Small disequilibria in the location of a free surface can cause instabilities
in the surface position and result in a "sloshing" instability.
This may be countered with a quasi-implicit free surface integration scheme
described in (Kaus, M&uuml;hlhaus, and May 2010). This scheme enters the
governing equations as a small stabilizing surface traction that prevents the
free surface advection from overshooting its true position at the next time
step. ASPECT implements this stabilization, the
details of which may be found in (Kaus, M&uuml;hlhaus, and May 2010).

An example of a simple model which uses a free surface may be found in Section
[5.2.6][].

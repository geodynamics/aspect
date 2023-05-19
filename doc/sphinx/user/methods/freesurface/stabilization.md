
# Free surface stabilization

Small disequilibria in the location of a free surface can cause instabilities
in the surface position and result in a "sloshing" instability.
This may be countered with a quasi-implicit free surface integration scheme
described in {cite:t}`kaus:etal:2010`. This scheme enters the
governing equations as a small stabilizing surface traction that prevents the
free surface advection from overshooting its true position at the next time
step. ASPECT implements this stabilization, the
details of which may be found in {cite:t}`kaus:etal:2010`.

An example of a simple model which uses a free surface may be found in
{ref}`sec:cookbooks:free-surface`.

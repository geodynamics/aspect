# Using lower order elements for the temperature/compositional discretization

The default settings of ASPECT use quadratic
finite elements for the velocity. Given that the temperature and compositional
fields essentially (up to material parameters) satisfy advection equations of
the kind $\partial_t T +
\mathbf u \cdot \nabla T = \ldots$, it seems appropriate to also use quadratic
finite element shape functions for the temperature and compositional fields.

However, this is not mandatory. If you do not care about high accuracy in
these fields and are mostly interested in the velocity or pressure field, you
can select lower-order finite elements in the input file. The polynomial
degrees are controlled with the parameters in the *discretization* section of
the input file, see {ref}`parameters:Discretization`, in
particular by `Temperature polynomial degree` and
`Composition polynomial degree`.

As with the other parameters discussed above and below, it is worthwhile
comparing the results you get with different values of these parameters when
making a decision whether you want to save on accuracy in order to reduce
compute time. An example of how this choice affects the accuracy you get is
discussed in {ref}`sec:cookbooks:convection-box`.


# Approximation of the buoyancy term

The buoyancy term (right-hand side of the momentum equation) always uses the density that is provided by the material model (see {ref}`sec:extending:plugin-types:material-models`).
Depending on the material model, this density could for example depend on temperature and pressure (such as in ALA), or on temperature and depth (as in TALA); and the model can also be set up in a way that it uses density deviations from a reference state instead of a full density (see {ref}`sec:methods:static-v-dynamic`).

:::{note}
In the current version of ASPECT, it is the responsibility of the user to select a material model that is consistent with the formulation they want to use in their model.
In the future, we plan to make it more obvious which approximations are supported by a particular material model.
:::

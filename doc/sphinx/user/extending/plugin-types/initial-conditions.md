# Initial conditions for temperature and composition

The initial temperature conditions plugins are responsible for describing the initial
temperature distribution throughout the domain. Corresponding classes describe the
initial values for compositional fields.
Both in essence only have to provide
a function that for each point can return the initial temperature or composition. (Note that
the model {math:numref}`eq:stokes-1`&ndash;{math:numref}`eq:temperature` does not require
initial values for the pressure or velocity, and so there are no corresponding
plugins for these variables. However, if coefficients are
nonlinear, one can significantly reduce the number of initial nonlinear
iterations if a good guess for them is available; consequently,
ASPECT initializes the pressure with the
adiabatically computed hydrostatic pressure, and a zero velocity.)

To implement a new initial temperature model, you need to overload the
`aspect::InitialTemperature::Interface` class and use the
`ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL` macro to register your new class. The
implementation of the new class should be in namespace
`aspect::InitialTemperature`. For initial compositional values, the class and
macro names have to be adjusted correspondingly.

The functions you need to overload are extensively
discussed in the documentation of this interface class at
[aspect::InitialTemperature::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1InitialTemperature.html) and
[aspect::InitialComposition::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1InitialComposition.html).

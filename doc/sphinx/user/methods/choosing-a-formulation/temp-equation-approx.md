
# Temperature equation approximation

The density occurs multiple times in the temperature equation.
Depending on the selected approximation it is computed in one of two different ways.
Which of these options is used can be chosen in the parameter file in the subsection `Formulation/Temperature equation` (see also {ref}`parameters:Formulation/Temperature_20equation`):

-   "real density": Use the full density $\rho(p,T)$ that equals the one also used in the buoyancy term of the force balance equation; this is also the value that is computed by the material models when asked for the density,

-   "reference density profile": Use the density as computed for the reference profile (which can be constant, an adiabatic profile, or an entirely different function, and is determined by the adiabatic conditions model).

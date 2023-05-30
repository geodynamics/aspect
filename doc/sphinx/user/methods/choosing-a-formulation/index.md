(sec:methods:choosing-a-formulation)=
# Choosing a formulation

After discussing different reasonable approximations for modeling compressible or incompressible mantle convection, we will now describe the different steps one has to take to use one of these approximations in a computation.
As noted towards the beginning of {ref}`sec:methods:approximate-equations`, one should choose a formulation that is adequate for the situation one wants to simulate, taking into account that other authors in previous studies may have used formulations not because they were suited to the situation but also because, possibly, they did not have software available that implemented the most suitable formulation.

The choices involved in selecting a formulation include:

1.  Choosing an approximation for the mass conservation equation;

2.  Choosing an approximation for the density in the energy balance, and
    deciding which heating terms should be included;

3.  Formulating the buoyancy term in the material model to be used on the
    right-hand side of the momentum equation;

4.  Prescribing a suitable reference state for the temperature, pressure, and
    density; i.e. the adiabatic profile, if necessary for the approximations
    chosen in the first three steps.

All of these choices can be made in the input file by selecting the corresponding parameters (see {ref}`parameters:Formulation` and {ref}`parameters:Adiabatic_20conditions_20model`).
A description of how to run ASPECT and the basic structure of the input file can be found in {ref}`cha:run-aspect`.

:::{toctree}
mass-conservation-approx.md
temp-equation-approx.md
buoyancy-approx.md
adiabatic-profile.md
combined-formulations.md
:::

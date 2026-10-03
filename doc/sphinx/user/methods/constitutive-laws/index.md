(sec:methods:constitutive-laws)=
# Constitutive Models (Rheology)

This section provides a detailed description of the various constitutive models available in ASPECT and how they are organized with Material models. Material models in ASPECT are responsible for calculating various coefficients in the equations that ASPECT solves, which includes terms related to the equation of state (e.g., density, thermal conductivity) and viscosity. Constitutive models (hereby refer to as rheological models) primarily handle the computation of viscosity, but also often require computing additional terms such as viscoelastic stress and finite strain, which are then used to modify the viscosity.

Historically, Material models handled the computation of equation of state and rheological coefficients within a single file. However, as the complexity and number of material models expanded distinct rheology and equation of state modules were developed to simplify the code structure and create reusable elements across material models. The sections below describe both the formulations of distinct rheology modules and how they are combined within material models in ASPECT.

:::{toctree}
constant-viscosity.md
temperature-dependent-viscosity.md
scaling-experimental-prefactors.md
overview-invariants.md
diffusion-creep.md
dislocation-creep.md
peierls-creep.md
simple-composite-viscous-creep.md
grain-size-evolution.md
plasticity.md
strain-dependent-rheology.md
linear-viscoelasticity.md
viscoelastic-plastic-deformation.md
iterative-composite-rheology.md
anisotropic-viscosity.md
:::

This benchmark is a technical benchmark for the entropy advection method.
The model setup is described in entropy_adiabat.prm and requires the shared
libraries in the plugins/ directory. In particular it requires a material
model (plugins/entropy_model) that is formulated using pressure and entropy as independent variables, a different assembler that assembles the entropy equation (plugins/entropy_advection) and a different way to compute the adiabatic background profile (plugins/compute_entropy_profile).

(sec:benchmarks:entropy_adiabat)=
# Entropy adiabat benchmark

[This folder](https://github.com/geodynamics/aspect/tree/main/benchmarks/entropy_adiabat)
contains a technical benchmark for the entropy advection method.
The model setup is described in `entropy_adiabat.prm` and requires the shared
libraries in the `plugins/` directory. In particular it requires a material
model (`plugins/entropy_model.cc/h`) that is formulated using pressure and entropy as independent variables, a different assembler that assembles the entropy equation (`plugins/entropy_advection.cc/h`) and a different way to compute the adiabatic background profile (`plugins/compute_entropy_profile.cc/h`).

The benchmark and its results are described in {cite:t}`dannberg:etal:2022`. Please remember to cite this paper when using the entropy advection method.

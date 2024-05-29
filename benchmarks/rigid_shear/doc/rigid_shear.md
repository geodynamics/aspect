# The rigid shear benchmark

*This section was contributed by R. Gassmoeller.*

This directory contains the "Rigid shear" benchmark for which a known solution
is available. The necessary code lives in the plugin/ subdirectory, while
instantaneous, steady-state and transient versions of the benchmark are in the
corresponding subdirectories. Each directory has a "run_all_models" shell
script that will run the benchmark setups and a "create_formatted_tables.sh"
script that formats the output into readable text files. For a more
detailed description of the benchmark see {cite:t}`gassmoller:etal:2019`, "Evaluating
the Accuracy of Hybrid Finite Element/Particle-In-Cell Methods for
Modeling Incompressible Stokes Flow".

Additionally, the "transient" directory contains an extension of the benchmark
to time-dependent flow. The benchmark and its results are described in
{cite:t}`gassmoeller:etal:2023`, "Benchmarking the accuracy of higher order
particle methods in geodynamic models of transient flow",
see there for a more detailed description.

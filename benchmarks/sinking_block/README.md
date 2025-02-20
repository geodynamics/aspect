# The sinking block benchmark files

[This folder](https://github.com/geodynamics/aspect/tree/main/benchmarks/sinking_block)
contains the sinking block benchmark.

In order to run the benchmark for many resolution, density and viscosity combinations use the
`run_benchmark` shell script. This makes use of the `blank.prm` file which serves as a basis for
all the input files generated inside the triple loop.
The `sinking_block.prm` input file is the standalone prm file for the sinking block benchmark.

For more information see [](doc/sinking_block).

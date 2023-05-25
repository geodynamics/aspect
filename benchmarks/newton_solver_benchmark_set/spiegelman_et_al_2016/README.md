# Newton Solver Benchmark Set - Spiegelman et at. (2016)

The files in [this directory](https://github.com/geodynamics/aspect/tree/main/benchmarks/newton_solver_benchmark_set/spiegelman_et_al_2016)
can be used to recreate the Spiegelman et al. (2016) flow figures from
{cite:t}`fraters:etal:2019`. Adjust both the `metabash.sh` and `bash.sh` files to
reflect the parameter search you are interested in and run the metabash script.
If you have run exactly the same runs as shown in the paper, you can use the
`cp_P-minLT-res-it-vel_BV.sh` script to recreate the figures in that paper.
This script requires that the drucker prager compositions plugin is installed.
It should be possible to replicate the same behavior with the visco plastic
plugin.

# Newton Solver Benchmark Set - Tosi et at. (2015)

The files in [this directory](https://github.com/geodynamics/aspect/tree/main/benchmarks/newton_solver_benchmark_set/tosi_et_al_2015)
can be used to recreate the Tosi et al. (2015) benchmark figures
from {cite:t}`fraters:etal:2019`. Adjust both `metabash.sh` and `bash.sh`
to reflect the parameter search you are interested in and run the metabash
script. If you have run exactly the same runs as shown in the paper, you can
use the gnuplot files in plotfiles to recreate the figures in that paper. This
benchmark requires the shared library from the `benchmarks/tosi_et_al_2015_gcubed`
folder in the ASPECT repository.

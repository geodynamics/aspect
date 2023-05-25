# Newton Solver Benchmark Set - Nonlinear Channel Flow

The files in [this directory](https://github.com/geodynamics/aspect/tree/main/benchmarks/newton_solver_benchmark_set/nonlinear_channel_flow)
can be used to recreate the nonlinear channel flow figures from
{cite:t}`fraters:etal:2019`. Adjust both `metabash.sh` and `bash.sh` to reflect
the parameter search you are interested in and run the metabash script. If you
have run exactly the same runs as shown in the paper, you can use the gnuplot
files in `plotfiles/` to recreate the figures in that paper. This benchmark
requires the normal nonlinear channel flow benchmark to be compiled.

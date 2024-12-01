# Newton Solver Benchmark Set - Nonlinear Channel Flow

The files in [this directory](https://github.com/geodynamics/aspect/tree/main/benchmarks/newton_solver_benchmark_set/nonlinear_channel_flow)
can be used to recreate the nonlinear channel flow figures from
{cite:t}`fraters:etal:2019`.
Start by compiling the `simple_nonlinear` plugin in this directory as described (at the example of another benchmark)
in {ref}`sec:benchmark-run` .
After compiling the plugin adjust `metabash.sh` and/or `bash.sh` to reflect
the parameter search you are interested in and run the script.
Running the `bash.sh` script without changes will recreate the model runs of Fig. 1 and Fig. 2 in {cite:t}`fraters:etal:2019`
(by default the script will assume your processor has 4 compute cores, but you can adjust the number in the script).
After running the model series you can use the gnuplot files in `plotfiles/` to recreate the following figures.

```{figure-md} fig:benchmark-newton-nonlinear-channel-flow-tractions
<img src="doc/figure_t.png" />

Convergence history for several methods for a rheology with n = 3 where in- and outflow are described by prescribing the traction. Top row: Computations in which we switch abruptly from Picard iterations to Newton iterations. Bottom row: Same as top, but using the residual scaling method in which we gradually introduce the Newton method. Left-hand column: Unmodified Newton iterations. Right-hand column: Results where we applied the SPD stabilizations to the Newton matrix. Horizontal axes: number of the non-linear (outer) iteration. Vertical axes: non-linear residual.
```

```{figure-md} fig:benchmark-newton-nonlinear-channel-flow-velocities
<img src="doc/figure_v.png" />

Convergence history for several methods for a rheology with n = 3 where in- and outflow are described by prescribing the velocity. Top row: Computations in which we switch abruptly from Picard iterations to Newton iterations. Bottom row: Same as top, but using the residual scaling method in which we gradually introduce the Newton method. Left-hand column: Unmodified Newton iterations. Right-hand column: Results where we applied the SPD stabilizations to the Newton matrix. Horizontal axes: number of the non-linear (outer) iteration. Vertical axes: non-linear residual.
```

The nonlinear convergence behavior is shown in {numref}`fig:benchmark-newton-nonlinear-channel-flow-tractions` and {numref}`fig:benchmark-newton-nonlinear-channel-flow-velocities`.

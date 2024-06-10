# Benchmarking run time

When changing code it can sometimes be useful to if and how much your code changes make the code faster or slower. To test this, one could compile and run aspect-release compiled from the main branch and from the branch with the changed code (e.g. feature branch) and time both. The problems is that is is not very reliable, given that your computer might be doing something on the background while you are running one, but not while you are running the other.

```{note}
When running timing benchmarks, always run in release mode.
```

One solution is to this problem is to run both executables after each other in a random order. This is what the program [cbdr](https://crates.io/crates/cbdr) can do for you. It can then also analyze and report back the result.

This program can be installed with [cargo](https://crates.io/) (see Install Cargo button), with `cargo install cbdr`. To generate the run data, program is run with `cbdr sample --timeout=900s main:/path/to/build-main/aspect-release /path/to/test.prm features:/path/to/build-feature/aspect-release /path/to/test.prm > results.csv`. This first run runs both the main and feature branch executable to "warm up" and then it keeps randomly selection one of the two executables and write the timings out to `results.csv` until the timer reaches 900 seconds.

The `results.csv` can then be summarized with `cbdr analyze <results.csv`, which produces the following output:

```
           main           feature        difference (99.9% CI)
sys_time   0.868 ± 0.071  0.916 ± 0.051  [  +1.5% ..   +9.7%]
user_time  7.216 ± 0.133  7.112 ± 0.100  [  -2.4% ..   -0.5%]
wall_time  6.556 ± 0.045  6.498 ± 0.046  [  -1.3% ..   -0.5%]
samples    66             74
```

You can ignore the sys_time in the output above. I just shows how much time the programs spends on kernel calls. The user time is how much the program actually spending on the cpu and the wall time is how much time the program takes if you measure it by have a stop-watch next to your computer. In this case it is nearly the same, but if you run it with many processors, the user/cpu time might be much larger then the wall time.

The `cbdr` program can also create a graphical output in a vega-lite format through the command `cbdr plot <results.csv > results.vl`. This can then be converted to a png through [vl-convert](https://crates.io/crates/vl-convert) with the command `vl-convert vl2png --input results.vl --output results.png`.

```{figure-md} fig:benchmark_run_test_graph
<img src="../../_static/images/timing_benchmark_cbdr_graph.*" alt="benchmark runtime graph"  width="100%"/>

Example graph created for the same run as the table above.
```


# Geometric Multigrid

### Geometric Multigrid

ASPECT can optionally use a Geometric Multigrid
Solver (GMG) for the efficient solution of the Stokes system (velocity and
pressure). When used correctly, this can reduce the compute time spent in the
solver by about a factor of 3, and decrease the memory requirements by a
factor of 10. For more details about the method see&nbsp;(Clevenger and
Heister 2021; Clevenger et al. 2020).

To take advantage of the GMG solver, you need to:

1.  Enable it in your parameter file, namely:

    ``` prmfile
    ```

    See&nbsp;{ref}`parameters:Solver_20parameters/Stokes_20solver_20parameters`45]
    for other parameters that influence the solver behavior.

2.  The GMG solver requires that the viscosity is averaged, either as a
    constant (for example by using harmonic averaging) or as a $Q_1$
    averaging. (See {ref}`5.2.8][] for more about averaging.)
    Averaging other properties is optional. You can use

    ``` prmfile
    ```

    for example. Note that $Q_1$ averaging is a bit slower than averaging to a
    constant per cell, but it might provide more accurate solutions.

3.  Run in release mode. The GMG solver depends on running optimized code, so
    using optimized mode is more important than for other parts of ASPECT to
    get good performance. (Of course, the GMG solver also runs in debug mode,
    and you should do so while setting up a model. You will just not get the
    same speedup from the non-GMG to the GMG solver in debug mode as you get
    in release mode.)

4.  Enable vectorization optimizations. The GMG solver takes advantage of
    special instructions (AVX2, AVX512) in modern CPUs and requires these do
    be enabled when compiling deal.II. This can
    be achieved by passing the compiler flag `CMAKE_CXX_FLAGS="-march=native"`
    to CMake or setting `NATIVE_OPTIMIZATIONS=true` in `candi.cfg` when using
    candi (see&nbsp;[3][] for more information). When you have vectorization
    enabled, ASPECT will report something like this:

    ``` ksh
    -----------------------------------------------------------------------------
    -- This is ASPECT, the Advanced Solver for Problems in Earth's ConvecTion.
    --     . version 2.3.0-pre
    --     . using deal.II 9.2.0
    --     .       with 64 bit indices and vectorization level 3 (512 bits)
    --     . using Trilinos 12.10.1
    --     . using p4est 2.2.0
    --     . running in OPTIMIZED mode
    --     . running with 114688 MPI processes
    -----------------------------------------------------------------------------

    Vectorization over 8 doubles = 512 bits (AVX512), VECTORIZATION_LEVEL=3
    ```

    Without optimizations enabled, the output will be "and vectorization
    level 1 (128 bits)" in the fourth line above.

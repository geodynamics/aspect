(sec:benchmarks:solkz)=
# The SolKz Stokes benchmark

The SolKz benchmark is another variation on the same theme as {ref}`sec:benchmarks:solcx`,
it solves a Stokes problem with a spatially variable
viscosity, but this time the viscosity is not a discontinuous function.
Instead, it grows exponentially with the vertical coordinate so that its
overall variation is again $10^6$. The forcing is again chosen by imposing a
spatially variable density variation. For details, refer again to {cite}`duretz:etal:2011`.

The following input file, only a small variation of the one in the previous
section, solves the benchmark (see [benchmarks/solkz/](https://github.com/geodynamics/aspect/tree/main/benchmarks/solkz)):

```{literalinclude} solkz.prm

```

The output when running ASPECT on this parameter file looks similar to the one shown
for the SolCx case. The solution when computed with one more level of global
refinement is visualized in {numref}`fig:solkz1` and {numref}`fig:solkz2`.
The velocity solution computed with three more levels of global refinement
and plotted over the viscosity field is shown in {numref}`fig:solkz3`.

```{figure-md} fig:solkz1
<img src="solkz-solution.*" alt="SolKz Stokes benchmark. The density perturbation field overlaid with velocity vectors. The viscosity grows exponentially in the vertical direction, leading to small velocities at the top despite the large density variations." />

SolKz Stokes benchmark. The density perturbation field overlaid with velocity vectors. The viscosity grows exponentially in the vertical direction, leading to small velocities at the top despite the large density variations.
```
```{figure-md} fig:solkz2
<img src="solkz-solution-pressure.*" title="fig:" id="fig:solkz2" alt="SolKz Stokes benchmark pressure." />

SolKz Stokes benchmark pressure.
```

```{figure-md} fig:solkz3
<img src="solkz-solution-viscosity.*" id="fig:solkz3" alt="SolKz Stokes benchmark. Another view of the velocity vectors, this time plotted over the viscosity field." />

SolKz Stokes benchmark. Another view of the velocity vectors, this time plotted over the viscosity field.
```

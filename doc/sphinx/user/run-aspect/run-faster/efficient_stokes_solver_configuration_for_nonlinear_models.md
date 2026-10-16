# Efficient Stokes solver configuration for nonlinear models

## Background

One of the major challenges in geodynamic modeling is that the constitutive equations can be highly nonlinear.
A common cause is the nonlinear rheological behavior of mantle rocks, which can make viscosity depend on a power of the strain-rate invariant, as in dislocation creep and plastic yielding.

ASPECT already provides nonlinear solver schemes to address this numerical challenge. Rather than diving into their scientific meaning or technical implementation, this cookbook focuses on a more specific and practical question that every user must deal with: How should the existing nonlinear solver be configured to achieve robust and efficient convergence in a simulation?

The primary consideration is, of course, choosing Stokes solver parameters that ensure convergence. However, even when convergence is achieved, the solver configuration may still have the problem of being overly expensive. For example, the convergence tolerance may be set unnecessarily tight, or the chosen solver scheme may be too rigid for the problem.

Ideally, the solver configuration should be just sufficient to compute a solution that is sufficiently accurate (i.e., within the desired numerical tolerance) while minimizing computational cost.


## A short introduction to ASPECT's nonlinear solver scheme

When a problem is nonlinear, each time step requires solving a nonlinear Stokes system.
Here we show an example with the set up in Section&nbsp;{ref}`sec:cookbooks:crustal-deformation` cookbook.
The nonlinear solver first defines the required solution accuracy through the parameter `Nonlinear solver tolerance`. It then achieves this target through a sequence of nonlinear iterations, with each iteration bringing the current numerical approximation closer to the exact solution, stopping when the desired accuracy has been reached.
During each nonlinear iteration, this is accomplished by solving a linearized Stokes system, as configured in the `Stokes solver parameters` section.
Consequently, the accuracy required for this linear solve affects both the computational cost of an individual nonlinear iteration and the total number of nonlinear iterations required to reach convergence at a timestep.

To make things more complicated, each linear solve in ASPECT is itself a mixture of cheap and expensive Stokes solver steps. As their names suggest, cheap Stokes solver steps are less robust but computationally inexpensive, whereas expensive Stokes solver steps are more robust but computationally more costly.

Ideally, the configuration of the "inner" linear solver should be well matched to that of the "outer" nonlinear solver, so that the linear solver is just rigid enough to converge the nonlinear problem to the desired tolerance.
Therefore the computational cost is minimized.

```text
set Nonlinear solver scheme                         = single Advection, iterated Stokes
set Nonlinear solver tolerance                      = 2e-6
set Max nonlinear iterations                        = 100
set CFL number                                      = 0.5

subsection Solver parameters
  subsection Stokes solver parameters
    set Number of cheap Stokes solver steps = 0
    set Number of expensive Stokes solver steps = 5000
    set Linear solver tolerance = 1e-5
  end
end
```

## Numerical hardness in the initial time steps

Unfortunately, a single solver configuration often does not provide optimal performance throughout an entire simulation.

The reason is that many nonlinear models are most difficult to solve during the first few time steps. The initial conditions are often far from mechanical equilibrium, leading to large changes in the velocity, pressure, and viscosity fields during the first nonlinear iterations.
On the other hand, as the simulation progresses, the solution evolves toward a more stable state, making the problem easier to solve.

We again use the setup in the&nbsp;{ref}`sec:cookbooks:crustal-deformation` cookbook as an example to demonstrate this. In a typical log output generated during the simulation, each time step contains a sequence of nonlinear iterations, with each iteration represented by a block such as the one below.

```text
   Solving Stokes system (AMG)... 0+177 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 3: 0.0259232
```

This block indicates that nonlinear iteration 3 uses the AMG solver, requiring 0 cheap Stokes iterations and 177 expensive Stokes iterations. After this nonlinear iteration, the relative nonlinear residual is reduced to 0.0259232.

As shown below, the number of nonlinear iterations is much larger during the timestep 0 (i.e. 59) than during the later timesteps 41 (i.e. 10). Furthermore, within each nonlinear iteration, the linear Stokes solve also requires more iterations, particularly during the first few nonlinear iterations of time step 0 (e.g. 159, 177, 125).


```text
Number of mesh deformation degrees of freedom: 2,754
   Solving mesh displacement system... 0 iterations.
*** Timestep 0:  t=0 years, dt=0 years
   Solving mesh displacement system... 0 iterations.
   Skipping temperature solve because RHS is zero.
   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+20 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 1: 1

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+159 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 2: 0.0893438

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+177 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 3: 0.0259232

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+125 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 4: 0.00630584

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+103 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 5: 0.00174308

   ...

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+15 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 58: 2.01973e-06

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+15 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 59: 1.98143e-06
```

```text
*** Timestep 41:  t=1e+06 years, dt=1.23784e-10 years
   Solving mesh displacement system... 15 iterations.
   Skipping temperature solve because RHS is zero.
   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+29 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 1: 5.82296e-05

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+28 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 2: 1.61619e-05

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+26 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 3: 8.10853e-06

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+26 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 4: 5.27581e-06

   ...

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+23 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 9: 2.12497e-06

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 0+23 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 10: 1.88002e-06
```

We record the wall-clock time required to run the model to 1 Myr. This job is launched with two cpus of 11th Gen Intel(R) Core(TM) i5-11600KF @ 3.90GHz.

```text
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |       185s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble Stokes system          |      1242 |      5.05s |       2.7% |
| Assemble temperature system     |        43 |     0.212s |      0.11% |
| Build Stokes preconditioner     |      1242 |      8.14s |       4.4% |
| Initialization                  |         1 |     0.112s |         0% |
| Mesh deformation                |        43 |     0.662s |      0.36% |
| Mesh deformation initialize     |         3 |    0.0394s |         0% |
| Postprocessing                  |        43 |      0.26s |      0.14% |
| Refine mesh structure, part 1   |        43 |    0.0994s |         0% |
| Refine mesh structure, part 2   |        43 |    0.0962s |         0% |
| Setup dof systems               |        44 |     0.156s |         0% |
| Setup initial conditions        |         2 |   0.00554s |         0% |
| Setup matrices                  |        43 |     0.242s |      0.13% |
| Solve Stokes system             |      1242 |       169s |        91% |
+---------------------------------+-----------+------------+------------+
```

In general, we want a robust solver configuration that is capable of handling the difficult initial stage of the simulation.
However, such a configuration can also be unnecessarily expensive because the same strict linear solver settings are retained throughout the simulation, even after the nonlinear problem has become much easier to solve. This is the case in the previous example, which uses a large number of expensive Stokes solver iterations throughout the simulation.

Moreover, this type of costly solver configuration often results from a typical workflow. When the first simulation attempts fail because of solver convergence issues, it is tempting to simply adopt a more robust solver configuration without a second thought on the transient nature of the non-convergence.

## Strategy A: Introducing a number of cheap Stokes iterations

One simple strategy is to always include a reasonable number of inexpensive Stokes iterations, even if they provide little benefit during the initial, difficult transient stage.
The rationale is that, once these transients have passed, many nonlinear iterations can already be resolved using the inexpensive Stokes solver.
This behavior can be exploited by including a number of inexpensive Stokes iterations to accelerate the later stages of the simulation.

For example, below is a slight modification of the previous example:


```text
subsection Solver parameters
  subsection Stokes solver parameters
    set Stokes solver type = block AMG
    set Linear solver tolerance = 1e-9
    set Number of cheap Stokes solver steps = 200
    set Maximum number of expensive Stokes solver steps = 5000
    set GMRES solver restart length = 200
  end
end
```

Notably, these values are not chosen with the goal of converging the solution during the initial, more challenging time steps using only the inexpensive Stokes solver.
As shown in the simulation log below, the nonlinear iterations during the initial time steps still require the expensive Stokes solver to achieve convergence.


```text
*** Timestep 0:  t=0 years, dt=0 years
   Solving mesh displacement system... 0 iterations.
   Skipping temperature solve because RHS is zero.
   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 39+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 1: 1

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 200+154 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 2: 0.0893438

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 200+168 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 3: 0.0259232

   ...

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 46+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 58: 2.01977e-06

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 49+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 59: 1.98146e-06
```

The motivation for including a reasonable number of inexpensive Stokes solver iterations is that these iterations can successfully converge the nonlinear solve during later time steps:

```text
*** Timestep 41:  t=1e+06 years, dt=1.23784e-10 years
   Solving mesh displacement system... 15 iterations.
   Skipping temperature solve because RHS is zero.
   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 98+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 1: 5.91113e-05

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 98+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 2: 1.64402e-05

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 95+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 3: 8.25114e-06

   ...

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 83+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 9: 2.05799e-06

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 81+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 10: 1.79089e-06
```

In effect, this strategy replaces the "expensive" Stokes step with "cheap" Stokes step in later timesteps.
Even from the name of the steps itself as "cheap" versus "expensive", one could expect this to be faster.
Indeed, compared with the original prm, this strategy saves more than 50 % of the computational cost.

```text
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      88.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble Stokes system          |      1277 |      5.19s |       5.9% |
| Assemble temperature system     |        43 |     0.211s |      0.24% |
| Build Stokes preconditioner     |      1277 |      8.31s |       9.4% |
| Initialization                  |         1 |    0.0896s |       0.1% |
| Mesh deformation                |        43 |     0.647s |      0.73% |
| Mesh deformation initialize     |         3 |    0.0385s |         0% |
| Postprocessing                  |        43 |     0.283s |      0.32% |
| Refine mesh structure, part 1   |        43 |    0.0994s |      0.11% |
| Refine mesh structure, part 2   |        43 |    0.0946s |      0.11% |
| Setup dof systems               |        44 |     0.155s |      0.18% |
| Setup initial conditions        |         2 |    0.0054s |         0% |
| Setup matrices                  |        43 |     0.238s |      0.27% |
| Solve Stokes system             |      1277 |      71.8s |        81% |
+---------------------------------+-----------+------------+------------+
```

In this example, one could indeed enlarge the "Number of cheap Stokes solver steps" to converge the nonlinear iteration solely with cheap Stokes solver, for instance, in the current example a number of 2000 would converge most nonlinear iterations.
However, this is not the purpose of this strategy, which is to reduce computation cost.
There are no substantial differences in the wall-clock time with a number of 200 vs 2000.


## Strategy B: Introducing a linear solver failure strategy to continue with the nonlinear solver

Strategy A demonstrates a philosophy for handling the challenging initial time steps by prioritizing any configuration that allows the nonlinear solver to converge, while choosing linear solver parameters that provide good performance during the later, easier stages of the simulation.

However, Strategy A still requires the linear solver to be sufficiently robust to achieve the required linear solver tolerance.
Otherwise, the simulation terminates when a linear solve fails.

This philosophy can be extended further by allowing the linear solver to fail occasionally during nonlinear iterations.
Instead of terminating the simulation, the nonlinear solver continues using the solution from the failed linear solve and proceeds to the next nonlinear iteration, where another linear solve is attempted. The underlying philosophy here is that in a nonlinear scheme, it is not necessary to solve each linear system *exactly*; all we care about is that the computed solution of the linear system helps us move *closer* to the solution of the nonlinear problem. That might be the case even if we can't reach the desired linear solver tolerance.
ASPECT provides the `Linear solver failure strategy` option to support this behavior:

```text
set Linear solver failure strategy = continue with nonlinear solver
subsection Solver parameters
  subsection Stokes solver parameters
    set Stokes solver type = block AMG
    set Linear solver tolerance = 1e-9
    set Number of cheap Stokes solver steps = 200
    set Maximum number of expensive Stokes solver steps = 0
    set GMRES solver restart length = 200
  end
end
```

This option allows the linear solver to be configured entirely based on what is needed for the later and easier time steps (i.e. 200 cheap Stokes steps), without requiring every linear solve during the initial difficult stage to satisfy the required linear tolerance.

The results for time step 0 are shown below. The linear solve fails during nonlinear iterations 2 and 3 after reaching the limit of 200 cheap Stokes solver steps. Compared with the corresponding results from Strategy A, the failure is caused by setting the number of expensive Stokes solver steps to zero.

However, these linear solver failures are not critical, as evidenced by the continued reduction of the nonlinear residual after each failed solve. Because the `Linear solver failure strategy` is set to "continue with the nonlinear solver", the nonlinear iterations proceed normally and eventually converge for the time step.

```text
*** Timestep 0:  t=0 years, dt=0 years
   Solving mesh displacement system... 0 iterations.
   Skipping temperature solve because RHS is zero.
   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 39+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 1: 1

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 200+0 iterations.
 linear solver failed, continuing
      Relative nonlinear residual (Stokes system) after nonlinear iteration 2: 0.0893438

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 200+0 iterations.
 linear solver failed, continuing
      Relative nonlinear residual (Stokes system) after nonlinear iteration 3: 0.0260592

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 196+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 4: 0.00628391

   ...

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 49+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 58: 2.02388e-06

   Rebuilding Stokes preconditioner...
   Solving Stokes system (AMG)... 49+0 iterations.
      Relative nonlinear residual (Stokes system) after nonlinear iteration 59: 1.98477e-06
```

During the later and easier time steps, Strategy B behaves very similarly to Strategy A. Most linearized Stokes solves are completed using only the cheap Stokes solver, without requiring any expensive solver steps. Consequently, the overall wall-clock time is nearly identical to that of Strategy A:

```text
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      77.3s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble Stokes system          |      1141 |      4.63s |         6% |
| Assemble temperature system     |        43 |     0.209s |      0.27% |
| Build Stokes preconditioner     |      1141 |       7.4s |       9.6% |
| Initialization                  |         1 |    0.0968s |      0.13% |
| Mesh deformation                |        43 |     0.655s |      0.85% |
| Mesh deformation initialize     |         3 |    0.0388s |         0% |
| Postprocessing                  |        43 |     0.242s |      0.31% |
| Refine mesh structure, part 1   |        43 |       0.1s |      0.13% |
| Refine mesh structure, part 2   |        43 |    0.0943s |      0.12% |
| Setup dof systems               |        44 |     0.156s |       0.2% |
| Setup initial conditions        |         2 |   0.00545s |         0% |
| Setup matrices                  |        43 |     0.238s |      0.31% |
| Solve Stokes system             |      1141 |      62.2s |        80% |
+---------------------------------+-----------+------------+------------+

WARNING: During this computation 5 linear solver failures occurred!
```

The total number of linear solver failures is also reported at the end of the log output. In this relatively small example with a low-resolution model, there are only five failed linear solves in total. For more realistic and higher-resolution models, one would expect this number to be significantly larger when using this strategy.

While this strategy has proven effective in a variety of applications, it still relies on the assumption that a valid solution exists for every time step, regardless of whether the nonlinear system is easy or difficult to solve. If the nonlinear system itself becomes ill-posed, this strategy alone is unlikely to improve convergence.

Another limitation is that one must have a reasonable estimate of the number of cheap Stokes solver steps required for the later and easier time steps. This value cannot be too far off. The underlying rationale is that if the selected number of cheap Stokes solver steps is sufficient to converge the later and easier linear solves, it should also be sufficient to substantially reduce the residual of the earlier and more difficult linear solves despite the failure.

Conceptually, the usage of this strategy is not limited to the difficult initial time steps. It can help when a simulation suddenly becomes more difficult later in the computation. For example, in a subduction model, events such as slab folding or slab detachment may suddenly make a time step significantly more challenging to solve. In such cases, allowing a few linear solver failures while continuing the nonlinear iterations can potentially converge the nonlinear solver, and release the pain of searching for a more rigid solver scheme.

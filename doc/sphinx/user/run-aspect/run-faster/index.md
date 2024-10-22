
# Making ASPECT run faster

When developing ASPECT, we are guided by the
principle that the default for all settings should be *safe*. In particular,
this means that you should get errors when something goes wrong, the program
should not let you choose an input file parameter so that it doesn't
make any sense, and we should solve the equations to best ability without
cutting corners. The goal is that when you start working with
ASPECT that we give you the best answer we can. The
downside is that this also makes ASPECT run
slower than may be possible. This section describes ways of making
ASPECT run faster - assuming that you know what
you are doing and are making conscious decisions.



:::{toctree}
debug-vs-optimized.md
solver-tolerance.md
preconditioner-tolerance.md
lower-order-elements.md
limiting-postprocessing.md
pressure-norm-off.md
regularize.md
multithreading.md
file-system-io.md
:::

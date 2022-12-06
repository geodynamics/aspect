# Using multithreading

In most cases using as many MPI processes as possible is the optimal
parallelization strategy for ASPECT models, but
if you are limited by the amount of MPI communication it can be beneficial to
use multiple threads per MPI process. While not utilized by our linear
solvers, this parallelization can speed up the assembly of the system
matrices, e.g. by around 10-15% if you utilize unused logical cores, or nearly
linearly if you use otherwise unused physical cores. This can also reduce the
performance cost if you are memory limited and need to run your model on less
than the available number of cores per node on a cluster to increase the
available memory per core. Running with for example two threads per process
will offset some of the performance loss you will see in these situations.

Multithreading is controlled by setting the command line parameter `-j` or
`--threads`. If the parameter is not set,
ASPECT will create exactly one thread per MPI
process, i.e. multithreading is disabled. Appending the parameter allows
ASPECT to spawn several threads per MPI process. Note
that the internally used TBB library will determine the number of threads
based on the number of available cores, i.e., if you start 2&nbsp;MPI
processes on a quadcore machine with hyperthreading (8 logical cores),
ASPECT will spawn 4 threads on each MPI process. Also
note that there is no guarantee that the final number of threads will exactly
match the number of available logical cores if you start with a number of
processes that is not a divisor of your logical cores (e.g. 3 MPI processes
for 8 logical cores).

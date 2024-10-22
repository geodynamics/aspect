(sec:run-aspect:run-faster:file-system-io)=
# File system I/O

Depending on how exactly you run ASPECT and how large your computations are,
ASPECT can create output that can be hundreds of gigabytes or more. This output
also includes individual files that can be many gigabytes or more, for example
for graphical output and, in particular, for checkpointing the simulation.
For computations that run on thousands of individual MPI processes, writing
all of this information to disk can be a bottleneck.

Most large clusters have extensive documentation on how to tune file storage
for optimal performance, and if you are doing large computations, it is worth
reading through this documentation. At least for Lustre file systems, you can
also employ the {ref}`parameters:Output_20directory_20LFS_20stripe_20count`.
parameter to set LFS stripe counts -- see the documentation of that parameter
for more information.

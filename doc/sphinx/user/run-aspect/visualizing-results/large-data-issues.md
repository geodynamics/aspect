# Large data issues for parallel computations

#### Large data issues for parallel computations

Among the challenges in visualizing the results of parallel computations is
dealing with the large amount of data. The first bottleneck this presents is
during run-time when ASPECT wants to write the
visualization data of a time step to disk. Using the compressed VTU format,
ASPECT generates on the order of 10 bytes of
output for each degree of freedom in 2d and more in 3d; thus, output of a
single time step can run into the range of gigabytes that somehow have to get
from compute nodes to disk. This stresses both the cluster interconnect as
well as the data storage array.

There are essentially two strategies supported by 
ASPECT for this scenario:

-   If your cluster has a fast interconnect, for example Infiniband, and if
    your cluster has a fast, distributed file system, then <span
    ASPECT can produce output files that are already
    located in the correct output directory (see the options in
    {ref}`parameters:global`17]) on the global file system. <span
    ASPECT uses MPI I/O calls to this end, ensuring
    that the local machines do not have to access these files using slow
    NFS-mounted global file systems.

-   If your cluster has a slow interconnect, e.g., if it is simply a
    collection of machines connected via Ethernet, then writing data to a
    central file server may block the rest of the program for a while. On the
    other hand, if your machines have fast local storage for temporary file
    systems, then ASPECT can write data first
    into such a file and then move it in the background to its final
    destination while already continuing computations. To select this mode,
    set the appropriate variables discussed in
    {ref}`parameters:Postprocess/Visualization`52]. Note,
    however, that this scheme only makes sense if every machine on which MPI
    processes run has fast local disk space for temporary storage.

<div class="center">

</div>

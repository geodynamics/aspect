# Checkpoint/restart support

### Checkpoint/restart support

If you do long runs, especially when using parallel computations, there are a
number of reasons to periodically save the state of the program:

-   If the program crashes for whatever reason, the entire computation may be
    lost. A typical reason is that a program has exceeded the requested
    wallclock time allocated by a batch scheduler on a cluster.

-   Most of the time, no realistic initial conditions for strongly convecting
    flow are available. Consequently, one typically starts with a somewhat
    artificial state and simply waits for a long while till the convective
    state enters the phase where it shows its long-term behavior. However,
    getting there may take a good amount of CPU time and it would be silly to
    always start from scratch for each different parameter setting. Rather,
    one would like to start such parameter studies with a saved state that has
    already passed this initial, unphysical, transient stage.

To this end, ASPECT creates a set of files in
the output directory (selected in the parameter file) every N time steps
(controlled by the number of steps or wall time as specified in
`subsection Checkpointing`, see
{ref}`parameters:Checkpointing`61]) in which the entire state of
the program is saved so that a simulation can later be continued at this
point. The previous checkpoint files will then be deleted. To resume
operations from the last saved state, you need to set the `Resume computation`
flag in the input parameter file to `true`, see
{ref}`parameters:Resume computation`62].

<div class="center">

</div>

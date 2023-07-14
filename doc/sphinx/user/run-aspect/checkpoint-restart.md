(sec:checkpoint-restart)=
# Checkpoint/restart support

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
{ref}`parameters:Checkpointing`) in which the entire state of
the program is saved so that a simulation can later be continued at this
point. The previous checkpoint files will then be deleted. To resume
operations from the last saved state, you need to set the `Resume computation`
flag in the input parameter file to `true`, see
{ref}`parameters:global`.

:::{note}
It is not imperative that the parameters selected in the input file are exactly the same
when resuming a program from a saved state than what they were at the time when this state
was saved. For example, one may want to choose a different parameterization of the material law,
or add or remove postprocessors that should be run at the end of each time step. Likewise, the
end time, the times at which some additional mesh refinement steps should happen, etc., can be
different.
Yet, it is clear that some other things can’t be changed: For example, the geometry model that
was used to generate the coarse mesh and describe the boundary must be the same before and
after resuming a computation. Likewise, you can not currently restart a computation with a
different number of processors than initially used to checkpoint the simulation. Not all invalid
combinations are easy to detect, and ASPECT may not always realize immediate what is going on
if you change a setting that can’t be changed. However, you will almost invariably get nonsensical
results after some time.
:::

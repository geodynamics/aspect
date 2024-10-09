# Criteria for terminating a simulation

ASPECT allows for different ways of terminating
a simulation. For example, the simulation may have reached a final time
specified in the input file. However, it also allows for ways to terminate a
simulation when it has reached a steady state (or, rather, some criterion
determines that it is close enough to steady state), or by an external action
such as placing a specially named file in the output directory. The criteria
determining termination of a simulation are all implemented in plugins. The
parameters describing these criteria are listed in {ref}`parameters:Termination_20criteria`.

To implement a termination criterion, you need to overload the
`aspect::TerminationCriteria::Interface` class and use the
`ASPECT_REGISTER_TERMINATION_CRITERION` macro to register your new class. The
implementation of the new class should be in namespace
`aspect::TerminationCriteria`.

The principal function that needs to be overloaded returns a value that indicates whether the
simulation should be terminated. Typical examples can be found in the existing
implementations in the `source/termination_criteria/` directory. As usual, your
termination criterion implementation will likely need to be derived from the
`SimulatorAccess` to get access to the current state of the simulation.

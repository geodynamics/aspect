# Criteria for terminating a simulation

ASPECT allows for different ways of terminating
a simulation. For example, the simulation may have reached a final time
specified in the input file. However, it also allows for ways to terminate a
simulation when it has reached a steady state (or, rather, some criterion
determines that it is close enough to steady state), or by an external action
such as placing a specially named file in the output directory. The criteria
determining termination of a simulation are all implemented in plugins. The
parameters describing these criteria are listed in {ref}`sec:3.195][].

To implement a termination criterion, you need to overload the
[aspect::TerminationCriteria::Interface][] class and use the
`ASPECT_REGISTER_TERMINATION_CRITERION` macro to register your new class. The
implementation of the new class should be in namespace
`aspect::TerminationCriteria`.

Specifically, your new class needs to implement the following basic interface:

``` c++
template <int dim>
    class aspect::TerminationCriteria::Interface
    {
      public:
        virtual
        bool
        execute () const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
```

The first of these functions returns a value that indicates whether the
simulation should be terminated. Typical examples can be found in the existing
implementations in the `source/termination_criteria` directory. As usual, your
termination criterion implementation will likely need to be derived from the
`SimulatorAccess` to get access to the current state of the simulation.

The remaining functions are obvious, and are also discussed in the
documentation of this interface class at
[aspect::TerminationCriteria::Interface][]. The purpose of the last two
functions has been discussed in the general overview of plugins above.

(sec:extending-signals)=
# Extending ASPECT through signals

Not all things you may want to do fit neatly into the list of plugins of the
previous sections. Rather, there are cases where you may want to change things
that are more of the one-off kind and that require code that is at a lower
level and requires more knowledge about
ASPECT's internal workings. For such changes,
we still want to stick with the general principle outlined at the beginning of
{ref}`cha:extending`: You should be able to make all of your changes and
extensions in your own files, without having to modify
ASPECT's own sources.

To support this, ASPECT uses a
"signals" mechanism. Signals are, in essence, objects that
represent *events*, for example the fact that the solver has finished a time
step. The core of ASPECT defines a number of
such signals, and *triggers* them at the appropriate points. The idea of
signals is now that you can *connect* to them: you can tell the signal that it
should call a particular function every time the signal is triggered. The
functions that are connected to a signal are called "slots" in
common diction. One, several, or no slots may be connected to each signal.

There are two kinds of signals that ASPECT
provides:

-   Signals that are triggered at startup of the program: These are, in
    essence, signals that live in some kind of global scope. Examples are
    signals that declare additional parameters for use in input files, or that
    read the values of these parameters from a `ParameterHandler` object.
    These signals are static member variables of the structure that contains
    them and consequently exist only once for the entire program.

-   Signals that reference specific events that happen inside a simulator
    object. These are regular member variables of the structure that contains
    them, and because each simulator object has such a structure, the signals
    exist once per simulator object. (Which in practice is only once per
    program, of course.)

For both of these kinds, a user-written plugin file can (but does not need) to
register functions that connect functions in this file (i.e., slots) to their
respective signals.

In the first case, code that registers slots with global signals would look
like this:

```{code-block} c++
// A function that will be called at the time when parameters are declared.
// It receives the dimension in which ASPECT will be run as the first argument,
// and the ParameterHandler object that holds the runtime parameter
// declarations as second argument.
void declare_parameters(const unsigned int dim,
                        ParameterHandler &prm)
{
  prm.declare_entry("My parameter", ...);
}


// The same for parsing parameters. 'my_parameter' is a parameter
// that stores something we want to read from the input file
// and use in other functions in this file (which we don't show here).
// For simplicity, we assume that it is an integer.
//
// The function also receives a first argument that contains all
// of the other (already parsed) arguments of the simulation, in
// case what you want to do here wants to refer to other parameters.
int my_parameter;

template <int dim>
void parse_parameters(const Parameters<dim> &parameters,
                      ParameterHandler &prm)
{
  my_parameter = prm.get_integer ("My parameter");
}


// Now have a function that connects slots (i.e., the two functions
// above) to the static signals. Do this for both the 2d and 3d
// case for generality.
void parameter_connector ()
{
  SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
  SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

  SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
  SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
}


// Finally register the connector function above to make sure it gets run
// whenever we load a user plugin that is mentioned among the additional
// shared libraries in the input file:
ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
```

The second kind of signal can be connected to once a simulator object has been
created. As above, one needs to define the slots, define a connector function,
and register the connector function. The following gives an example:

```{code-block} c++
// A function that is called at the end of creating the current constraints
// on degrees of freedom (i.e., the constraints that describe, for example,
// hanging nodes, boundary conditions, etc).
template <int dim>
void post_constraints_creation (const SimulatorAccess<dim> &simulator_access,
                                ConstraintMatrix &current_constraints)
{
  ...; // do whatever you want to do here
}


// A function that is called from the simulator object and that can connect
// a slot (such as the function above) to any of the signals declared in the
// structure passed as argument:
template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.post_constraints_creation.connect (&post_constraints_creation<dim>);
}


// Finally register the connector function so that it is called whenever
// a simulator object has been set up. For technical reasons, we need to
// register both 2d and 3d versions of this function:
ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
```

As mentioned above, each signal may be connected to zero, one, or many slots.
Consequently, you could have multiple plugins each of which connect to the
same slot, or the connector function above may just connect multiple slots
(i.e., functions in your program) to the same signal.

So what could one do in a place like this? One option would be to just monitor
what is going on, e.g., in code like this that simply outputs into the
statistics file (see {ref}`sec:run-aspect:visualize:stat-data`):

```{code-block} c++
template <int dim>
void post_constraints_creation (const SimulatorAccess<dim> &simulator_access,
                                ConstraintMatrix &current_constraints)
{
  simulator_access.get_statistics_object()
    .add_value ("number of constraints",
                current_constraints.n_constraints());
}
```

This will produce, for every time step (because this is how often the signal
is called) an entry in a new column in the statistics file that records the
number of constraints. On the other hand, it is equally possible to also
modify the constraints object at this point. An application would be if you
wanted to run a simulation where you prescribe the velocity in a part of the
domain, e.g., for a subducting slab (see {ref}`sec:cookbooks:prescribed_velocity`).

Signals exist for various waypoints in a simulation and you can consequently
monitor and change what is happening inside a simulation by connecting your
own functions to these signals. It would be pointless to list here what
signals actually exist &ndash; simply refer to the documentation of the
[SimulatorSignals](https://aspect.geodynamics.org/doc/doxygen/classes.html)
class for a complete list of signals you can connect to.

As a final note, it is generally true that writing functions that can connect
to signals require significantly more internal knowledge of the workings of
ASPECT than writing plugins through the
mechanisms outlined above. It also allows to affect the course of a simulation
by working on the internal data structures of
ASPECT in ways that are not available to the largely
passive and reactive plugins discussed in previous sections. With this
obviously also comes the potential for trouble. On the other hand, it also
allows to do things with ASPECT that were not
initially intended by the authors, and that would be hard or impossible to
implement through plugins. An example would be to couple different codes by
exchanging details of the internal data structures, or even update the
solution vectors using information received from another code.

:::{note}
Chances are that if you think about using the signal mechanism, there is not yet a signal
that is triggered at exactly the point where you need it. Consequently, you will be tempted to
just put your code into the place where it fits inside ASPECT where it fits best. This is poor
practice: it prevents you from upgrading to a newer version of ASPECT at a later time because
this would overwrite the code you inserted.
Rather, a more productive approach would be to either define a new signal that is triggered
where you need it, and connect a function (slot) in your own plugin file to this signal using
the mechanisms outlined above. Then send the code that defines and triggers the signal to the
developers of ASPECT to make sure that it is also included in the next release. Alternatively,
you can also simply ask on the forum for someone to add such a signal in the place where you
want it. Either way, adding signals is something that is easy to do, and we would much rather
add signals than have people who modify the ASPECT source files for their own needs and are
then stuck on a particular version.
:::

# Postprocessors: Evaluating the solution after each time step

Postprocessors are arguably the most complex and powerful of the plugins
available in ASPECT since they do not only
passively provide any information but can actually compute quantities derived
from the solution. They are executed once at the end of each time step and,
unlike all the other plugins discussed above, there can be an arbitrary number
of active postprocessors in the same program (for the plugins discussed in
previous sections it was clear that there is always exactly one material
model, geometry model, etc.).

##### Motivation.

The original motivation for postprocessors is that the goal of a simulation is
of course not the simulation itself, but that we want to do something with the
solution. Examples for already existing postprocessors are:

-   Generating output in file formats that are understood by visualization
    programs. This is facilitated by the
    [aspect::Postprocess::Visualization][] class and a separate class of
    visualization postprocessors, see {ref}`sec:1.4.9][].

-   Computing statistics about the velocity field (e.g., computing minimal,
    maximal, and average velocities), temperature field (minimal, maximal, and
    average temperatures), or about the heat fluxes across boundaries of the
    domain. This is provided by the
    [aspect::Postprocess::VelocityStatistics][],
    [aspect::Postprocess::TemperatureStatistics][],
    [aspect::Postprocess::HeatFluxStatistics][] classes, respectively.

Since writing this text, there may have been other additions as well.

However, postprocessors can be more powerful than this. For example, while the
ones listed above are by and large stateless, i.e., they do not carry
information from one invocation at one timestep to the next invocation,[6]
there is nothing that prohibits postprocessors from doing so. For example, the
following ideas would fit nicely into the postprocessor framework:

-   *Passive particles:* If one would like to follow the trajectory of
    material as it is advected along with the flow field, one technique is to
    use particles. To implement this, one would start with an initial
    population of particles distributed in a certain way, for example close to
    the core-mantle boundary. At the end of each time step, one would then
    need to move them forward with the flow field by one time increment. As
    long as these particles do not affect the flow field (i.e., they do not
    carry any information that feeds into material properties; in other words,
    they are *passive*), their location could well be stored in a
    postprocessor object and then be output in periodic intervals for
    visualization. In fact, such a passive particle postprocessor is already
    available.

-   *Surface or crustal processes:* Another possibility would be to keep track
    of surface or crustal processes induced by mantle flow. An example would
    be to keep track of the thermal history of a piece of crust by updating it
    every time step with the heat flux from the mantle below. One could also
    imagine integrating changes in the surface topography by considering the
    surface divergence of the surface velocity computed in the previous time
    step: if the surface divergence is positive, the topography is lowered,
    eventually forming a trench; if the divergence is negative, a mountain
    belt eventually forms.

In all of these cases, the essential limitation is that postprocessors are
*passive*, i.e., that they do not affect the simulation but only observe it.

##### The statistics file.

Postprocessors fall into two categories: ones that produce lots of output
every time they run (e.g., the visualization postprocessor), and ones that
only produce one, two, or in any case a small and fixed number of often
numerical results (e.g., the postprocessors computing velocity, temperature,
or heat flux statistics). While the former are on their own in implementing
how they want to store their data to disk, there is a mechanism in place that
allows the latter class of postprocessors to store their data into a central
file that is updated at the end of each time step, after all postprocessors
are run.

To this end, the function that executes each of the postprocessors is given a
reference to a `dealii::TableHandler` object that allows to store data in
named columns, with one row for each time step. This table is then stored in
the `statistics` file in the directory designated for output in the input
parameter file. It allows for easy visualization of trends over all time
steps. To see how to put data into this statistics object, take a look at the
existing postprocessor objects.

Note that the data deposited into the statistics object need not be numeric in
type, though it often is. An example of text-based entries in this table is
the visualization class that stores the name of the graphical output file
written in a particular time step.

##### Implementing a postprocessor.

Ultimately, implementing a new postprocessor is no different than any of the
other plugins. Specifically, you&rsquo;ll have to write a class that overloads
the [aspect::Postprocess::Interface][] base class and use the
`ASPECT_REGISTER_POSTPROCESSOR` macro to register your new class. The
implementation of the new class should be in namespace `aspect::Postprocess`.

In reality, however, implementing new postprocessors is often more difficult.
Primarily, this difficulty results from two facts:

-   Postprocessors are not self-contained (only providing information) but in
    fact need to access the solution of the model at each time step. That is,
    of course, the purpose of postprocessors, but it requires that the writer
    of a plugin has a certain amount of knowledge of how the solution is
    computed by the main `Simulator` class, and how it is represented in data
    structures. To alleviate this somewhat, and to insulate the two worlds
    from each other, postprocessors do not directly access the data structures
    of the simulator class. Rather, in addition to deriving from the
    [aspect::Postprocess::Interface][] base class, postprocessors also derive
    from the [SimulatorAccess][aspect::SimulatorAccess class] class that has a
    number of member functions postprocessors can call to obtain read-only
    access to some of the information stored in the main class of <span
    ASPECT. See [the documentation of this
    class][aspect::SimulatorAccess class] to see what kind of information is
    available to postprocessors. See also {ref}`sec:1.1][] for more
    information about the `SimulatorAccess` class.

-   Writing a new postprocessor typically requires a fair amount of knowledge
    how to leverage the DEAL.II library to
    extract information from the solution. The existing postprocessors are
    certainly good examples to start from in trying to understand how to do
    this.

Given these comments, the interface a postprocessor class has to implement is
rather basic:

``` c++
template <int dim>
    class aspect::Postprocess::Interface
    {
      public:
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) = 0;

        virtual
        void
        save (std::map<std::string, std::string> &status_strings) const;

        virtual
        void
        load (const std::map<std::string, std::string> &status_strings);

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
```

The purpose of these functions is described in detail in the documentation of
the [aspect::Postprocess::Interface][] class. While the first one is
responsible for evaluating the solution at the end of a time step, the
`save/load` functions are used in checkpointing the program and restarting it
at a previously saved point during the simulation. The first of these
functions therefore needs to store the status of the object as a string under
a unique key in the database described by the argument, while the latter
function restores the same state as before by looking up the status string
under the same key. The default implementation of these functions is to do
nothing; postprocessors that do have non-static member variables that contain
a state need to overload these functions.

There are numerous postprocessors already implemented. If you want to
implement a new one, it would be helpful to look at the existing ones to see
how they implement their functionality.

##### Postprocessors and checkpoint/restart.

Postprocessors have `save()` and `load()` functions that are used to write the
data a postprocessor has into a checkpoint file, and to load it again upon
restart. This is important since many postprocessors store some state &ndash;
say, a temporal average over all the time steps seen so far, or the number of
the last graphical output file generated so that we know how the next one
needs to be numbered.

The typical case is that this state is the same across all processors of a
parallel computation. Consequently, what ASPECT
writes into the checkpoint file is only the state obtained from the
postprocessors on processor 0 of a parallel computation. On restart, all
processors read from the same file and the postprocessors on *all* processors
will be initialized by what the same postprocessor on processor 0 wrote.

There are situations where postprocessors do in fact store complementary
information on different processors. At the time of writing this, one example
is the postprocessor that supports advecting passive particles along the
velocity field: on every processor, it handles only those particles that lie
inside the part of the domain that is owned by this MPI rank. The
serialization approach outlined above can not work in this case, for obvious
reasons. In cases like this, one needs to implement the `save()` and `load()`
differently than usual: one needs to put all variables that are common across
processors into the maps of string as usual, but one then also needs to save
all state that is different across processors, from all processors. There are
two ways: If the amount of data is small, you can use MPI communications to
send the state of all processors to processor zero, and have processor zero
store it in the result so that it gets written into the checkpoint file; in
the `load()` function, you will then have to identify which part of the text
written by processor 0 is relevant to the current processor. Or, if your
postprocessor stores a large amount of data, you may want to open a restart
file specifically for this postprocessor, use MPI I/O or other ways to write
into it, and do the reverse operation in `load()`.

Note that this approach requires that ASPECT
actually calls the `save()` function on all processors. This in fact happens
&ndash; though ASPECT also discards the result
on all but processor zero.

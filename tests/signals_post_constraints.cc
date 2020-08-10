#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

using namespace aspect;


template <int dim>
void post_constraints_creation (const SimulatorAccess<dim> &simulator_access,
                                AffineConstraints<double> &current_constraints)
{
  simulator_access.get_statistics_object()
  .add_value ("number of constraints",
              current_constraints.n_constraints());
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.post_constraints_creation.connect (&post_constraints_creation<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)

#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

using namespace aspect;


template <int dim>
void post_stokes_solver (const SimulatorAccess<dim> &sim,
                         const bool success,
                         const std::vector<double> &history)
{
  std::cout << "\npost_stokes_solver:\n";
  for (unsigned int i=0; i<history.size(); ++i)
    std::cout << history[i] << std::endl;
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  std::cout << "Connecting signals" << std::endl;
  signals.post_stokes_solver.connect (&post_stokes_solver<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)

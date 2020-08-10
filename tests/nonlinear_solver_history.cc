#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

namespace aspect
{
  template <int dim>
  void post_nonlinear_solver (const SolverControl &nonlinear_solver_control)
  {
    const bool success = nonlinear_solver_control.last_check() == SolverControl::success;
    std::cout << "\nnumber of nonlinear iterations: " << nonlinear_solver_control.last_step() << ". State: " << success << ".\n";
  }


  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    std::cout << "Connecting signals" << std::endl;
    signals.post_nonlinear_solver.connect (&post_nonlinear_solver<dim>);
  }


  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                    signal_connector<3>)
}

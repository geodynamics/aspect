#include <deal.II/base/parameter_handler.h>
#include <aspect/global.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  using namespace dealii;

  // Global variables (to be set by parameters)
  int switch_step;
  bool switched;

  /**
   * Declare additional parameters.
   */
  void declare_parameters(const unsigned int,
                          ParameterHandler &prm)
  {
    prm.declare_entry("Switch step", "0",
                       Patterns::Integer(0),
                       "Switch the bottom boundary condition from fixed temperature to no-flux"
                       "at the timestep given.");
  }

  template <int dim>
  void parse_parameters(const Parameters<dim>,
                        ParameterHandler &prm)
  {
    switch_step = prm.get_integer("Switch step");
    switched = false;
  }

  template <int dim>
  void change_boundary_condition (const SimulatorAccess<dim> &simulator_access,
                                  Parameters<dim> &parameters)
  {
    types::boundary_id bottom = 2; //bottom boundary indicator for 2D box

    if ( simulator_access.get_timestep_number() >= switch_step && !switched )
      {
        simulator_access.get_pcout()<<"Switching boundary condition!"<<std::endl;
        parameters.fixed_temperature_boundary_indicators.erase(bottom);
        switched = true;
      }
  }

  // Connect declare_parameters and parse_parameters to appropriate signals.
  void parameter_connector ()
  {
    SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
    SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

    SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
    SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
  }

  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.edit_parameters_pre_setup_dofs.connect (&change_boundary_condition<dim>);
  }

  // Tell Aspect to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}

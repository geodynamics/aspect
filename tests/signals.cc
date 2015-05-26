// make sure we can include deal.II and aspect files
#include <aspect/simulator_signals.h>

#include <iostream>

using namespace aspect;


// create a function that is run upon loading the plugin
// when declaring parameters, and that produces some output
void declare_parameters(const unsigned int dim,
                        ParameterHandler &prm)
{
  std::cout << "declaring parameters" << std::endl;
  prm.declare_entry("abc", "42", Patterns::Integer(21,43));
}


// same for parsing parameters
template <int dim>
void parse_parameters(const Parameters<dim> &parameters,
                      ParameterHandler &prm)
{
  std::cout << "parsing parameters: abc="
            << prm.get("abc") << std::endl;
  Assert (prm.get_integer ("abc") == 21, ExcInternalError());
}


void parameter_connector ()
{
  SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
  SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

  SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
  SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
}


ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)

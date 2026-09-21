// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/simulator_signals.h>

#include <list>


namespace aspect
{
  // provide locations for the static variables of this class
  template <int dim>
  boost::signals2::signal<void (const unsigned int,
                                ParameterHandler &)>  SimulatorSignals<dim>::declare_additional_parameters;

  template <int dim>
  boost::signals2::signal<void (const Parameters<dim> &,
                                ParameterHandler &)>  SimulatorSignals<dim>::parse_additional_parameters;


  namespace internals
  {
    namespace SimulatorSignals
    {
      std::list<std::function<void (aspect::SimulatorSignals<2> &)>> connector_functions_2d;
      std::list<std::function<void (aspect::SimulatorSignals<3> &)>> connector_functions_3d;

      static bool connector_functions_have_been_called = false;

      // add a user-provided connector to the list of connectors we keep
      void register_connector_function_2d (const std::function<void (aspect::SimulatorSignals<2> &)> &connector)
      {
        Assert(!connector_functions_have_been_called,
               ExcMessage("Registration of signal connector happened after connection has already been called!"));
        connector_functions_2d.push_back (connector);
      }

      void register_connector_function_3d (const std::function<void (aspect::SimulatorSignals<3> &)> &connector)
      {
        Assert(!connector_functions_have_been_called,
               ExcMessage("Registration of signal connector happened after connection has already been called!"));
        connector_functions_3d.push_back (connector);
      }


      // call connectors to ensure that plugins get a change to register their slots
      template <>
      void call_connector_functions (aspect::SimulatorSignals<2> &signals)
      {
        Assert(!connector_functions_have_been_called,
               ExcInternalError());

        for (const auto &p : connector_functions_2d)
          p(signals);

        connector_functions_have_been_called = true;
      }

      template <>
      void call_connector_functions (aspect::SimulatorSignals<3> &signals)
      {
        Assert(!connector_functions_have_been_called,
               ExcInternalError());

        for (const auto &p : connector_functions_3d)
          p(signals);

        connector_functions_have_been_called = true;
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct SimulatorSignals<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}

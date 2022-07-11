/*
  Copyright (C) 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/interface.h>
#include <aspect/simulator.h>
#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

using namespace aspect;


template <int dim>
void post_constraints_creation (const SimulatorAccess<dim> &sim,
                                AffineConstraints<double> &)
{
  static bool first = true;
  if (first)
    {
      std::stringstream ss;
      sim.get_simulator().write_plugin_graph(ss);
      std::string content = ss.str();

      const char *checks[] =
      {
        "digraph Plugins",
        "label=\"TestMaterial\"",
        "N6aspect12TestMaterialILi2EEE -> N6aspect13MaterialModel9InterfaceILi2EEE",
        "label=\"simpler\"",
        "Simulator -> SimulatorAccess",
        "N6aspect19BoundaryComposition9AsciiDataILi2EEE [label=\"ascii data\""
      };

      for (unsigned int i = 0; i < sizeof(checks)/sizeof(checks[0]); ++i)
        {
          bool found = content.find(checks[i])!=std::string::npos;

          if (found)
            std::cout << "found " << checks[i] << std::endl;
          Assert(found, ExcInternalError());
        }

    }
  first = false;
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.post_constraints_creation.connect (&post_constraints_creation<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)


// make a dummy material model and instantiate and register it
namespace aspect
{
  template <int dim>
  class TestMaterial : public MaterialModel::Interface<dim>
  {
    public:
      virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                            MaterialModel::MaterialModelOutputs<dim> &out) const
      {}

      virtual bool is_compressible () const
      {
        return false;
      }
  };


  ASPECT_REGISTER_MATERIAL_MODEL(TestMaterial,
                                 "TestMaterial",
                                 "A material model.")
}

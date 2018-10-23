/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include "lithosphere.h"
#include <aspect/initial_temperature/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    double
    Lithosphere<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double temperature = this->get_initial_temperature_manager().initial_temperature(position);

      // Crustal composition
      if (depth < moho && n_comp == 0)
        return 1.;
      else if (depth >= moho && temperature < LAB_isotherm && n_comp == 1)
        return 1.;
      else
        return 0.;
    }


    template <int dim>
    void
    Lithosphere<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Lithosphere");
        {
          prm.declare_entry ("Moho", "30000.0",
                             Patterns::Double (0),
                             "Moho depth. Units: $m.");
          prm.declare_entry ("LAB isotherm", "1673",
                             Patterns::Double (0),
                             "Temperature at the base of the lithosphere. Units: $K.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Lithosphere<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Lithosphere");
        {
          moho                            = prm.get_double ("Moho");
          LAB_isotherm                    = prm.get_double ("LAB isotherm");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Lithosphere,
                                              "lithosphere",
                                              "Specify the composition in terms of an explicit formula. The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}

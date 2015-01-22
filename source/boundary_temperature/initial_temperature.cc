/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/boundary_temperature/initial_temperature.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ InitialTemperature -------------------

    template <int dim>
    double
    InitialTemperature<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location) const
    {
      return this->get_initial_conditions().initial_temperature(location);
    }


    template <int dim>
    double
    InitialTemperature<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return min_temperature;
    }



    template <int dim>
    double
    InitialTemperature<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return max_temperature;
    }



    template <int dim>
    void
    InitialTemperature<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Initial temperature");
        {
          prm.declare_entry ("Minimal temperature", "0",
                             Patterns::Double (),
                             "Minimal temperature. Units: K.");
          prm.declare_entry ("Maximal temperature", "3773",
                             Patterns::Double (),
                             "Maximal temperature. Units: K.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    InitialTemperature<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Initial temperature");
        {
          min_temperature = prm.get_double ("Minimal temperature");
          max_temperature = prm.get_double ("Maximal temperature");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(InitialTemperature,
                                               "initial temperature",
                                               "A model in which the temperature at the boundary "
                                               "is chosen to be the same as given in the initial "
                                               "conditions."
                                               "\n\n"
                                               "Because this class simply takes what the initial "
                                               "temperature had described, this class can not "
                                               "know certain pieces of information such as the "
                                               "minimal and maximal temperature on the boundary. "
                                               "For operations that require this, for example in "
                                               "postprocessing, this boundary temperature model "
                                               "must therefore be told what the minimal and "
                                               "maximal values on the boundary are. This is done "
                                               "using parameters set in section ``Boundary temperature model/Initial temperature''.")
  }
}

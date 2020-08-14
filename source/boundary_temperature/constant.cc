/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/boundary_temperature/constant.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/signaling_nan.h>

#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ Constant -------------------

    template <int dim>
    double
    Constant<dim>::
    boundary_temperature (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/) const
    {
      const std::map<types::boundary_id, double>::const_iterator it = boundary_temperatures.find(boundary_indicator);
      if (it != boundary_temperatures.end())
        return it->second;
      else
        {
          Assert (false,
                  ExcMessage ("Unknown boundary indicator with number <" + Utilities::int_to_string(boundary_indicator) + ">. "
                              "You may not have specified the temperature for this boundary indicator "
                              "in the input file."));
          return numbers::signaling_nan<double>();
        }
    }


    template <int dim>
    double
    Constant<dim>::
    minimal_temperature (const std::set<types::boundary_id> &) const
    {
      std::map<types::boundary_id, double>::const_iterator it = boundary_temperatures.begin();
      double min = it->second;
      ++it;

      for ( ; it != boundary_temperatures.end(); ++it)
        if ( it->second < min )
          min = it->second;

      return min;
    }



    template <int dim>
    double
    Constant<dim>::
    maximal_temperature (const std::set<types::boundary_id> &) const
    {
      std::map<types::boundary_id, double>::const_iterator it = boundary_temperatures.begin();
      double max = it->second;
      ++it;

      for ( ; it != boundary_temperatures.end(); ++it)
        if ( it->second > max )
          max = it->second;

      return max;
    }



    template <int dim>
    void
    Constant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Constant");
        {
          prm.declare_entry ("Boundary indicator to temperature mappings", "",
                             Patterns::Map (Patterns::Anything(),
                                            Patterns::Double()),
                             "A comma separated list of mappings between boundary "
                             "indicators and the temperature associated with the "
                             "boundary indicators. The format for this list is "
                             "``indicator1 : value1, indicator2 : value2, ...'', "
                             "where each indicator is a valid boundary indicator "
                             "(either a number or the symbolic name of a boundary as provided "
                             "by the geometry model) "
                             "and each value is the temperature of that boundary." );
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Constant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Constant");
        {
          // get the list of mappings
          const std::vector<std::string> x_boundary_temperatures
            = Utilities::split_string_list(prm.get ("Boundary indicator to temperature mappings"));


          for (const auto &boundary_id_string : x_boundary_temperatures)
            {
              // each entry has the format (white space is optional):
              // <id> : <value (might have spaces)>
              const std::vector<std::string> parts = Utilities::split_string_list (boundary_id_string, ':');

              AssertThrow (parts.size() == 2,
                           ExcMessage (std::string("Invalid entry trying to describe boundary "
                                                   "temperatures. Each entry needs to have the form "
                                                   "<boundary_id : name>, "
                                                   "but there is an entry of the form <") + boundary_id_string + ">"));

              types::boundary_id boundary_id;
              try
                {
                  boundary_id
                    = this->get_geometry_model().translate_symbolic_boundary_name_to_id (parts[0]);
                }
              catch (const std::string &error)
                {
                  AssertThrow (false, ExcMessage ("While parsing the entry <Boundary temperature model/Constant>, "
                                                  "there was an error. Specifically, "
                                                  "the conversion function complained as follows: "
                                                  + error));
                }

              AssertThrow (boundary_temperatures.find(boundary_id) == boundary_temperatures.end(),
                           ExcMessage ("Boundary indicator <" + Utilities::int_to_string(boundary_id) +
                                       "> appears more than once in the list of indicators "
                                       "for constant temperature boundary conditions."));

              boundary_temperatures[boundary_id] = Utilities::string_to_double (parts[1]);
            }
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
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Constant,
                                               "constant",
                                               "A model in which the temperature is chosen constant on "
                                               "a given boundary indicator.  Parameters are read from the "
                                               "subsection 'Constant'.")
  }
}

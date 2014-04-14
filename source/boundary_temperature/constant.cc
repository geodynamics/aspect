/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#include <aspect/boundary_temperature/constant.h>
#include <deal.II/base/utilities.h>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ Constant -------------------

    template <int dim>
    double
    Constant<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location) const
    {
      types::boundary_id bid = boundary_indicator;
      std::map< types::boundary_id, double>::const_iterator it = boundary_temperatures.find(bid);
      if (it != boundary_temperatures.end())
        return it->second;
      else
        {
          Assert (false, ExcMessage ("Unknown boundary indicator."));
          return std::numeric_limits<double>::quiet_NaN();
        }
    }


    template <int dim>
    double
    Constant<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
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
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
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
                             Patterns::Map (Patterns::Integer(0, std::numeric_limits<types::boundary_id>::max()),
                                            Patterns::Double()),
                             "A comma separated list of mappings between boundary "
                             "indicators and the temperature associated with the "
                             "boundary indicators. The format for this list is "
                             "``indicator1 : value1, indicator2 : value2, ...'', "
                             "where each indicator is a valid boundary indicator "
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
          //get the list of mappings
          const std::vector<std::string> x_boundary_temperatures
            = Utilities::split_string_list(prm.get ("Boundary indicator to temperature mappings"));


          for (std::vector<std::string>::const_iterator it = x_boundary_temperatures.begin();
               it != x_boundary_temperatures.end(); ++it)
            {
              // each entry has the format (white space is optional):
              // <id> [x][y][z] : <value (might have spaces)>
              std::string comp = "";
              std::string value = "";

              std::stringstream ss(*it);
              int b_id;
              ss >> b_id; // need to read as int, not char
              types::boundary_id boundary_id = b_id;

              char c;
              while (ss.peek()==' ') ss.get(c); // eat spaces

              if (ss.peek() != ':')
                {
                  Assert(false, ExcMessage("Cannot parse boundary temperature list, format"
                                           "``boundary_id : value'' appears to be missing"));
                }
              else
                ss.get(c); // read the ':'

              while (ss.peek()==' ') ss.get(c); // eat spaces
              std::getline(ss,value); // read until the end of the string

              AssertThrow (boundary_temperatures.find(boundary_id) == boundary_temperatures.end(),
                           ExcMessage ("Boundary indicator <" + Utilities::int_to_string(boundary_id) +
                                       "> appears more than once in the list of indicators "
                                       "for constant temperature boundary conditions."));

              boundary_temperatures[boundary_id] = Utilities::string_to_double(value);
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

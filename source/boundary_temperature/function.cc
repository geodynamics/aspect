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


#include <aspect/boundary_temperature/function.h>
#include <aspect/global.h>

namespace aspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    Function<dim>::Function ()
      :
      boundary_temperature_function (1)
    {}



    template <int dim>
    double
    Function<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &position) const
    {
      return boundary_temperature_function.value(position);
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_temperature_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_temperature_function.set_time (this->get_time());
    }



    template <int dim>
    double
    Function<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return min_temperature;
    }



    template <int dim>
    double
    Function<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return max_temperature;
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);

          prm.declare_entry ("Minimal temperature", "273",
                             Patterns::Double (),
                             "Minimal temperature. Units: K.");
          prm.declare_entry ("Maximal temperature", "3773",
                             Patterns::Double (),
                             "Maximal temperature. Units: K.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Function");
        try
          {
            boundary_temperature_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary temperature model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
        min_temperature = prm.get_double ("Minimal temperature");
        max_temperature = prm.get_double ("Maximal temperature");
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Function,
                                               "function",
                                               "Implementation of a model in which the boundary "
                                               "temperature is given in terms of an explicit formula "
                                               "that is elaborated in the parameters in section "
                                               "``Boundary temperature model|Function''. "
                                               "\n\n"
                                               "Since the symbol $t$ indicating time "
                                               "may appear in the formulas for the prescribed "
                                               "temperatures, it is interpreted as having units "
                                               "seconds unless the global input parameter ``Use "
                                               "years in output instead of seconds'' is set, in "
                                               "which case we interpret the formula expressions "
                                               "as having units year."
                                               "\n\n"
                                               "Because this class simply takes what the "
                                               "function calculates, this class can not "
                                               "know certain pieces of information such as the "
                                               "minimal and maximal temperature on the boundary. "
                                               "For operations that require this, for example in "
                                               "postprocessing, this boundary temperature model "
                                               "must therefore be told what the minimal and "
                                               "maximal values on the boundary are. This is done "
                                               "using parameters set in section ``Boundary temperature model/Initial temperature''."
                                               "\n\n"
                                               "The format of these "
                                               "functions follows the syntax understood by the "
                                               "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}

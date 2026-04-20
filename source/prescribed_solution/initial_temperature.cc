/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE. If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/initial_temperature/interface.h>
#include <aspect/prescribed_solution/initial_temperature.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace PrescribedSolution
  {
    template <int dim>
    InitialTemperature<dim>::InitialTemperature ()
      :
      prescribed_temperature_indicator_function(1)
    {}



    template <int dim>
    void
    InitialTemperature<dim>::initialize ()
    {
      initial_temperature_manager =
        this->get_initial_temperature_manager_pointer();
    }



    template <int dim>
    void
    InitialTemperature<dim>::update ()
    {
      if (this->convert_output_to_years())
        prescribed_temperature_indicator_function.set_time(this->get_time() / year_in_seconds);
      else
        prescribed_temperature_indicator_function.set_time(this->get_time());
    }



    template <int dim>
    void
    InitialTemperature<dim>::constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &/*cell*/,
                                                 const std::vector<Point<dim>> &positions,
                                                 const std::vector<unsigned int> &component_indices,
                                                 std::vector<bool> &should_be_constrained,
                                                 std::vector<double> &solution)
    {

      const unsigned int temperature_index =
        this->introspection().component_indices.temperature;

      for (unsigned int q=0; q<positions.size(); ++q)
        {
          if (component_indices[q] != temperature_index)
            continue;

          const auto point =
            this->get_geometry_model().cartesian_to_other_coordinates(positions[q], coordinate_system);

          const double indicator =
            prescribed_temperature_indicator_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), 0);

          if (indicator > 0.5)
            {
              should_be_constrained[q] = true;
              solution[q] = initial_temperature_manager->initial_temperature(positions[q]);
            }
        }
    }



    template <int dim>
    void
    InitialTemperature<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed solution");
      {
        prm.enter_subsection("Initial temperature");
        {
          prm.declare_entry("Coordinate system",
                            "cartesian",
                            Patterns::Selection("cartesian|spherical|depth"),
                            "A selection that determines the assumed coordinate "
                            "system for the indicator function variables. "
                            "Allowed values are `cartesian', `spherical', and `depth'.");

          prm.enter_subsection("Indicator function");
          {
            Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    InitialTemperature<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed solution");
      {
        prm.enter_subsection("Initial temperature");
        {
          coordinate_system =
            Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));

          prm.enter_subsection("Indicator function");
          {
            try
              {
                prescribed_temperature_indicator_function.parse_parameters(prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed solution.Initial temperature.Indicator function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'\n";
                throw;
              }
          }
          prm.leave_subsection();
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
  namespace PrescribedSolution
  {
    ASPECT_REGISTER_PRESCRIBED_SOLUTION(InitialTemperature,
                                        "initial temperature",
                                        "Prescribe the temperature in a selected region using the active "
                                        "initial temperature model. The selected region is defined through "
                                        "an indicator function. At locations where the indicator value is greater than 0.5, "
                                        "the temperature is constrained to the initial temperature evaluated "
                                        "at that position.")
  }
}

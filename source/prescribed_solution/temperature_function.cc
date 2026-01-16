/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#include <aspect/prescribed_solution/temperature_function.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace PrescribedSolution
  {
    template <int dim>
    TemperatureFunction<dim>::TemperatureFunction ()
      :
      prescribed_temperature_indicator_function (1),
      prescribed_temperature_function (1)
    {}



    template <int dim>
    void TemperatureFunction<dim>::update()
    {
      if (this->convert_output_to_years())
        {
          prescribed_temperature_function.set_time(this->get_time() / year_in_seconds);
          prescribed_temperature_indicator_function.set_time(this->get_time() / year_in_seconds);
        }
      else
        {
          prescribed_temperature_function.set_time(this->get_time());
          prescribed_temperature_indicator_function.set_time(this->get_time());
        }
    }



    template <int dim>
    void
    TemperatureFunction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Prescribed solution");
      {
        prm.enter_subsection ("Temperature function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'.");

          prm.enter_subsection ("Function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();

          prm.enter_subsection ("Indicator function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    TemperatureFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Prescribed solution");
      {
        prm.enter_subsection ("Temperature function");
        {
          coordinate_system =
            Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));

          prm.enter_subsection ("Function");
          {
            try
              {
                prescribed_temperature_function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed solution.Temperature function.Function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'\n";
                throw;
              }
          }
          prm.leave_subsection();

          prm.enter_subsection ("Indicator function");
          {
            try
              {
                prescribed_temperature_indicator_function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed solution.Temperature function.Indicator function'\n"
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



    template <int dim>
    void
    TemperatureFunction<dim>::constrain_solution (
      const typename DoFHandler<dim>::active_cell_iterator &/*cell*/,
      const std::vector<Point<dim>> &positions,
      const std::vector<unsigned int> &component_indices,
      std::vector<bool> &should_be_constrained,
      std::vector<double> &solution)
    {
      const unsigned int temperature_index =
        this->introspection().component_indices.temperature;

      for (unsigned int q=0; q<positions.size(); ++q)
        {
          if (component_indices[q] == temperature_index)
            {
              const Utilities::NaturalCoordinate<dim> point =
                this->get_geometry_model()
                .cartesian_to_other_coordinates(positions[q], coordinate_system);

              const double indicator =
                prescribed_temperature_indicator_function.value(
                  Utilities::convert_array_to_point<dim>(point.get_coordinates()), 0);

              if (indicator > 0.5)
                {
                  should_be_constrained[q] = true;

                  solution[q] =
                    prescribed_temperature_function.value(
                      Utilities::convert_array_to_point<dim>(point.get_coordinates()), 0);
                }
            }
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace PrescribedSolution
  {
    ASPECT_REGISTER_PRESCRIBED_SOLUTION(TemperatureFunction,
                                        "temperature function",
                                        "Prescribe the temperature in terms of an explicit formula. "
                                        "The format of these functions follows the syntax understood by the "
                                        "muparser library.")
  }
}

/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/initial_composition/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/prescribed_solution/initial_composition.h>

namespace aspect
{
  namespace PrescribedSolution
  {

    template <int dim>
    InitialComposition<dim>::InitialComposition ()
      :
      indicator_function(1)
    {}



    template <int dim>
    void
    InitialComposition<dim>::initialize ()
    {
      initial_composition_manager =
        this->get_initial_composition_manager_pointer();
    }



    template <int dim>
    void
    InitialComposition<dim>::update ()
    {
      if (this->convert_output_to_years())
        indicator_function.set_time(this->get_time() / year_in_seconds);
      else
        indicator_function.set_time(this->get_time());
    }



    template <int dim>
    void
    InitialComposition<dim>::constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &,
                                                 const std::vector<Point<dim>> &positions,
                                                 const std::vector<unsigned int> &component_indices,
                                                 std::vector<bool> &should_be_constrained,
                                                 std::vector<double> &solution)
    {
      // Determine the component range corresponding to compositional fields.
      // These component indices are originally mapped from local DoFs and
      // determine which DoFs in the system belong to compositional fields.
      const unsigned int first_comp =
        this->introspection().component_indices.compositional_fields[0];

      const unsigned int last_comp =
        this->introspection().component_indices.compositional_fields[this->introspection().n_compositional_fields-1];

      const auto &geometry_model = this->get_geometry_model();

      for (unsigned int q=0; q<positions.size(); ++q)
        {
          const unsigned int component = component_indices[q];

          if (component < first_comp || component > last_comp)
            continue;

          const auto point = geometry_model.cartesian_to_other_coordinates(
                               positions[q], coordinate_system);

          const double indicator = indicator_function.value(
                                     Utilities::convert_array_to_point<dim>(point.get_coordinates()));

          if (indicator > 0.5)
            {
              const unsigned int field_index = component - first_comp;

              solution[q] =
                initial_composition_manager->initial_composition(
                  positions[q], field_index);

              should_be_constrained[q] = true;
            }
        }
    }



    template <int dim>
    void
    InitialComposition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed solution");
      {
        prm.enter_subsection("Initial composition");
        {

          prm.declare_entry("Coordinate system",
                            "cartesian",
                            Patterns::Selection("cartesian|spherical|depth"),
                            "Coordinate system used for evaluating the indicator function.");

          prm.enter_subsection("Indicator function");
          {
            Functions::ParsedFunction<dim>::declare_parameters(prm,1);
          }
          prm.leave_subsection();

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    InitialComposition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed solution");
      {
        prm.enter_subsection("Initial composition");
        {

          coordinate_system =
            Utilities::Coordinates::string_to_coordinate_system(
              prm.get("Coordinate system"));

          prm.enter_subsection("Indicator function");
          {
            indicator_function.parse_parameters(prm);
          }
          prm.leave_subsection();

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}



namespace aspect
{
  namespace PrescribedSolution
  {

    ASPECT_REGISTER_PRESCRIBED_SOLUTION(InitialComposition,
                                        "initial composition",
                                        "Prescribe compositional fields from the initial composition model "
                                        "within a region defined by an indicator function.")

  }
}

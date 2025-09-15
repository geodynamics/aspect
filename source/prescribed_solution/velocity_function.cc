/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/prescribed_solution/velocity_function.h>


namespace aspect
{
  namespace PrescribedSolution
  {
    template <int dim>
    VelocityFunction<dim>::VelocityFunction ()
      :
      prescribed_velocity_indicator_function (dim),
      prescribed_velocity_function (dim)
    {}



    template <int dim>
    void VelocityFunction<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        prescribed_velocity_function.set_time (this->get_time() / year_in_seconds);
      else
        prescribed_velocity_function.set_time (this->get_time());

      if (this->convert_output_to_years())
        prescribed_velocity_indicator_function.set_time (this->get_time() / year_in_seconds);
      else
        prescribed_velocity_indicator_function.set_time (this->get_time());
    }



    template <int dim>
    void
    VelocityFunction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Prescribed solution");
      {
        prm.enter_subsection("Velocity function");
        {
          /**
           * Choose the coordinates to evaluate the maximum refinement level
           * function. The function can be declared in dependence of depth,
           * cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2d/3d "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          prm.enter_subsection("Function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
          }
          prm.leave_subsection();

          prm.enter_subsection("Indicator function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    VelocityFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Prescribed solution");
      {
        prm.enter_subsection("Velocity function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));

          prm.enter_subsection("Function");
          {
            try
              {
                prescribed_velocity_function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed fields model.Function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'"
                          << "More information about the cause of the parse error \n"
                          << "is shown below.\n";
                throw;
              }
          }
          prm.leave_subsection();

          prm.enter_subsection("Indicator function");
          {
            try
              {
                prescribed_velocity_indicator_function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed fields model.Indicator function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'"
                          << "More information about the cause of the parse error \n"
                          << "is shown below.\n";
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
    VelocityFunction<dim>::constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &/*cell*/,
                                               const std::vector<Point<dim>> &positions,
                                               const std::vector<unsigned int> component_indices,
                                               std::vector<bool> &should_be_constrained,
                                               std::vector<double> &solution)
    {
      for (unsigned int q=0; q<positions.size(); ++q)
        {
          const unsigned int c_idx = component_indices[q];

          // If we're on one of the velocity DOFs
          if ((c_idx >= this->introspection().component_indices.velocities[0])
              &&
              (c_idx <= this->introspection().component_indices.velocities[dim-1]))
            {
              // Which velocity component is this DOF associated with?
              const unsigned int component_direction
                = (c_idx - this->introspection().component_indices.velocities[0]);

              // TODO: the position needs to be converted into the appropriate coordinate system see boundary_velocity/function.cc

              const double indicator = prescribed_velocity_indicator_function.value(positions[q],
                                                                                    component_direction);

              if (indicator > 0.5)
                {
                  should_be_constrained[q] = true;

                  // TODO: the velocity needs to account for "Use years in seconds" see boundary_velocity/function.cc
                  solution[q] = prescribed_velocity_function.value(positions[q],
                                                                   component_direction);
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
    ASPECT_REGISTER_PRESCRIBED_SOLUTION(VelocityFunction,
                                        "velocity function",
                                        "Specify the prescribed fields in terms of an "
                                        "explicit formula. The format of these "
                                        "functions follows the syntax understood by the "
                                        "muparser library, see {ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}

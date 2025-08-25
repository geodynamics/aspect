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


#include <aspect/prescribed_fields/velocity_function.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_signals.h>

#include <iostream>

namespace aspect
{
  namespace PrescribedFields
  {
    /**
    * A set of helper functions that either return the point passed to it (if
    * the current dimension is the same) or return a dummy value (otherwise).
    */
    template <int dim>
    VelocityFunction<dim>::VelocityFunction ()
      :
      prescribed_velocity_indicator_function (dim),
      prescribed_velocity_function (dim)
    {}


    template <int dim>
    void VelocityFunction<dim>::update(const SimulatorAccess<dim> &simulator_access)
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (simulator_access.convert_output_to_years())
        prescribed_velocity_function.set_time (simulator_access.get_time() / year_in_seconds);
      else
        prescribed_velocity_function.set_time (simulator_access.get_time());

      if (simulator_access.convert_output_to_years())
        prescribed_velocity_indicator_function.set_time (simulator_access.get_time() / year_in_seconds);
      else
        prescribed_velocity_indicator_function.set_time (simulator_access.get_time());
    }


    template <int dim>
    void
    VelocityFunction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Prescribed fields");
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
      prm.enter_subsection ("Prescribed fields");
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
    VelocityFunction<dim>::constrain_internal_fields (const SimulatorAccess<dim> &simulator_access,
                                                      AffineConstraints<double> &current_constraints)
    {
      const std::vector<Point<dim>> points = aspect::Utilities::get_unit_support_points(simulator_access);
      const Quadrature<dim> quadrature (points);
      FEValues<dim> fe_values (simulator_access.get_fe(), quadrature, update_quadrature_points);
      typename DoFHandler<dim>::active_cell_iterator cell;

      // Loop over all cells
      for (cell = simulator_access.get_dof_handler().begin_active();
           cell != simulator_access.get_dof_handler().end();
           ++cell)
        if (! cell->is_artificial())
          {
            fe_values.reinit (cell);
            std::vector<types::global_dof_index> local_dof_indices(simulator_access.get_fe().dofs_per_cell);
            cell->get_dof_indices (local_dof_indices);

            for (unsigned int q=0; q<quadrature.size(); q++)
              // If it's okay to constrain this DOF
              if (current_constraints.can_store_line(local_dof_indices[q]) &&
                  !current_constraints.is_constrained(local_dof_indices[q]))
                {
                  // Get the velocity component index
                  const unsigned int c_idx =
                    simulator_access.get_fe().system_to_component_index(q).first;

                  // If we're on one of the velocity DOFs
                  if ((c_idx >=
                       simulator_access.introspection().component_indices.velocities[0])
                      &&
                      (c_idx <=
                       simulator_access.introspection().component_indices.velocities[dim-1]))
                    {
                      // Which velocity component is this DOF associated with?
                      const unsigned int component_direction
                        = (c_idx
                           - simulator_access.introspection().component_indices.velocities[0]);

                      const Point<dim> p = fe_values.quadrature_point(q);

                      // Because we defined and parsed our parameter
                      // file differently for 2d and 3d we need to
                      // be sure to query the correct object for
                      // function values. The function parser
                      // objects expect points of a certain
                      // dimension, but Point p will be compiled for
                      // both 2d and 3d, so we need to do some trickery
                      // to make this compile.
                      const double indicator = prescribed_velocity_indicator_function.value
                                               (p,
                                                component_direction);
                      const double u_i = prescribed_velocity_function.value
                                         (p,
                                          component_direction);

                      if (indicator > 0.5)
                        {
                          // Add a constraint of the form dof[q] = u_i
                          // to the list of constraints.
                          current_constraints.add_line (local_dof_indices[q]);
                          current_constraints.set_inhomogeneity (local_dof_indices[q], u_i);
                        }
                    }
                }
          }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedFields
  {
    ASPECT_REGISTER_PRESCRIBED_FIELDS(VelocityFunction,
                                      "velocity function",
                                      "Specify the prescribed fields in terms of an "
                                      "explicit formula. The format of these "
                                      "functions follows the syntax understood by the "
                                      "muparser library, see {ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}

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
#include <aspect/geometry_model/interface.h>


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
      // to reinterpret it in years if the global flag
      // for using years instead of seconds is given
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
           * Choose the coordinates to evaluate the velocity. The function can be
           * declared in dependence of depth, cartesian coordinates or spherical
           * coordinates. Note that the order of spherical coordinates is
           * r,phi,theta and not r,theta,phi, since this allows for dimension
           * independent expressions.
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
          prm.declare_entry ("Use spherical unit vectors", "false",
                             Patterns::Bool (),
                             "Specify velocity as $r$, $\\phi$, and $\\theta$ components "
                             "instead of $x$, $y$, and $z$. Positive velocities point up, east, "
                             "and north (in 3d) or out and clockwise (in 2d). "
                             "This setting only makes sense for spherical geometries.");

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
          use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
          if (use_spherical_unit_vectors)
            AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                         ExcMessage ("Spherical unit vectors should not be used "
                                     "when geometry model is not spherical."));

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
                                               const std::vector<unsigned int> &component_indices,
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


              // The position is converted into the appropriate coordinate system.
              const Utilities::NaturalCoordinate<dim> point =
                this->get_geometry_model().cartesian_to_other_coordinates(positions[q], coordinate_system);

              double indicator = 0.0;
              if (use_spherical_unit_vectors)
                {
                  // If we use spherical unit vectors a velocity in the direction of any
                  // spherical basis vector will constrain all Cartesian velocity components
                  for (unsigned int d=0; d<dim; ++d)
                    {
                      indicator = std::max(indicator, prescribed_velocity_indicator_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),
                                                                                                   d));
                    }
                }
              else
                {
                  indicator = prescribed_velocity_indicator_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),
                                                                           component_direction);
                }

              if (indicator > 0.5)
                {
                  should_be_constrained[q] = true;


                  Tensor<1,dim> velocity;

                  for (unsigned int d=0; d<dim; ++d)
                    velocity[d] = prescribed_velocity_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),
                                                                     d);
                  if (use_spherical_unit_vectors)
                    velocity = Utilities::Coordinates::spherical_to_cartesian_vector(velocity, positions[q]);

                  // ASPECT always wants things in MKS system. however, as described
                  // in the documentation of this class, we interpret the formulas
                  // given to this plugin as meters per year if the global flag
                  // for using years instead of seconds is given. so if someone
                  // write "5" in their parameter file and sets the flag, then this
                  // means "5 meters/year" and we need to convert it to the ASPECT
                  // time system by dividing by the number of seconds per year
                  // to get MKS units
                  if (this->convert_output_to_years())
                    solution[q] = velocity[component_direction] / year_in_seconds;
                  else
                    solution[q] = velocity[component_direction];
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
                                        "explicit formula in a region that is "
                                        "functions follows the syntax understood by the "
                                        "muparser library, see {ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}

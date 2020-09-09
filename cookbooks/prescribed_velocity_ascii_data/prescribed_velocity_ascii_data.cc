/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_signals.h>
#include <aspect/initial_composition/ascii_data.h>

namespace aspect
{
  using namespace dealii;

  // Global variables (to be set by parameters)
  bool prescribe_internal_velocities;
  bool prescribed_velocity_ascii_file_loaded;

  // Because we do not initially know what dimension we're in, we need
  // function parser objects for both 2d and 3d.

  aspect::Utilities::AsciiDataBase<2> prescribed_velocity_field_params_2d;
  aspect::Utilities::AsciiDataBase<3> prescribed_velocity_field_params_3d;

  aspect::Utilities::AsciiDataLookup<2> prescribed_velocity_field_data_2d(1, 1.);
  aspect::Utilities::AsciiDataLookup<3> prescribed_velocity_field_data_3d(1, 1.);

  Functions::ParsedFunction<2> prescribed_velocity_function_2d (2);
  Functions::ParsedFunction<3> prescribed_velocity_function_3d (3);

  /**
   * Declare additional parameters.
   */
  void declare_parameters(const unsigned int dim,
                          ParameterHandler &prm)
  {
    prm.declare_entry ("Prescribe internal velocities", "false",
                       Patterns::Bool (),
                       "Whether or not to use any prescribed internal velocities. "
                       "Locations in which to prescribe velocities are defined "
                       "in section ``Prescribed velocities/Indicator function'' "
                       "and the velocities are defined in section ``Prescribed "
                       "velocities/Velocity function''. Indicators are evaluated "
                       "at the center of each cell, and all DOFs associated with "
                       "the specified velocity component at the indicated cells "
                       "are constrained."
                      );

    prm.enter_subsection ("Prescribed velocities");
    {
      {
        if (dim == 2)
          {
            prescribed_velocity_field_params_2d.declare_parameters(prm, "blob", "blob", "Ascii data model");
            //prescribed_velocity_field_data_2d.initialize(1);
          }
        else
          {
            prescribed_velocity_field_params_3d.declare_parameters(prm, "blob", "blob", "Ascii data model");
            //prescribed_velocity_field_data_3d.initialize(1);
          }
      }

      prm.enter_subsection ("Velocity function");
      {
        if (dim == 2)
          Functions::ParsedFunction<2>::declare_parameters (prm, 2);
        else
          Functions::ParsedFunction<3>::declare_parameters (prm, 3);
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void parse_parameters(const Parameters<dim> &,
                        ParameterHandler &prm)
  {
    prescribe_internal_velocities = prm.get_bool ("Prescribe internal velocities");
    prm.enter_subsection ("Prescribed velocities");
    {
      {
        if (dim == 2)
          {
            prescribed_velocity_field_params_2d.parse_parameters(prm, "Ascii data model");
          }
        else
          {
            prescribed_velocity_field_params_3d.parse_parameters(prm, "Ascii data model");
          }
      }

      prm.enter_subsection("Velocity function");
      {
        try
          {
            if (dim == 2)
              prescribed_velocity_function_2d.parse_parameters (prm);
            else
              prescribed_velocity_function_3d.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed velocities.Velocity function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
      }
      prm.leave_subsection();
    }
    prm.leave_subsection ();

    prescribed_velocity_ascii_file_loaded = false;
  }



  /**
   * This function retrieves the unit support points (in the unit cell) for the current element.
   * The DGP element used when 'set Use locally conservative discretization = true' does not
   * have support points. If these elements are in use, a fictitious support point at the cell
   * center is returned for each shape function that corresponds to the pressure variable,
   * whereas the support points for the velocity are correct; the fictitious points don't matter
   * because we only use this function when interpolating the velocity variable, and ignore the
   * evaluation at the pressure support points.
   */
  template <int dim>
  std::vector< Point<dim> >
  get_unit_support_points_for_velocity(const SimulatorAccess<dim> &simulator_access)
  {
    std::vector< Point<dim> > unit_support_points;
    if ( !simulator_access.get_parameters().use_locally_conservative_discretization )
      {
        return simulator_access.get_fe().get_unit_support_points();
      }
    else
      {
        //special case for discontinuous pressure elements, which lack unit support points
        std::vector< Point<dim> > unit_support_points;
        const unsigned int dofs_per_cell = simulator_access.get_fe().dofs_per_cell;
        for (unsigned int dof=0; dof < dofs_per_cell; ++dof)
          {
            // base will hold element, base_index holds node/shape function within that element
            const unsigned int
            base       = simulator_access.get_fe().system_to_base_index(dof).first.first,
            base_index = simulator_access.get_fe().system_to_base_index(dof).second;
            // get the unit support points for the relevant element
            std::vector< Point<dim> > my_support_points = simulator_access.get_fe().base_element(base).get_unit_support_points();
            if ( my_support_points.size() == 0 )
              {
                //manufacture a support point, arbitrarily at cell center
                if ( dim==2 )
                  unit_support_points.push_back( Point<dim> (0.5,0.5) );
                if ( dim==3 )
                  unit_support_points.push_back( Point<dim> (0.5,0.5,0.5) );
              }
            else
              {
                unit_support_points.push_back(my_support_points[base_index]);
              }
          }
        return unit_support_points;
      }
  }

  /**
   * A set of helper functions that either return the point passed to it (if
   * the current dimension is the same) or return a dummy value (otherwise).
   */
  namespace
  {
    const Point<2> as_2d(const Point<3> &/*p*/)
    {
      return Point<2>();
    }

    const Point<2> &as_2d(const Point<2> &p)
    {
      return p;
    }

    const Point<3> as_3d(const Point<2> &/*p*/)
    {
      return Point<3>();
    }

    const Point<3> &as_3d(const Point<3> &p)
    {
      return p;
    }

  }


  /**
   * This function is called by a signal which is triggered after the other constraints
   * have been calculated. This enables us to define additional constraints in the mass
   * matrix on any arbitrary degree of freedom in the model space.
   */
  template <int dim>
  void constrain_internal_velocities (const SimulatorAccess<dim> &simulator_access,
                                      AffineConstraints<double> &current_constraints)
  {
    if (!prescribed_velocity_ascii_file_loaded)
      {
        if (dim == 2)
          {
            const std::string filename = prescribed_velocity_field_params_2d.data_directory + prescribed_velocity_field_params_2d.data_file_name;
            prescribed_velocity_field_data_2d.load_file(filename, simulator_access.get_mpi_communicator());

            std::cout << "Loaded 2D prescribed velocity ascii file" << std::endl;
          }
        else
          {
            const std::string filename = prescribed_velocity_field_params_3d.data_directory + prescribed_velocity_field_params_3d.data_file_name;
            prescribed_velocity_field_data_3d.load_file(filename, simulator_access.get_mpi_communicator());

            std::cout << "Loaded 3D prescribed velocity ascii file" << std::endl;
          }
        prescribed_velocity_ascii_file_loaded = true;
      }

    if (prescribe_internal_velocities)
      {
        const std::vector< Point<dim> > points = get_unit_support_points_for_velocity(simulator_access);
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

                        // we get time passed as seconds (always) but may want
                        // to reinterpret it in years
                        double time = simulator_access.get_time();
                        if (simulator_access.convert_output_to_years())
                          time /= year_in_seconds;

                        prescribed_velocity_function_2d.set_time (time);
                        prescribed_velocity_function_3d.set_time (time);

                        const Point<dim> p = fe_values.quadrature_point(q);

                        // Because we defined and parsed our parameter
                        // file differently for 2d and 3d we need to
                        // be sure to query the correct object for
                        // function values. The function parser
                        // objects expect points of a certain
                        // dimension, but Point p will be compiled for
                        // both 2d and 3d, so we need to do some trickery
                        // to make this compile.
                        double indicator, u_i;
                        if (dim == 2)
                          {
                            indicator = prescribed_velocity_field_data_2d.get_data(as_2d(p),0);
                            u_i       = prescribed_velocity_function_2d.value
                                        (as_2d(p),
                                         component_direction);
                          }
                        else
                          {
                            indicator = prescribed_velocity_field_data_3d.get_data(as_3d(p),0);
                            u_i       = prescribed_velocity_function_3d.value
                                        (as_3d(p),
                                         component_direction);
                          }

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

  // Connect declare_parameters and parse_parameters to appropriate signals.
  void parameter_connector ()
  {
    SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
    SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

    SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
    SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
  }

  // Connect constraints function to correct signal.
  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.post_constraints_creation.connect (&constrain_internal_velocities<dim>);
  }

  // Tell ASPECT to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}

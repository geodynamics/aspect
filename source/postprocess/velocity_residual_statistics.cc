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



#include <aspect/postprocess/velocity_residual_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    VelocityResidualStatistics<dim>::initialize ()
    {
      // The input ascii table contains dim data columns (velocity components) in addition to the coordinate columns.
      data_lookup = std::make_unique<Utilities::StructuredDataLookup<dim>>(dim, scale_factor);
      data_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
    }



    template <int dim>
    Tensor<1,dim>
    VelocityResidualStatistics<dim>::get_data_velocity (const Point<dim> &p) const
    {
      Tensor<1,dim> data_velocity;
      Point<dim> internal_position = p;

      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
        {
          const std::array<double,dim> spherical_position = this->get_geometry_model().
                                                            cartesian_to_natural_coordinates(p);

          for (unsigned int d = 0; d < dim; ++d)
            internal_position[d] = spherical_position[d];
        }

      for (unsigned int d = 0; d < dim; ++d)
        data_velocity[d] = data_lookup->get_data(internal_position, d);
      if (use_spherical_unit_vectors == true)
        data_velocity = Utilities::Coordinates::spherical_to_cartesian_vector(data_velocity, p);

      return data_velocity;
    }



    template <int dim>
    std::pair<std::string,std::string>
    VelocityResidualStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula for the velocity.
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.velocities;
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<Tensor<1,dim>> velocity_values(n_q_points);

      // compute the maximum, minimum, and squared*area velocity residual
      // magnitude and the face area.
      double local_vel_residual_square_integral = 0;
      double local_max_vel_residual             = std::numeric_limits<double>::lowest();
      double local_min_vel_residual             = std::numeric_limits<double>::max();
      double local_fe_volume                    = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const Point<dim> data_point = fe_values.quadrature_point(q);
                Tensor<1,dim> data_velocity = get_data_velocity(data_point);

                const double vel_residual_mag = (velocity_values[q] - data_velocity).norm();

                local_max_vel_residual = std::max(vel_residual_mag, local_max_vel_residual);
                local_min_vel_residual = std::min(vel_residual_mag, local_min_vel_residual);

                local_vel_residual_square_integral += ((vel_residual_mag * vel_residual_mag) * fe_values.JxW(q));
                local_fe_volume += fe_values.JxW(q);
              }
          }

      const double global_vel_res_square_integral
        = Utilities::MPI::sum (local_vel_residual_square_integral, this->get_mpi_communicator());
      const double global_max_vel_residual
        = Utilities::MPI::max (local_max_vel_residual, this->get_mpi_communicator());
      const double global_min_vel_residual
        = Utilities::MPI::min (local_min_vel_residual, this->get_mpi_communicator());

      const double vrms_residual = std::sqrt(global_vel_res_square_integral) /
                                   std::sqrt(this->get_volume());

      // now add the computed max, min, and rms velocities to the statistics object
      // and create a single string that can be output to the screen
      const std::string units = (this->convert_output_to_years() == true) ? "m/year" : "m/s";
      const double unit_scale_factor = (this->convert_output_to_years() == true) ? year_in_seconds : 1.0;

      const std::vector<std::string> column_names = {"RMS velocity residual  (" + units + ")",
                                                     "Max. velocity residual (" + units + ")",
                                                     "Min. velocity residual (" + units + ")"
                                                    };

      statistics.add_value (column_names[0],
                            vrms_residual * unit_scale_factor);
      statistics.add_value (column_names[1],
                            global_max_vel_residual * unit_scale_factor);
      statistics.add_value (column_names[2],
                            global_min_vel_residual * unit_scale_factor);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      for (auto &column : column_names)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      std::ostringstream screen_text;
      screen_text.precision(3);
      screen_text << vrms_residual *unit_scale_factor
                  << ' ' << units << ", "
                  << global_max_vel_residual *unit_scale_factor
                  << ' ' << units << ", "
                  << global_min_vel_residual *unit_scale_factor
                  << ' ' << units;

      return std::pair<std::string, std::string> ("RMS, max, and min velocity residual velocity in the model:",
                                                  screen_text.str());
    }



    template <int dim>
    void
    VelocityResidualStatistics<dim>::declare_parameters (ParameterHandler  &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Velocity residual statistics");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/postprocess/velocity-residual/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the input ascii data "
                             "against which the velocity residual in the model is computed. "
                             "This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name", "box_2d_velocity.txt",
                             Patterns::Anything (),
                             "The file name of the input ascii velocity data. "
                             "The file is provided in the same format as described in "
                             " 'ascii data' initial composition plugin." );
          prm.declare_entry ("Scale factor", "1.",
                             Patterns::Double (),
                             "Scalar factor, which is applied to the model data. "
                             "You might want to use this to scale the input to a "
                             "reference model. Another way to use this factor is to "
                             "convert units of the input files. For instance, if you "
                             "provide velocities in cm/year set this factor to 0.01.");
          prm.declare_entry ("Use spherical unit vectors", "false",
                             Patterns::Bool (),
                             "Specify velocity as r, phi, and theta components "
                             "instead of x, y, and z. Positive velocities point up, east, "
                             "and north (in 3d) or out and clockwise (in 2d). "
                             "This setting only makes sense for spherical geometries.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    VelocityResidualStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Velocity residual statistics");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          data_file_name    = prm.get ("Data file name");
          scale_factor      = prm.get_double ("Scale factor");

          use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");

          if (use_spherical_unit_vectors)
            AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                         ExcMessage ("Spherical unit vectors should not be used "
                                     "when geometry model is not spherical."));
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VelocityResidualStatistics,
                                  "velocity residual statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the velocity residual in the model. The velocity residual "
                                  "is the difference between the model solution velocities and the input "
                                  "ascii data velocities.")
  }
}

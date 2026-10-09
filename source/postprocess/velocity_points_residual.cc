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



#include <aspect/postprocess/velocity_points_residual.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    VelocityPointsResidual<dim>::VelocityPointsResidual ()
      :
      output_file_number (numbers::invalid_unsigned_int)
    {}

    template <int dim>
    void
    VelocityPointsResidual<dim>::initialize ()
    {
      // The input ascii table contains one column that just represents points id.
      // The scale factor is set to 1 because the columns contain the observed point
      // coordinates and velocities, which we do not want to scale.
      data_lookup = std::make_unique<Utilities::StructuredDataLookup<1>>(2*dim, 1.0);
      data_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
    }



    template <int dim>
    std::pair <Point<dim>, Tensor<1,dim>>
    VelocityPointsResidual<dim>::get_observed_data (const unsigned int p) const
    {
      Tensor<1,dim> data_velocity;
      Point<dim> data_position;
      Point<1> point;
      point[0] = p;

      for (unsigned int d = 0; d < dim; ++d)
        {
          data_position[d] = data_lookup->get_data(point, d);
          data_velocity[d] = data_lookup->get_data(point, d+dim);
        }

      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
          use_spherical_unit_vectors == true)
        {
          std::array<double,dim> internal_position;

          for (unsigned int d = 0; d < dim; ++d)
            internal_position[d] = data_position[d];

          Point<dim> cartesian_position = this->get_geometry_model().
                                          natural_to_cartesian_coordinates(internal_position);

          data_position = cartesian_position;
          data_velocity = Utilities::Coordinates::spherical_to_cartesian_vector(data_velocity, cartesian_position);
        }

      return std::make_pair (data_position, data_velocity);
    }



    template <int dim>
    std::pair<std::string,std::string>
    VelocityPointsResidual<dim>::execute (TableHandler &statistics)
    {
      // Add the rms residual in the statistics file.
      const std::string units = (this->convert_output_to_years() == true) ? "m/year" : "m/s";
      const double unit_scale_factor = (this->convert_output_to_years() == true) ? year_in_seconds : 1.0;

      const bool increase_file_number = (this->get_nonlinear_iteration() == 0) ||
                                        (!this->get_parameters().run_postprocessors_on_nonlinear_iterations);
      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else if (increase_file_number)
        ++output_file_number;

      const std::string file_prefix = "velocity-residual-" + Utilities::int_to_string (output_file_number, 5);
      const std::string filename = (this->get_output_directory()
                                    + file_prefix);

      const unsigned int n_rows = data_lookup->get_interpolation_point_coordinates(0).size();

      // evaluate the solution at all of our evaluation points
      std::vector<Vector<double>>
      current_point_values (n_rows, Vector<double> (this->introspection().n_components));

      // velocity residual at the evaluation points
      std::vector<Tensor<1,dim>> velocity_residual (n_rows);

      double rms_velocity_residual = 0.;

      for (unsigned int p=0; p<n_rows; ++p)
        {
          // try to evaluate the solution at this point. in parallel, the point
          // will be on only one processor's owned cells, so the others are
          // going to throw an exception. make sure at least one processor
          // finds the given point
          bool point_found = false;

          try
            {
              VectorTools::point_value(this->get_mapping(),
                                       this->get_dof_handler(),
                                       this->get_solution(),
                                       get_observed_data(p).first,
                                       current_point_values[p]);
              point_found = true;
            }
          catch (const VectorTools::ExcPointNotAvailableHere &)
            {
              // ignore
            }

          // ensure that at least one processor found things
          const int n_procs = Utilities::MPI::sum (point_found ? 1 : 0, this->get_mpi_communicator());
          AssertThrow (n_procs > 0,
                       ExcMessage ("While trying to evaluate the solution at point " +
                                   Utilities::to_string(get_observed_data(p).first[0]) + ", " +
                                   Utilities::to_string(get_observed_data(p).first[1]) +
                                   (dim == 3
                                    ?
                                    ", " + Utilities::to_string(get_observed_data(p).first[2])
                                    :
                                    "") + "), " +
                                   "no processors reported that the point lies inside the " +
                                   "set of cells they own. Are you trying to evaluate the " +
                                   "solution at a point that lies outside of the domain?"
                                  ));

          // Reduce all collected values into local Vector
          Utilities::MPI::sum (current_point_values[p], this->get_mpi_communicator(),
                               current_point_values[p]);

          Tensor<1,dim> modeled_velocity_solution;

          // Normalize in cases where points are claimed by multiple processors
          for (unsigned int d=0; d<dim; ++d)
            modeled_velocity_solution[d] = current_point_values[p][d];

          Tensor<1,dim> observed_velocity;
          if (this->convert_output_to_years() == true)
            observed_velocity = get_observed_data(p).second/year_in_seconds;

          velocity_residual[p] = modeled_velocity_solution - observed_velocity;

          rms_velocity_residual += velocity_residual[p].norm();
        }

      const double global_rms_velocity_residual = rms_velocity_residual/n_rows;

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {

          std::ofstream f (filename);
          f << ("# <point_id> "
                "<evaluation_point_x> "
                "<evaluation_point_y> ")
            << (dim == 3 ? "<evaluation_point_z> " : "")
            << ("<velocity_residual_x> "
                "<velocity_residual_y> ")
            << (dim == 3 ? "<velocity_residual_z> " : "");

          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            f << " <" << this->introspection().name_for_compositional_index(c) << '>';
          f << '\n';
          f << "# POINTS: " << n_rows << '\n';

          for (unsigned int i=0; i<n_rows; ++i)
            {
              f << /* time = */ i
                << ' '
                << /* evaluation point = */ get_observed_data(i).first << ' ';

              for (unsigned int c=0; c<dim; ++c)
                {
                  // output a data element. internally, we store all point
                  // values in the same format in which they were computed,
                  // but we convert velocities to meters per year if so
                  // requested
                  if ((this->introspection().component_masks.velocities[c] == false)
                      ||
                      (this->convert_output_to_years() == false))
                    f << velocity_residual[i][c];
                  else
                    f << velocity_residual[i][c] * year_in_seconds;

                  f << ' ';
                }
              // have an empty line between points
              f << '\n';
            }

          AssertThrow (f, ExcMessage("Writing data to <" + filename +
                                     "> did not succeed in the `point values' "
                                     "postprocessor."));
        }

      // return what should be printed to the screen. note that we had
      // just incremented the number, so use the previous value

      const std::vector<std::string> column_names = {"RMS velocity residual  (" + units + ")"};

      statistics.add_value (column_names[0],
                            global_rms_velocity_residual * unit_scale_factor);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      for (auto &column : column_names)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      return std::make_pair (std::string ("Velocity residual file and RMS:"),
                             filename + " " + Utilities::to_string(global_rms_velocity_residual * unit_scale_factor));
    }



    template <int dim>
    void
    VelocityPointsResidual<dim>::declare_parameters (ParameterHandler  &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Velocity points residual");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/postprocess/velocity-points-residual/",
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
          prm.declare_entry ("Data file name", "point_values_velocity.txt",
                             Patterns::Anything (),
                             "The file name of the input ascii velocity data. "
                             "The file is provided in a structured data format in one dimension "
                             "such that the first column represents point ids from zero to n-1 "
                             " (for n number of points) and the subsequent columns are data coordinates "
                             "and the observed velocities at those points.");
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
    VelocityPointsResidual<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Velocity points residual");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          data_file_name    = prm.get ("Data file name");

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
    ASPECT_REGISTER_POSTPROCESSOR(VelocityPointsResidual,
                                  "velocity points residual",
                                  "A postprocessor that computes root mean square velocity "
                                  "residual as \\frac{1}{n} \\sum_n \\lVert \\mathbf u_i^{\\text{solution}} -"
                                  "\\mathbf u_i^{\\text{input}} \\rVert_2, where \\lVert\\cdot\\rVert "
                                  "denotes the L2-norm of the difference between the model solution "
                                  "velocities and the input ascii data velocities at a specified "
                                  "evaluation point i in the input data for n numer of data points.")
  }
}

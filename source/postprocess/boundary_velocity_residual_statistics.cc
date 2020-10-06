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



#include <aspect/postprocess/boundary_velocity_residual_statistics.h>
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
    BoundaryVelocityResidualStatistics<dim>::initialize ()
    {
      if (use_ascii_data)
        {
          // The input ascii table contains dim data columns (velocity components) in addition to the coordinate columns.
          ascii_data_lookup = std_cxx14::make_unique<Utilities::AsciiDataLookup<dim> >(dim, scale_factor);
          ascii_data_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
        }
      else
        {
          // The two points are used in GPlates to find the 2D plane in which
          // the model lies.  These values are not used for 3D geometries and
          // are thus set to the default.
          Point<2> point_one(numbers::PI/2., 0.);
          Point<2> point_two(numbers::PI/2., numbers::PI/2.);

          gplates_lookup =  std_cxx14::make_unique<BoundaryVelocity::internal::GPlatesLookup<dim> > (
                              point_one, point_two);

          gplates_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
        }
    }



    template <int dim>
    Tensor<1,dim>
    BoundaryVelocityResidualStatistics<dim>::get_data_velocity (const Point<dim> &p) const
    {
      Tensor<1,dim> data_velocity;
      if (use_ascii_data)
        {
          Point<dim> internal_position = p;

          if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
            {
              const std::array<double,dim> spherical_position = this->get_geometry_model().
                                                                cartesian_to_natural_coordinates(p);

              for (unsigned int d = 0; d < dim; ++d)
                internal_position[d] = spherical_position[d];
            }

          for (unsigned int d = 0; d < dim; ++d)
            data_velocity[d] = ascii_data_lookup->get_data(internal_position, d);
          if (use_spherical_unit_vectors == true)
            data_velocity = Utilities::Coordinates::spherical_to_cartesian_vector(data_velocity, p);
        }
      else
        {
          data_velocity =  gplates_lookup->surface_velocity(p);
        }

      return data_velocity;
    }



    template <int dim>
    std::pair<std::string,std::string>
    BoundaryVelocityResidualStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula for the velocity.
      const QGauss<dim-1> quadrature_formula (this->introspection().polynomial_degree.velocities+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_values |
                                        update_JxW_values |
                                        update_quadrature_points);

      std::vector<Tensor<1,dim> > velocities (fe_face_values.n_quadrature_points);

      std::map<types::boundary_id, double> local_max_vel;
      std::map<types::boundary_id, double> local_min_vel;
      std::map<types::boundary_id, double> local_velocity_square_integral;
      std::map<types::boundary_id, double> local_boundary_area;

      std::set<types::boundary_id> boundary_indicators;
      boundary_indicators.insert(this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"));

      for (const auto p : boundary_indicators)
        {
          local_max_vel[p] = -std::numeric_limits<double>::max();
          local_min_vel[p] = std::numeric_limits<double>::max();
        }

      // for every face that is part of the selected geometry boundary
      // and that is owned by this processor,
      // compute the maximum, minimum, and squared*area velocity residual
      // magnitude and the face area.

      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_face_values.reinit (cell, f);

                fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                    velocities);

                // determine the max, min, and squared velocity residual on the face
                // also determine the face area
                double local_max = -std::numeric_limits<double>::max();
                double local_min = std::numeric_limits<double>::max();
                double local_sqvel = 0.0;
                double local_fe_face_area = 0.0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const Point<dim> point_at_surface = fe_face_values.quadrature_point(q);
                    // Extract data velocity.
                    Tensor<1,dim> data_velocity = get_data_velocity(point_at_surface);

                    if (this->convert_output_to_years() == true)
                      data_velocity = data_velocity/year_in_seconds;

                    // The velocity residual is calculated here.
                    const double vel_residual_mag = (velocities[q] - data_velocity).norm();

                    local_max = std::max(vel_residual_mag,
                                         local_max);
                    local_min = std::min(vel_residual_mag,
                                         local_min);
                    local_sqvel += ((vel_residual_mag * vel_residual_mag) * fe_face_values.JxW(q));
                    local_fe_face_area += fe_face_values.JxW(q);
                  }

                // then merge them with the min/max/squared velocities
                // and face areas we found for other faces with the same boundary indicator
                local_max_vel[top_boundary_id] = std::max(local_max,
                                                          local_max_vel[top_boundary_id]);
                local_min_vel[top_boundary_id] = std::min(local_min,
                                                          local_min_vel[top_boundary_id]);
                local_velocity_square_integral[top_boundary_id] += local_sqvel;
                local_boundary_area[top_boundary_id] += local_fe_face_area;
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_max_vel;
      std::map<types::boundary_id, double> global_min_vel;
      std::map<types::boundary_id, double> global_rms_vel;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        std::vector<double> local_max_values;
        std::vector<double> local_min_values;
        std::vector<double> local_velocity_square_integral_values;
        std::vector<double> local_boundary_area_values;

        for (const auto p : boundary_indicators)
          {
            local_max_values.push_back (local_max_vel[p]);
            local_min_values.push_back (local_min_vel[p]);

            local_velocity_square_integral_values.push_back (local_velocity_square_integral[p]);
            local_boundary_area_values.push_back (local_boundary_area[p]);
          }
        // then collect contributions from all processors
        std::vector<double> global_max_values (local_max_values.size());
        Utilities::MPI::max (local_max_values, this->get_mpi_communicator(), global_max_values);
        std::vector<double> global_min_values (local_min_values.size());
        Utilities::MPI::min (local_min_values, this->get_mpi_communicator(), global_min_values);

        std::vector<double> global_velocity_square_integral_values (local_velocity_square_integral_values.size());
        Utilities::MPI::sum(local_velocity_square_integral_values,this->get_mpi_communicator(),global_velocity_square_integral_values);
        std::vector<double> global_boundary_area_values (local_boundary_area_values.size());
        Utilities::MPI::sum(local_boundary_area_values,this->get_mpi_communicator(),global_boundary_area_values);

        // and now take them apart into the global map again
        // At the same time, calculate the rms velocity for each boundary
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          {
            global_max_vel[*p] = global_max_values[index];
            global_min_vel[*p] = global_min_values[index];
            global_rms_vel[*p] = std::sqrt(global_velocity_square_integral_values[index] / global_boundary_area_values[index]);
          }
      }

      // now add the computed max, min, and rms velocities to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_max_vel.begin(), a = global_min_vel.begin(), rms = global_rms_vel.begin();
           p != global_max_vel.end() && a != global_min_vel.end() && rms != global_rms_vel.end();
           ++p, ++a, ++rms, ++index)
        {
          if (this->convert_output_to_years() == true)
            {
              const std::string name_max = "Maximum velocity residual magnitude on boundary with indicator "
                                           + Utilities::int_to_string(p->first)
                                           + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                         .translate_id_to_symbol_name (p->first))
                                           + " (m/yr)";
              statistics.add_value (name_max, p->second*year_in_seconds);
              const std::string name_min = "Minimum velocity residual magnitude on boundary with indicator "
                                           + Utilities::int_to_string(a->first)
                                           + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                         .translate_id_to_symbol_name (a->first))
                                           + " (m/yr)";
              statistics.add_value (name_min, a->second*year_in_seconds);
              const std::string name_rms = "RMS velocity residual on boundary with indicator "
                                           + Utilities::int_to_string(rms->first)
                                           + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                         .translate_id_to_symbol_name (rms->first))
                                           + " (m/yr)";
              statistics.add_value (name_rms, rms->second*year_in_seconds);
              // also make sure that the other columns filled by this object
              // all show up with sufficient accuracy and in scientific notation
              statistics.set_precision (name_max, 8);
              statistics.set_scientific (name_max, true);
              statistics.set_precision (name_min, 8);
              statistics.set_scientific (name_min, true);
              statistics.set_precision (name_rms, 8);
              statistics.set_scientific (name_rms, true);
            }
          else
            {
              const std::string name_max = "Maximum velocity residual magnitude on boundary with indicator "
                                           + Utilities::int_to_string(p->first)
                                           + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                         .translate_id_to_symbol_name (p->first))
                                           + " (m/s)";
              statistics.add_value (name_max, p->second);
              const std::string name_min = "Minimum velocity residual magnitude on boundary with indicator "
                                           + Utilities::int_to_string(a->first)
                                           + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                         .translate_id_to_symbol_name (a->first))
                                           + " (m/s)";
              statistics.add_value (name_min, a->second);
              const std::string name_rms = "RMS velocity residual on boundary with indicator "
                                           + Utilities::int_to_string(rms->first)
                                           + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                         .translate_id_to_symbol_name (rms->first))
                                           + " (m/s)";
              statistics.add_value (name_rms, rms->second);
              // also make sure that the other columns filled by the this object
              // all show up with sufficient accuracy and in scientific notation
              statistics.set_precision (name_max, 8);
              statistics.set_scientific (name_max, true);
              statistics.set_precision (name_min, 8);
              statistics.set_scientific (name_min, true);
              statistics.set_precision (name_rms, 8);
              statistics.set_scientific (name_rms, true);
            }

          // finally have something for the screen
          screen_text.precision(4);
          if (this->convert_output_to_years() == true)
            {
              screen_text << p->second *year_in_seconds << " m/yr, "
                          << a->second *year_in_seconds << " m/yr, "
                          << rms->second *year_in_seconds << " m/yr"
                          << (index == global_max_vel.size()-1 ? "" : ", ");
            }
          else
            {
              screen_text << p->second << " m/s, "
                          << a->second << " m/s, "
                          << rms->second << " m/s"
                          << (index == global_max_vel.size()-1 ? "" : ", ");
            }

        }

      return std::pair<std::string, std::string> ("Max, min, and RMS residual velocity along boundary parts:",
                                                  screen_text.str());
    }



    template <int dim>
    void
    BoundaryVelocityResidualStatistics<dim>::declare_parameters (ParameterHandler  &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Boundary velocity residual statistics");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the GPlates model or the "
                             "ascii data. This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name", "current_day.gpml",
                             Patterns::Anything (),
                             "The file name of the input velocity as a GPlates model or an ascii data. "
                             "For the GPlates model, provide file in the same format as described "
                             "in the 'gplates' boundary velocity plugin. "
                             "For the ascii data, provide file in the same format as described in "
                             " 'ascii data' initial composition plugin." );
          prm.declare_entry ("Scale factor", "1.",
                             Patterns::Double (),
                             "Scalar factor, which is applied to the model data. "
                             "You might want to use this to scale the input to a "
                             "reference model. Another way to use this factor is to "
                             "convert units of the input files. For instance, if you "
                             "provide velocities in cm/yr set this factor to 0.01.");
          prm.declare_entry ("Use spherical unit vectors", "false",
                             Patterns::Bool (),
                             "Specify velocity as r, phi, and theta components "
                             "instead of x, y, and z. Positive velocities point up, east, "
                             "and north (in 3D) or out and clockwise (in 2D). "
                             "This setting only makes sense for spherical geometries."
                             "GPlates data is always interpreted to be in east and north directions "
                             "and is not affected by this parameter.");
          prm.declare_entry ("Use ascii data", "false",
                             Patterns::Bool (),
                             "Use ascii data files (e.g., GPS) for computing residual velocities "
                             "instead of GPlates data.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    BoundaryVelocityResidualStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Boundary velocity residual statistics");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          data_file_name    = prm.get ("Data file name");
          scale_factor      = prm.get_double ("Scale factor");

          use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
          use_ascii_data = prm.get_bool("Use ascii data");

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
    ASPECT_REGISTER_POSTPROCESSOR(BoundaryVelocityResidualStatistics,
                                  "boundary velocity residual statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the velocity residual along the top boundary. The velocity residual "
                                  "is the difference between the model solution velocities and the input "
                                  "velocities (GPlates model or ascii data). Currently, the velocity residual "
                                  "statistics, i.e., min, max and the rms magnitude, is computed at the top "
                                  "suface.")
  }
}

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



#include <aspect/postprocess/velocity_boundary_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VelocityBoundaryStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula for the velocity.
      const Quadrature<dim-1> &quadrature_formula = this->introspection().face_quadratures.velocities;

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_values |
                                        update_JxW_values |
                                        update_quadrature_points);

      std::vector<Tensor<1,dim>> velocities (fe_face_values.n_quadrature_points);

      std::map<types::boundary_id, double> local_max_vel;
      std::map<types::boundary_id, double> local_min_vel;
      std::map<types::boundary_id, double> local_velocity_square_integral;
      std::map<types::boundary_id, double> local_boundary_area;

      const std::set<types::boundary_id>
      boundary_indicators
        = this->get_geometry_model().get_used_boundary_indicators ();
      for (const auto p : boundary_indicators)
        {
          local_max_vel[p] = std::numeric_limits<double>::lowest();
          local_min_vel[p] = std::numeric_limits<double>::max();
        }

      // for every surface face that is part of a geometry boundary
      // and that is owned by this processor,
      // compute the maximum, minimum, and squared*area velocity magnitude,
      // and the face area.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->face(f)->at_boundary())
              {
                fe_face_values.reinit (cell, f);
                // extract velocities
                fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                    velocities);
                // determine the max, min, and squared velocity on the face
                // also determine the face area
                double local_max = std::numeric_limits<double>::lowest();
                double local_min = std::numeric_limits<double>::max();
                double local_sqvel = 0.0;
                double local_fe_face_area = 0.0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const double vel_mag = velocities[q].norm();
                    local_max = std::max(vel_mag,
                                         local_max);
                    local_min = std::min(vel_mag,
                                         local_min);
                    local_sqvel += ((vel_mag * vel_mag) * fe_face_values.JxW(q));
                    local_fe_face_area += fe_face_values.JxW(q);
                  }

                // then merge them with the min/max/squared velocities
                // and face areas we found for other faces with the same boundary indicator
                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();

                local_max_vel[boundary_indicator] = std::max(local_max,
                                                             local_max_vel[boundary_indicator]);
                local_min_vel[boundary_indicator] = std::min(local_min,
                                                             local_min_vel[boundary_indicator]);
                local_velocity_square_integral[boundary_indicator] += local_sqvel;
                local_boundary_area[boundary_indicator] += local_fe_face_area;
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

        for (const auto boundary_id : boundary_indicators)
          {
            local_max_values.push_back (local_max_vel[boundary_id]);
            local_min_values.push_back (local_min_vel[boundary_id]);

            local_velocity_square_integral_values.push_back (local_velocity_square_integral[boundary_id]);
            local_boundary_area_values.push_back (local_boundary_area[boundary_id]);
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
        for (const auto boundary_id: boundary_indicators)
          {
            global_max_vel[boundary_id] = global_max_values[index];
            global_min_vel[boundary_id] = global_min_values[index];
            global_rms_vel[boundary_id] = std::sqrt(global_velocity_square_integral_values[index] / global_boundary_area_values[index]);
            ++index;
          }
      }

      // now add the computed max, min, and rms velocities to the statistics object
      // and create a single string that can be output to the screen
      const std::string units = (this->convert_output_to_years() == true) ? "m/year" : "m/s";
      const double unit_scale_factor = (this->convert_output_to_years() == true) ? year_in_seconds : 1.0;
      std::ostringstream screen_text;
      unsigned int index = 0;

      for (std::map<types::boundary_id, double>::const_iterator
           max_velocity = global_max_vel.begin(), min_velocity = global_min_vel.begin(), rms = global_rms_vel.begin();
           max_velocity != global_max_vel.end() && min_velocity != global_min_vel.end() && rms != global_rms_vel.end();
           ++max_velocity, ++min_velocity, ++rms, ++index)
        {
          const std::string name_max = "Maximum velocity magnitude on boundary with indicator "
                                       + Utilities::int_to_string(max_velocity->first)
                                       + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                     .translate_id_to_symbol_name (max_velocity->first))
                                       + " (" + units + ")";
          statistics.add_value (name_max, max_velocity->second*unit_scale_factor);
          const std::string name_min = "Minimum velocity magnitude on boundary with indicator "
                                       + Utilities::int_to_string(min_velocity->first)
                                       + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                     .translate_id_to_symbol_name (min_velocity->first))
                                       + " (" + units + ")";
          statistics.add_value (name_min, min_velocity->second*unit_scale_factor);
          const std::string name_rms = "RMS velocity on boundary with indicator "
                                       + Utilities::int_to_string(rms->first)
                                       + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                     .translate_id_to_symbol_name (rms->first))
                                       + " (" + units + ")";
          statistics.add_value (name_rms, rms->second*unit_scale_factor);
          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name_max, 8);
          statistics.set_scientific (name_max, true);
          statistics.set_precision (name_min, 8);
          statistics.set_scientific (name_min, true);
          statistics.set_precision (name_rms, 8);
          statistics.set_scientific (name_rms, true);

          // finally have something for the screen
          screen_text.precision(4);
          screen_text << max_velocity->second *unit_scale_factor << ' ' << units << ", "
                      << min_velocity->second *unit_scale_factor << ' ' << units << ", "
                      << rms->second *unit_scale_factor << ' ' << units
                      << (index == global_max_vel.size()-1 ? "" : ", ");
        }

      return std::pair<std::string, std::string> ("Max, min, and rms velocity along boundary parts:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VelocityBoundaryStatistics,
                                  "velocity boundary statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the velocity along the boundaries. For each boundary "
                                  "indicator (see your geometry description for which boundary "
                                  "indicators are used), the min and max velocity magnitude is computed.")
  }
}

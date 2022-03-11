/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/postprocess/mobility_statistics.h>
#include <aspect/material_model/simple.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MobilityStatistics<dim>::execute (TableHandler &statistics)
    {
      // If this is the first time we get here, set the next output time
      // to the current time. This makes sure we always produce data during
      // the first time step.
      if (std::isnan(last_output_time))
        {  
          last_output_time = this->get_time() - output_interval;
          last_average_time = this->get_time() - average_interval;
         
          //Initialise combined_mobility
           combined_mobility = 0;
        }
      else
        {
          // Get a pointer to the combined_mobility
          const Postprocess::MobilityStatistics<dim> &mobility_statistics =
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::MobilityStatistics<dim>>();
          double combined_mobility = mobility_statistics.get_combined_mobility();
        }    
      
      // see if output is requested at this time      
      if (this->get_time() < last_output_time + output_interval)
        return std::pair<std::string,std::string>();

      // create a quadrature formula for the velocity for the volume of cells
      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<Tensor<1,dim>> velocity_values(n_q_points);

      // create a quadrature formula for the velocity for the surface of cells
      const QGauss<dim-1> quadrature_formula_face (this->introspection().polynomial_degree.velocities+1);
      const unsigned int n_q_points_face = quadrature_formula_face.size();

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_JxW_values |
                                        update_quadrature_points);

      std::vector<Tensor<1,dim>> surface_velocity_values (n_q_points_face);

      double local_velocity_square_integral = 0;
      double local_velocity_square_integral_top_boundary = 0;
      double local_top_boundary_area = 0;

      const std::set<types::boundary_id>
      boundary_indicators
        = this->get_geometry_model().get_used_boundary_indicators ();
      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            // extract velocities
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                local_velocity_square_integral += velocity_values[q].norm_square() *
                                                  fe_values.JxW(q);
              }

            // now for the surface cells only
            for (const unsigned int f : cell->face_indices())
              if (cell->face(f)->boundary_id() == top_boundary_id)
                {
                  fe_face_values.reinit (cell, f);
                  // extract velocities
                  fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                      surface_velocity_values);

                  // determine the squared velocity on the face
                  for (unsigned int q = 0; q < n_q_points_face; ++q)
                    {
                      const double JxW = fe_face_values.JxW(q);
                      local_velocity_square_integral_top_boundary += surface_velocity_values[q].norm_square() * JxW;
                      local_top_boundary_area += JxW;
                    }
                }
          }

      // now communicate to get the global values and convert to rms
      const double global_velocity_square_integral
        = Utilities::MPI::sum (local_velocity_square_integral, this->get_mpi_communicator());

      const double global_rms_vel = std::sqrt(global_velocity_square_integral) /
                                    std::sqrt(this->get_volume());

      const double global_top_velocity_square_integral = Utilities::MPI::sum(local_velocity_square_integral_top_boundary, this->get_mpi_communicator());
      const double global_top_boundary_area_integral = Utilities::MPI::sum(local_top_boundary_area, this->get_mpi_communicator());

      const double global_top_rms_vel = std::sqrt(global_top_velocity_square_integral) /
                                        std::sqrt(global_top_boundary_area_integral);

      // now add the computed rms velocities to the statistics object
      // and create a single string that can be output to the screen
      const double mobility  = global_top_rms_vel / global_rms_vel;
      
      //Add mobility every time to be averaged
      combined_mobility = combined_mobility + mobility;
      
      // see if averaging is requested at this time      
      if (this->get_time() >= last_average_time + average_interval)
        {  
          average_mobility = combined_mobility/(average_interval/output_interval);
          combined_mobility = 0;
        } 

      const std::string name_mobility = "Mobility";

      statistics.add_value (name_mobility,
                            mobility);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name_mobility, 8);
      statistics.set_scientific (name_mobility, true);

      std::ostringstream output;
      output.precision(3);
      output << mobility;

      // Update time
      set_last_output_time (this->get_time());
      set_last_average_time (this->get_time());
      
      return std::pair<std::string, std::string> ("Mobility:",
                                                  output.str());
    }


    template <int dim>
    void
    MobilityStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Mobility statistics");
        { 
          prm.declare_entry ("Time between mobility output", "0.",  
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "mobility output. A value of zero indicates "
                              "that output should be generated in each time step. "
                              "Units: years if the "
                              "'Use years in output instead of seconds' parameter is set; "
                              "seconds otherwise.");
          prm.declare_entry ("Time between averaging mobility", "0.",
                             Patterns::Double (0.),
                             "The time interval between averaging the "
                             "mobility output. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MobilityStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Mobility statistics");
        {
          output_interval = prm.get_double ("Time between mobility output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
          average_interval = prm.get_double ("Time between averaging mobility");
          if (this->convert_output_to_years())
            average_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    template <class Archive>
    void MobilityStatistics<dim>::serialize (Archive &ar, const unsigned int)
    {
      & last_output_time;
      & last_average_time;
    }

    
    template <int dim>
    void
    MobilityStatistics<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.      
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }   

    template <int dim>
    void
    MobilityStatistics<dim>::set_last_average_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.       if (output_interval > 0)
      if (average_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done           // by forcing std::floor to round 1.0-eps to 1.0.      
          // by forcing std::floor to round 1.0-eps to 1.0.      
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_average_time = last_average_time + std::floor((current_time-last_average_time)/average_interval*magic) * average_interval/magic;
        }
    }
    
    template <int dim>
    double
    MobilityStatistics<dim>::get_combined_mobility() const
    {
      //Elodie Feb 2022
      return combined_mobility;
    }

    template <int dim>
    double
    MobilityStatistics<dim>::get_average_mobility() const
    {
      //Elodie Feb 2022
      return average_mobility;
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MobilityStatistics,
                                  "mobility statistics",
                                  "A postprocessor that computes some statistics about mobility "
                                  "following Tackley (2000) and Lourenco et al. (2020).")
  }
}

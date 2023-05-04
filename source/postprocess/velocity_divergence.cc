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



#include <aspect/postprocess/velocity_divergence.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/tensor.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VelocityDivergence<dim>::execute (TableHandler &statistics)
    {
      Assert (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()) ||
              Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()),
              ExcMessage ("This postprocessor can only be used if the geometry "
                          "is a sphere or spherical shell."));

      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int n_q_points = quadrature_formula.size(); 

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values |
                               update_gradients);
      std::vector<Tensor<1,dim>> velocity_values(n_q_points);
      std::vector<Tensor<2,dim>> velocity_gradients(n_q_points);

     // Create vector of flags for each cell to indicate whether it is a subduction or rift zone
      std::vector<bool> is_subduction(this->get_dof_handler().n_locally_owned_dofs(), false);
     // Create vector of flags for each cell to indicate whether it is a subduction or rift zone
      std::vector<bool> is_rift(this->get_dof_handler().n_locally_owned_dofs(), false);


      double local_rad_velocity_square_integral = 0;
      double local_tan_velocity_square_integral = 0;
      double local_divergence_integral = 0;
      const double threshold_multiplier_subduction = 2;
      const double threshold_multiplier_rift = 2;
      const double radius_top = 6371000 - 100000;


      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {

        if (cell->is_locally_owned())
          {

            // bool subduction = false;
            // bool rifts = false;       
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            const std::vector<Point<dim>> &position_point = fe_values.get_quadrature_points();
            
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double radius = position_point[q].norm();
              if (radius < radius_top)
              {
                

                  // create unit vector in radial direction
                  const Tensor<1,dim> radial_unit_vector = position_point[q] / position_point[q].norm();

                  // compute the radial velocity by multiplying with the radial unit vector
                  const double radial_vel = (velocity_values[q] * radial_unit_vector);
                  local_rad_velocity_square_integral += (radial_vel * radial_vel) * fe_values.JxW(q);

                  // compute the tangential velocity by subtracting the radial velocity from the velocity
                  const Tensor<1,dim> tangential_vel = velocity_values[q] -  radial_vel * radial_unit_vector;
                  local_tan_velocity_square_integral += (tangential_vel * tangential_vel) * fe_values.JxW(q);
              }
            }
          }
        }


      // compute the global sums
      const double global_rad_velocity_square_integral
        = Utilities::MPI::sum (local_rad_velocity_square_integral, this->get_mpi_communicator());
      const double global_tan_velocity_square_integral
        = Utilities::MPI::sum (local_tan_velocity_square_integral, this->get_mpi_communicator());

      // compute the final output by dividing by the volume over which
      // we integrated and taking the sqrt
      const double rad_vrms = std::sqrt(global_rad_velocity_square_integral / this->get_volume());
      const double tan_vrms = std::sqrt(global_tan_velocity_square_integral / this->get_volume());
      const double vrms = std::sqrt(rad_vrms*rad_vrms + tan_vrms*tan_vrms);



      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        // Get cell index
        const unsigned int cell_index = cell->active_cell_index();

        if (cell->is_locally_owned())
          {

            // bool subduction = false;
            // bool rifts = false;       
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_solution(),
                                                                                         velocity_gradients);
            const std::vector<Point<dim>> &position_point = fe_values.get_quadrature_points();
            

            for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double radius = position_point[q].norm();
              if (radius > radius_top)
              {
                  // compute the RMS velocity
                  // const double local_vrms_cell = std::sqrt(local_rad_velocity_square_integral_cell + local_tan_velocity_square_integral_cell) / n_q_points;    

                  Tensor<2,dim> grad_vel = velocity_gradients[q];
                  double divergence = 0.0;
                  for (unsigned int d = 0; d < dim; ++d)
                      divergence += grad_vel[d][d];
                  local_divergence_integral += divergence * fe_values.JxW(q);
              

                // Check if divergence is less than or greater than threshold
                // Check if divergence is less than or greater than threshold and set flag
              if (divergence < -threshold_multiplier_subduction* vrms)
                  is_subduction[cell_index] = true;
              else if (divergence > threshold_multiplier_rift* vrms)
                  is_rift[cell_index] = true;   
              }             

            }
          }

          // If flag is set for current cell, set flag for all neighboring cells to true
          if (is_subduction[cell_index])
          {
              for (unsigned int face_index = 0; face_index < GeometryInfo<dim>::faces_per_cell; ++face_index)
              {              
                  const auto &neighbor = cell->neighbor(face_index);
                  if (neighbor->is_locally_owned())
                  {
                    const unsigned int neighbor_index = neighbor->active_cell_index();
                    is_subduction[neighbor_index] = true;
                  }  
              }
          }

          // If flag is set for current cell, set flag for all neighboring cells to true
          if (is_rift[cell_index])
          {
              for (unsigned int face_index = 0; face_index < GeometryInfo<dim>::faces_per_cell; ++face_index)
              {              
                  const auto &neighbor = cell->neighbor(face_index);
                  if (neighbor->is_locally_owned())
                  {
                    const unsigned int neighbor_index = neighbor->active_cell_index();
                    is_rift[neighbor_index] = true;
                  }  
              }
          }           

      }

      const double global_divergence_integral = Utilities::MPI::sum(local_divergence_integral, this->get_mpi_communicator());
      const double average_divergence = global_divergence_integral / this->get_volume();

      // Count number of subduction zones 
      unsigned int num_subduction_zones = 0;
      unsigned int num_rifts = 0;
      for (const auto &flag : is_subduction)
      {
        if (flag)
          flag > 0 ? ++num_subduction_zones : 0;
      }      

      // Count number of rifts
      for (const auto &flag : is_rift)
      {
        if (flag)
          flag > 0 ? ++num_rifts : 0;
      }      




      if (this->convert_output_to_years() == true)
        {
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const char *columns[] = { "Radial RMS velocity (m/yr)",
                                    "Tangential RMS velocity (m/yr)",
                                    "Total RMS velocity (m/yr)",
                                    "Number of subduction zones", 
                                    "Number of rifts",
                                    "Divergence"
                                  };
          statistics.add_value (columns[0],
                                rad_vrms * year_in_seconds);
          statistics.add_value (columns[1],
                                tan_vrms * year_in_seconds);
          statistics.add_value (columns[2],
                                vrms * year_in_seconds);
          statistics.add_value(columns[3], 
                                num_subduction_zones);
          statistics.add_value(columns[4], 
                                num_rifts);   
          statistics.add_value(columns[5], 
                                average_divergence);                                                             
          for (auto &column : columns)
            {
              statistics.set_precision (column, 8);
              statistics.set_scientific (column, true);
            }
        }
      else
        {
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const char *columns[] = { "Radial RMS velocity (m/s)",
                                    "Tangential RMS velocity (m/s)",
                                    "Total RMS velocity (m/s)",
                                    "Number of subduction zones", 
                                    "Number of rifts",
                                    "Divergence"                                    
                                  };
          statistics.add_value (columns[0], rad_vrms);
          statistics.add_value (columns[1], tan_vrms);
          statistics.add_value (columns[2], vrms);
          statistics.add_value(columns[3], num_subduction_zones);
          statistics.add_value(columns[4], num_rifts); 
          statistics.add_value(columns[5], average_divergence);                
          for (auto &column : columns)
            {
              statistics.set_precision (column, 8);
              statistics.set_scientific (column, true);
            }
        }

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << rad_vrms *year_in_seconds
               << " m/yr, "
               << tan_vrms *year_in_seconds
               << " m/yr, "
               << vrms *year_in_seconds
               << " m/yr"
               << num_subduction_zones 
               << " subduction zone(s) "
               << num_rifts 
               << " rift(s)"
               << average_divergence
               << " divergence 1/m";
      else
        output << rad_vrms
               << " m/s, "
               << tan_vrms
               << " m/s, "
               << vrms
               << " m/s"
               << num_subduction_zones 
               << " subduction zone(s) "
               << num_rifts 
               << " rift(s)"
               << average_divergence
               << " divergence 1/m";


      return std::pair<std::string, std::string> ("Radial RMS, tangential RMS, total RMS velocity, subduction zones, rifts, divergence:",
                                                  output.str());               

      // return std::pair<std::string, std::string> ("Radial RMS, tangential RMS, total RMS velocity:",
      //                                             output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VelocityDivergence,
                                  "spherical velocity divergence",
                                  "A postprocessor that computes radial, tangential and total RMS "
                                  "velocity.")
  }
}

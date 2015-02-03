/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/



#include <aspect/postprocess/spherical_velocity_statistics.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    SphericalVelocityStatistics<dim>::execute (TableHandler &statistics)
    {
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0 ||
              dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != 0,
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
                               update_JxW_values);
      std::vector<Tensor<1,dim> > velocity_values(n_q_points);

      double local_rad_velocity_square_integral = 0;
      double local_tan_velocity_square_integral = 0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            const std::vector<Point<dim> > &position_point = fe_values.get_quadrature_points();
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                // create unit vector in radial direction
                const Tensor<1,dim> radial_unit_vector = position_point[q] / position_point[q].norm();

                // compute the radial velocity by multiplying with the radial unit vector
                const double radial_vel = (velocity_values[q] * radial_unit_vector);
                local_rad_velocity_square_integral += (radial_vel * radial_vel) *
                                                      fe_values.JxW(q);
                // compute the tangential velocity by subtracting the radial velocity from the velocity
                const Tensor<1,dim> tangential_vel = velocity_values[q] -  radial_vel * radial_unit_vector;
                local_tan_velocity_square_integral += (tangential_vel * tangential_vel) *
                                                      fe_values.JxW(q);
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

      if (this->convert_output_to_years() == true)
        {
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const char *columns[] = { "Radial RMS velocity (m/yr)",
                                    "Tangential RMS velocity (m/yr)",
                                    "Total RMS velocity (m/yr)"
                                  };
          statistics.add_value (columns[0],
                                rad_vrms * year_in_seconds);
          statistics.add_value (columns[1],
                                tan_vrms * year_in_seconds);
          statistics.add_value (columns[2],
                                vrms * year_in_seconds);
          for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
            {
              statistics.set_precision (columns[i], 8);
              statistics.set_scientific (columns[i], true);
            }
        }
      else
        {
          // make sure that the columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          const char *columns[] = { "Radial RMS velocity (m/s)",
                                    "Tangential RMS velocity (m/s)",
                                    "Total RMS velocity (m/s)"
                                  };
          statistics.add_value (columns[0], rad_vrms);
          statistics.add_value (columns[1], tan_vrms);
          statistics.add_value (columns[2], vrms);
          for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
            {
              statistics.set_precision (columns[i], 8);
              statistics.set_scientific (columns[i], true);
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
               << " m/yr";
      else
        output << rad_vrms
               << " m/s, "
               << tan_vrms
               << " m/s, "
               << vrms
               << " m/s";

      return std::pair<std::string, std::string> ("Radial RMS, tangential RMS, total RMS velocity:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(SphericalVelocityStatistics,
                                  "spherical velocity statistics",
                                  "A postprocessor that computes radial, tangential and total RMS "
                                  "velocity.")
  }
}

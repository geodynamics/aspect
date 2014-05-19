/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: spherical_velocity_statistics.cc 1969 2013-10-16 19:32:59Z heien $  */


#include <aspect/postprocess/spherical_velocity_statistics.h>
#include <aspect/geometry_model/spherical_shell.h>
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
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0,
              ExcMessage ("This postprocessor can only be used if the geometry "
                          "is a spherical shell."));

      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(0).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<Tensor<1,dim> > velocity_values(n_q_points);
      std::vector<Point<dim> > position_point(n_q_points);
//      std::vector<double > position;

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
            position_point = fe_values.get_quadrature_points();
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                Tensor<1,dim> position;
                for (unsigned int i = 0; i<dim; ++i)
                position[i] = position_point[q][i];
                double length = (std::sqrt(position*position));
                position /= length;
                // compute the radial velocity by multiplying with the unit vector in the radial direction twice
                Tensor<1,dim> vel = (velocity_values[q] * position) * position;
                local_rad_velocity_square_integral += (vel * vel) *
                                                   fe_values.JxW(q);
                // compute the tangential velocity by subtractin the radial velocity from the velocity
                vel /= -1;
                vel += velocity_values[q];
                local_tan_velocity_square_integral += (vel * vel) *
                                                   fe_values.JxW(q);
              }
          }

      const double global_rad_velocity_square_integral
        = Utilities::MPI::sum (local_rad_velocity_square_integral, this->get_mpi_communicator());
      const double global_tan_velocity_square_integral
        = Utilities::MPI::sum (local_tan_velocity_square_integral, this->get_mpi_communicator());

      const double rad_vrms = std::sqrt(global_rad_velocity_square_integral) /
                          std::sqrt(this->get_volume());
      const double tan_vrms = std::sqrt(global_tan_velocity_square_integral) /
                          std::sqrt(this->get_volume());
      const double vrms = std::sqrt(rad_vrms*rad_vrms + tan_vrms*tan_vrms);

      if (this->convert_output_to_years() == true)
        {
          statistics.add_value ("Radial RMS velocity (m/year)",
                                rad_vrms * year_in_seconds);
          statistics.add_value ("Tangential RMS velocity (m/year)",
                                tan_vrms * year_in_seconds);
          statistics.add_value ("Total RMS velocity (m/year)",
                                vrms * year_in_seconds);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "Radial RMS velocity (m/year)",
                                      "Tangential RMS velocity (m/year)",
                                      "Total RMS velocity (m/year)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }
      else
        {
          statistics.add_value ("Radial RMS velocity (m/s)", rad_vrms);
          statistics.add_value ("Tangential RMS velocity (m/s)", tan_vrms);
          statistics.add_value ("Total RMS velocity (m/s)", vrms);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "Radial RMS velocity (m/s)",
                                      "Tangential velocity (m/s)"
                                      "Total RMS velocity (m/s)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << rad_vrms *year_in_seconds
               << " m/year, "
               << tan_vrms *year_in_seconds
               << " m/year, "
               << vrms *year_in_seconds
               << " m/year";
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

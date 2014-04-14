/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/postprocess/velocity_statistics.h>
#include <aspect/material_model/simple.h>
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
    VelocityStatistics<dim>::execute (TableHandler &statistics)
    {
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

      double local_velocity_square_integral = 0;
      double local_max_velocity = 0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                local_velocity_square_integral += ((velocity_values[q] * velocity_values[q]) *
                                                   fe_values.JxW(q));
                local_max_velocity = std::max (std::sqrt(velocity_values[q]*velocity_values[q]),
                                               local_max_velocity);
              }
          }

      const double global_velocity_square_integral
        = Utilities::MPI::sum (local_velocity_square_integral, this->get_mpi_communicator());
      const double global_max_velocity
        = Utilities::MPI::max (local_max_velocity, this->get_mpi_communicator());

      const double vrms = std::sqrt(global_velocity_square_integral) /
                          std::sqrt(this->get_volume());

      if (this->convert_output_to_years() == true)
        {
          statistics.add_value ("RMS velocity (m/year)",
                                vrms * year_in_seconds);
          statistics.add_value ("Max. velocity (m/year)",
                                global_max_velocity * year_in_seconds);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "RMS velocity (m/year)",
                                      "Max. velocity (m/year)"
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
          statistics.add_value ("RMS velocity (m/s)", vrms);
          statistics.add_value ("Max. velocity (m/s)", global_max_velocity);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "RMS velocity (m/s)",
                                      "Max. velocity (m/s)"
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
        output << vrms *year_in_seconds
               << " m/year, "
               << global_max_velocity *year_in_seconds
               << " m/year";
      else
        output << vrms
               << " m/s, "
               << global_max_velocity
               << " m/s";

// TODO: This really doesn't belong here
      if (this->get_time() == 0e0)
        {
          if (dynamic_cast<const MaterialModel::Simple<dim> *>(&this->get_material_model()) != 0)
            {

              const MaterialModel::Simple<dim> &material_model
                = dynamic_cast<const MaterialModel::Simple<dim> &>(this->get_material_model());

              const double h = this->get_geometry_model().maximal_depth();

              const double dT = this->get_boundary_temperature().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
                                - this->get_boundary_temperature().minimal_temperature(this->get_fixed_temperature_boundary_indicators());
              // we do not compute the compositions but give the functions below the value 0.0 instead
              std::vector<double> composition_values(this->n_compositional_fields(),0.0);

              Point<dim> representative_point = this->get_geometry_model().representative_point(0);
              const double gravity = this->get_gravity_model().gravity_vector(representative_point).norm();
              const double Ra = material_model.reference_density()*
                                gravity*
                                material_model.reference_thermal_expansion_coefficient()*
                                dT*std::pow(h,3)/
                                (material_model.reference_thermal_diffusivity()*
                                 material_model.reference_viscosity());

              this->get_pcout()<<  std::endl;
              this->get_pcout()<< "     Reference density (kg/m^3):                    "
                               << material_model.reference_density()
                               << std::endl;
              this->get_pcout()<< "     Reference gravity (m/s^2):                     "
                               << gravity
                               << std::endl;
              this->get_pcout()<< "     Reference thermal expansion (1/K):             "
                               << material_model.reference_thermal_expansion_coefficient()
                               << std::endl;
              this->get_pcout()<< "     Temperature contrast across model domain (K): "
                               << dT
                               << std::endl;
              this->get_pcout()<< "     Model domain depth (m):                        "
                               << h
                               << std::endl;
              this->get_pcout()<< "     Reference thermal diffusivity (m^2/s):         "
                               << material_model.reference_thermal_diffusivity()
                               << std::endl;
              this->get_pcout()<< "     Reference viscosity (Pas):                     "
                               << material_model.reference_viscosity()
                               << std::endl;
              this->get_pcout()<< "     Ra number:                                     "
                               << Ra
                               << std::endl;
              this->get_pcout()<< "     k_value:                                       "
                               << material_model.thermal_conductivity(dT, dT, composition_values, representative_point) //TODO: dT for the pressure is wrong
                               << std::endl;
              this->get_pcout()<< "     reference_cp:                                  "
                               << material_model.reference_cp()
                               << std::endl;
              this->get_pcout()<< "     reference_thermal_diffusivity:                 "
                               << material_model.thermal_conductivity(dT, dT, composition_values, representative_point)/(material_model.reference_density()*material_model.reference_cp()) //TODO: dT for the pressure is wrong
                               << std::endl;
              this->get_pcout()<<  std::endl;
            }
        }
      return std::pair<std::string, std::string> ("RMS, max velocity:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VelocityStatistics,
                                  "velocity statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "velocity field.")
  }
}

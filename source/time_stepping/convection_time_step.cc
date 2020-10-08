/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/time_stepping/convection_time_step.h>

namespace aspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    ConvectionTimeStep<dim>::execute()
    {
      const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                               this->get_parameters().stokes_velocity_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values);

      const unsigned int n_q_points = quadrature_formula.size();

      std::vector<Tensor<1,dim> > velocity_values(n_q_points);
      std::vector<Tensor<1,dim> > fluid_velocity_values(n_q_points);

      double max_local_speed_over_meshsize = 0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);

            double max_local_velocity = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
              max_local_velocity = std::max (max_local_velocity,
                                             velocity_values[q].norm());

            if (this->get_parameters().include_melt_transport)
              {
                const FEValuesExtractors::Vector ex_u_f = this->introspection().variable("fluid velocity").extractor_vector();
                fe_values[ex_u_f].get_function_values (this->get_solution(), fluid_velocity_values);

                for (unsigned int q=0; q<n_q_points; ++q)
                  max_local_velocity = std::max (max_local_velocity,
                                                 fluid_velocity_values[q].norm());
              }

            max_local_speed_over_meshsize = std::max(max_local_speed_over_meshsize,
                                                     max_local_velocity
                                                     /
                                                     cell->minimum_vertex_distance());

          }

      const double max_global_speed_over_meshsize
        = Utilities::MPI::max (max_local_speed_over_meshsize, this->get_mpi_communicator());

      double min_convection_timestep = std::numeric_limits<double>::max();

      if (max_global_speed_over_meshsize != 0.0)
        min_convection_timestep = this->get_parameters().CFL_number / (this->get_parameters().temperature_degree * max_global_speed_over_meshsize);

      AssertThrow (min_convection_timestep > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(min_convection_timestep) + ". "
                              "Please check for non-positive material properties."));

      return min_convection_timestep;
    }


  }
}

// explicit instantiations
namespace aspect
{
  namespace TimeStepping
  {
    ASPECT_REGISTER_TIME_STEPPING_MODEL(ConvectionTimeStep,
                                        "convection time step",
                                        "This model computes the convection time step as the minimum "
                                        "over all cells of $ CFL h^2 \\cdot \\rho C_p / k$, "
                                        "where k is the thermal conductivity. This plugin will always "
                                        "request advancing to the next time step.")
  }
}

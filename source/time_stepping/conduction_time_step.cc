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
#include <aspect/time_stepping/conduction_time_step.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    ConductionTimeStep<dim>::execute()
    {
      double min_local_conduction_timestep = std::numeric_limits<double>::max();

      const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                               this->get_parameters().stokes_velocity_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();

      MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                 this->introspection().n_compositional_fields);
      MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                   this->introspection().n_compositional_fields);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            in.reinit(fe_values,
                      cell,
                      this->introspection(),
                      this->get_solution());

            this->get_material_model().evaluate(in, out);

            if (this->get_parameters().formulation_temperature_equation
                == Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
              {
                // Overwrite the density by the reference density coming from the
                // adiabatic conditions as required by the formulation
                for (unsigned int q=0; q<n_q_points; ++q)
                  out.densities[q] = this->get_adiabatic_conditions().density(in.position[q]);
              }
            else if (this->get_parameters().formulation_temperature_equation
                     == Parameters<dim>::Formulation::TemperatureEquation::real_density)
              {
                // use real density
              }
            else
              AssertThrow(false, ExcNotImplemented());

            // Evaluate thermal diffusivity at each quadrature point and
            // calculate the corresponding conduction timestep, if applicable
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double k = out.thermal_conductivities[q];
                const double rho = out.densities[q];
                const double c_p = out.specific_heat[q];

                Assert(rho * c_p > 0,
                       ExcMessage ("The product of density and c_P needs to be a "
                                   "non-negative quantity."));

                const double thermal_diffusivity = k/(rho*c_p);

                if (thermal_diffusivity > 0)
                  {
                    min_local_conduction_timestep = std::min(min_local_conduction_timestep,
                                                             this->get_parameters().CFL_number*pow(cell->minimum_vertex_distance(),2.)
                                                             / thermal_diffusivity);
                  }
              }
          }

      const double min_conduction_timestep = Utilities::MPI::min (min_local_conduction_timestep, this->get_mpi_communicator());

      AssertThrow (min_conduction_timestep > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(min_conduction_timestep) + ". "
                              "Please check for non-positive material properties."));

      return min_conduction_timestep;
    }


  }
}

// explicit instantiations
namespace aspect
{
  namespace TimeStepping
  {
    ASPECT_REGISTER_TIME_STEPPING_MODEL(ConductionTimeStep,
                                        "conduction time step",
                                        "This model computes the conduction time step as the minimum "
                                        "over all cells of $ CFL h^2 \\cdot \\rho C_p / k$, "
                                        "where k is the thermal conductivity. This plugin will always "
                                        "request advancing to the next time step.")
  }
}

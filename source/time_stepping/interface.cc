/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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


#include <aspect/time_stepping/interface.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace TimeStepping
  {

    template <int dim>
    double
    Manager<dim>::
    compute_time_step_size() const
    {
      const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                               this->get_parameters().stokes_velocity_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               ((this->get_parameters().use_conduction_timestep || this->get_parameters().include_melt_transport)
                                ?
                                update_quadrature_points
                                :
                                update_default));

      const unsigned int n_q_points = quadrature_formula.size();


      std::vector<Tensor<1,dim> > velocity_values(n_q_points);
      std::vector<Tensor<1,dim> > fluid_velocity_values(n_q_points);
      std::vector<std::vector<double> > composition_values (this->introspection().n_compositional_fields,std::vector<double> (n_q_points));

      double max_local_speed_over_meshsize = 0;
      double min_local_conduction_timestep = std::numeric_limits<double>::max();

      MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                 this->introspection().n_compositional_fields);
      MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                   this->introspection().n_compositional_fields);

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

            if (this->get_parameters().use_conduction_timestep)
              {
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
          }

      const double max_global_speed_over_meshsize
        = Utilities::MPI::max (max_local_speed_over_meshsize, this->get_mpi_communicator());

      double min_convection_timestep = std::numeric_limits<double>::max();
      double min_conduction_timestep = std::numeric_limits<double>::max();

      if (max_global_speed_over_meshsize != 0.0)
        min_convection_timestep = this->get_parameters().CFL_number / (this->get_parameters().temperature_degree * max_global_speed_over_meshsize);

      if (this->get_parameters().use_conduction_timestep)
        min_conduction_timestep = - Utilities::MPI::max (-min_local_conduction_timestep, this->get_mpi_communicator());

      double new_time_step = std::min(min_convection_timestep,
                                      min_conduction_timestep);

      AssertThrow (new_time_step > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(new_time_step) + ". "
                              "Please check for non-positive material properties."));

      // Make sure we do not exceed the maximum time step length. This can happen
      // if velocities get too small or even zero in models that are stably stratified
      // or use prescribed velocities.
      new_time_step = std::min(new_time_step, this->get_parameters().maximum_time_step);

      // Make sure that the time step doesn't increase too fast
      if (this->get_timestep() != 0)
        new_time_step = std::min(new_time_step, this->get_timestep() + this->get_timestep() * this->get_parameters().maximum_relative_increase_time_step);

      // Make sure we do not exceed the maximum length for the first time step
      if (this->get_timestep_number() == 0)
        new_time_step = std::min(new_time_step, this->get_parameters().maximum_first_time_step);

      // Make sure we reduce the time step length appropriately if we terminate after this step
      return termination_manager.check_for_last_time_step(new_time_step);
    }



    template <int dim>
    bool
    Manager<dim>::
    need_checkpoint_on_terminate() const
    {
      return do_checkpoint_on_terminate;
    }



    template <int dim>
    bool
    Manager<dim>::should_simulation_terminate_now() const
    {
      return termination_manager.execute();
    }



    template <int dim>
    void
    Manager<dim>::initialize_simulator (const Simulator<dim> &simulator_object)
    {
      SimulatorAccess<dim>::initialize_simulator(simulator_object);
      termination_manager.initialize_simulator(simulator_object);
    }



    template <int dim>
    void
    Manager<dim>:: declare_parameters (ParameterHandler &prm)
    {
      TerminationCriteria::Manager<dim>::declare_parameters(prm);
      prm.enter_subsection("Termination criteria");
      {
        // Whether to checkpoint the simulation right before termination
        prm.declare_entry("Checkpoint on termination", "false",
                          Patterns::Bool (),
                          "Whether to checkpoint the simulation right before termination.");
      }
      prm.leave_subsection();

    }



    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      termination_manager.parse_parameters(prm);
      prm.enter_subsection("Termination criteria");
      {
        do_checkpoint_on_terminate = prm.get_bool("Checkpoint on termination");
      }
      prm.leave_subsection();
    }

  }
}




// explicit instantiation of the functions we implement in this file
namespace aspect
{

  namespace TimeStepping
  {
#define INSTANTIATE(dim) \
  \
  template \
  class \
  Manager<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

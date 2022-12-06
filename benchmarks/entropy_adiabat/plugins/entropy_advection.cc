/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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

#include "entropy_advection.h"
#include "entropy_model.h"

#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <aspect/plugins.h>

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    EntropyAdvectionSystem<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                          internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;
      if (!advection_field.is_temperature()
          && introspection.name_for_compositional_index(advection_field.compositional_variable) != "entropy")
        return;

      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const bool   use_bdf2_scheme = (this->get_timestep_number() > 1);
      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();

      const double bdf2_factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                     (time_step + old_time_step)) : 1.0;

      const unsigned int solution_component = advection_field.component_index(introspection);
      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      scratch.finite_element_values[introspection.extractors.temperature].get_function_values (this->get_old_solution(),
          scratch.old_temperature_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          // precompute the values of shape functions and their gradients.
          // We only need to look up values of shape functions if they
          // belong to 'our' component. They are zero otherwise anyway.
          // Note that we later only look at the values that we do set here.
          for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
            {
              if (fe.system_to_component_index(i).first == solution_component)
                {
                  scratch.grad_phi_field[i_advection] = scratch.finite_element_values[solution_field].gradient (i,q);
                  scratch.phi_field[i_advection]      = scratch.finite_element_values[solution_field].value (i,q);
                  ++i_advection;
                }
              ++i;
            }

          const double JxW = scratch.finite_element_values.JxW(q);

          // solve the diffusion equation for the temperature
          if (advection_field.is_temperature())
            {
              const double density_c_P = scratch.material_model_outputs.densities[q] *
                                         scratch.material_model_outputs.specific_heat[q];

              const double field_term_for_rhs = scratch.old_field_values[q] * density_c_P;

              // do the actual assembly. note that we only need to loop over the advection
              // shape functions because these are the only contributions we compute here
              for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                {
                  data.local_rhs(i)
                  += field_term_for_rhs * scratch.phi_field[i] * JxW;

                  for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                    {
                      data.local_matrix(i,j)
                      += (
                           (time_step * scratch.material_model_outputs.thermal_conductivities[q]
                            * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                           + (scratch.phi_field[i] * scratch.phi_field[j]) * density_c_P
                         )
                         * JxW;
                    }
                }
            }
          else
            // solve the entropy equation
            {
              const double rho_T              =
                scratch.material_model_outputs.densities[q] *
                scratch.material_model_inputs.temperature[q];

              Assert (rho_T >= 0,
                      ExcMessage ("The product of density and temperature needs to be a "
                                  "non-negative quantity."));

              const double gamma =
                scratch.heating_model_outputs.heating_source_terms[q];

              const double field_term_for_rhs
                = (use_bdf2_scheme ?
                   (scratch.old_field_values[q] *
                    (1 + time_step/old_time_step)
                    -
                    scratch.old_old_field_values[q] *
                    (time_step * time_step) /
                    (old_time_step * (time_step + old_time_step)))
                   :
                   scratch.old_field_values[q])
                  *
                  (rho_T);

              Tensor<1,dim> current_u = scratch.current_velocity_values[q];
              // Subtract off the mesh velocity for ALE corrections if necessary
              if (this->get_parameters().mesh_deformation_enabled)
                current_u -= scratch.mesh_velocity_values[q];

              // We compute the amount of diffusion based on the solution of the temperature equation.
              const double diffusion_term = (scratch.material_model_inputs.temperature[q] - scratch.old_temperature_values[q])
                                            * scratch.material_model_outputs.densities[q] * scratch.material_model_outputs.specific_heat[q];

              // do the actual assembly. note that we only need to loop over the advection
              // shape functions because these are the only contributions we compute here
              for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                {
                  data.local_rhs(i)
                  += (field_term_for_rhs * scratch.phi_field[i]
                      + diffusion_term * scratch.phi_field[i]
                      + time_step * scratch.phi_field[i] * gamma)
                     *
                     JxW;

                  for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                    {
                      data.local_matrix(i,j)
                      += (
                           (time_step * scratch.artificial_viscosity
                            * scratch.grad_phi_field[i] * scratch.grad_phi_field[j])
                           + ((time_step * (scratch.phi_field[i] * (current_u * scratch.grad_phi_field[j])))
                              + (bdf2_factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                           (rho_T)
                         )
                         * JxW;
                    }
                }
            }
        }
    }



    template <int dim>
    std::vector<double>
    EntropyAdvectionSystem<dim>::compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      std::vector<double> residuals(n_q_points,0.0);

      if (advection_field.is_temperature() ||
          this->introspection().name_for_compositional_index(advection_field.compositional_variable) != "entropy")
        return residuals;

      this->get_heating_model_manager().evaluate(scratch.material_model_inputs,
                                                 scratch.material_model_outputs,
                                                 scratch.heating_model_outputs);

      for (unsigned int q=0; q < n_q_points; ++q)
        {
          const Tensor<1,dim> u = (scratch.old_velocity_values[q] +
                                   scratch.old_old_velocity_values[q]) / 2;

          const double dField_dt = (this->get_old_timestep() == 0.0) ? 0.0 :
                                   (
                                     ((scratch.old_field_values)[q] - (scratch.old_old_field_values)[q])
                                     / this->get_old_timestep());
          const double u_grad_field = u * (scratch.old_field_grads[q] +
                                           scratch.old_old_field_grads[q]) / 2;

          const double density      = scratch.material_model_outputs.densities[q];
          const double gamma        = scratch.heating_model_outputs.heating_source_terms[q];

          // Because we solve the diffusion equation for the temperature before the advection equation, we can use
          // the current and old temperature here, together with the current time step.
          const double diffusion_term = (this->get_timestep() == 0.0) ? 0.0
                                        :
                                        (scratch.material_model_inputs.temperature[q] - scratch.old_temperature_values[q]) / this->get_timestep()
                                        * scratch.material_model_outputs.densities[q] * scratch.material_model_outputs.specific_heat[q];


          residuals[q]
            = std::abs((density * scratch.material_model_inputs.temperature[q]) * (dField_dt + u_grad_field) - gamma - diffusion_term);
        }
      return residuals;
    }



    template <int dim>
    std::vector<double>
    EntropyAdvectionSystem<dim>::advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      std::vector<double> prefactors(scratch.material_model_inputs.n_evaluation_points(), 0.0);

      for (unsigned int i=0; i<prefactors.size(); ++i)
        prefactors[i] = scratch.material_model_outputs.densities[i] * scratch.material_model_inputs.temperature[i];

      return prefactors;
    }



    template <int dim>
    std::vector<double>
    EntropyAdvectionSystem<dim>::diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      std::vector<double> prefactors(scratch.material_model_inputs.n_evaluation_points(), 0.0);

      for (unsigned int i=0; i<prefactors.size(); ++i)
        prefactors[i] = scratch.material_model_outputs.thermal_conductivities[i] *
                        scratch.material_model_inputs.temperature[i] /
                        scratch.material_model_outputs.specific_heat[i];

      return prefactors;
    }
  }

// The function below replaces the normal temperature advection assembler with the entropy advection assembler of this file
  template <int dim>
  void set_assemblers_entropy_advection(const SimulatorAccess<dim> &simulator_access,
                                        Assemblers::Manager<dim> &assemblers)
  {
    AssertThrow (Plugins::plugin_type_matches<MaterialModel::EntropyModel<dim>>
                 (simulator_access.get_material_model()),
                 ExcMessage ("The entropy advection assembler can only be used with the "
                             "material model 'entropy model'!"));

    AssertThrow (simulator_access.get_heating_model_manager().adiabatic_heating_enabled() == false,
                 ExcMessage("The entropy advection assembler requires "
                            "that adiabatic heating is disabled."));

    // Replace all existing assemblers by the one for the entropy equation.
    assemblers.advection_system.resize(1);
    assemblers.advection_system[0] = std::make_unique<Assemblers::EntropyAdvectionSystem<dim>>();

    assemblers.advection_system_assembler_properties[0].needed_update_flags = update_hessians;
  }
} // namespace aspect

template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  signals.set_assemblers.connect (&aspect::set_assemblers_entropy_advection<dim>);
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)

// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class EntropyAdvectionSystem<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}

/*
  Copyright (C) 2016 - 2023 by the authors of the ASPECT code.

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

#include <aspect/simulator/assemblers/advection.h>

#include <aspect/melt.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Assemblers
  {
    namespace
    {
      /* These functions implement a reduced form of the code from deal.II's TriaAccessor::measure().
       * In the 3d dG case, a call to face->measure() is not implemented for non-planar faces.
       * Since we only care about the scaling here, it is enough to have an approximation instead.
       * The 2d case remains unchanged.
       */
      double
      approximate_face_measure(const DoFHandler<2>::face_iterator &face)
      {
        return (face->vertex(0)-face->vertex(1)).norm();
      }

      double
      approximate_face_measure(const DoFHandler<3>::face_iterator &face)
      {
        const Tensor<1,3> v03 = face->vertex(3) - face->vertex(0);
        const Tensor<1,3> v12 = face->vertex(2) - face->vertex(1);
        const Tensor<1,3> twice_area = cross_product_3d(v03, v12);
        return 0.5 * twice_area.norm();
      }
    }

    template <int dim>
    void
    AdvectionSystem<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                   internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const bool use_supg = (this->get_parameters().advection_stabilization_method
                             == Parameters<dim>::AdvectionStabilizationMethod::supg);
      const bool   use_bdf2_scheme = (this->get_timestep_number() > 1);
      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();
      const double bdf2_factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                     (time_step + old_time_step)) : 1.0;

      const bool advection_field_is_temperature = advection_field.is_temperature();
      const unsigned int solution_component = advection_field.component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      Assert(advection_field.advection_method(introspection)
             == Parameters<dim>::AdvectionFieldMethod::fem_field,
             ExcMessage("The 'AdvectionSystem' assembler can only be executed for fields "
                        "that use the advection method 'field'."));

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
                  if (use_supg)
                    scratch.laplacian_phi_field[i_advection] = trace(scratch.finite_element_values[solution_field].hessian (i,q));

                  scratch.grad_phi_field[i_advection] = scratch.finite_element_values[solution_field].gradient (i,q);
                  scratch.phi_field[i_advection]      = scratch.finite_element_values[solution_field].value (i,q);
                  ++i_advection;
                }
              ++i;
            }

          const double density_c_P              =
            ((advection_field_is_temperature)
             ?
             scratch.material_model_outputs.densities[q] *
             scratch.material_model_outputs.specific_heat[q]
             :
             1.0);

          AssertThrow (density_c_P >= 0,
                       ExcMessage ("The product of density and c_P needs to be a "
                                   "non-negative quantity."));

          const double latent_heat_LHS =
            ((advection_field_is_temperature)
             ?
             scratch.heating_model_outputs.lhs_latent_heat_terms[q]
             :
             0.0);
          AssertThrow (density_c_P + latent_heat_LHS >= 0,
                       ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                   "to the left hand side needs to be a non-negative quantity."));

          const double gamma =
            ((advection_field_is_temperature)
             ?
             scratch.heating_model_outputs.heating_source_terms[q]
             :
             0.0);

          const double reaction_term =
            ((advection_field_is_temperature)
             ?
             0.0
             :
             scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]);

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
              (density_c_P + latent_heat_LHS);

          Tensor<1,dim> current_u = scratch.current_velocity_values[q];
          // Subtract off the mesh velocity for ALE corrections if necessary
          if (this->get_parameters().mesh_deformation_enabled)
            current_u -= scratch.mesh_velocity_values[q];

          const double JxW = scratch.finite_element_values.JxW(q);

          // For the diffusion constant, use the larger of the physical
          // and the artificial viscosity/conductivity/diffusion constant.
          // One could also choose the sum of the two, but if the
          // physical diffusion is larger than the artificial one,
          // then (because the latter is chosen sufficiently large to
          // make the problem stable) one may as well stick with the
          // physical one. And if the physical diffusion is too small to
          // make the problem stable, then we ought to choose the smallest
          // diffusivity value that makes the problem stable -- which is
          // exactly the artificial viscosity.
          const double conductivity = (advection_field_is_temperature
                                       ?
                                       scratch.material_model_outputs.thermal_conductivities[q]
                                       :
                                       0.0);

          const double diffusion_constant = (use_supg)
                                            ?
                                            conductivity
                                            :
                                            std::max (conductivity, scratch.artificial_viscosity);

          const double tau = (use_supg) ? scratch.artificial_viscosity : 0.0;

          // do the actual assembly. note that we only need to loop over the advection
          // shape functions because these are the only contributions we compute here
          for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
            {
              data.local_rhs(i)
              += (field_term_for_rhs * scratch.phi_field[i]
                  + time_step *
                  scratch.phi_field[i]
                  * gamma
                  + scratch.phi_field[i]
                  * reaction_term)
                 *
                 JxW;

              if (use_supg)
                data.local_rhs(i)
                += tau *
                   (
                     (current_u * (density_c_P + latent_heat_LHS)) *
                     scratch.grad_phi_field[i] *
                     (
                       field_term_for_rhs
                       +
                       time_step * gamma
                       +
                       reaction_term
                     )
                   ) * JxW;


              for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                {
                  data.local_matrix(i,j)
                  += (
                       (time_step * diffusion_constant
                        * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                       + ((time_step * (scratch.phi_field[i] * (current_u * scratch.grad_phi_field[j])))
                          + (bdf2_factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                       (density_c_P + latent_heat_LHS)
                     )
                     * JxW;
                }

              if (use_supg)
                {
                  for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                    {
                      // Note that we assume that the conductivity is constant, otherwise we would need to
                      // compute div (kappa grad T), which we don't have access to.
                      data.local_matrix(i,j)
                      += tau *
                         (
                           (current_u * (density_c_P + latent_heat_LHS)) *
                           scratch.grad_phi_field[i] *
                           (
                             -time_step * conductivity * scratch.laplacian_phi_field[j]
                             +
                             (
                               (time_step * current_u * scratch.grad_phi_field[j])
                               +
                               (bdf2_factor * scratch.phi_field[j])
                             ) *
                             (density_c_P + latent_heat_LHS)
                           )
                         ) * JxW;
                    }
                }
            }
        }
    }



    template <int dim>
    std::vector<double>
    AdvectionSystem<dim>::compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      std::vector<double> residuals(n_q_points);

      if (advection_field.is_temperature())
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

          if (advection_field.is_temperature())
            {
              const double density       = scratch.material_model_outputs.densities[q];
              const double conductivity  = scratch.material_model_outputs.thermal_conductivities[q];
              const double c_P           = scratch.material_model_outputs.specific_heat[q];
              // Note that we assume that the conductivity is constant, otherwise we would need to
              // compute div (kappa grad T), which we don't have access to.
              const double k_Delta_field = conductivity * (scratch.old_field_laplacians[q] +
                                                           scratch.old_old_field_laplacians[q]) / 2;

              const double gamma           = scratch.heating_model_outputs.heating_source_terms[q];
              const double latent_heat_LHS = scratch.heating_model_outputs.lhs_latent_heat_terms[q];

              residuals[q]
                = std::abs((density * c_P + latent_heat_LHS) * (dField_dt + u_grad_field) - k_Delta_field - gamma);
            }
          else
            {
              const double dreaction_term_dt = (this->get_old_timestep() == 0) ? 0.0 :
                                               scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]
                                               / this->get_old_timestep();

              residuals[q] = std::abs(dField_dt + u_grad_field - dreaction_term_dt);
            }
        }
      return residuals;
    }



    template <int dim>
    void
    DiffusionSystem<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                   internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>&> (data_base);

      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;

      Assert(advection_field.advection_method(introspection)
             == Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion,
             ExcMessage("The 'DiffusionSystem' assembler can only be executed for fields "
                        "that use the advection method 'prescribed field with diffusion'."));

      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const unsigned int solution_component = advection_field.component_index(introspection);
      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

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

          // do the actual assembly. note that we only need to loop over the advection
          // shape functions because these are the only contributions we compute here
          for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
            {
              data.local_rhs(i)
              += scratch.old_field_values[q] * scratch.phi_field[i]
                 *
                 JxW;

              for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                {
                  data.local_matrix(i,j)
                  += (parameters.diffusion_length_scale * parameters.diffusion_length_scale *
                      (scratch.grad_phi_field[i] * scratch.grad_phi_field[j])
                      + (scratch.phi_field[i] * scratch.phi_field[j])
                     )
                     * JxW;
                }
            }
        }
    }



    template <int dim>
    std::vector<double>
    DiffusionSystem<dim>::compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      std::vector<double> residuals(scratch.finite_element_values.n_quadrature_points, 0.0);

      return residuals;
    }



    template <int dim>
    std::vector<double>
    DiffusionSystem<dim>::advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      return std::vector<double> (scratch.material_model_inputs.n_evaluation_points(), 0.0);
    }



    template <int dim>
    std::vector<double>
    DiffusionSystem<dim>::diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      const double prefactor = this->get_timestep() > 0
                               ?
                               this->get_parameters().diffusion_length_scale
                               * this->get_parameters().diffusion_length_scale
                               / this->get_timestep()
                               :
                               0.0;

      return std::vector<double> (scratch.material_model_inputs.n_evaluation_points(), prefactor);
    }



    template <int dim>
    void
    AdvectionSystemBoundaryHeatFlux<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                                  internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;

      if (!advection_field.is_temperature())
        return;

      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(face_no);

      const unsigned int n_face_q_points    = scratch.face_finite_element_values->n_quadrature_points;
      const double time_step = this->get_timestep();

      // also have the number of dofs that correspond just to the element for
      // the system we are currently trying to assemble
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
      Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());

      const unsigned int solution_component = advection_field.component_index(introspection);
      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      if (this->get_fixed_heat_flux_boundary_indicators().find(face->boundary_id())
          != this->get_fixed_heat_flux_boundary_indicators().end())
        {
          // We are in the case of a Neumann temperature boundary.
          // Impose the Neumann value weakly using a RHS term.

          const std::vector<Tensor<1,dim>> heat_flux
            = this->get_boundary_heat_flux().heat_flux(
                face->boundary_id(),
                scratch.face_material_model_inputs,
                scratch.face_material_model_outputs,
                scratch.face_finite_element_values->get_normal_vectors());

          for (unsigned int q=0; q<n_face_q_points; ++q)
            {
              // precompute the values of shape functions.
              // We only need to look up values of shape functions if they
              // belong to 'our' component. They are zero otherwise anyway.
              // Note that we later only look at the values that we do set here.
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      scratch.face_phi_field[i_advection] = (*scratch.face_finite_element_values)[solution_field].value (i, q);
                      ++i_advection;
                    }
                  ++i;
                }

              for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                {
                  data.local_rhs(i)
                  -= time_step * scratch.face_phi_field[i] *
                     (heat_flux[q] * scratch.face_finite_element_values->normal_vector(q))
                     *
                     scratch.face_finite_element_values->JxW(q);
                }
            }
        }
    }



    template <int dim>
    void
    DarcySystem<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      MeltHandler<dim>::create_material_model_outputs(outputs);
    }

    template <int dim>
    void
    DarcySystem<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                               internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>& > (data_base);
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;

      Assert(advection_field.advection_method(introspection)
             == Parameters<dim>::AdvectionFieldMethod::fem_darcy_field,
             ExcMessage("The 'DarcySystem' assembler can only be executed for fields "
                        "that use the advection method 'darcy field'."));

      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const bool use_supg = (this->get_parameters().advection_stabilization_method
                             == Parameters<dim>::AdvectionStabilizationMethod::supg);
      AssertThrow(use_supg == false,
                  ExcMessage("The Darcy field advection method does not support the use of SUPG"));

      const bool   use_bdf2_scheme = (this->get_timestep_number() > 1);
      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();

      const double bdf2_factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                     (time_step + old_time_step)) : 1.0;
      const unsigned int solution_component = advection_field.component_index(introspection);
      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);
      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim>>();

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

          const double reaction_term = scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable];

          const double field_term_for_rhs
            = (use_bdf2_scheme ?
               (scratch.old_field_values[q] *
                (1 + time_step/old_time_step)
                -
                scratch.old_old_field_values[q] *
                (time_step * time_step) /
                (old_time_step * (time_step + old_time_step)))
               :
               scratch.old_field_values[q]);


          // We calculate the fluid velocity current_u_f using an approximation of Darcy's Law:
          // u_f = u_s - K_D / phi * (rho_s * g - rho_f * g)
          // u_f = fluid velocity
          // u_s = solid velocity
          // K_D = Darcy Coefficient
          // phi = porosity
          // rho_f = fluid density
          // rhos_s = solid density
          // g = gravity
          // The second term on the rhs of Darcy's Law only contributes to the component of the
          // fluid velocity parallel to the gravity vector, which is proportional to the
          // buoyancy of the fluid phase.
          const double rho_s = scratch.material_model_outputs.densities[q];
          const double rho_f = melt_outputs->fluid_densities[q];
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));
          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity         = std::max(scratch.material_model_inputs.composition[q][porosity_index],1e-10);
          const double K_D = melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q] / porosity;
          const Tensor<1,dim> current_u = scratch.current_velocity_values[q];
          Tensor<1,dim> current_u_f = current_u - K_D * (rho_s - rho_f) * gravity;

          // Subtract off the mesh velocity for ALE corrections if necessary
          if (this->get_parameters().mesh_deformation_enabled)
            current_u_f -= scratch.mesh_velocity_values[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          // do the actual assembly. note that we only need to loop over the advection
          // shape functions because these are the only contributions we compute here
          for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
            {
              data.local_rhs(i)
              += (field_term_for_rhs * scratch.phi_field[i]
                  + scratch.phi_field[i]
                  * reaction_term)
                 * JxW;

              for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                {
                  data.local_matrix(i,j)
                  += (
                       (time_step * scratch.artificial_viscosity
                        * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                       + ((time_step * (scratch.phi_field[i] * (current_u_f * scratch.grad_phi_field[j])))
                          + (bdf2_factor * scratch.phi_field[i] * scratch.phi_field[j]))
                     )
                     * JxW;
                }
            }
        }
    }



    template <int dim>
    std::vector<double>
    DarcySystem<dim>::compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>& > (scratch_base);
      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const Introspection<dim> &introspection = this->introspection();
      std::vector<double> residuals(n_q_points);
      if (advection_field.is_temperature())
        {
          return residuals;
        }


      const MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
      for (unsigned int q=0; q < n_q_points; ++q)
        {

          const double K_D = melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q];
          const double rho_s = scratch.material_model_outputs.densities[q];
          const double rho_f = melt_outputs->fluid_densities[q];
          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity         = std::max(scratch.material_model_inputs.composition[q][porosity_index], 1e-10);
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const Tensor<1,dim> u = (scratch.old_velocity_values[q] +
                                   scratch.old_old_velocity_values[q]) / 2
                                  - K_D*(rho_s - rho_f)*gravity/porosity;

          const double dField_dt = (this->get_old_timestep() == 0.0) ? 0.0 :
                                   (
                                     ((scratch.old_field_values)[q] - (scratch.old_old_field_values)[q])
                                     / this->get_old_timestep());
          const double u_grad_field = u * (scratch.old_field_grads[q] +
                                           scratch.old_old_field_grads[q]) / 2;

          const double dreaction_term_dt = (this->get_old_timestep() == 0) ? 0.0 :
                                           scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]
                                           / this->get_old_timestep();

          residuals[q] = std::abs(dField_dt + u_grad_field - dreaction_term_dt);
        }
      return residuals;
    }



    template <int dim>
    void
    AdvectionSystemBoundaryFace<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                              internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>&> (data_base);

      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;

      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(face_no);

      const unsigned int n_q_points    = scratch.face_finite_element_values->n_quadrature_points;

      const double time_step = this->get_timestep();

      Assert(advection_field.is_discontinuous(introspection)
             && advection_field.advection_method(introspection)
             == Parameters<dim>::AdvectionFieldMethod::fem_field,
             ExcMessage("The 'AdvectionSystemBoundaryFace' assembler can only be executed for fields "
                        "that use the advection method 'field' and a discontinuous discretization."));

      // also have the number of dofs that correspond just to the element for
      // the system we are currently trying to assemble
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
      Assert (scratch.face_grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
      Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());

      const unsigned int solution_component = advection_field.component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);


      if (((this->get_fixed_temperature_boundary_indicators().find(
              face->boundary_id()
            )
            != this->get_fixed_temperature_boundary_indicators().end())
           && (advection_field.is_temperature()))
          ||
          (( this->get_fixed_composition_boundary_indicators().find(
               face->boundary_id()
             )
             != this->get_fixed_composition_boundary_indicators().end())
           && (!advection_field.is_temperature())))
        {
          /*
           * We are in the case of a Dirichlet temperature or composition boundary.
           * In the temperature case, impose the Dirichlet value weakly using a matrix term
           * and RHS term. In the composition case, Dirichlet conditions can only be imposed
           * on inflow boundaries, and we only have the flow-dependent terms, so we only
           * assemble the corresponding flow-dependent and matrix and RHS terms
           * if we are on an inflow boundary.
           */

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              // precompute the values of shape functions and their gradients.
              // We only need to look up values of shape functions if they
              // belong to 'our' component. They are zero otherwise anyway.
              // Note that we later only look at the values that we do set here.
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      scratch.face_grad_phi_field[i_advection] = (*scratch.face_finite_element_values)[solution_field].gradient (i, q);
                      scratch.face_phi_field[i_advection]      = (*scratch.face_finite_element_values)[solution_field].value (i, q);
                      ++i_advection;
                    }
                  ++i;
                }

              const double density_c_P              =
                ((advection_field.is_temperature())
                 ?
                 scratch.face_material_model_outputs.densities[q] *
                 scratch.face_material_model_outputs.specific_heat[q]
                 :
                 1.0);

              AssertThrow (density_c_P >= 0,
                           ExcMessage ("The product of density and c_P needs to be a "
                                       "non-negative quantity."));

              const double conductivity =
                ((advection_field.is_temperature())
                 ?
                 scratch.face_material_model_outputs.thermal_conductivities[q]
                 :
                 0.0);
              const double latent_heat_LHS =
                ((advection_field.is_temperature())
                 ?
                 scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                 :
                 0.0);
              AssertThrow (density_c_P + latent_heat_LHS >= 0,
                           ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                       "to the left hand side needs to be a non-negative quantity."));

              const double penalty = (advection_field.is_temperature()
                                      ?
                                      parameters.discontinuous_penalty
                                      * parameters.temperature_degree
                                      * parameters.temperature_degree
                                      / approximate_face_measure(face)
                                      * conductivity
                                      / (density_c_P + latent_heat_LHS)
                                      :
                                      0.0);

              const double dirichlet_value = (advection_field.is_temperature()
                                              ?
                                              this->get_boundary_temperature_manager().boundary_temperature(
                                                face->boundary_id(),
                                                scratch.face_finite_element_values->quadrature_point(q))
                                              :
                                              this->get_boundary_composition_manager().boundary_composition(
                                                face->boundary_id(),
                                                scratch.face_finite_element_values->quadrature_point(q),
                                                advection_field.compositional_variable));

              Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
              // Subtract off the mesh velocity for ALE corrections if necessary
              if (parameters.mesh_deformation_enabled)
                current_u -= scratch.face_mesh_velocity_values[q];

              /**
               * The discontinuous Galerkin method uses 2 types of jumps over edges:
               * undirected and directed jumps. Undirected jumps are dependent only
               * on the order of the numbering of cells. Directed jumps are dependent
               * on the direction of the flow. Thus the flow-dependent terms below are
               * only calculated if the edge is an inflow edge.
               */
              const bool inflow = ((current_u * scratch.face_finite_element_values->normal_vector(q)) < 0.);

              for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                {
                  data.local_rhs(i)
                  += (- time_step *  conductivity
                      * scratch.face_grad_phi_field[i]
                      * scratch.face_finite_element_values->normal_vector(q)
                      * dirichlet_value

                      + time_step
                      * (density_c_P + latent_heat_LHS)
                      * penalty
                      * scratch.face_phi_field[i]
                      * dirichlet_value

                      + (inflow
                         ?
                         - (density_c_P + latent_heat_LHS)
                         * time_step
                         * (current_u
                            * scratch.face_finite_element_values->normal_vector(q))
                         * dirichlet_value
                         * scratch.face_phi_field[i]
                         :
                         0.)
                     )
                     *
                     scratch.face_finite_element_values->JxW(q);

                  // local_matrix terms
                  for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                    {
                      data.local_matrix(i,j)
                      += (- time_step *  conductivity
                          * scratch.face_grad_phi_field[i]
                          * scratch.face_finite_element_values->normal_vector(q)
                          * scratch.face_phi_field[j]

                          - time_step *  conductivity
                          * scratch.face_grad_phi_field[j]
                          * scratch.face_finite_element_values->normal_vector(q)
                          * scratch.face_phi_field[i]

                          + time_step
                          * (density_c_P + latent_heat_LHS)
                          * penalty
                          * scratch.face_phi_field[i]
                          * scratch.face_phi_field[j]

                          + (inflow
                             ?
                             - (density_c_P + latent_heat_LHS)
                             * time_step
                             * (current_u
                                * scratch.face_finite_element_values->normal_vector(q))
                             * scratch.face_phi_field[i]
                             * scratch.face_phi_field[j]
                             :
                             0.)
                         )
                         * scratch.face_finite_element_values->JxW(q);
                    }
                }
            }
        }
      else
        {
          // Neumann temperature term - no non-zero contribution as only homogeneous Neumann boundary conditions are implemented elsewhere for temperature
        }
    }



    template <int dim>
    void
    AdvectionSystemInteriorFace<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                              internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>&> (data_base);

      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename DoFHandler<dim>::active_cell_iterator cell = scratch.cell;
      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);

      const typename Simulator<dim>::AdvectionField advection_field = *scratch.advection_field;

      const unsigned int n_q_points    = scratch.face_finite_element_values->n_quadrature_points;

      const double time_step = this->get_timestep();

      if (!advection_field.is_discontinuous(introspection))
        return;

      // also have the number of dofs that correspond just to the element for
      // the system we are currently trying to assemble
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int dofs_per_cell = fe.dofs_per_cell;

      Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
      Assert (scratch.face_grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
      Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
      Assert (n_q_points == scratch.subface_finite_element_values->n_quadrature_points, ExcInternalError());
      Assert (n_q_points == scratch.neighbor_face_finite_element_values->n_quadrature_points, ExcInternalError());

      const unsigned int solution_component = advection_field.component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      // interior face or periodic face - no contribution on RHS

      const typename DoFHandler<dim>::cell_iterator
      neighbor = cell->neighbor_or_periodic_neighbor (face_no);
      // note: "neighbor" defined above is NOT active_cell_iterator, so this includes cells that are refined
      // for example: cell with periodic boundary.
      Assert (neighbor.state() == IteratorState::valid,
              ExcInternalError());
      const bool cell_has_periodic_neighbor = cell->has_periodic_neighbor (face_no);

      if (!neighbor->has_children())
        {
          if (neighbor->level () == cell->level () &&
              neighbor->is_active() &&
              (((neighbor->is_locally_owned()) && (cell->index() < neighbor->index()))
               ||
               ((!neighbor->is_locally_owned()) && (cell->subdomain_id() < neighbor->subdomain_id()))))
            {
              Assert (cell->is_locally_owned(), ExcInternalError());
              // cell and neighbor are equal-sized, and cell has been chosen to assemble this face, so calculate from cell

              const unsigned int neighbor2 =
                (cell->has_periodic_neighbor(face_no)
                 ?
                 // how does the periodic neighbor talk about this cell?
                 cell->periodic_neighbor_of_periodic_neighbor( face_no )
                 :
                 // how does the neighbor talk about this cell?
                 cell->neighbor_of_neighbor(face_no));

              // set up neighbor values
              scratch.neighbor_face_finite_element_values->reinit (neighbor, neighbor2);

              scratch.neighbor_face_material_model_inputs.reinit  (*scratch.neighbor_face_finite_element_values,
                                                                   neighbor,
                                                                   this->introspection(),
                                                                   this->get_current_linearization_point());

              this->create_additional_material_model_outputs(scratch.neighbor_face_material_model_outputs);
              this->get_heating_model_manager().create_additional_material_model_inputs_and_outputs(scratch.neighbor_face_material_model_inputs,
                  scratch.neighbor_face_material_model_outputs);

              this->get_material_model().fill_additional_material_model_inputs(scratch.neighbor_face_material_model_inputs,
                                                                               this->get_current_linearization_point(),
                                                                               *scratch.neighbor_face_finite_element_values,
                                                                               this->introspection());

              this->get_material_model().evaluate(scratch.neighbor_face_material_model_inputs,
                                                  scratch.neighbor_face_material_model_outputs);

              if (parameters.formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      scratch.neighbor_face_material_model_outputs.densities[q] = this->get_adiabatic_conditions().density(scratch.neighbor_face_material_model_inputs.position[q]);
                    }
                }

              this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                         scratch.neighbor_face_material_model_outputs,
                                                         scratch.neighbor_face_heating_model_outputs);

              std::vector<types::global_dof_index> neighbor_dof_indices (dofs_per_cell);
              // get all dof indices on the neighbor, then extract those
              // that correspond to the solution_field we are interested in
              neighbor->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      data.neighbor_dof_indices[nth_interface_matrix(fe.reference_cell(), face_no)][i_advection] = neighbor_dof_indices[i];
                      ++i_advection;
                    }
                  ++i;
                }
              data.assembled_matrices[nth_interface_matrix(fe.reference_cell(), face_no)] = true;

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
                          scratch.face_grad_phi_field[i_advection]          = (*scratch.face_finite_element_values)[solution_field].gradient (i, q);
                          scratch.face_phi_field[i_advection]               = (*scratch.face_finite_element_values)[solution_field].value (i, q);
                          scratch.neighbor_face_grad_phi_field[i_advection] = (*scratch.neighbor_face_finite_element_values)[solution_field].gradient (i, q);
                          scratch.neighbor_face_phi_field[i_advection]      = (*scratch.neighbor_face_finite_element_values)[solution_field].value (i, q);
                          ++i_advection;
                        }
                      ++i;
                    }

                  const double density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.densities[q] *
                     scratch.face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  AssertThrow (density_c_P >= 0,
                               ExcMessage ("The product of density and c_P needs to be a "
                                           "non-negative quantity."));

                  const double conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  AssertThrow (density_c_P + latent_heat_LHS >= 0,
                               ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                           "to the left hand side needs to be a non-negative quantity."));

                  const double penalty = (advection_field.is_temperature()
                                          ?
                                          parameters.discontinuous_penalty
                                          * parameters.temperature_degree
                                          * parameters.temperature_degree
                                          / approximate_face_measure(face)
                                          * conductivity
                                          / (density_c_P + latent_heat_LHS)
                                          :
                                          0.0);

                  Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                  // Subtract off the mesh velocity for ALE corrections if necessary
                  if (parameters.mesh_deformation_enabled)
                    current_u -= scratch.face_mesh_velocity_values[q];

                  const double neighbor_density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.densities[q] *
                     scratch.neighbor_face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  AssertThrow (neighbor_density_c_P >= 0,
                               ExcMessage ("The product of density and c_P on the neighbor needs to be a "
                                           "non-negative quantity."));

                  const double neighbor_conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double neighbor_latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  AssertThrow (neighbor_density_c_P + neighbor_latent_heat_LHS >= 0,
                               ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                           "to the left hand side on the neighbor needs to be a non-negative quantity."));

                  const double neighbor_penalty = (advection_field.is_temperature()
                                                   ?
                                                   parameters.discontinuous_penalty
                                                   * parameters.temperature_degree
                                                   * parameters.temperature_degree
                                                   / approximate_face_measure(neighbor->face(neighbor2))
                                                   * neighbor_conductivity
                                                   / (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                                   :
                                                   0.0);

                  const double max_penalty = std::max(penalty, neighbor_penalty);

                  const double max_density_c_P_and_latent_heat =
                    std::max(density_c_P + latent_heat_LHS,
                             neighbor_density_c_P + neighbor_latent_heat_LHS);

                  AssertThrow (numbers::is_finite(max_density_c_P_and_latent_heat),
                               ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a finite quantity."));
                  AssertThrow (max_density_c_P_and_latent_heat >= 0,
                               ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a "
                                           "non-negative quantity."));

                  /**
                   * The discontinuous Galerkin method uses 2 types of jumps over edges:
                   * undirected and directed jumps. Undirected jumps are dependent only
                   * on the order of the numbering of cells. Directed jumps are dependent
                   * on the direction of the flow. Thus the flow-dependent terms below are
                   * only calculated if the edge is an inflow edge.
                   */
                  const bool inflow = ((current_u * scratch.face_finite_element_values->normal_vector(q)) < 0.);

                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    {
                      for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j)
                          += (- 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[i]
                              * scratch.face_phi_field[j]

                              - (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);

                          data.local_matrices_int_ext[nth_interface_matrix(fe.reference_cell(), face_no)](i,j)
                          += (- 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[j]
                              * scratch.face_phi_field[i]

                              + (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);

                          data.local_matrices_ext_int[nth_interface_matrix(fe.reference_cell(), face_no)](i,j)
                          += (+ 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              - 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[j]
                              * scratch.neighbor_face_phi_field[i]

                              - (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);

                          data.local_matrices_ext_ext[nth_interface_matrix(fe.reference_cell(), face_no)](i,j)
                          += (+ 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              + 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[i]
                              * scratch.neighbor_face_phi_field[j]

                              + (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);
                        }
                    }
                }
            }
          else
            {
              /* neighbor is taking responsibility for assembly of this face, because
               * either (1) neighbor is coarser, or
               *        (2) neighbor is equally-sized and
               *           (a) neighbor is on a different subdomain, with lower subdmain_id(), or
               *           (b) neighbor is on the same subdomain and has lower index().
              */
            }
        }
      // neighbor has children, so always assemble from here.
      else
        {
          const unsigned int neighbor2 =
            (cell_has_periodic_neighbor
             ?
             cell->periodic_neighbor_face_no(face_no)
             :
             cell->neighbor_face_no(face_no));

          // Loop over subfaces. We know that the neighbor is finer, so we could loop over the subfaces of the current
          // face. but if we are at a periodic boundary, then the face of the current cell has no children, so instead use
          // the children of the periodic neighbor's corresponding face since we know that the letter does indeed have
          // children (because we know that the neighbor is refined).
          typename DoFHandler<dim>::face_iterator neighbor_face=neighbor->face(neighbor2);
          for (unsigned int subface_no=0; subface_no<neighbor_face->n_children(); ++subface_no)
            {
              const typename DoFHandler<dim>::active_cell_iterator neighbor_child
                = ( cell_has_periodic_neighbor
                    ?
                    cell->periodic_neighbor_child_on_subface(face_no,subface_no)
                    :
                    cell->neighbor_child_on_subface (face_no, subface_no));

              // set up subface values
              scratch.subface_finite_element_values->reinit (cell, face_no, subface_no);

              // subface->face
              (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_current_linearization_point(),
                  scratch.face_current_velocity_values);

              // get the mesh velocity, as we need to subtract it off of the advection systems
              if (parameters.mesh_deformation_enabled)
                (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_mesh_velocity(),
                    scratch.face_mesh_velocity_values);

              scratch.face_material_model_inputs.reinit  (*scratch.subface_finite_element_values,
                                                          cell,
                                                          this->introspection(),
                                                          this->get_current_linearization_point());
              this->get_material_model().evaluate(scratch.face_material_model_inputs,
                                                  scratch.face_material_model_outputs);

              if (parameters.formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      scratch.face_material_model_outputs.densities[q] = this->get_adiabatic_conditions().density(scratch.face_material_model_inputs.position[q]);
                    }
                }

              this->get_heating_model_manager().evaluate(scratch.face_material_model_inputs,
                                                         scratch.face_material_model_outputs,
                                                         scratch.face_heating_model_outputs);

              // set up neighbor values
              scratch.neighbor_face_finite_element_values->reinit (neighbor_child, neighbor2);

              scratch.neighbor_face_material_model_inputs.reinit  (*scratch.neighbor_face_finite_element_values,
                                                                   neighbor_child,
                                                                   this->introspection(),
                                                                   this->get_current_linearization_point());
              this->get_material_model().evaluate(scratch.neighbor_face_material_model_inputs,
                                                  scratch.neighbor_face_material_model_outputs);

              if (parameters.formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      scratch.neighbor_face_material_model_outputs.densities[q] = this->get_adiabatic_conditions().density(scratch.neighbor_face_material_model_inputs.position[q]);
                    }
                }

              this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                         scratch.neighbor_face_material_model_outputs,
                                                         scratch.neighbor_face_heating_model_outputs);

              std::vector<types::global_dof_index> neighbor_dof_indices (fe.dofs_per_cell);
              // get all dof indices on the neighbor, then extract those
              // that correspond to the solution_field we are interested in
              neighbor_child->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      data.neighbor_dof_indices[nth_interface_matrix(fe.reference_cell(), face_no, subface_no)][i_advection] = neighbor_dof_indices[i];
                      ++i_advection;
                    }
                  ++i;
                }
              data.assembled_matrices[nth_interface_matrix(fe.reference_cell(), face_no, subface_no)] = true;

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  // precompute the values of shape functions and their gradients.
                  // We only need to look up values of shape functions if they
                  // belong to 'our' component. They are zero otherwise anyway.
                  // Note that we later only look at the values that we do set here.
                  for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
                    {
                      if (fe.system_to_component_index(i).first == solution_component)
                        {
                          scratch.face_grad_phi_field[i_advection]          = (*scratch.subface_finite_element_values)[solution_field].gradient (i, q);
                          scratch.face_phi_field[i_advection]               = (*scratch.subface_finite_element_values)[solution_field].value (i, q);
                          scratch.neighbor_face_grad_phi_field[i_advection] = (*scratch.neighbor_face_finite_element_values)[solution_field].gradient (i, q);
                          scratch.neighbor_face_phi_field[i_advection]      = (*scratch.neighbor_face_finite_element_values)[solution_field].value (i, q);
                          ++i_advection;
                        }
                      ++i;
                    }

                  const double density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.densities[q] *
                     scratch.face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  AssertThrow (density_c_P >= 0,
                               ExcMessage ("The product of density and c_P needs to be a "
                                           "non-negative quantity."));

                  const double conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  AssertThrow (density_c_P + latent_heat_LHS >= 0,
                               ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                           "to the left hand side needs to be a non-negative quantity."));

                  const double penalty = (advection_field.is_temperature()
                                          ?
                                          parameters.discontinuous_penalty
                                          * parameters.temperature_degree
                                          * parameters.temperature_degree
                                          / approximate_face_measure(face)
                                          * conductivity
                                          / (density_c_P + latent_heat_LHS)
                                          :
                                          0.0);

                  Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                  // Subtract off the mesh velocity for ALE corrections if necessary
                  if (parameters.mesh_deformation_enabled)
                    current_u -= scratch.face_mesh_velocity_values[q];

                  const double neighbor_density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.densities[q] *
                     scratch.neighbor_face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  AssertThrow (neighbor_density_c_P >= 0,
                               ExcMessage ("The product of density and c_P on the neighbor needs to be a "
                                           "non-negative quantity."));

                  const double neighbor_conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double neighbor_latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  AssertThrow (neighbor_density_c_P + neighbor_latent_heat_LHS >= 0,
                               ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                           "to the left hand side on the neighbor needs to be a non-negative quantity."));

                  const double neighbor_penalty = (advection_field.is_temperature()
                                                   ?
                                                   parameters.discontinuous_penalty
                                                   * parameters.temperature_degree
                                                   * parameters.temperature_degree
                                                   / approximate_face_measure(neighbor_child->face(neighbor2))
                                                   * neighbor_conductivity
                                                   / (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                                   :
                                                   0.0);

                  const double max_penalty = std::max(penalty, neighbor_penalty);

                  const double max_density_c_P_and_latent_heat =
                    std::max(density_c_P + latent_heat_LHS,
                             neighbor_density_c_P + neighbor_latent_heat_LHS);

                  AssertThrow (numbers::is_finite(max_density_c_P_and_latent_heat),
                               ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a finite quantity."));
                  AssertThrow (max_density_c_P_and_latent_heat >= 0,
                               ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a "
                                           "non-negative quantity."));

                  /**
                   * The discontinuous Galerkin method uses 2 types of jumps over edges:
                   * undirected and directed jumps. Undirected jumps are dependent only
                   * on the order of the numbering of cells. Directed jumps are dependent
                   * on the direction of the flow. Thus the flow-dependent terms below are
                   * only calculated if the edge is an inflow edge.
                   */
                  const bool inflow = ((current_u * scratch.subface_finite_element_values->normal_vector(q)) < 0.);

                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    {
                      for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j)
                          += (- 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[i]
                              * scratch.face_phi_field[j]

                              - (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);

                          data.local_matrices_int_ext[nth_interface_matrix(fe.reference_cell(), face_no, subface_no)](i,j)
                          += (- 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[j]
                              * scratch.face_phi_field[i]

                              + (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);

                          data.local_matrices_ext_int[nth_interface_matrix(fe.reference_cell(), face_no, subface_no)](i,j)
                          += (+ 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              - 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[j]
                              * scratch.neighbor_face_phi_field[i]

                              - (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);

                          data.local_matrices_ext_ext[nth_interface_matrix(fe.reference_cell(), face_no, subface_no)](i,j)
                          += (+ 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              + 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[i]
                              * scratch.neighbor_face_phi_field[j]

                              + (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);
                        }
                    }
                }
            }
        }
    }
  }
} // namespace aspect

// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class AdvectionSystem<dim>; \
  template class DiffusionSystem<dim>; \
  template class DarcySystem<dim>; \
  template class AdvectionSystemBoundaryFace<dim>; \
  template class AdvectionSystemInteriorFace<dim>; \
  template class AdvectionSystemBoundaryHeatFlux<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

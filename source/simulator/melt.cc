/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#include <aspect/melt.h>
#include <aspect/simulator.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

using namespace dealii;


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void MeltOutputs<dim>::average (const MaterialAveraging::AveragingOperation operation,
                                    const FullMatrix<double>  &projection_matrix,
                                    const FullMatrix<double>  &expansion_matrix)
    {
      average_property (operation, projection_matrix, expansion_matrix,
                        compaction_viscosities);
      average_property (operation, projection_matrix, expansion_matrix,
                        fluid_viscosities);
      average_property (operation, projection_matrix, expansion_matrix,
                        permeabilities);
      average_property (operation, projection_matrix, expansion_matrix,
                        fluid_densities);

      // The fluid density gradients are unfortunately stored in reverse
      // indexing and averaging is not implemented for tensors (only for
      // doubles). It's also not quite clear whether these should
      // really be averaged, so avoid this for now
    }
  }

  namespace Assemblers
  {
    template <int dim>
    void
    MeltEquations<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      MeltHandler<dim>::create_material_model_outputs(outputs);
    }


    template <int dim>
    void
    MeltEquations<dim>::
    local_assemble_stokes_preconditioner_melt (const double                                             pressure_scaling,
                                               internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                                               internal::Assembly::CopyData::StokesPreconditioner<dim> &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();
      const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
      const unsigned int   n_q_points      = scratch.finite_element_values.n_quadrature_points;

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();
      const FEValuesExtractors::Scalar ex_p_c = introspection.variable("compaction pressure").extractor_scalar();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              scratch.grads_phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(k,q);
              scratch.phi_p[k]       = scratch.finite_element_values[ex_p_f].value (k, q);
              scratch.phi_p_c[k] = scratch.finite_element_values[ex_p_c].value (k, q);
              scratch.grad_phi_p[k] = scratch.finite_element_values[ex_p_f].gradient (k, q);
            }

          const double eta = scratch.material_model_outputs.viscosities[q];

          /*
            - R = 1/eta M_p + K_D L_p for p
            S = - (1/eta + 1/viscosity_c)  M_p  for p_c
          */
          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.material_model_inputs.composition[q][porosity_index], 0.0);

          const double K_D = (porosity > this->get_melt_handler().melt_transport_threshold
                              ?
                              melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]
                              :
                              0.0);
          const double viscosity_c = melt_outputs->compaction_viscosities[q];

          const SymmetricTensor<4,dim> &stress_strain_director =
            scratch.material_model_outputs.stress_strain_directors[q];
          const bool use_tensor = (stress_strain_director != dealii::identity_tensor<dim> ());

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              if (fe.system_to_component_index(i).first
                  ==
                  fe.system_to_component_index(j).first)
                data.local_matrix(i,j) += ((use_tensor ?
                                            eta * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                            :
                                            eta * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                           +
                                           (1./eta *
                                            pressure_scaling *
                                            pressure_scaling)
                                           * scratch.phi_p[i] * scratch.phi_p[j]
                                           +
                                           (K_D *
                                            pressure_scaling *
                                            pressure_scaling) *
                                           scratch.grad_phi_p[i] *
                                           scratch.grad_phi_p[j]
                                           +
                                           (1./eta + 1./viscosity_c) *
                                           pressure_scaling *
                                           pressure_scaling *
                                           (scratch.phi_p_c[i] * scratch.phi_p_c[j])
                                          )
                                          * scratch.finite_element_values.JxW(q);
        }
    }



    template <int dim>
    void
    MeltEquations<dim>::
    local_assemble_stokes_system_melt (const typename DoFHandler<dim>::active_cell_iterator &/*cell*/,
                                       const double                                     pressure_scaling,
                                       const bool                                       rebuild_stokes_matrix,
                                       internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                       internal::Assembly::CopyData::StokesSystem<dim> &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      const FEValuesExtractors::Scalar extractor_pressure = introspection.variable("fluid pressure").extractor_scalar();
      const FEValuesExtractors::Scalar ex_p_c = introspection.variable("compaction pressure").extractor_scalar();
      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              scratch.phi_u[k]   = scratch.finite_element_values[introspection.extractors.velocities].value (k,q);
              scratch.phi_p[k]   = scratch.finite_element_values[extractor_pressure].value (k, q);
              scratch.phi_p_c[k] = scratch.finite_element_values[ex_p_c].value (k, q);
              scratch.grad_phi_p[k] = scratch.finite_element_values[extractor_pressure].gradient (k, q);

              if (rebuild_stokes_matrix)
                {
                  scratch.grads_phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(k,q);
                  scratch.div_phi_u[k]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (k, q);
                }
            }

          // Viscosity scalar
          const double eta = (rebuild_stokes_matrix
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              std::numeric_limits<double>::quiet_NaN());

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));
          const SymmetricTensor<4,dim> &stress_strain_director =
            scratch.material_model_outputs.stress_strain_directors[q];
          const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

          const double compressibility
            = (this->get_material_model().is_compressible()
               ?
               scratch.material_model_outputs.compressibilities[q]
               :
               std::numeric_limits<double>::quiet_NaN() );

          const double density_s = scratch.material_model_outputs.densities[q]; // density of the solid

          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.material_model_inputs.composition[q][porosity_index],0.000);
          const double K_D = (porosity > this->get_melt_handler().melt_transport_threshold
                              ?
                              melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]
                              :
                              0.0);
          const double viscosity_c = melt_outputs->compaction_viscosities[q];
          const Tensor<1,dim> density_gradient_f = melt_outputs->fluid_density_gradients[q];
          const double density_f = melt_outputs->fluid_densities[q];
          const double p_f_RHS = compute_fluid_pressure_rhs(scratch,
                                                            scratch.material_model_inputs,
                                                            scratch.material_model_outputs,
                                                            q);
          const double bulk_density = (1.0 - porosity) * density_s + porosity * density_f;


          if (rebuild_stokes_matrix)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                data.local_matrix(i,j) += ( (use_tensor ?
                                             eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                             :
                                             eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                            - (use_tensor ?
                                               eta * 2.0/3.0 * (scratch.div_phi_u[i] * trace(stress_strain_director * scratch.grads_phi_u[j]))
                                               :
                                               eta * 2.0/3.0 * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                              )
                                            - (pressure_scaling *
                                               scratch.div_phi_u[i] * scratch.phi_p[j])
                                            // finally the term -div(u). note the negative sign to make this
                                            // operator adjoint to the grad(p) term
                                            - (pressure_scaling *
                                               scratch.phi_p[i] * scratch.div_phi_u[j])
                                            +
                                            (- pressure_scaling * pressure_scaling / viscosity_c
                                             * scratch.phi_p_c[i] * scratch.phi_p_c[j])
                                            - pressure_scaling * scratch.div_phi_u[i] * scratch.phi_p_c[j]
                                            - pressure_scaling * scratch.phi_p_c[i] * scratch.div_phi_u[j]
                                            - K_D * pressure_scaling * pressure_scaling *
                                            (scratch.grad_phi_p[i] * scratch.grad_phi_p[j])
                                            + (this->get_material_model().is_compressible()
                                               ?
                                               K_D * pressure_scaling * pressure_scaling / density_f
                                               * scratch.phi_p[i] * (scratch.grad_phi_p[j] * density_gradient_f)
                                               :
                                               0.0))
                                          * scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            data.local_rhs(i) += (
                                   (bulk_density * gravity * scratch.phi_u[i])
                                   +
                                   // add the term that results from the compressibility. compared
                                   // to the manual, this term seems to have the wrong sign, but this
                                   // is because we negate the entire equation to make sure we get
                                   // -div(u) as the adjoint operator of grad(p) (see above where
                                   // we assemble the matrix)
                                   (this->get_material_model().is_compressible()
                                    ?
                                    (pressure_scaling *
                                     compressibility * density_s *
                                     (scratch.velocity_values[q] * gravity) *
                                     scratch.phi_p[i])
                                    :
                                    0)
                                   + pressure_scaling *
                                   p_f_RHS * scratch.phi_p[i]
                                   - pressure_scaling *
                                   K_D * density_f *
                                   (scratch.grad_phi_p[i] * gravity)
                                 )
                                 * scratch.finite_element_values.JxW(q);
        }
    }



    template <int dim>
    void
    MeltEquations<dim>::
    local_assemble_stokes_system_melt_boundary (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                const unsigned int                                    face_no,
                                                const double                                          pressure_scaling,
                                                internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                                internal::Assembly::CopyData::StokesSystem<dim>      &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const unsigned int n_face_q_points = scratch.face_finite_element_values.n_quadrature_points;
      const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
      const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.face_material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      std::vector<double> grad_p_f(n_face_q_points);
      this->get_melt_handler().fluid_pressure_boundary_conditions->fluid_pressure_gradient(
        cell->face(face_no)->boundary_id(),
        scratch.face_material_model_inputs,
        scratch.face_material_model_outputs,
        scratch.face_finite_element_values.get_all_normal_vectors(),
        grad_p_f);

      for (unsigned int q=0; q<n_face_q_points; ++q)
        {
          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.face_finite_element_values.quadrature_point(q));
          const double density_f = melt_outputs->fluid_densities[q];

          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.face_material_model_inputs.composition[q][porosity_index],0.0);

          const double K_D = (porosity > this->get_melt_handler().melt_transport_threshold
                              ?
                              melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]
                              :
                              0.0);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              // apply the fluid pressure boundary condition
              data.local_rhs(i) += (scratch.face_finite_element_values[ex_p_f].value(i, q)
                                    * pressure_scaling * K_D *
                                    (density_f
                                     * (scratch.face_finite_element_values.normal_vector(q) * gravity)
                                     - grad_p_f[q])
                                    * scratch.face_finite_element_values.JxW(q));
            }
        }
    }


    template <int dim>
    void
    MeltEquations<dim>::
    local_assemble_advection_system_melt (const typename Simulator<dim>::AdvectionField &advection_field,
                                          const double artificial_viscosity,
                                          internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                          internal::Assembly::CopyData::AdvectionSystem<dim> &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const bool use_bdf2_scheme = (this->get_timestep_number() > 1);
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();

      const unsigned int solution_component
        = (advection_field.is_temperature()
           ?
           introspection.component_indices.temperature
           :
           introspection.component_indices.compositional_fields[advection_field.compositional_variable]
          );

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      std::vector<Tensor<1,dim> > fluid_velocity_values(n_q_points);
      const FEValuesExtractors::Vector ex_u_f = introspection.variable("fluid velocity").extractor_vector();
      scratch.finite_element_values[ex_u_f].get_function_values (this->get_solution(),fluid_velocity_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          // precompute the values of shape functions and their gradients.
          // We only need to look up values of shape functions if they
          // belong to 'our' component. They are zero otherwise anyway.
          // Note that we later only look at the values that we do set here.
          for (unsigned int k=0; k<advection_dofs_per_cell; ++k)
            {
              scratch.grad_phi_field[k] = scratch.finite_element_values[solution_field].gradient (scratch.finite_element_values.get_fe().component_to_system_index(solution_component, k),q);
              scratch.phi_field[k]      = scratch.finite_element_values[solution_field].value (scratch.finite_element_values.get_fe().component_to_system_index(solution_component, k), q);
            }

          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.material_model_inputs.composition[q][porosity_index],0.0);
          const double bulk_density = (1.0 - porosity) * scratch.material_model_outputs.densities[q] + porosity * melt_outputs->fluid_densities[q];

          const double density_c_P              =
            ((advection_field.is_temperature())
             ?
             bulk_density *
             scratch.material_model_outputs.specific_heat[q]
             :
             1.0);

          Assert (density_c_P >= 0,
                  ExcMessage ("The product of density and c_P needs to be a "
                              "non-negative quantity."));

          const double conductivity =
            (advection_field.is_temperature()
             ?
             scratch.material_model_outputs.thermal_conductivities[q]
             :
             0.0);
          const double latent_heat_LHS =
            ((advection_field.is_temperature())
             ?
             scratch.heating_model_outputs.lhs_latent_heat_terms[q]
             :
             0.0);
          Assert (density_c_P + latent_heat_LHS >= 0,
                  ExcMessage ("The sum of density times c_P and the latent heat contribution "
                              "to the left hand side needs to be a non-negative quantity."));

          const double gamma =
            ((advection_field.is_temperature())
             ?
             scratch.heating_model_outputs.heating_source_terms[q]
             :
             0.0);

          const double reaction_term =
            ((advection_field.is_temperature() || this->get_melt_handler().is_porosity(advection_field))
             ?
             0.0
             :
             scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]);

          const double melt_transport_RHS = compute_melting_RHS (scratch,
                                                                 scratch.material_model_inputs,
                                                                 scratch.material_model_outputs,
                                                                 advection_field,
                                                                 q);

          const double field_term_for_rhs
            = (use_bdf2_scheme ?
               ((scratch.old_field_values)[q] *
                (1 + time_step/old_time_step)
                -
                (scratch.old_old_field_values)[q] *
                (time_step * time_step) /
                (old_time_step * (time_step + old_time_step)))
               :
               (scratch.old_field_values)[q])
              *
              (density_c_P + latent_heat_LHS);

          Tensor<1,dim> current_u = scratch.current_velocity_values[q];
          //Subtract off the mesh velocity for ALE corrections if necessary
          if (this->get_parameters().free_surface_enabled)
            current_u -= scratch.mesh_velocity_values[q];

          const double melt_transport_LHS =
            (this->get_melt_handler().is_porosity(advection_field)
             ?
             scratch.current_velocity_divergences[q]
             + (this->get_material_model().is_compressible()
                ?
                scratch.material_model_outputs.compressibilities[q]
                * scratch.material_model_outputs.densities[q]
                * current_u
                * this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q))
                :
                0.0)
             :
             0.0);

          const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                    (time_step + old_time_step)) : 1.0;

          // for advective transport of heat by melt, we have to consider melt velocities
          const Tensor<1,dim> current_u_f = fluid_velocity_values[q];
          double density_c_P_solid = density_c_P;
          double density_c_P_melt = 0.0;

          if (advection_field.is_temperature() && porosity >= this->get_melt_handler().melt_transport_threshold
              && this->get_melt_handler().heat_advection_by_melt)
            {
              density_c_P_solid = (1.0 - porosity) * scratch.material_model_outputs.densities[q] * scratch.material_model_outputs.specific_heat[q];
              density_c_P_melt = porosity * melt_outputs->fluid_densities[q] * scratch.material_model_outputs.specific_heat[q];
            }

          // do the actual assembly. note that we only need to loop over the advection
          // shape functions because these are the only contributions we compute here
          for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
            {
              data.local_rhs(i)
              += (field_term_for_rhs * scratch.phi_field[i]
                  + time_step *
                  scratch.phi_field[i]
                  * (gamma + melt_transport_RHS)
                  + scratch.phi_field[i]
                  * reaction_term)
                 *
                 scratch.finite_element_values.JxW(q);

              for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                {
                  data.local_matrix(i,j)
                  += (
                       (time_step * (conductivity + artificial_viscosity)
                        * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                       + ((time_step * (scratch.phi_field[i] * (current_u * scratch.grad_phi_field[j])))
                          + (factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                       (density_c_P_solid + latent_heat_LHS)
                       + ((time_step * (scratch.phi_field[i] * (current_u_f * scratch.grad_phi_field[j])))
                          + (factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                       (density_c_P_melt)
                       + time_step * scratch.phi_field[i] * scratch.phi_field[j] * melt_transport_LHS
                     )
                     * scratch.finite_element_values.JxW(q);
                }
            }
        }
    }

    template <int dim>
    std::vector<double>
    MeltEquations<dim>::
    compute_advection_system_residual_melt(const typename Simulator<dim>::AdvectionField     &advection_field,
                                           internal::Assembly::Scratch::AdvectionSystem<dim> &scratch) const
    {
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      std::vector<double> residuals(n_q_points);

      HeatingModel::HeatingModelOutputs heating_model_outputs(n_q_points, this->get_parameters().n_compositional_fields);
      this->get_heating_model_manager().evaluate(scratch.material_model_inputs,
                                                 scratch.material_model_outputs,
                                                 heating_model_outputs);

      for (unsigned int q=0; q < n_q_points; ++q)
        {
          const Tensor<1,dim> u = (scratch.old_velocity_values[q] +
                                   scratch.old_old_velocity_values[q]) / 2;

          const double dField_dt = (this->get_old_timestep() == 0.0) ? 0 :
                                   (
                                     ((scratch.old_field_values)[q] - (scratch.old_old_field_values)[q])
                                     / this->get_old_timestep());
          const double u_grad_field = u * (scratch.old_field_grads[q] +
                                           scratch.old_old_field_grads[q]) / 2;

          const double density              = ((advection_field.is_temperature())
                                               ? scratch.material_model_outputs.densities[q] : 1.0);
          const double conductivity = ((advection_field.is_temperature()) ? scratch.material_model_outputs.thermal_conductivities[q] : 0.0);
          const double c_P                  = ((advection_field.is_temperature()) ? scratch.material_model_outputs.specific_heat[q] : 1.0);
          const double k_Delta_field = conductivity
                                       * (scratch.old_field_laplacians[q] +
                                          scratch.old_old_field_laplacians[q]) / 2;

          const double field = ((scratch.old_field_values)[q] + (scratch.old_old_field_values)[q]) / 2;

          const double gamma =
            ((advection_field.is_temperature())
             ?
             heating_model_outputs.heating_source_terms[q]
             :
             0.0);

          const double latent_heat_LHS =
            ((advection_field.is_temperature())
             ?
             heating_model_outputs.lhs_latent_heat_terms[q]
             :
             0.0);

          const double dreaction_term_dt =
            (advection_field.is_temperature() || this->get_old_timestep() == 0 || (this->get_melt_handler().is_porosity(advection_field)
                                                                                   && this->include_melt_transport()))
            ?
            0.0
            :
            (scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]
             / this->get_old_timestep());

          const double melt_transport_RHS = compute_melting_RHS (scratch,
                                                                 scratch.material_model_inputs,
                                                                 scratch.material_model_outputs,
                                                                 advection_field,
                                                                 q);

          const double melt_transport_LHS =
            (this->get_melt_handler().is_porosity(advection_field)
             ?
             scratch.current_velocity_divergences[q]
             + (this->get_material_model().is_compressible()
                ?
                scratch.material_model_outputs.compressibilities[q]
                * scratch.material_model_outputs.densities[q]
                * u
                * this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q))
                :
                0.0)
             :
             0.0);

          residuals[q]
            = std::abs((density * c_P + latent_heat_LHS) * (dField_dt + u_grad_field) - k_Delta_field  + melt_transport_LHS * field
                       - gamma  - melt_transport_RHS - dreaction_term_dt);
        }
      return residuals;
    }



    template <int dim>
    double
    MeltEquations<dim>::
    compute_fluid_pressure_rhs(const internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                               MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                               MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                               const unsigned int q_point) const
    {
      if (!this->get_parameters().include_melt_transport)
        return 0.0;

      const Introspection<dim> &introspection = this->introspection();
      MaterialModel::MeltOutputs<dim> *melt_out = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      Assert(melt_out != NULL, ExcInternalError());

      const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
      const unsigned int is_compressible = this->get_material_model().is_compressible();

      const double melting_rate     = material_model_outputs.reaction_terms[q_point][porosity_index];
      const double solid_density    = material_model_outputs.densities[q_point];
      const double fluid_density    = melt_out->fluid_densities[q_point];
      const double solid_compressibility = material_model_outputs.compressibilities[q_point];
      const Tensor<1,dim> fluid_density_gradient = melt_out->fluid_density_gradients[q_point];
      const Tensor<1,dim> current_u = scratch.velocity_values[q_point];
      const double porosity         = std::max(material_model_inputs.composition[q_point][porosity_index],0.0);
      const double K_D = (porosity > this->get_melt_handler().melt_transport_threshold
                          ?
                          melt_out->permeabilities[q_point] / melt_out->fluid_viscosities[q_point]
                          :
                          0.0);

      const Tensor<1,dim>
      gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q_point));

      double fluid_pressure_RHS = 0.0;

      // melting term
      fluid_pressure_RHS -= melting_rate * (1.0/fluid_density - 1.0/solid_density);

      // compression term
      // The whole expression for the first term on the RHS would be
      // (u_s \cdot g) (\phi \rho_f \kappa_f + (1 - \phi) \rho_s \kappa_s).
      // However, we already have the term (u_s \cdot g) \rho_s \kappa_s in the
      // assembly of the stokes system without melt. Because of that, we only
      // need to have -\phi \rho_s \kappa_s here.
      fluid_pressure_RHS += is_compressible
                            ?
                            (current_u * fluid_density_gradient) * porosity / fluid_density
                            - (current_u * gravity) * porosity * solid_density * solid_compressibility
                            + K_D * (fluid_density_gradient * gravity)
                            :
                            0.0;

      return fluid_pressure_RHS;
    }


    template <int dim>
    double
    MeltEquations<dim>::
    compute_melting_RHS(const internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                        const typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
                        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
                        const typename Simulator<dim>::AdvectionField     &advection_field,
                        const unsigned int q_point) const
    {
      if ((!this->get_melt_handler().is_porosity(advection_field)) || (!this->include_melt_transport()))
        return 0.0;

      Assert (material_model_outputs.densities[q_point] > 0,
              ExcMessage ("The density needs to be a positive quantity "
                          "when melt transport is included in the simulation."));

      const double melting_rate         = material_model_outputs.reaction_terms[q_point][advection_field.compositional_variable];
      const double density              = material_model_outputs.densities[q_point];
      const double current_phi          = material_model_inputs.composition[q_point][advection_field.compositional_variable];
      const double divergence_u         = scratch.current_velocity_divergences[q_point];
      const double compressibility      = (this->get_material_model().is_compressible()
                                           ?
                                           material_model_outputs.compressibilities[q_point]
                                           :
                                           0.0);
      const Tensor<1,dim> current_u     = scratch.current_velocity_values[q_point];
      const Tensor<1,dim>
      gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q_point));

      double melt_transport_RHS = melting_rate / density
                                  + divergence_u + compressibility * density * (current_u * gravity);

      if (current_phi < this->get_melt_handler().melt_transport_threshold
          && melting_rate < this->get_melt_handler().melt_transport_threshold)
        melt_transport_RHS = melting_rate / density;

      return melt_transport_RHS;
    }


    namespace OtherTerms
    {
      template <int dim>
      void
      pressure_rhs_compatibility_modification_melt (const SimulatorAccess<dim>                      &simulator_access,
                                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                    internal::Assembly::CopyData::StokesSystem<dim> &data)
      {
        const Introspection<dim> &introspection = simulator_access.introspection();

        const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

        const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();

        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              scratch.phi_p[i] = scratch.finite_element_values[ex_p_f].value (i, q);
              data.local_pressure_shape_function_integrals(i) += scratch.phi_p[i] * scratch.finite_element_values.JxW(q);
            }
      }
    }
  }


  template <int dim>
  void
  MeltHandler<dim>::
  compute_melt_variables(LinearAlgebra::BlockVector &solution)
  {
    if (!this->include_melt_transport())
      return;

    Assert (this->include_melt_transport(), ExcMessage ("'Include melt transport' has to be on to "
                                                        "compute melt variables"));

    LinearAlgebra::BlockVector distributed_vector (this->introspection().index_sets.system_partitioning,
                                                   this->get_mpi_communicator());

    const unsigned int por_idx = this->introspection().compositional_index_for_name("porosity");


    {
      // compute fluid_velocity
      // u_f =  u_s - K_D (nabla p_f - rho_f g) / phi  or = 0
      // by solving a mass matrix problem

      // TODO: lots of cleanup/optimization opportunities here:
      // - store matrix between timesteps
      // - can we maybe reuse system matrix/sparsity pattern?
      // - only construct matrix for the u_f block instead of a block matrix?

      LinearAlgebra::BlockSparseMatrix matrix;
      LinearAlgebra::BlockDynamicSparsityPattern sp;
#ifdef ASPECT_USE_PETSC
      sp.reinit (this->introspection().index_sets.system_relevant_partitioning);
#else
      sp.reinit (this->introspection().index_sets.system_partitioning,
                 this->introspection().index_sets.system_partitioning,
                 this->introspection().index_sets.system_relevant_partitioning,
                 this->get_mpi_communicator());
#endif

      Table<2,DoFTools::Coupling> coupling (this->introspection().n_components,
                                            this->introspection().n_components);
      const unsigned int first_fluid_c_i = this->introspection().variable("fluid velocity").first_component_index;
      for (unsigned int c=0; c<dim; ++c)
        for (unsigned int d=0; d<dim; ++d)
          coupling[first_fluid_c_i+c][first_fluid_c_i+d] = DoFTools::always;

      DoFTools::make_sparsity_pattern (this->get_dof_handler(),
                                       coupling, sp,
                                       this->get_current_constraints(), false,
                                       Utilities::MPI::
                                       this_mpi_process(this->get_mpi_communicator()));

#ifdef ASPECT_USE_PETSC
      SparsityTools::distribute_sparsity_pattern(sp,
                                                 this->get_dof_handler().locally_owned_dofs_per_processor(),
                                                 this->get_mpi_communicator(), this->introspection().index_sets.system_relevant_set);

      sp.compress();
      matrix.reinit (this->introspection().index_sets.system_partitioning,
                     this->introspection().index_sets.system_partitioning,
                     sp, this->get_mpi_communicator());
#else
      sp.compress();
      matrix.reinit (sp);
#endif

      LinearAlgebra::BlockVector rhs, distributed_solution;
      rhs.reinit(this->introspection().index_sets.system_partitioning, this->get_mpi_communicator());
      distributed_solution.reinit(this->introspection().index_sets.system_partitioning, this->get_mpi_communicator());

      const QGauss<dim> quadrature(this->get_parameters().stokes_velocity_degree+1);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         n_q_points    = fe_values.n_quadrature_points;

      std::vector<unsigned int> cell_dof_indices (dofs_per_cell);
      Vector<double> cell_vector (dofs_per_cell);
      FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

      std::vector<double> porosity_values(quadrature.size());
      std::vector<Tensor<1,dim> > grad_p_f_values(quadrature.size());
      std::vector<Tensor<1,dim> > u_s_values(quadrature.size());

      MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
      create_material_model_outputs(out);

      typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active(),
                                                     endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            cell_vector = 0;
            cell_matrix = 0;
            cell->get_dof_indices (cell_dof_indices);
            fe_values.reinit (cell);

            fe_values[this->introspection().extractors.compositional_fields[por_idx]].get_function_values (
              solution, porosity_values);
            fe_values[this->introspection().extractors.velocities].get_function_values (
              solution, u_s_values);

            fe_values[this->introspection().variable("fluid pressure").extractor_scalar()].get_function_gradients (
              solution, grad_p_f_values);
            this->compute_material_model_input_values (solution,
                                                       fe_values,
                                                       cell,
                                                       true,
                                                       in);

            this->get_material_model().evaluate(in, out);

            MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
            Assert(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for computing the melt variables."));
            const FEValuesExtractors::Vector fluid_velocity_extractor = this->introspection().variable("fluid velocity").extractor_vector();

            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += fe_values[fluid_velocity_extractor].value(j,q) *
                                        fe_values[fluid_velocity_extractor].value(i,q) *
                                        fe_values.JxW(q);

                  const double phi = std::max(0.0, porosity_values[q]);

                  // u_f =  u_s - K_D (nabla p_f - rho_f g) / phi  or = 0
                  if (phi > this->get_melt_handler().melt_transport_threshold)
                    {
                      const double K_D = melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q];
                      const Tensor<1,dim>  gravity = this->get_gravity_model().gravity_vector(in.position[q]);
                      cell_vector(i) += (u_s_values[q] - K_D * (grad_p_f_values[q] - melt_outputs->fluid_densities[q]*gravity) / phi)
                                        * fe_values[fluid_velocity_extractor].value(i,q)
                                        * fe_values.JxW(q);
                    }

                }

            this->get_current_constraints().distribute_local_to_global (cell_matrix, cell_vector,
                                                                        cell_dof_indices, matrix, rhs, false);
          }

      rhs.compress (VectorOperation::add);
      matrix.compress (VectorOperation::add);



      LinearAlgebra::PreconditionAMG preconditioner;
      LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
#ifdef ASPECT_USE_PETSC
      Amg_data.symmetric_operator = false;
#else
      //Amg_data.constant_modes = constant_modes;
      Amg_data.elliptic = true;
      Amg_data.higher_order_elements = false;
      Amg_data.smoother_sweeps = 2;
      Amg_data.aggregation_threshold = 0.02;
#endif
      const unsigned int block_idx = this->introspection().variable("fluid velocity").block_index;
      preconditioner.initialize(matrix.block(block_idx, block_idx));

      SolverControl solver_control(5*rhs.size(), 1e-8*rhs.block(block_idx).l2_norm());
      SolverCG<LinearAlgebra::Vector> cg(solver_control);

      cg.solve (matrix.block(block_idx, block_idx), distributed_solution.block(block_idx), rhs.block(block_idx), preconditioner);
      this->get_pcout() << "   Solving for u_f in " << solver_control.last_step() <<" iterations."<< std::endl;

      this->get_current_constraints().distribute (distributed_solution);
      solution.block(block_idx) = distributed_solution.block(block_idx);
    }

    //compute solid pressure
    {
      const unsigned int block_p = this->introspection().block_indices.pressure;

      // Think what we need to do if the pressure is not an FE_Q...
      Assert(this->get_parameters().use_locally_conservative_discretization == false, ExcNotImplemented());
      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.pressure).get_unit_support_points());
      std::vector<double> porosity_values(quadrature.size());
      std::vector<double> p_c_values(quadrature.size());
      std::vector<double> p_f_values(quadrature.size());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values);

      std::vector<types::global_dof_index> local_dof_indices (this->get_fe().dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            cell->get_dof_indices (local_dof_indices);
            fe_values[this->introspection().extractors.compositional_fields[por_idx]].get_function_values (
              solution, porosity_values);
            fe_values[this->introspection().variable("compaction pressure").extractor_scalar()].get_function_values (
              solution, p_c_values);
            fe_values[this->introspection().variable("fluid pressure").extractor_scalar()].get_function_values (
              solution, p_f_values);

            for (unsigned int j=0; j<this->get_fe().base_element(this->introspection().base_elements.pressure).dofs_per_cell; ++j)
              {
                const unsigned int pressure_idx
                  = this->get_fe().component_to_system_index(this->introspection().component_indices.pressure,
                                                             /*dof index within component=*/ j);

                // skip entries that are not locally owned:
                if (!this->get_dof_handler().locally_owned_dofs().is_element(local_dof_indices[pressure_idx]))
                  continue;

                const double phi = std::max(0.0, porosity_values[j]);

                double p = p_f_values[j];
                if (phi < 1.0-this->get_melt_handler().melt_transport_threshold)
                  p = (p_c_values[j] - (phi-1.0) * p_f_values[j]) / (1.0-phi);

                distributed_vector(local_dof_indices[pressure_idx]) = p;
              }
          }
      distributed_vector.block(block_p).compress(VectorOperation::insert);
      solution.block(block_p) = distributed_vector.block(block_p);
    }
  }

  template <int dim>
  bool
  MeltHandler<dim>::
  is_porosity(const typename Simulator<dim>::AdvectionField &advection_field) const
  {
    if (advection_field.field_type != Simulator<dim>::AdvectionField::compositional_field)
      return false;
    else
      return (this->introspection().name_for_compositional_index(advection_field.compositional_variable) == "porosity");
  }


  template <int dim>
  void
  MeltHandler<dim>::
  edit_finite_element_variables(const Parameters<dim> &parameters,
                                std::vector<VariableDeclaration<dim> > &variables)
  {

    if (!parameters.include_melt_transport)
      return;

    variables.insert(variables.begin()+1,
                     VariableDeclaration<dim>(
                       "fluid pressure",
                       (parameters.use_locally_conservative_discretization)
                       ?
                       std_cxx11::shared_ptr<FiniteElement<dim> >(new FE_DGP<dim>(parameters.stokes_velocity_degree-1))
                       :
                       std_cxx11::shared_ptr<FiniteElement<dim> >(new FE_Q<dim>(parameters.stokes_velocity_degree-1)),
                       1,
                       0)); // same block as p_c even without a direct solver!

    variables.insert(variables.begin()+2,
                     VariableDeclaration<dim>(
                       "compaction pressure",
                       (parameters.use_locally_conservative_discretization)
                       ?
                       std_cxx11::shared_ptr<FiniteElement<dim> >(new FE_DGP<dim>(parameters.stokes_velocity_degree-1))
                       :
                       std_cxx11::shared_ptr<FiniteElement<dim> >(new FE_Q<dim>(parameters.stokes_velocity_degree-1)),
                       1,
                       1));

    variables.insert(variables.begin()+3,
                     VariableDeclaration<dim>("fluid velocity",
                                              std_cxx11::shared_ptr<FiniteElement<dim> >(
                                                new FE_Q<dim>(parameters.stokes_velocity_degree)),
                                              dim,
                                              1));

  }


  template <int dim>
  void
  MeltHandler<dim>::
  declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Melt settings");
    {
      prm.declare_entry ("Melt transport threshold", "1e-3",
                         Patterns::Double (),
                         "The porosity limit for melt migration. For smaller porosities, the equations "
                         "reduce to the Stokes equations and do neglect melt transport. Only used "
                         "if Include melt transport is true. ");
      prm.declare_entry ("Heat advection by melt", "false",
                         Patterns::Bool (),
                         "Whether to use a porosity weighted average of the melt and solid velocity "
                         "to advect heat in the temperature equation or not. If this is set to true, "
                         "additional terms are assembled on the left-hand side of the temperature "
                         "advection equation. Only used if Include melt transport is true. "
                         "If this is set to false, only the solid velocity is used (as in models "
                         "without melt migration).");
    }
    prm.leave_subsection();

    FluidPressureBoundaryConditions::declare_parameters<dim> (prm);
  }

  template <int dim>
  MeltHandler<dim>::MeltHandler (ParameterHandler &prm)
  {
    parse_parameters(prm);
  }

  template <int dim>
  void
  MeltHandler<dim>::initialize_simulator (const Simulator<dim> &simulator_object)
  {
    this->SimulatorAccess<dim>::initialize_simulator(simulator_object);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(fluid_pressure_boundary_conditions.get()))
      sim->initialize_simulator (simulator_object);
    fluid_pressure_boundary_conditions->initialize ();
  }

  template <int dim>
  void
  MeltHandler<dim>::
  parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Melt settings");
    {
      melt_transport_threshold = prm.get_double("Melt transport threshold");
      heat_advection_by_melt = prm.get_bool("Heat advection by melt");
    }
    prm.leave_subsection();

    fluid_pressure_boundary_conditions.reset(FluidPressureBoundaryConditions::create_fluid_pressure_boundary<dim>(prm));
    fluid_pressure_boundary_conditions->parse_parameters (prm);
  }


  template <int dim>
  void
  MeltHandler<dim>::
  create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output)
  {
    if (output.template get_additional_output<MaterialModel::MeltOutputs<dim> >() != NULL)
      return;

    const unsigned int n_points = output.viscosities.size();
    const unsigned int n_comp = output.reaction_terms[0].size();
    output.additional_outputs.push_back(
      std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
      (new MaterialModel::MeltOutputs<dim> (n_points, n_comp)));
  }
}




// explicit instantiation of the functions we implement in this file
namespace aspect
{

#define INSTANTIATE(dim) \
  \
  template \
  class \
  MeltHandler<dim>; \
  \
  namespace Assemblers \
  { \
    template class MeltEquations<dim>; \
    \
    \
    namespace OtherTerms \
    { \
      template void pressure_rhs_compatibility_modification_melt<dim> (const SimulatorAccess<dim>                      &simulator_access, \
                                                                       internal::Assembly::Scratch::StokesSystem<dim>  &scratch, \
                                                                       internal::Assembly::CopyData::StokesSystem<dim> &data); \
    } \
  } \
   
  ASPECT_INSTANTIATE(INSTANTIATE)

}

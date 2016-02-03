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

using namespace dealii;


namespace aspect
{
  namespace Assemblers
  {
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

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              scratch.grads_phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(k,q);
              scratch.phi_p[k]       = scratch.finite_element_values[introspection.extractors.fluid_pressure].value (k, q);
              scratch.phi_p_c[k] = scratch.finite_element_values[introspection.extractors.compaction_pressure].value (k, q);
              scratch.grad_phi_p[k] = scratch.finite_element_values[introspection.extractors.fluid_pressure].gradient (k, q);
            }

          const double eta = scratch.material_model_outputs.viscosities[q];

          /*
            - R = 1/eta M_p + K_D L_p for p
            S = - (1/eta + 1/viscosity_c)  M_p  for p_c
          */
          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          double porosity = std::max(scratch.material_model_inputs.composition[q][porosity_index], 0.0);

          double K_D = (porosity > this->get_parameters().melt_transport_threshold
                        ?
                        melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]
                        :
                        0.0);
          double viscosity_c = melt_outputs->compaction_viscosities[q];

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
    local_assemble_stokes_system_melt (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const double                                     pressure_scaling,
                                       const bool                                       rebuild_stokes_matrix,
                                       internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                       internal::Assembly::CopyData::StokesSystem<dim> &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      const FEValuesExtractors::Scalar &extractor_pressure = introspection.extractors.fluid_pressure;
      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              scratch.phi_u[k]   = scratch.finite_element_values[introspection.extractors.velocities].value (k,q);
              scratch.phi_p[k]   = scratch.finite_element_values[extractor_pressure].value (k, q);
              scratch.phi_p_c[k] = scratch.finite_element_values[introspection.extractors.compaction_pressure].value (k, q);
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
          const double K_D = (porosity > this->get_parameters().melt_transport_threshold
                              ?
                              melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]
                              :
                              0.0);
          const double viscosity_c = melt_outputs->compaction_viscosities[q];
          const double compressibility_f = melt_outputs->fluid_compressibilities[q];
          const double density_f = melt_outputs->fluid_densities[q];
          const double p_f_RHS = compute_fluid_pressure_RHS(scratch,
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
                                               K_D * pressure_scaling * pressure_scaling *
                                               compressibility_f * density_f
                                               * scratch.phi_p[i] * (scratch.grad_phi_p[j] * gravity)
                                               :
                                               0.0))
                                          * scratch.finite_element_values.JxW(q);

          Tensor<1,dim> force_u;
          for (unsigned int d=0; d<dim; ++d)
            force_u[d] = scratch.material_model_outputs.force_vector[q][d];
          const double force_p = scratch.material_model_outputs.force_vector[q][dim];

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            data.local_rhs(i) += (
                                   ((bulk_density * gravity + force_u) * scratch.phi_u[i])
                                   + (pressure_scaling * force_p * scratch.phi_p[i])
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

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.face_material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      std::vector<Tensor<1,dim> > grad_p_f(n_face_q_points);
      this->get_fluid_pressure_boundary_conditions().fluid_pressure_gradient(scratch.face_material_model_inputs,
                                                                             scratch.face_material_model_outputs,
                                                                             grad_p_f);

      for (unsigned int q=0; q<n_face_q_points; ++q)
        {
          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.face_finite_element_values.quadrature_point(q));
          const double density_f = melt_outputs->fluid_densities[q];
          const double density_s = scratch.face_material_model_outputs.densities[q];

          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.face_material_model_inputs.composition[q][porosity_index],0.000);

          const double K_D = (porosity > this->get_parameters().melt_transport_threshold
                              ?
                              melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]
                              :
                              0.0);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              // apply the fluid pressure boundary condition
              data.local_rhs(i) += (scratch.face_finite_element_values[introspection.extractors.fluid_pressure].value(i, q)
                                    * pressure_scaling * K_D *
                                    (density_f
                                     * (scratch.face_finite_element_values.get_normal_vectors()[q] * gravity)
                                     - (scratch.face_finite_element_values.get_normal_vectors()[q] * grad_p_f[q]))
                                    * scratch.face_finite_element_values.JxW(q));
            }
        }
    }



    template <int dim>
    double
    MeltEquations<dim>::
    compute_fluid_pressure_RHS(const internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
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
      const double fluid_compressibility = melt_out->fluid_compressibilities[q_point];
      const Tensor<1,dim> current_u = scratch.velocity_values[q_point];
      const double porosity         = std::max(material_model_inputs.composition[q_point][porosity_index],0.0);
      const double K_D = (porosity > this->get_parameters().melt_transport_threshold
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
                            (current_u * gravity) * (porosity * fluid_density * fluid_compressibility
                                                     - porosity * solid_density * solid_compressibility)
                            + K_D * fluid_compressibility * fluid_density * fluid_density * (gravity * gravity)
                            :
                            0.0;

      return fluid_pressure_RHS;
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

        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              scratch.phi_p[i] = scratch.finite_element_values[introspection.extractors.fluid_pressure].value (i, q);
              data.local_pressure_shape_function_integrals(i) += scratch.phi_p[i] * scratch.finite_element_values.JxW(q);
            }
      }
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class MeltEquations<dim>; \
  namespace OtherTerms \
  { \
    template void pressure_rhs_compatibility_modification_melt<dim> (const SimulatorAccess<dim>                      &simulator_access, \
                                                                     internal::Assembly::Scratch::StokesSystem<dim>  &scratch, \
                                                                     internal::Assembly::CopyData::StokesSystem<dim> &data); \
  } \
   
    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}

/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/newton.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/material_model/rheology/elasticity.h>

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    NewtonInterface<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      if (this->get_newton_handler().parameters.newton_derivative_scaling_factor == 0)
        return;

      NewtonHandler<dim>::create_material_model_outputs(outputs);
    }



    template <int dim>
    void
    NewtonStokesPreconditioner<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;
      const double derivative_scaling_factor = this->get_newton_handler().parameters.newton_derivative_scaling_factor;
      const double pressure_scaling = this->get_pressure_scaling();

      const MaterialModel::MaterialAveraging::AveragingOperation
      material_averaging = this->get_parameters().material_averaging;

      // If elasticity is enabled, then we need ElasticOutputs to get the viscoelastic
      // strain rate.
      const MaterialModel::ElasticOutputs<dim> *elastic_out =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();
      if (this->get_parameters().enable_elasticity)
        AssertThrow(elastic_out != nullptr,
                    ExcMessage("Error: The Newton method requires ElasticOutputs when elasticity is enabled."));

      // First loop over all dofs and find those that are in the Stokes system
      // save the component (pressure and dim velocities) each belongs to.
      for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
        {
          if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
            {
              scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
              ++i_stokes;
            }
          ++i;
        }

      const MaterialModel::MaterialModelDerivatives<dim> *derivatives
        = scratch.material_model_outputs.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim>>();

      std::vector<double> deta_deps_times_grads_phi_u(stokes_dofs_per_cell);
      std::vector<double> eps_times_grads_phi_u(stokes_dofs_per_cell);

      // Calculate the weighted average of viscosity derivatives if
      // material averaging is applied.
      if (derivative_scaling_factor > 0 &&
          material_averaging != MaterialModel::MaterialAveraging::none)
        {
          AssertThrow(derivatives != nullptr,
                      ExcMessage ("Error: The newton method requires derivatives from the material model."));

          for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  double avg = 0;
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    avg += derivatives->viscosity_derivative_averaging_weights[q] *
                           (derivatives->viscosity_derivative_wrt_strain_rate[q] *
                            scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i, q));

                  deta_deps_times_grads_phi_u[i_stokes] = avg;

                  ++i_stokes;
                }
              ++i;
            }
        }

      // Loop over all quadrature points and assemble their contributions to
      // the preconditioner matrix
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.grads_phi_u[i_stokes] =
                    scratch.finite_element_values[introspection.extractors
                                                  .velocities].symmetric_gradient(i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection
                                                                          .extractors.pressure].value(i, q);

#if DEBUG
                  // This is needed to test the velocity part of the matrix for
                  // being symmetric positive-definite.
                  scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
#endif
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];
          const double one_over_eta = 1. / eta;
          const double JxW = scratch.finite_element_values.JxW(q);

          // TODO: Find out why in this version of ASPECT adding the derivative to the preconditioning
          // is way worse than the normal preconditioning
          if (derivative_scaling_factor == 0)
            {
              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                  if (scratch.dof_component_indices[i] ==
                      scratch.dof_component_indices[j])
                    data.local_matrix(i, j) += (
                                                 // top left block: for the current case with
                                                 // derivative_scaling_factor==0 the top left block
                                                 // of the system matrix only contains the usual
                                                 // Stokes term. So approximate this block in the
                                                 // same way as we do in the preconditioner when
                                                 // we solve the regular Stokes problem, i.e.,
                                                 // by replacing the symmetric gradient by the
                                                 // regular gradient to make the block more sparse
                                                 (2.0 * eta * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                 +
                                                 // bottom right block: approximate the pressure
                                                 // Schur complement by the pressure mass matrix.
                                                 // if the derivative scaling factor is zero,
                                                 // then we only use the viscosity in the top
                                                 // left block, and so clearly the mass matrix
                                                 // approximation in the bottom right needs to
                                                 // only be scaled by 1/eta, without considering
                                                 // the derivatives
                                                 one_over_eta
                                                 * pressure_scaling
                                                 * pressure_scaling
                                                 * (scratch.phi_p[i] * scratch.phi_p[j]))
                                               * JxW;
            }
          else
            {
              const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
              const typename Newton::Parameters::Stabilization
              preconditioner_stabilization = this->get_newton_handler().parameters.preconditioner_stabilization;

              SymmetricTensor<2,dim> effective_strain_rate = scratch.material_model_inputs.strain_rate[q];
              if (elastic_out != nullptr)
                effective_strain_rate = elastic_out->viscoelastic_strain_rate[q];
              else if ((preconditioner_stabilization & Newton::Parameters::Stabilization::PD) != Newton::Parameters::Stabilization::none)
                effective_strain_rate = deviator(effective_strain_rate);

              // use the spd factor when the stabilization is PD or SPD
              const double alpha = (preconditioner_stabilization & Newton::Parameters::Stabilization::PD) != Newton::Parameters::Stabilization::none ?
                                   Utilities::compute_spd_factor<dim>(eta, effective_strain_rate, viscosity_derivative_wrt_strain_rate,
                                                                      this->get_newton_handler().parameters.SPD_safety_factor)
                                   :
                                   1;

              // pre-compute the tensor contractions
              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                {
                  eps_times_grads_phi_u[i] = effective_strain_rate * scratch.grads_phi_u[i];
                  if (material_averaging == MaterialModel::MaterialAveraging::none)
                    deta_deps_times_grads_phi_u[i] = viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[i];
                }

              // symmetrize when the stabilization is symmetric or SPD
              const bool symmetrize = ((preconditioner_stabilization & Newton::Parameters::Stabilization::symmetric)
                                       != Newton::Parameters::Stabilization::none);

              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                  if (scratch.dof_component_indices[i] ==
                      scratch.dof_component_indices[j])
                    {
                      data.local_matrix(i, j)
                      += (
                           // top left block: approximate J^{uu}
                           (2.0 * eta * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                           +
                           derivative_scaling_factor * alpha *
                           ( symmetrize ?
                             (eps_times_grads_phi_u[i] * deta_deps_times_grads_phi_u[j] +
                              eps_times_grads_phi_u[j] * deta_deps_times_grads_phi_u[i] ) :
                             2.0 * eps_times_grads_phi_u[i] * deta_deps_times_grads_phi_u[j] )
                           +
                           // bottom right block: approximate the
                           // pressure Schur complement by the
                           // pressure mass matrix.  strictly
                           // speaking, we probably ought to also
                           // consider the derivatives deta/deps
                           // here, but we leave this as a TODO
                           one_over_eta
                           * pressure_scaling
                           * pressure_scaling
                           * (scratch.phi_p[i] * scratch.phi_p[j])
                         )
                         * JxW;
                    }
            }
        }
#if DEBUG
      {
        if (this->get_newton_handler().parameters.preconditioner_stabilization == Newton::Parameters::Stabilization::SPD)
          {
            // regardless of whether we do or do not add the Newton
            // linearization terms, we ought to test whether the top-left
            // block of the matrix is Symmetric Positive Definite (SPD).
            //
            // the reason why this is not entirely obvious is described in
            // the paper that discusses the Newton implementation
            Vector<double> tmp (stokes_dofs_per_cell);
            for (unsigned int sample = 0; sample < 100; ++sample)
              {
                // fill a vector with random numbers for the Stokes DoFs
                for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                  if (scratch.dof_component_indices[i] < dim)
                    tmp[i] = Utilities::generate_normal_random_number (0, 1);
                  else
                    tmp[i] = 0;

                // then verify that the product tmp*(A*tmp) -- which by construction
                // of the vector tmp only covers the top-left block of the matrix A
                // -- is indeed positive (or nearly so)
                Assert (data.local_matrix.matrix_norm_square(tmp)/(tmp*tmp)
                        >=
                        -1e-12*data.local_matrix.frobenius_norm(),
                        ExcMessage ("The top left block of the local Newton-Stokes "
                                    "matrix is not positive definite but has an "
                                    "eigenvalue less than "
                                    + Utilities::to_string (data.local_matrix.matrix_norm_square(tmp)/(tmp*tmp))
                                    + " < 0. This should not happen."));
              }
          }
      }
#endif
    }


    template <int dim>
    void
    NewtonStokesPreconditioner<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      if (this->get_parameters().enable_elasticity &&
          outputs.template get_additional_output<MaterialModel::ElasticOutputs<dim>>() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::ElasticOutputs<dim>> (outputs.n_evaluation_points()));
        }

      if (this->get_newton_handler().parameters.newton_derivative_scaling_factor != 0)
        NewtonHandler<dim>::create_material_model_outputs(outputs);
    }


    template <int dim>
    void
    NewtonStokesIncompressibleTerms<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double derivative_scaling_factor = this->get_newton_handler().parameters.newton_derivative_scaling_factor;

      const MaterialModel::MaterialAveraging::AveragingOperation
      material_averaging = this->get_parameters().material_averaging;

      const bool enable_additional_stokes_rhs = this->get_parameters().enable_additional_stokes_rhs;

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> *force = enable_additional_stokes_rhs ?
                                                                            scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>()
                                                                            : nullptr;

      // If elasticity is enabled, then we need ElasticOutputs to get the viscoelastic
      // strain rate and the elastic force.
      const bool enable_elasticity = this->get_parameters().enable_elasticity;
      const MaterialModel::ElasticOutputs<dim> *elastic_out =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();
      if (enable_elasticity)
        AssertThrow(elastic_out != nullptr,
                    ExcMessage("Error: The Newton method requires ElasticOutputs when elasticity is enabled."));

      const bool enable_prescribed_dilation = this->get_parameters().enable_prescribed_dilation;

      const MaterialModel::PrescribedPlasticDilation<dim>
      *prescribed_dilation = enable_prescribed_dilation ?
                             scratch.material_model_outputs.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim>>()
                             : nullptr;

      const bool material_model_is_compressible = (this->get_material_model().is_compressible());

      const MaterialModel::MaterialModelDerivatives<dim> *derivatives
        = scratch.material_model_outputs.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim>>();



      std::vector<double> deta_deps_times_grads_phi_u(stokes_dofs_per_cell);
      std::vector<double> deta_dp_times_phi_p(stokes_dofs_per_cell);
      std::vector<double> eps_times_grads_phi_u(stokes_dofs_per_cell);

      // Calculate the weighted average of viscosity derivatives if
      // material averaging is applied.
      if (derivative_scaling_factor > 0 &&
          material_averaging != MaterialModel::MaterialAveraging::none)
        {
          // This one is only available in debug mode, because normally
          // the AssertThrow in the preconditioner should already have
          // caught the problem.
          Assert(derivatives != nullptr,
                 ExcMessage ("Error: The Newton method requires the material model to "
                             "compute derivatives."));

          for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  double avg_wrt_eps = 0, avg_wrt_p = 0;
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      avg_wrt_eps += derivatives->viscosity_derivative_averaging_weights[q] *
                                     (derivatives->viscosity_derivative_wrt_strain_rate[q] *
                                      scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i, q));
                      avg_wrt_p   += derivatives->viscosity_derivative_averaging_weights[q] *
                                     (derivatives->viscosity_derivative_wrt_pressure[q] *
                                      scratch.finite_element_values[introspection.extractors.pressure].value(i, q));
                    }

                  deta_deps_times_grads_phi_u[i_stokes] = avg_wrt_eps;
                  deta_dp_times_phi_p[i_stokes] = avg_wrt_p;

                  ++i_stokes;
                }
              ++i;
            }
        }

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);

#if DEBUG
                  // This is needed to test the velocity part of the matrix for
                  // being symmetric positive-definite.
                  scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
#endif

                  ++i_stokes;
                }
              ++i;
            }

          // Viscosity scalar
          const double eta = scratch.material_model_outputs.viscosities[q];
          const double pressure = scratch.material_model_inputs.pressure[q];
          const double velocity_divergence = scratch.velocity_divergence[q];

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];

          const double JxW = scratch.finite_element_values.JxW(q);
          const double pressure_scaling = this->get_pressure_scaling();

          // first assemble the rhs
          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) -= (eta * 2.0 * (scratch.grads_phi_u[i] *
                                                 scratch.material_model_inputs.strain_rate[q])
                                    - (scratch.div_phi_u[i] * pressure)
                                    - (pressure_scaling * scratch.phi_p[i] * velocity_divergence)
                                    -(density * gravity * scratch.phi_u[i]))
                                   * JxW;

              if (enable_elasticity)
                data.local_rhs(i) += ( deviator(elastic_out->elastic_force[q])
                                       * scratch.grads_phi_u[i]
                                     ) * JxW;

              if (enable_additional_stokes_rhs)
                data.local_rhs(i) += (force->rhs_u[q] * scratch.phi_u[i]
                                      + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i])
                                     * JxW;

              if (enable_prescribed_dilation)
                data.local_rhs(i) += (
                                       // RHS of - (div u,q) = - (R,q)
                                       - pressure_scaling
                                       * prescribed_dilation->dilation[q]
                                       * scratch.phi_p[i]
                                     ) * JxW;

              // Only assemble this term if we are running incompressible, otherwise this term
              // is already included on the LHS of the equation.
              if (enable_prescribed_dilation && !material_model_is_compressible)
                data.local_rhs(i) += (
                                       // RHS of momentum eqn: - \int 2/3 eta R, div v
                                       - 2.0 / 3.0 * eta
                                       * prescribed_dilation->dilation[q]
                                       * scratch.div_phi_u[i]
                                     ) * JxW;
            }

          // and then the matrix, if necessary
          if (scratch.rebuild_newton_stokes_matrix)
            {
              // always compute the common terms in the Newton matrix
              for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += (
                                                eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j])
                                                // assemble \nabla p as -(p, div v):
                                                - (pressure_scaling *
                                                   (scratch.div_phi_u[i] * scratch.phi_p[j]))
                                                // assemble the term -div(u) as -(div u, q).
                                                // Note the negative sign to make this
                                                // operator adjoint to the grad p term:
                                                - (pressure_scaling *
                                                   (scratch.phi_p[i] * scratch.div_phi_u[j])))
                                              * JxW;
                  }

              // then also see whether we have to add terms due to the
              // Newton linearization
              if (derivative_scaling_factor != 0)
                {
                  const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
                  const double viscosity_derivative_wrt_pressure = derivatives->viscosity_derivative_wrt_pressure[q];
                  const Newton::Parameters::Stabilization velocity_block_stabilization
                    = this->get_newton_handler().parameters.velocity_block_stabilization;

                  SymmetricTensor<2,dim> effective_strain_rate = scratch.material_model_inputs.strain_rate[q];
                  if (elastic_out != nullptr)
                    effective_strain_rate = elastic_out->viscoelastic_strain_rate[q];
                  else if ((velocity_block_stabilization & Newton::Parameters::Stabilization::PD) != Newton::Parameters::Stabilization::none)
                    effective_strain_rate = deviator(effective_strain_rate);

                  // use the spd factor when the stabilization is PD or SPD
                  const double alpha =  (velocity_block_stabilization & Newton::Parameters::Stabilization::PD)
                                        != Newton::Parameters::Stabilization::none
                                        ?
                                        Utilities::compute_spd_factor<dim>(eta, effective_strain_rate, viscosity_derivative_wrt_strain_rate,
                                                                           this->get_newton_handler().parameters.SPD_safety_factor)
                                        :
                                        1;

                  // pre-compute the tensor contractions
                  for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                    {
                      eps_times_grads_phi_u[i] = effective_strain_rate * scratch.grads_phi_u[i];
                      if (material_averaging == MaterialModel::MaterialAveraging::none)
                        {
                          deta_deps_times_grads_phi_u[i] = viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[i];
                          deta_dp_times_phi_p[i]         = viscosity_derivative_wrt_pressure * scratch.phi_p[i];
                        }
                    }

                  // symmetrize when the stabilization is symmetric or SPD
                  const bool symmetrize = ((velocity_block_stabilization & Newton::Parameters::Stabilization::symmetric)
                                           != Newton::Parameters::Stabilization::none);

                  for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                    {
                      for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j) += ( alpha *
                                                      ( symmetrize ?
                                                        ( eps_times_grads_phi_u[i] * deta_deps_times_grads_phi_u[j] +
                                                          eps_times_grads_phi_u[j] * deta_deps_times_grads_phi_u[i] ) :
                                                        2.0 * eps_times_grads_phi_u[i] * deta_deps_times_grads_phi_u[j]
                                                      )
                                                      + pressure_scaling * 2.0 * eps_times_grads_phi_u[i] * deta_dp_times_phi_p[j]
                                                    )
                                                    * derivative_scaling_factor
                                                    * JxW;

                          Assert(dealii::numbers::is_finite(data.local_matrix(i,j)),
                                 ExcMessage ("Error: Assembly matrix is not finite." +
                                             Utilities::to_string(data.local_matrix(i,j)) +
                                             " = " + Utilities::to_string(eta)));
                        }
                    }
                }
            }
        }

#if DEBUG
      if (scratch.rebuild_newton_stokes_matrix)
        {
          if (this->get_newton_handler().parameters.velocity_block_stabilization
              == Newton::Parameters::Stabilization::SPD)
            {
              // regardless of whether we do or do not add the Newton
              // linearization terms, we ought to test whether the top-left
              // block of the matrix is Symmetric Positive Definite (SPD).
              //
              // the reason why this is not entirely obvious is described in
              // the paper that discusses the Newton implementation
              Vector<double> tmp (stokes_dofs_per_cell);
              for (unsigned int sample = 0; sample < 100; ++sample)
                {
                  // fill a vector with random numbers for the Stokes DoFs
                  for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                    if (scratch.dof_component_indices[i] < dim)
                      tmp[i] = Utilities::generate_normal_random_number (0, 1);
                    else
                      tmp[i] = 0;

                  // then verify that the product tmp*(A*tmp) -- which by construction
                  // of the vector tmp only covers the top-left block of the matrix A
                  // -- is indeed positive (or nearly so)
                  Assert (data.local_matrix.matrix_norm_square(tmp)/(tmp*tmp)
                          >=
                          -1e-12*data.local_matrix.frobenius_norm(),
                          ExcMessage ("The top left block of the local Newton-Stokes "
                                      "matrix is not positive definite but has an "
                                      "eigenvalue less than "
                                      + Utilities::to_string (data.local_matrix.matrix_norm_square(tmp)/(tmp*tmp))
                                      + " < 0. This should not happen."));
                }
            }
        }
#endif
    }

    template <int dim>
    void
    NewtonStokesIncompressibleTerms<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.n_evaluation_points();

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));
        }

      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>()->rhs_u.size()
             == n_points, ExcInternalError());

      if ((this->get_parameters().enable_elasticity) &&
          outputs.template get_additional_output<MaterialModel::ElasticOutputs<dim>>() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::ElasticOutputs<dim>> (n_points));
        }

      Assert(!this->get_parameters().enable_elasticity
             ||
             outputs.template get_additional_output<MaterialModel::ElasticOutputs<dim>>()->elastic_force.size()
             == n_points, ExcInternalError());

      // prescribed dilation:
      if (this->get_parameters().enable_prescribed_dilation
          && outputs.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim>>() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedPlasticDilation<dim>> (n_points));
        }

      Assert(!this->get_parameters().enable_prescribed_dilation
             ||
             outputs.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim>>()->dilation.size()
             == n_points, ExcInternalError());

      if (this->get_newton_handler().parameters.newton_derivative_scaling_factor != 0)
        NewtonHandler<dim>::create_material_model_outputs(outputs);
    }


    template <int dim>
    void
    NewtonStokesCompressibleStrainRateViscosityTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const double derivative_scaling_factor = this->get_newton_handler().parameters.newton_derivative_scaling_factor;

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);

                  ++i_stokes;
                }
              ++i;
            }

          // Viscosity scalar
          const double two_thirds = 2.0 / 3.0;
          const double eta_two_thirds = scratch.material_model_outputs.viscosities[q] * two_thirds;
          const double velocity_divergence = scratch.velocity_divergence[q];

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) -= (- eta_two_thirds * (scratch.div_phi_u[i]*velocity_divergence)) * JxW;
            }

          if (scratch.rebuild_stokes_matrix)
            {
              if (derivative_scaling_factor == 0)
                {
                  if (scratch.rebuild_newton_stokes_matrix)
                    for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                      for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j) += (- eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])) * JxW;
                        }
                }
              else
                {
                  const MaterialModel::MaterialModelDerivatives<dim> *derivatives = scratch.material_model_outputs.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim>>();

                  // This one is only available in debug mode, because normally
                  // the AssertThrow in the preconditioner should already have
                  // caught the problem.
                  Assert(derivatives != nullptr, ExcMessage ("Error: The newton method requires the derivatives"));

                  const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
                  const double viscosity_derivative_wrt_pressure = derivatives->viscosity_derivative_wrt_pressure[q];

                  for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                    for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                      {
                        data.local_matrix(i,j) += (-eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                                   - derivative_scaling_factor * two_thirds * scratch.div_phi_u[i] * ( (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[j]) * velocity_divergence)
                                                   - derivative_scaling_factor * two_thirds * (scratch.div_phi_u[i] * viscosity_derivative_wrt_pressure * scratch.phi_p[j]) * velocity_divergence
                                                  )
                                                  * JxW;
                      }
                }
            }
        }
    }



    template <int dim>
    void
    NewtonStokesReferenceDensityCompressibilityTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      // assemble RHS of:
      //  - div u = 1/rho * drho/dz g/||g||* u
      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::reference_density_profile,
             ExcInternalError());

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));
          const double drho_dz_u = scratch.reference_densities_depth_derivative[q]
                                   * (gravity * scratch.velocity_values[q]) / gravity.norm();
          const double one_over_rho = 1.0/scratch.reference_densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            data.local_rhs(i) += (pressure_scaling *
                                  one_over_rho * drho_dz_u * scratch.phi_p[i])
                                 * JxW;
        }
    }



    template <int dim>
    void
    NewtonStokesImplicitReferenceDensityCompressibilityTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      // assemble compressibility term of:
      //  - div u - 1/rho * drho/dz g/||g||* u = 0
      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile,
             ExcInternalError());

      if (!scratch.rebuild_stokes_matrix)
        return;

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i,q);
                  ++i_stokes;
                }
              ++i;
            }

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));
          const Tensor<1,dim> drho_dz = scratch.reference_densities_depth_derivative[q]
                                        * gravity / gravity.norm();
          const double one_over_rho = 1.0/scratch.reference_densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
              data.local_matrix(i,j) += (pressure_scaling *
                                         one_over_rho * drho_dz * scratch.phi_u[j] * scratch.phi_p[i])
                                        * JxW;
        }
    }



    template <int dim>
    void
    NewtonStokesIsentropicCompressionTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      // assemble RHS of:
      //  - div u = 1/rho * drho/dp rho * g * u
      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::isentropic_compression,
             ExcInternalError());

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double compressibility
            = scratch.material_model_outputs.compressibilities[q];

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            data.local_rhs(i) += (
                                   // add the term that results from the compressibility. compared
                                   // to the manual, this term seems to have the wrong sign, but this
                                   // is because we negate the entire equation to make sure we get
                                   // -div(u) as the adjoint operator of grad(p)
                                   (pressure_scaling *
                                    compressibility * density *
                                    (scratch.velocity_values[q] * gravity) *
                                    scratch.phi_p[i])
                                 )
                                 * JxW;
        }
    }



    template <int dim>
    void
    NewtonStokesProjectedDensityFieldTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      // assemble RHS of:
      // $ - \nabla \cdot \mathbf{u} = \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \frac{1}{\rho} \nabla \rho \cdot \mathbf{u}$

      // Compared to the manual, this term seems to have the wrong sign, but
      // this is because we negate the entire equation to make sure we get
      // -div(u) as the adjoint operator of grad(p)

      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::projected_density_field,
             ExcInternalError());

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();
      const unsigned int density_idx = this->introspection().find_composition_type(CompositionalFieldDescription::density);

      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();

      std::vector<Tensor<1,dim>> density_gradients(n_q_points);
      std::vector<double> density(n_q_points);
      std::vector<double> density_old(n_q_points);
      std::vector<double> density_old_old(n_q_points);

      scratch.finite_element_values[introspection.extractors.compositional_fields[density_idx]].get_function_gradients (this->get_current_linearization_point(),
          density_gradients);
      scratch.finite_element_values[introspection.extractors.compositional_fields[density_idx]].get_function_values (this->get_current_linearization_point(),
          density);
      scratch.finite_element_values[introspection.extractors.compositional_fields[density_idx]].get_function_values (this->get_old_solution(),
          density_old);
      scratch.finite_element_values[introspection.extractors.compositional_fields[density_idx]].get_function_values (this->get_old_old_solution(),
          density_old_old);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  ++i_stokes;
                }
              ++i;
            }

          double drho_dt;

          if (this->get_timestep_number() > 1)
            drho_dt = (1.0/time_step) *
                      (density[q] *
                       (2*time_step + old_time_step) / (time_step + old_time_step)
                       -
                       density_old[q] *
                       (1 + time_step/old_time_step)
                       +
                       density_old_old[q] *
                       (time_step * time_step) / (old_time_step * (time_step + old_time_step)));
          else if (this->get_timestep_number() == 1)
            drho_dt =
              (density[q] - density_old[q]) / time_step;
          else
            drho_dt = 0.0;

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            data.local_rhs(i) += (
                                   (pressure_scaling *
                                    (1.0 / density[q]) *
                                    (density_gradients[q] *
                                     scratch.velocity_values[q]) *
                                    scratch.phi_p[i])
                                   + pressure_scaling * (1.0 / density[q]) * drho_dt * scratch.phi_p[i]
                                 )
                                 * JxW;
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
  template class NewtonInterface<dim>; \
  template class NewtonStokesPreconditioner<dim>; \
  template class NewtonStokesIncompressibleTerms<dim>; \
  template class NewtonStokesCompressibleStrainRateViscosityTerm<dim>; \
  template class NewtonStokesReferenceDensityCompressibilityTerm<dim>; \
  template class NewtonStokesImplicitReferenceDensityCompressibilityTerm<dim>; \
  template class NewtonStokesIsentropicCompressionTerm<dim>; \
  template class NewtonStokesProjectedDensityFieldTerm<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

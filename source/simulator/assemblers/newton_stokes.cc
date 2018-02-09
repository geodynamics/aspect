/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    NewtonInterface<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      NewtonHandler<dim>::create_material_model_outputs(outputs);
    }

    template <int dim>
    void
    NewtonStokesPreconditioner<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;
      const double derivative_scaling_factor = this->get_newton_handler().get_newton_derivative_scaling_factor();

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

          const SymmetricTensor<4, dim> &stress_strain_director = scratch
                                                                  .material_model_outputs.stress_strain_directors[q];

          AssertThrow(stress_strain_director == dealii::identity_tensor<dim>(),
                      ExcMessage ("Error: The newton method currently doesn't support anisotropic viscosities."));

          const double JxW = scratch.finite_element_values.JxW(q);

          // TODO: Find out why in this version of ASPECT adding the derivative to the preconditioning
          // is way worse than the normal preconditioning
          if (derivative_scaling_factor == 0)
            {
              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                  if (scratch.dof_component_indices[i] ==
                      scratch.dof_component_indices[j])
                    data.local_matrix(i, j) += ((2.0 * eta * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                + one_over_eta * this->get_pressure_scaling() * this->get_pressure_scaling()
                                                * (scratch.phi_p[i] * scratch.phi_p[j]))
                                               * JxW;
            }
          else
            {

              const MaterialModel::MaterialModelDerivatives<dim> *derivatives = scratch.material_model_outputs.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

              AssertThrow(derivatives != NULL, ExcMessage ("Error: The newton method requires the derivatives"));

              const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
              const SymmetricTensor<2,dim> strain_rate = scratch.material_model_inputs.strain_rate[q];

              typename NewtonHandler<dim>::NewtonStabilization preconditioner_stabilization = this->get_newton_handler().get_preconditioner_stabilization();
              // todo: make this 0.9 into a global input parameter
              // use the spd factor when the stabilization is PD or SPD
              const double alpha = (preconditioner_stabilization | NewtonHandler<dim>::NewtonStabilization::PD) != NewtonHandler<dim>::NewtonStabilization::none ?
                                   Utilities::compute_spd_factor<dim>(eta, strain_rate, viscosity_derivative_wrt_strain_rate, 0.9)
                                   :
                                   1;

              // symmetrize when the stabilization is symmetric or SPD
              if ((preconditioner_stabilization | NewtonHandler<dim>::NewtonStabilization::symmetric) != NewtonHandler<dim>::NewtonStabilization::none)
                {
                  for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                      if (scratch.dof_component_indices[i] ==
                          scratch.dof_component_indices[j])
                        {
                          data.local_matrix(i, j) += ((2.0 * eta * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                      + derivative_scaling_factor * alpha * (scratch.grads_phi_u[i] * (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[j]) * strain_rate
                                                                                             + scratch.grads_phi_u[j] * (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[i]) * strain_rate)
                                                      + one_over_eta * this->get_pressure_scaling()
                                                      * this->get_pressure_scaling()
                                                      * (scratch.phi_p[i] * scratch
                                                         .phi_p[j]))
                                                     * JxW;
                        }
                }
              else
                {
                  for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                      if (scratch.dof_component_indices[i] ==
                          scratch.dof_component_indices[j])
                        {
                          data.local_matrix(i, j) += ((2.0 * eta * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                      + derivative_scaling_factor * alpha * 2.0 * (scratch.grads_phi_u[i] * (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[j]) * strain_rate)
                                                      + one_over_eta * this->get_pressure_scaling()
                                                      * this->get_pressure_scaling()
                                                      * (scratch.phi_p[i] * scratch
                                                         .phi_p[j]))
                                                     * JxW;
                        }
                }


            }
        }
#if DEBUG
      {
        if (this->get_newton_handler().get_preconditioner_stabilization() == NewtonHandler<dim>::NewtonStabilization::SPD)
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
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double derivative_scaling_factor = this->get_newton_handler().get_newton_derivative_scaling_factor();

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

          const SymmetricTensor<4,dim> &stress_strain_director =
            scratch.material_model_outputs.stress_strain_directors[q];

          AssertThrow(stress_strain_director == dealii::identity_tensor<dim>(),
                      ExcMessage ("Error: The newton method currently doesn't support anisotropic viscosities."));

          const SymmetricTensor<2,dim> strain_rate = scratch.material_model_inputs.strain_rate[q];
          const double pressure = scratch.material_model_inputs.pressure[q];
          const double velocity_divergence = scratch.velocity_divergence[q];

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];

          const double JxW = scratch.finite_element_values.JxW(q);

          // first assemble the rhs
          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            data.local_rhs(i) -= (eta * 2.0 * (scratch.grads_phi_u[i] * strain_rate)
                                  - (scratch.div_phi_u[i] * pressure)
                                  - (this->get_pressure_scaling() * scratch.phi_p[i] * velocity_divergence)
                                  -(density * gravity * scratch.phi_u[i]))
                                 * JxW;

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
                                                - (this->get_pressure_scaling() *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // assemble the term -div(u) as -(div u, q).
                                                // Note the negative sign to make this
                                                // operator adjoint to the grad p term:
                                                - (this->get_pressure_scaling() *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * JxW;
                  }

              // then also see whether we have to add terms due to the
              // Newton linearization
              if (derivative_scaling_factor != 0)
                {
                  const MaterialModel::MaterialModelDerivatives<dim> *derivatives
                    = scratch.material_model_outputs.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

                  // This one is only available in debug mode, because normally
                  // the AssertTrow in the preconditioner should already have
                  // caught the problem.
                  Assert(derivatives != NULL,
                         ExcMessage ("Error: The Newton method requires the material to "
                                     "compute derivatives."));

                  const SymmetricTensor<2,dim> viscosity_derivative_wrt_strain_rate = derivatives->viscosity_derivative_wrt_strain_rate[q];
                  const double viscosity_derivative_wrt_pressure = derivatives->viscosity_derivative_wrt_pressure[q];
                  typename NewtonHandler<dim>::NewtonStabilization velocity_block_stabilization = this->get_newton_handler().get_velocity_block_stabilization();
                  // todo: make this 0.9 into a global input parameter
                  // use the spd factor when the stabilization is PD or SPD
                  const double alpha =  (velocity_block_stabilization | NewtonHandler<dim>::NewtonStabilization::PD) != NewtonHandler<dim>::NewtonStabilization::none ?
                                        Utilities::compute_spd_factor<dim>(eta, strain_rate, viscosity_derivative_wrt_strain_rate, 0.9)
                                        :
                                        1;

                  // symmetrize when the stabilization is symmetric or SPD
                  if ((velocity_block_stabilization | NewtonHandler<dim>::NewtonStabilization::symmetric) != NewtonHandler<dim>::NewtonStabilization::none)
                    {
                      for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                        for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                          {
                            data.local_matrix(i,j) += ( derivative_scaling_factor * alpha * (scratch.grads_phi_u[i] * (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[j]) * strain_rate
                                                                                             + scratch.grads_phi_u[j] * (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[i]) * strain_rate)
                                                        + derivative_scaling_factor * this->get_pressure_scaling() * scratch.grads_phi_u[i] * 2.0 * viscosity_derivative_wrt_pressure * scratch.phi_p[j] * strain_rate )
                                                      * JxW;

                            Assert(dealii::numbers::is_finite(data.local_matrix(i,j)),
                                   ExcMessage ("Error: Assembly matrix is not finite." +
                                               Utilities::to_string(data.local_matrix(i,j)) +
                                               " = " + Utilities::to_string(eta)));
                          }
                    }
                  else
                    {
                      for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                        for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                          {
                            data.local_matrix(i,j) += ( derivative_scaling_factor * alpha * 2.0 * (scratch.grads_phi_u[i] * (viscosity_derivative_wrt_strain_rate * scratch.grads_phi_u[j]) * strain_rate)
                                                        + derivative_scaling_factor * this->get_pressure_scaling() * scratch.grads_phi_u[i] * 2.0 * viscosity_derivative_wrt_pressure * scratch.phi_p[j] * strain_rate )
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
          if (this->get_newton_handler().get_velocity_block_stabilization() == NewtonHandler<dim>::NewtonStabilization::SPD)
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
    NewtonStokesCompressibleStrainRateViscosityTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      if (!scratch.rebuild_stokes_matrix)
        return;



      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const double derivative_scaling_factor = this->get_newton_handler().get_newton_derivative_scaling_factor();

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

          const SymmetricTensor<4,dim> &stress_strain_director =
            scratch.material_model_outputs.stress_strain_directors[q];

          AssertThrow(stress_strain_director == dealii::identity_tensor<dim>(),
                      ExcMessage ("Error: The newton method currently doesn't support anisotropic viscosities."));

          const double velocity_divergence = trace(scratch.grads_phi_u[q]);

          const double JxW = scratch.finite_element_values.JxW(q);


          if (derivative_scaling_factor == 0)
            {
              for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += (- eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])) * JxW;
                  }
            }
          else
            {
              const MaterialModel::MaterialModelDerivatives<dim> *derivatives = scratch.material_model_outputs.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

              // This one is only avaiable in debug mode, because normally
              // the AssertTrow in the preconditioner should already have
              // caught the problem.
              Assert(derivatives != NULL, ExcMessage ("Error: The newton method requires the derivatives"));

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



    template <int dim>
    void
    NewtonStokesReferenceDensityCompressibilityTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      // assemble RHS of:
      //  - div u = 1/rho * drho/dz g/||g||* u
      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::reference_density_profile,
             ExcInternalError());

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

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
            data.local_rhs(i) += (this->get_pressure_scaling() *
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
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

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
              data.local_matrix(i,j) += (this->get_pressure_scaling() *
                                         one_over_rho * drho_dz * scratch.phi_u[j] * scratch.phi_p[i])
                                        * JxW;
        }
    }



    template <int dim>
    void
    NewtonStokesIsothermalCompressionTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      // assemble RHS of:
      //  - div u = 1/rho * drho/dp rho * g * u
      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::isothermal_compression,
             ExcInternalError());

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

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
                                   (this->get_pressure_scaling() *
                                    compressibility * density *
                                    (scratch.velocity_values[q] * gravity) *
                                    scratch.phi_p[i])
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
  template class NewtonStokesIsothermalCompressionTerm<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}

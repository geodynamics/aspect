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

#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    StokesPreconditioner<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

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

      // When using the Q1-Q1 equal order element, we need to compute the
      // projection of the Q1 pressure shape functions onto the constants
      // and use this projection in the computation of matrix terms.
      // Do this here by computing the integral of the shape functions
      // over the cell and then dividing by the area of the cell.
      std::vector<double> average_pressure_shape_function (this->get_parameters().use_equal_order_interpolation_for_stokes
                                                           ?
                                                           stokes_dofs_per_cell
                                                           :
                                                           0,
                                                           numbers::signaling_nan<double>());
      if (this->get_parameters().use_equal_order_interpolation_for_stokes)
        {
          // Check that we are really only using a Q1-Q1 element and
          // not a Q2-Q2 element. This is because in the latter case, the
          // projection isn't just on the piecewise constants, but onto
          // the piecewise (bi,tri)linears, and this is going to be a bit
          // more involved than just computing a single number per shape
          // function.
          Assert (this->get_parameters().stokes_velocity_degree==1,
                  ExcNotImplemented());

          double area       = 0;
          for (unsigned int q=0; q<n_q_points; ++q)
            area += scratch.finite_element_values.JxW(q);

          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  double int_over_p = 0;

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      int_over_p += scratch.finite_element_values[introspection.extractors.pressure].value(i,q)
                                    *
                                    scratch.finite_element_values.JxW(q);
                    }

                  average_pressure_shape_function[i_stokes] = int_over_p/area;
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
                  if (this->get_parameters().use_full_A_block_preconditioner == false)
                    scratch.grads_phi_u[i_stokes] =
                      scratch.finite_element_values[introspection.extractors
                                                    .velocities].symmetric_gradient(i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection
                                                                          .extractors.pressure].value(i, q);
                  if (this->get_parameters().use_bfbt == true)
                    {
                      scratch.grad_phi_p[i_stokes]=scratch.finite_element_values[introspection.extractors.pressure].gradient(i,q);
                      scratch.phi_u[i_stokes]=scratch.finite_element_values[introspection.extractors.velocities].value(i,q);
                    }
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];
          const double one_over_eta = 1. / eta;

          const double JxW = scratch.finite_element_values.JxW(q);

          if (this->get_parameters().use_full_A_block_preconditioner == false)
            {
              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                  if (scratch.dof_component_indices[i] ==
                      scratch.dof_component_indices[j])
                    {
                      data.local_matrix(i, j) += ((2.0 * eta * (scratch.grads_phi_u[i]
                                                                * scratch.grads_phi_u[j]))
                                                 )
                                                 * JxW;
                    }


            }
          if (this->get_parameters().use_bfbt == true)
            {
              const double sqrt_eta = std::sqrt(eta);
              const unsigned int pressure_component_index = this->introspection().component_indices.pressure;

              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                    {


                      // i and j are not pressures
                      if (scratch.dof_component_indices[i] != pressure_component_index && scratch.dof_component_indices[j] != pressure_component_index)
                        data.local_inverse_lumped_mass_matrix[i] += sqrt_eta*scalar_product(scratch.phi_u[i],scratch.phi_u[j])*JxW;


                      // i and j are pressures
                      if (scratch.dof_component_indices[i] == pressure_component_index && scratch.dof_component_indices[j] == pressure_component_index)
                        data.local_matrix(i, j) += (
                                                     1.0/sqrt_eta * pressure_scaling
                                                     * (scratch.grad_phi_p[i]
                                                        * scratch.grad_phi_p[j] + 1e-6*scratch.phi_p[i]*scratch.phi_p[j] ))
                                                   * JxW;
                    }
                }
            }
          else
            {
              for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
                  if (scratch.dof_component_indices[i] ==
                      scratch.dof_component_indices[j])
                    {
                      data.local_matrix(i, j) += (
                                                   one_over_eta * pressure_scaling
                                                   * pressure_scaling
                                                   * (scratch.phi_p[i]
                                                      * scratch.phi_p[j]))
                                                 * JxW;
                    }
            }

          // If we are using the equal order Q1-Q1 element, then we also need
          // to add the stabilization term to the (P,P) block of the matrix.
          // Note the change in sign from the one in the assembly of the
          // system matrix, which is due to the fact that the Schur complement
          // of the matrix
          //    [ A   B ]
          //    [ B^T C ]
          // is actually
          //    S  =  B^T A^{-1} B - C
          // with the minus sign in front of C. Because C is defined in the
          // method of Dohrmann and Bochev as a *negative* definite operator,
          // we here need to add the *positive* operator.
          if (this->get_parameters().use_equal_order_interpolation_for_stokes)
            {
              for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( one_over_eta * pressure_scaling * pressure_scaling *
                                                (scratch.phi_p[i] - average_pressure_shape_function[i]) *
                                                (scratch.phi_p[j] - average_pressure_shape_function[j]))
                                              * JxW;
                  }
            }
        }
    }



    template <int dim>
    void
    StokesCompressiblePreconditioner<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      Assert (this->get_parameters().use_equal_order_interpolation_for_stokes == false,
              ExcNotImplemented());

      Assert (this->get_parameters().use_full_A_block_preconditioner == false,
              ExcMessage("This assembler should only be called if the simplified A block "
                         "preconditioner is used."));

      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;

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
                  scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);

                  ++i_stokes;
                }
              ++i;
            }

          const double eta_two_thirds = scratch.material_model_outputs.viscosities[q] * 2.0 / 3.0;

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += (- eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                           )
                                           * JxW;
        }
    }



    template <int dim>
    void
    StokesIncompressibleTerms<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>();

      const MaterialModel::ElasticOutputs<dim>
      *elastic_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();

      const MaterialModel::PrescribedPlasticDilation<dim>
      *prescribed_dilation =
        (this->get_parameters().enable_prescribed_dilation)
        ? scratch.material_model_outputs.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim>>()
        : nullptr;

      const bool material_model_is_compressible = (this->get_material_model().is_compressible());

      // When using the Q1-Q1 equal order element, we need to compute the
      // projection of the Q1 pressure shape functions onto the constants
      // and use this projection in the computation of matrix terms.
      // Do this here by computing the integral of the shape functions
      // over the cell and then dividing by the area of the cell.
      std::vector<double> average_pressure_shape_function (this->get_parameters().use_equal_order_interpolation_for_stokes
                                                           ?
                                                           stokes_dofs_per_cell
                                                           :
                                                           0,
                                                           numbers::signaling_nan<double>());
      if (this->get_parameters().use_equal_order_interpolation_for_stokes)
        {
          // Check that we are really only using a Q1-Q1 element and
          // not a Q2-Q2 element. This is because in the latter case, the
          // projection isn't just on the piecewise constants, but onto
          // the piecewise (bi,tri)linears, and this is going to be a bit
          // more involved than just computing a single number per shape
          // function.
          Assert (this->get_parameters().stokes_velocity_degree==1,
                  ExcNotImplemented());

          double area       = 0;
          for (unsigned int q=0; q<n_q_points; ++q)
            area += scratch.finite_element_values.JxW(q);

          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  double int_over_p = 0;

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      int_over_p += scratch.finite_element_values[introspection.extractors.pressure].value(i,q)
                                    *
                                    scratch.finite_element_values.JxW(q);
                    }

                  average_pressure_shape_function[i_stokes] = int_over_p/area;
                  ++i_stokes;
                }
              ++i;
            }
        }

      // Next, do the integration of matrix and right hand side terms.
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  if (scratch.rebuild_stokes_matrix)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  else if (this->get_parameters().enable_elasticity)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                    }
                  else if (prescribed_dilation && !material_model_is_compressible)
                    {
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  ++i_stokes;
                }
              ++i;
            }


          // Viscosity scalar
          const double eta = ((scratch.rebuild_stokes_matrix || prescribed_dilation)
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              numbers::signaling_nan<double>());
          const double one_over_eta = (scratch.rebuild_stokes_matrix
                                       &&
                                       this->get_parameters().use_equal_order_interpolation_for_stokes
                                       ?
                                       1./eta
                                       :
                                       numbers::signaling_nan<double>());

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) += (density * gravity * scratch.phi_u[i])
                                   * JxW;

              if (force != nullptr && this->get_parameters().enable_additional_stokes_rhs)
                data.local_rhs(i) += (force->rhs_u[q] * scratch.phi_u[i]
                                      + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i])
                                     * JxW;

              if (elastic_outputs != nullptr && this->get_parameters().enable_elasticity)
                data.local_rhs(i) += (scalar_product(elastic_outputs->elastic_force[q],Tensor<2,dim>(scratch.grads_phi_u[i])))
                                     * JxW;

              if (prescribed_dilation != nullptr)
                data.local_rhs(i) += (
                                       // RHS of - (div u,q) = - (R,q)
                                       - pressure_scaling
                                       * prescribed_dilation->dilation[q]
                                       * scratch.phi_p[i]
                                     ) * JxW;

              // Only assemble this term if we are running incompressible, otherwise this term
              // is already included on the LHS of the equation.
              if (prescribed_dilation != nullptr && !material_model_is_compressible)
                data.local_rhs(i) += (
                                       // RHS of momentum eqn: - \int 2/3 eta R, div v
                                       - 2.0 / 3.0 * eta
                                       * prescribed_dilation->dilation[q]
                                       * scratch.div_phi_u[i]
                                     ) * JxW;

              if (scratch.rebuild_stokes_matrix)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( (eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                // assemble \nabla p as -(p, div v):
                                                - (pressure_scaling *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // assemble the term -div(u) as -(div u, q).
                                                // Note the negative sign to make this
                                                // operator adjoint to the grad p term:
                                                - (pressure_scaling *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * JxW;
                  }
            }

          // If we are using the equal order Q1-Q1 element, then we also need
          // to put the stabilization term into the (P,P) block of the matrix:
          if (scratch.rebuild_stokes_matrix
              &&
              this->get_parameters().use_equal_order_interpolation_for_stokes)
            {
              for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( - (one_over_eta * pressure_scaling * pressure_scaling *
                                                   (scratch.phi_p[i] - average_pressure_shape_function[i]) *
                                                   (scratch.phi_p[j] - average_pressure_shape_function[j])))
                                              * JxW;
                  }
            }
        }
    }



    template <int dim>
    void
    StokesIncompressibleTerms<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.n_evaluation_points();

      // Stokes RHS:
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

      // Elasticity:
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
    }



    template <int dim>
    void
    StokesCompressibleStrainRateViscosityTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

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
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);

                  ++i_stokes;
                }
              ++i;
            }

          // Viscosity scalar
          const double eta_two_thirds = scratch.material_model_outputs.viscosities[q] * 2.0 / 3.0;

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
              {
                data.local_matrix(i,j) += (- (eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                             ))
                                          * JxW;
              }
        }
    }



    template <int dim>
    void
    StokesReferenceDensityCompressibilityTerm<dim>::
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
    StokesImplicitReferenceDensityCompressibilityTerm<dim>::
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
              data.local_matrix(i,j) += -(pressure_scaling *
                                          one_over_rho * drho_dz * scratch.phi_u[j] * scratch.phi_p[i])
                                        * JxW;
        }
    }



    template <int dim>
    void
    StokesIsentropicCompressionTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      // assemble RHS of:
      //  - div \mathbf{u} = \frac{1}{\rho} \frac{\partial rho}{\partial p} \rho \mathbf{g} \cdot \mathbf{u}

      // Compared to the manual, this term seems to have the wrong sign, but
      // this is because we negate the entire equation to make sure we get
      // -div(u) as the adjoint operator of grad(p)

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
    StokesProjectedDensityFieldTerm<dim>::
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



    template <int dim>
    void
    StokesHydrostaticCompressionTerm<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      // assemble RHS of:
      // $ -\nabla \cdot \mathbf{u} = \left( \kappa \rho \textbf{g} - \alpha \nabla T \right) \cdot \textbf{u}$
      //
      // where $\frac{1}{\rho} \frac{\partial \rho}{\partial p} = \kappa$ is the compressibility,
      // $- \frac{1}{\rho}\frac{\partial \rho}{\partial T} = \alpha$ is the thermal expansion coefficient,
      // and both are defined in the material model.

      // Compared to the manual, this term seems to have the wrong sign, but
      // this is because we negate the entire equation to make sure we get
      // -div(u) as the adjoint operator of grad(p)

      Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::hydrostatic_compression,
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

          const double thermal_alpha
            = scratch.material_model_outputs.thermal_expansion_coefficients[q];

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            data.local_rhs(i) += (
                                   (pressure_scaling *
                                    (
                                      // pressure part:
                                      compressibility * density *
                                      (scratch.velocity_values[q] * gravity)
                                      // temperature part:
                                      - thermal_alpha *
                                      (scratch.velocity_values[q] * scratch.temperature_gradients[q])
                                    ) * scratch.phi_p[i])
                                 )
                                 * JxW;

        }
    }


    template <int dim>
    void
    StokesPressureRHSCompatibilityModification<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                                              internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();

      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          const double JxW = scratch.finite_element_values.JxW(q);
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  data.local_pressure_shape_function_integrals(i_stokes) += scratch.phi_p[i_stokes] * JxW;
                  ++i_stokes;
                }
              ++i;
            }
        }
    }



    template <int dim>
    void
    StokesBoundaryTraction<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                                          internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();

      // see if any of the faces are traction boundaries for which
      // we need to assemble force terms for the right hand side
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();

      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(scratch.face_number);

      const auto &traction_bis = this->get_boundary_traction_manager().get_prescribed_boundary_traction_indicators();

      if (traction_bis.find(face->boundary_id()) != traction_bis.end())
        {
          for (unsigned int q=0; q<scratch.face_finite_element_values.n_quadrature_points; ++q)
            {
              const Tensor<1,dim> traction
                = this->get_boundary_traction_manager().
                  boundary_traction (face->boundary_id(),
                                     scratch.face_finite_element_values.quadrature_point(q),
                                     scratch.face_finite_element_values.normal_vector(q));

              const double JxW = scratch.face_finite_element_values.JxW(q);

              for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
                {
                  if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                    {
                      data.local_rhs(i_stokes) += scratch.face_finite_element_values[introspection.extractors.velocities].value(i,q) *
                                                  traction * JxW;
                      ++i_stokes;
                    }
                  ++i;
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
  template class StokesPreconditioner<dim>; \
  template class StokesCompressiblePreconditioner<dim>; \
  template class StokesIncompressibleTerms<dim>; \
  template class StokesCompressibleStrainRateViscosityTerm<dim>; \
  template class StokesReferenceDensityCompressibilityTerm<dim>; \
  template class StokesImplicitReferenceDensityCompressibilityTerm<dim>; \
  template class StokesIsentropicCompressionTerm<dim>; \
  template class StokesHydrostaticCompressionTerm<dim>; \
  template class StokesProjectedDensityFieldTerm<dim>; \
  template class StokesPressureRHSCompatibilityModification<dim>; \
  template class StokesBoundaryTraction<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

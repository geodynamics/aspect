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
#include <aspect/utilities.h>
#include <aspect/assembly.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    StokesAssembler<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >() == NULL)
        {
          outputs.additional_outputs.push_back(
            std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >
            (new MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> (n_points)));
        }
      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
             == n_points, ExcInternalError());
    }

    template <int dim>
    void
    StokesAssembler<dim>::
    preconditioner (const double                                             pressure_scaling,
                    internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                    internal::Assembly::CopyData::StokesPreconditioner<dim> &data) const
    {
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
                  scratch.grads_phi_u[i_stokes] =
                    scratch.finite_element_values[introspection.extractors
                                                  .velocities].symmetric_gradient(i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection
                                                                          .extractors.pressure].value(i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];
          const double one_over_eta = 1. / eta;

          const SymmetricTensor<4, dim> &stress_strain_director = scratch
                                                                  .material_model_outputs.stress_strain_directors[q];
          const bool use_tensor = (stress_strain_director
                                   != dealii::identity_tensor<dim>());

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += ((
                                              use_tensor ?
                                              eta * (scratch.grads_phi_u[i]
                                                     * stress_strain_director
                                                     * scratch.grads_phi_u[j]) :
                                              eta * (scratch.grads_phi_u[i]
                                                     * scratch.grads_phi_u[j]))
                                            + one_over_eta * pressure_scaling
                                            * pressure_scaling
                                            * (scratch.phi_p[i] * scratch
                                               .phi_p[j]))
                                           * JxW;
        }
    }



    template <int dim>
    void
    StokesAssembler<dim>::
    incompressible_terms (const double                                     pressure_scaling,
                          const bool                                       rebuild_stokes_matrix,
                          internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                          internal::Assembly::CopyData::StokesSystem<dim> &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  if (rebuild_stokes_matrix)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  ++i_stokes;
                }
              ++i;
            }


          // Viscosity scalar
          const double eta = (rebuild_stokes_matrix
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              numbers::signaling_nan<double>());

          const SymmetricTensor<4,dim> &stress_strain_director =
            scratch.material_model_outputs.stress_strain_directors[q];
          const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) += (density * gravity * scratch.phi_u[i])
                                   * JxW;

              if (force != NULL)
                data.local_rhs(i) += (force->rhs_u[q] * scratch.phi_u[i]
                                      + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i])
                                     * JxW;

              if (rebuild_stokes_matrix)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( (use_tensor ?
                                                 eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                                 :
                                                 eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
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
        }
    }



    template <int dim>
    void
    StokesAssembler<dim>::
    compressible_strain_rate_viscosity_term (const double                                     /*pressure_scaling*/,
                                             const bool                                       rebuild_stokes_matrix,
                                             internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                             internal::Assembly::CopyData::StokesSystem<dim> &data) const
    {
      if (!rebuild_stokes_matrix)
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
                  scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);

                  ++i_stokes;
                }
              ++i;
            }

          // Viscosity scalar
          const double eta_two_thirds = scratch.material_model_outputs.viscosities[q] * 2.0 / 3.0;

          const SymmetricTensor<4,dim> &stress_strain_director =
            scratch.material_model_outputs.stress_strain_directors[q];
          const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
              {
                data.local_matrix(i,j) += (- (use_tensor ?
                                              eta_two_thirds * (scratch.div_phi_u[i] * trace(stress_strain_director * scratch.grads_phi_u[j]))
                                              :
                                              eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                             ))
                                          * JxW;
              }
        }
    }



    template <int dim>
    void
    StokesAssembler<dim>::
    reference_density_compressibility_term (const double                                     pressure_scaling,
                                            const bool                                       /*rebuild_stokes_matrix*/,
                                            internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                            internal::Assembly::CopyData::StokesSystem<dim> &data,
                                            const Parameters<dim>                           &parameters) const
    {
      // assemble RHS of:
      //  - div u = 1/rho * drho/dz g/||g||* u
      (void)parameters;
      Assert(parameters.formulation_mass_conservation ==
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
            data.local_rhs(i) += (pressure_scaling *
                                  one_over_rho * drho_dz_u * scratch.phi_p[i])
                                 * JxW;
        }
    }



    template <int dim>
    void
    StokesAssembler<dim>::
    implicit_reference_density_compressibility_term (const double                                     pressure_scaling,
                                                     const bool                                       rebuild_stokes_matrix,
                                                     internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                     internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                     const Parameters<dim>                           &parameters) const
    {
      // assemble compressibility term of:
      //  - div u - 1/rho * drho/dz g/||g||* u = 0
      (void)parameters;
      Assert(parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile,
             ExcInternalError());

      if (!rebuild_stokes_matrix)
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
              data.local_matrix(i,j) += -(pressure_scaling *
                                          one_over_rho * drho_dz * scratch.phi_u[j] * scratch.phi_p[i])
                                        * JxW;
        }
    }



    template <int dim>
    void
    StokesAssembler<dim>::
    isothermal_compression_term (const double                                     pressure_scaling,
                                 const bool                                       /*rebuild_stokes_matrix*/,
                                 internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                 internal::Assembly::CopyData::StokesSystem<dim> &data,
                                 const Parameters<dim>                           &parameters) const
    {
      // assemble RHS of:
      //  - div u = 1/rho * drho/dp rho * g * u
      (void)parameters;
      Assert(parameters.formulation_mass_conservation ==
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
                                   (pressure_scaling *
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
  template class \
  StokesAssembler<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}

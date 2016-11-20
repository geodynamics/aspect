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

#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <aspect/assembly.h>
#include <aspect/assemblers.h>
#include <aspect/simulator_access.h>

namespace aspect
{

  template <int dim>
  void
  StokesAssembler<dim>::
  local_assemble_stokes_preconditioner (const double                                             pressure_scaling,
                                        internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                                        internal::Assembly::CopyData::StokesPreconditioner<dim> &data) const
  {
    const Introspection<dim> &introspection = this->introspection();
    const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();
    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = scratch.finite_element_values.n_quadrature_points;

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          {
            scratch.grads_phi_u[k] =
              scratch.finite_element_values[introspection.extractors
                                            .velocities].symmetric_gradient(k, q);
            scratch.phi_p[k] = scratch.finite_element_values[introspection
                                                             .extractors.pressure].value(k, q);
          }
        const double eta = scratch.material_model_outputs.viscosities[q];
        const SymmetricTensor<4, dim> &stress_strain_director = scratch
                                                                .material_model_outputs.stress_strain_directors[q];
        const bool use_tensor = (stress_strain_director
                                 != dealii::identity_tensor<dim>());
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            if (fe.system_to_component_index(i).first ==
                fe.system_to_component_index(j).first)
              data.local_matrix(i, j) += ((
                                            use_tensor ?
                                            eta * (scratch.grads_phi_u[i]
                                                   * stress_strain_director
                                                   * scratch.grads_phi_u[j]) :
                                            eta * (scratch.grads_phi_u[i]
                                                   * scratch.grads_phi_u[j]))
                                          + (1. / eta) * pressure_scaling
                                          * pressure_scaling
                                          * (scratch.phi_p[i] * scratch
                                             .phi_p[j]))
                                         * scratch.finite_element_values.JxW(
                                           q);
      }
  }


  template <int dim>
  void
  StokesAssembler<dim>::
  local_assemble_stokes_incompressible (const double                                     pressure_scaling,
                                        const bool                                       rebuild_stokes_matrix,
                                        internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                        internal::Assembly::CopyData::StokesSystem<dim> &data) const
  {
    const Introspection<dim> &introspection = this->introspection();
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].value (k,q);
            if (rebuild_stokes_matrix)
              {
                scratch.phi_p[k] = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
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

        const SymmetricTensor<4,dim> &stress_strain_director =
          scratch.material_model_outputs.stress_strain_directors[q];
        const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

        const Tensor<1,dim>
        gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

        if (rebuild_stokes_matrix)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              data.local_matrix(i,j) += (
                                          // first assemble the symmetric e(u) : e(v) term:
                                          (use_tensor ?
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
                                             scratch.phi_p[i] * scratch.div_phi_u[j])
                                        )
                                        * scratch.finite_element_values.JxW(q);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs(i) += (
                                 // buoyancy term:
                                 (scratch.material_model_outputs.densities[q] * gravity * scratch.phi_u[i])
                               )
                               * scratch.finite_element_values.JxW(q);
      }
  }

  template <int dim>
  void
  StokesAssembler<dim>::
  local_assemble_stokes_compressible_diffusion (const double                                     /*pressure_scaling*/,
                                                const bool                                       rebuild_stokes_matrix,
                                                internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                internal::Assembly::CopyData::StokesSystem<dim> &data) const
  {
    if (!rebuild_stokes_matrix)
      return;

    const Introspection<dim> &introspection = this->introspection();
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grads_phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(k,q);
            scratch.div_phi_u[k]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (k, q);
          }

        // Viscosity scalar
        const double eta = scratch.material_model_outputs.viscosities[q];

        const SymmetricTensor<4,dim> &stress_strain_director =
          scratch.material_model_outputs.stress_strain_directors[q];
        const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            data.local_matrix(i,j) += (  - (use_tensor ?
                                            eta * 2.0/3.0 * (scratch.div_phi_u[i] * trace(stress_strain_director * scratch.grads_phi_u[j]))
                                            :
                                            eta * 2.0/3.0 * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                           )
                                      )
                                      * scratch.finite_element_values.JxW(q);
      }
  }


  template <int dim>
  void
  StokesAssembler<dim>::
  local_assemble_stokes_mass_density_gradient (const double                                     pressure_scaling,
                                               const bool                                       /*rebuild_stokes_matrix*/,
                                               internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                               internal::Assembly::CopyData::StokesSystem<dim> &data,
                                               const Parameters<dim> &parameters) const
  {
    // assemble RHS of:
    //  - div u = 1/rho * drho/dz g/||g||* u

    Assert(parameters.formulation_mass == Parameters<dim>::FormulationType::adiabatic,
           ExcInternalError());

    const Introspection<dim> &introspection = this->introspection();
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_p[k] = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
          }

        const Tensor<1,dim>
        gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));
        const double drho_dz_u = scratch.adiabatic_density_gradients[q]
                                 * (gravity * scratch.velocity_values[q]) / gravity.norm();
        const double one_over_rho = 1.0/scratch.mass_densities[q];

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs(i) += (pressure_scaling *
                                one_over_rho * drho_dz_u * scratch.phi_p[i])
                               * scratch.finite_element_values.JxW(q);
      }
  }

  template <int dim>
  void
  StokesAssembler<dim>::
  local_assemble_stokes_mass_density_implicit (const double                                     pressure_scaling,
                                               const bool                                       rebuild_stokes_matrix,
                                               internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                               internal::Assembly::CopyData::StokesSystem<dim> &data,
                                               const Parameters<dim> &parameters) const
  {
    Assert(parameters.formulation_mass == Parameters<dim>::FormulationType::implicit_adiabatic,
           ExcInternalError());

    const Introspection<dim> &introspection = this->introspection();
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].value (k,q);
            scratch.phi_p[k] = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
          }

        const Tensor<1,dim>
        gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

        const double compressibility
          = scratch.material_model_outputs.compressibilities[q];
        const double mass_density = scratch.mass_densities[q];

        if (rebuild_stokes_matrix)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              data.local_matrix(i,j) += (
                                          - (pressure_scaling * compressibility * mass_density
                                             *(scratch.phi_u[j] * gravity)
                                             * scratch.phi_p[i])
                                        )
                                        * scratch.finite_element_values.JxW(q);
      }

  }

  template <int dim>
  void
  StokesAssembler<dim>::
  local_assemble_stokes_mass_density_explicit (const double                                     pressure_scaling,
                                               const bool                                       /*rebuild_stokes_matrix*/,
                                               internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                               internal::Assembly::CopyData::StokesSystem<dim> &data,
                                               const Parameters<dim> &parameters) const
  {
    Assert(parameters.formulation_mass != Parameters<dim>::FormulationType::implicit_adiabatic,
           ExcInternalError());

    const Introspection<dim> &introspection = this->introspection();
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_p[k] = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
          }

        const Tensor<1,dim>
        gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

        const double compressibility
          = scratch.material_model_outputs.compressibilities[q];

        const double density = scratch.mass_densities[q];

        for (unsigned int i=0; i<dofs_per_cell; ++i)
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
                               * scratch.finite_element_values.JxW(q);
      }
  }

} // namespace aspect

// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class \
  StokesAssembler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}

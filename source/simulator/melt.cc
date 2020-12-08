/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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


#include <aspect/melt.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <aspect/citation_info.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator/assemblers/advection.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    MeltInputs<dim>::MeltInputs (const unsigned int n_points)
      :
      compaction_pressures(n_points, numbers::signaling_nan<double>()),
      fluid_velocities(n_points, numbers::signaling_nan<Tensor<1,dim> >())
    {}


    template <int dim>
    void MeltInputs<dim>::fill (const LinearAlgebra::BlockVector &solution,
                                const FEValuesBase<dim>          &fe_values,
                                const Introspection<dim>         &introspection)
    {
      compaction_pressures.resize(fe_values.n_quadrature_points);
      fluid_velocities.resize(fe_values.n_quadrature_points, Tensor<1,dim>());

      const FEValuesExtractors::Vector ex_u_f = introspection.variable("fluid velocity").extractor_vector();
      fe_values[ex_u_f].get_function_values (solution, fluid_velocities);

      const FEValuesExtractors::Scalar ex_p_c = introspection.variable("compaction pressure").extractor_scalar();
      fe_values[ex_p_c].get_function_values (solution, compaction_pressures);
    }

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

    template <int dim>
    double
    MeltInterface<dim>::p_c_scale (const MaterialModel::MaterialModelInputs<dim> &inputs,
                                   const MaterialModel::MaterialModelOutputs<dim> &outputs,
                                   const MeltHandler<dim> &melt_handler,
                                   const bool consider_is_melt_cell) const
    {
      Assert(inputs.current_cell.state() == IteratorState::valid,
             ExcMessage("You are trying to access the property of a cell (if it is "
                        "a melt cell or not), but the cell iterator in the provided material "
                        "model inputs does not point to a cell."));

      if (consider_is_melt_cell && !melt_handler.is_melt_cell(inputs.current_cell))
        return 0.0;

      const MaterialModel::MeltOutputs<dim> *melt_outputs = outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
      const double ref_K_D = this->reference_darcy_coefficient();

      double K_D = 0.0;
      double max_K_D = 0.0;
      const unsigned int N = melt_outputs->permeabilities.size();
      for (unsigned int q=0; q<N; ++q)
        {
          K_D += std::max(melt_outputs->permeabilities[q], 0.0) / (melt_outputs->fluid_viscosities[q] * N);
          max_K_D = std::max(max_K_D, melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q]);
        }

      // For melt cells, we return the average Darcy coefficient of the cell,
      // but always a value larger than a threshold depending on the reference Darcy coefficient,
      // so that K_D can never be zero (or close to zero) in melt cells.
      // The same threshold is used when computing which cells are melt cells (the else branch),
      // with the difference that we return a p_c_scale of zero to indicate that the cell is not
      // a melt cell if the maximum Darcy coefficient of the cell is below the threshold.
      const double K_D_threshold = melt_handler.melt_parameters.melt_scaling_factor_threshold;
      if (consider_is_melt_cell)
        K_D = std::max(K_D, K_D_threshold * ref_K_D);
      else
        K_D = (max_K_D < K_D_threshold * ref_K_D) ? 0.0 : K_D;

      // If the reference permeability is set to zero, there is no melt transport in the whole model and we return zero.
      return (ref_K_D > 0 ? std::sqrt(K_D / ref_K_D) : 0.0);
    }
  }

  namespace Assemblers
  {
    namespace
    {
      template <int dim>
      bool
      is_velocity_or_pressures (const Introspection<dim> &introspection,
                                const unsigned int p_c_component_index,
                                const unsigned int p_f_component_index,
                                const unsigned int component_index)
      {
        if (component_index == p_c_component_index)
          return true;

        if (component_index == p_f_component_index)
          return true;

        for (unsigned int i=0; i<dim; ++i)
          if (component_index == introspection.component_indices.velocities[i])
            return true;

        return false;
      }
    }

    template <int dim>
    void
    MeltInterface<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      MeltHandler<dim>::create_material_model_outputs(outputs);

      const unsigned int n_points = outputs.viscosities.size();

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >(n_points));
        }

      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
             == n_points, ExcInternalError());
    }


    template <int dim>
    void
    MeltStokesPreconditioner<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int   stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int   n_q_points      = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      Assert(Plugins::plugin_type_matches<const MaterialModel::MeltInterface<dim>>(this->get_material_model()),
             ExcMessage("Error: The current material model needs to be derived from MeltInterface to use melt transport."));

      const double p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(
                                 this->get_material_model()).p_c_scale(scratch.material_model_inputs,
                                                                       scratch.material_model_outputs,
                                                                       this->get_melt_handler(),
                                                                       true);

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();
      const FEValuesExtractors::Scalar ex_p_c = introspection.variable("compaction pressure").extractor_scalar();

      const unsigned int p_f_component_index = introspection.variable("fluid pressure").first_component_index;
      const unsigned int p_c_component_index = introspection.variable("compaction pressure").first_component_index;

      for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
        {
          const unsigned int component_index_i = fe.system_to_component_index(i).first;

          if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
            {
              data.local_dof_indices[i_stokes] = scratch.local_dof_indices[i];
              scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
              ++i_stokes;
            }
          ++i;
        }

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              const unsigned int component_index_i = fe.system_to_component_index(i).first;

              if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
                {
                  scratch.phi_p[i_stokes]       = scratch.finite_element_values[ex_p_f].value (i, q);
                  scratch.phi_p_c[i_stokes]     = scratch.finite_element_values[ex_p_c].value (i, q);
                  scratch.grad_phi_p[i_stokes]  = scratch.finite_element_values[ex_p_f].gradient (i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];
          const double one_over_eta = 1. / eta;
          const double K_D = this->get_melt_handler().limited_darcy_coefficient(melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q],
                                                                                p_c_scale > 0);

          /*
            - R = 1/eta M_p + K_D L_p for p
            S = - (1/eta + 1/viscosity_c)  M_p  for p_c
          */
          const double viscosity_c = melt_outputs->compaction_viscosities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
              {
                if (scratch.dof_component_indices[i] == scratch.dof_component_indices[j])
                  data.local_matrix(i,j) += ((one_over_eta *
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
                                             (1./eta + 1./viscosity_c) * p_c_scale * p_c_scale *
                                             pressure_scaling *
                                             pressure_scaling *
                                             (scratch.phi_p_c[i] * scratch.phi_p_c[j])
                                            )
                                            * JxW;

                // add S between p_c and p_f
                data.local_matrix(i,j) +=
                  (
                    (p_c_scale * one_over_eta *
                     pressure_scaling *
                     pressure_scaling)
                    * scratch.phi_p[i] * scratch.phi_p_c[j]
                    +
                    (p_c_scale * one_over_eta *
                     pressure_scaling *
                     pressure_scaling)
                    * scratch.phi_p_c[i] * scratch.phi_p[j]
                  )
                  * JxW;
              }
        }
    }

    namespace
    {
      /**
       * Compute the right-hand side of the fluid pressure equation of the Stokes
       * system in case the simulation uses melt transport. This includes a term
       * derived from Darcy's law, a term including the melting rate and a term dependent
       * on the densities and velocities of fluid and solid.
       */
      template <int dim>
      double
      compute_fluid_pressure_rhs(const SimulatorAccess<dim> *simulator_access,
                                 const internal::Assembly::Scratch::StokesSystem<dim> &scratch,
                                 const unsigned int q_point,
                                 const double K_D,
                                 const double operator_split_reaction)
      {
        if (!simulator_access->get_parameters().include_melt_transport)
          return 0.0;

        const Introspection<dim> &introspection = simulator_access->introspection();
        const MaterialModel::MeltOutputs<dim> *melt_out = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

        Assert(melt_out != nullptr, ExcInternalError());

        const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
        const unsigned int is_compressible = simulator_access->get_material_model().is_compressible();

        const double solid_density    = scratch.material_model_outputs.densities[q_point];
        const double fluid_density    = melt_out->fluid_densities[q_point];
        double melting_rate           = scratch.material_model_outputs.reaction_terms[q_point][porosity_index];

        if (simulator_access->get_parameters().use_operator_splitting)
          melting_rate = (simulator_access->get_timestep() > 0
                          ?
                          operator_split_reaction * solid_density / simulator_access->get_timestep()
                          :
                          0.0);

        const double solid_compressibility = scratch.material_model_outputs.compressibilities[q_point];
        const Tensor<1,dim> fluid_density_gradient = melt_out->fluid_density_gradients[q_point];
        const Tensor<1,dim> current_u = scratch.velocity_values[q_point];
        const double porosity         = std::max(scratch.material_model_inputs.composition[q_point][porosity_index],0.0);

        const Tensor<1,dim>
        gravity = simulator_access->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q_point));

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
    }



    template <int dim>
    void
    MeltStokesSystem<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      Assert(Plugins::plugin_type_matches<const MaterialModel::MeltInterface<dim>>(this->get_material_model()),
             ExcMessage("Error: The current material model needs to be derived from MeltInterface to use melt transport."));

      const double p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>
                               (this->get_material_model()).p_c_scale(scratch.material_model_inputs,
                                                                      scratch.material_model_outputs,
                                                                      this->get_melt_handler(),
                                                                      true);

      const FEValuesExtractors::Scalar extractor_pressure = introspection.variable("fluid pressure").extractor_scalar();
      const FEValuesExtractors::Scalar ex_p_c = introspection.variable("compaction pressure").extractor_scalar();

      const unsigned int p_f_component_index = introspection.variable("fluid pressure").first_component_index;
      const unsigned int p_c_component_index = introspection.variable("compaction pressure").first_component_index;

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

      const double pressure_scaling = this->get_pressure_scaling();

      std::vector<double> reactions(n_q_points, numbers::signaling_nan<double>());
      const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
      if (this->get_parameters().use_operator_splitting)
        scratch.finite_element_values[introspection.extractors.compositional_fields[porosity_index]].get_function_values(this->get_reaction_vector(),
            reactions);

      for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
        {
          const unsigned int component_index_i = fe.system_to_component_index(i).first;

          if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
            {
              data.local_dof_indices[i_stokes] = scratch.local_dof_indices[i];
              ++i_stokes;
            }
          ++i;
        }

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              const unsigned int component_index_i = fe.system_to_component_index(i).first;

              if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
                {
                  scratch.phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes]   = scratch.finite_element_values[extractor_pressure].value (i, q);
                  scratch.phi_p_c[i_stokes] = scratch.finite_element_values[ex_p_c].value (i, q);
                  scratch.grad_phi_p[i_stokes] = scratch.finite_element_values[extractor_pressure].gradient (i, q);

                  if (scratch.rebuild_stokes_matrix)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  ++i_stokes;
                }
              ++i;
            }

          // Viscosity scalar
          const double eta = (scratch.rebuild_stokes_matrix
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              numbers::signaling_nan<double>());
          const double eta_two_thirds = (scratch.rebuild_stokes_matrix
                                         ?
                                         scratch.material_model_outputs.viscosities[q] * 2.0 / 3.0
                                         :
                                         numbers::signaling_nan<double>());
          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double compressibility
            = (this->get_material_model().is_compressible()
               ?
               scratch.material_model_outputs.compressibilities[q]
               :
               numbers::signaling_nan<double>() );

          const double density_s = scratch.material_model_outputs.densities[q]; // density of the solid

          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.material_model_inputs.composition[q][porosity_index],0.000);

          const double K_D = this->get_melt_handler().limited_darcy_coefficient(melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q],
                                                                                p_c_scale > 0);
          const double viscosity_c = melt_outputs->compaction_viscosities[q];
          const Tensor<1,dim> density_gradient_f = melt_outputs->fluid_density_gradients[q];
          const double density_f = melt_outputs->fluid_densities[q];
          const double p_f_RHS = compute_fluid_pressure_rhs(this,
                                                            scratch,
                                                            q,
                                                            K_D,
                                                            reactions[q]);
          const double bulk_density = (1.0 - porosity) * density_s + porosity * density_f;

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) += (
                                     (bulk_density * gravity * scratch.phi_u[i])
                                     +
                                     // add a force to the RHS if present
                                     (force!=nullptr ?
                                      (force->rhs_u[q] * scratch.phi_u[i]
                                       + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i]
                                       + pressure_scaling * force->rhs_melt_pc[q] * p_c_scale * scratch.phi_p_c[i])
                                      : 0.0)
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
                                   * JxW;

              if (scratch.rebuild_stokes_matrix)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( (eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                - (eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                                  )
                                                - (pressure_scaling *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // finally the term -div(u). note the negative sign to make this
                                                // operator adjoint to the grad(p) term
                                                - (pressure_scaling *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j])
                                                +
                                                (- pressure_scaling * pressure_scaling / viscosity_c * p_c_scale * p_c_scale
                                                 * scratch.phi_p_c[i] * scratch.phi_p_c[j])
                                                - pressure_scaling * scratch.div_phi_u[i] * scratch.phi_p_c[j] * p_c_scale
                                                - pressure_scaling * scratch.phi_p_c[i] * scratch.div_phi_u[j] * p_c_scale
                                                - K_D * pressure_scaling * pressure_scaling *
                                                (scratch.grad_phi_p[i] * scratch.grad_phi_p[j])
                                                + (this->get_material_model().is_compressible()
                                                   ?
                                                   K_D * pressure_scaling * pressure_scaling / density_f
                                                   * scratch.phi_p[i] * (scratch.grad_phi_p[j] * density_gradient_f)
                                                   :
                                                   0.0))
                                              * JxW;
                  }
            }
        }

    }



    template <int dim>
    void
    MeltStokesSystemBoundary<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int n_face_q_points = scratch.face_finite_element_values.n_quadrature_points;
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();
      const unsigned int p_f_component_index = introspection.variable("fluid pressure").first_component_index;
      const unsigned int p_c_component_index = introspection.variable("compaction pressure").first_component_index;
      const bool is_melt_cell = this->get_melt_handler().is_melt_cell(scratch.material_model_inputs.current_cell);

      const typename DoFHandler<dim>::face_iterator face = scratch.face_material_model_inputs.current_cell->face(scratch.face_number);

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.face_material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      std::vector<double> grad_p_f(n_face_q_points);
      this->get_melt_handler().get_boundary_fluid_pressure().fluid_pressure_gradient(
        face->boundary_id(),
        scratch.face_material_model_inputs,
        scratch.face_material_model_outputs,
        scratch.face_finite_element_values.get_normal_vectors(),
        grad_p_f);

      for (unsigned int q=0; q<n_face_q_points; ++q)
        {
          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.face_finite_element_values.quadrature_point(q));
          const double density_f = melt_outputs->fluid_densities[q];
          const double K_D = this->get_melt_handler().limited_darcy_coefficient(melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q],
                                                                                is_melt_cell);

          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              const unsigned int component_index_i = fe.system_to_component_index(i).first;

              if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
                {
                  // apply the fluid pressure boundary condition
                  data.local_rhs(i_stokes) += (scratch.face_finite_element_values[ex_p_f].value(i, q)
                                               * this->get_pressure_scaling() * K_D *
                                               (density_f
                                                * (scratch.face_finite_element_values.normal_vector(q) * gravity)
                                                - grad_p_f[q])
                                               * scratch.face_finite_element_values.JxW(q));
                  ++i_stokes;
                }
              ++i;
            }
        }
    }

    namespace
    {
      /**
       * Compute the right-hand side for the porosity system.
       * It includes the melting rate and a term dependent
       * on the density and velocity.
       *
       * This function is implemented in
       * <code>source/simulator/melt.cc</code>.
       */
      template <int dim>
      double
      compute_melting_RHS(const SimulatorAccess<dim> *simulator_access,
                          const internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                          const unsigned int q_point,
                          const double divergence_u)
      {
        Assert (scratch.material_model_outputs.densities[q_point] > 0,
                ExcMessage ("The density needs to be a positive quantity "
                            "when melt transport is included in the simulation."));

        const double melting_rate         = scratch.material_model_outputs.reaction_terms[q_point][scratch.advection_field->compositional_variable];
        const double density              = scratch.material_model_outputs.densities[q_point];

        const double compressibility      = (simulator_access->get_material_model().is_compressible()
                                             ?
                                             scratch.material_model_outputs.compressibilities[q_point]
                                             :
                                             0.0);
        const Tensor<1,dim> current_u     = scratch.current_velocity_values[q_point];
        const Tensor<1,dim>
        gravity = simulator_access->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q_point));

        // Note that the full advection term would be (1 - phi) * (div u + kappa rho u g),
        // but we expanded the term and move the part (- phi) * (div u + kappa rho u g) over to the LHS.
        double melt_transport_RHS = melting_rate / density
                                    + divergence_u + compressibility * density * (current_u * gravity);

        if (!simulator_access->get_melt_handler().is_melt_cell(scratch.material_model_inputs.current_cell))
          melt_transport_RHS = melting_rate / density;

        return melt_transport_RHS;
      }
    }


    template <int dim>
    void
    MeltAdvectionSystem<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::AdvectionSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::AdvectionSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const bool use_bdf2_scheme = (this->get_timestep_number() > 1);
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();

      const unsigned int solution_component = scratch.advection_field->component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = scratch.advection_field->scalar_extractor(introspection);

      MaterialModel::MeltOutputs<dim> *melt_outputs = scratch.material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();

      Assert(melt_outputs->compaction_viscosities[0] > 0.0,
             ExcMessage ("MeltOutputs have to be filled for models with melt transport. "
                         "At the moment, these outputs are not filled, or they do not have "
                         "reasonable values."));

      std::vector<Tensor<1,dim> > fluid_velocity_values(n_q_points);
      const FEValuesExtractors::Vector ex_u_f = introspection.variable("fluid velocity").extractor_vector();
      scratch.finite_element_values[ex_u_f].get_function_values (this->get_solution(),fluid_velocity_values);

      // average divergence u over the cell (needed for porosity advection)
      double divergence_u = 0.0;
      if (this->get_melt_handler().is_porosity(*scratch.advection_field))
        for (unsigned int q=0; q<n_q_points; ++q)
          divergence_u += scratch.current_velocity_divergences[q] * 1./n_q_points;

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
                  scratch.grad_phi_field[i_advection] = scratch.finite_element_values[solution_field].gradient (i,q);
                  scratch.phi_field[i_advection]      = scratch.finite_element_values[solution_field].value (i,q);
                  ++i_advection;
                }
              ++i;
            }

          const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
          const double porosity = std::max(scratch.material_model_inputs.composition[q][porosity_index],0.0);
          const double bulk_density = (1.0 - porosity) * scratch.material_model_outputs.densities[q] + porosity * melt_outputs->fluid_densities[q];

          const double density_c_P =
            ((scratch.advection_field->is_temperature())
             ?
             bulk_density *
             scratch.material_model_outputs.specific_heat[q]
             :
             1.0);

          Assert (density_c_P >= 0,
                  ExcMessage ("The product of density and c_P needs to be a "
                              "non-negative quantity."));

          const double conductivity =
            (scratch.advection_field->is_temperature()
             ?
             scratch.material_model_outputs.thermal_conductivities[q]
             :
             0.0);
          const double latent_heat_LHS =
            ((scratch.advection_field->is_temperature())
             ?
             scratch.heating_model_outputs.lhs_latent_heat_terms[q]
             :
             0.0);
          Assert (density_c_P + latent_heat_LHS >= 0,
                  ExcMessage ("The sum of density times c_P and the latent heat contribution "
                              "to the left hand side needs to be a non-negative quantity."));

          const double gamma =
            ((scratch.advection_field->is_temperature())
             ?
             scratch.heating_model_outputs.heating_source_terms[q]
             :
             0.0);

          const double reaction_term =
            ((scratch.advection_field->is_temperature() || this->get_melt_handler().is_porosity(*scratch.advection_field))
             ?
             0.0
             :
             scratch.material_model_outputs.reaction_terms[q][scratch.advection_field->compositional_variable]);

          const double melt_transport_RHS = (this->get_melt_handler().is_porosity(*scratch.advection_field) ?
                                             compute_melting_RHS (this,
                                                                  scratch,
                                                                  q,
                                                                  divergence_u)
                                             :
                                             0.0);

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
          // Subtract off the mesh velocity for ALE corrections if necessary
          if (this->get_parameters().mesh_deformation_enabled)
            current_u -= scratch.mesh_velocity_values[q];

          const double melt_transport_LHS =
            (this->get_melt_handler().is_porosity(*scratch.advection_field)
             && this->get_melt_handler().is_melt_cell(scratch.material_model_inputs.current_cell)
             ?
             divergence_u
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

          // for advective transport of heat by melt, or advecting compositional fields
          // with the melt velocity, we have to consider melt velocities
          const Tensor<1,dim> current_u_f = fluid_velocity_values[q];
          double density_c_P_solid = density_c_P;
          double density_c_P_melt = 0.0;

          if (scratch.advection_field->is_temperature() && this->get_melt_handler().is_melt_cell(scratch.cell)
              && this->get_melt_handler().melt_parameters.heat_advection_by_melt)
            {
              density_c_P_solid = (1.0 - porosity) * scratch.material_model_outputs.densities[q] * scratch.material_model_outputs.specific_heat[q];
              density_c_P_melt = porosity * melt_outputs->fluid_densities[q] * scratch.material_model_outputs.specific_heat[q];
            }
          else if (!scratch.advection_field->is_temperature()
                   && !this->get_melt_handler().is_porosity(*scratch.advection_field)
                   && scratch.advection_field->advection_method(this->introspection()) == Parameters<dim>::AdvectionFieldMethod::fem_melt_field)
            {
              // if the field is advected with the fem_melt_field method,
              // we want to use the melt (fluid) velocity instead of the solid velocity for advecting it
              density_c_P_solid = 0.0;
              density_c_P_melt = 1.0;
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
                       (time_step * (conductivity + scratch.artificial_viscosity)
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
    MeltAdvectionSystem<dim>::
    compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch =
        dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>& > (scratch_base);

      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      std::vector<double> residuals(n_q_points);

      this->get_heating_model_manager().evaluate(scratch.material_model_inputs,
                                                 scratch.material_model_outputs,
                                                 scratch.heating_model_outputs);

      // average divergence u over the cell (needed for porosity advection)
      double divergence_u = 0.0;
      if (this->get_melt_handler().is_porosity(*scratch.advection_field))
        for (unsigned int q=0; q<n_q_points; ++q)
          divergence_u += scratch.current_velocity_divergences[q] * 1./n_q_points;

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

          if (scratch.advection_field->is_temperature())
            {
              const double density       = scratch.material_model_outputs.densities[q];
              const double conductivity  = scratch.material_model_outputs.thermal_conductivities[q];
              const double c_P           = scratch.material_model_outputs.specific_heat[q];
              const double k_Delta_field = conductivity * (scratch.old_field_laplacians[q] +
                                                           scratch.old_old_field_laplacians[q]) / 2;

              const double gamma           = scratch.heating_model_outputs.heating_source_terms[q];
              const double latent_heat_LHS = scratch.heating_model_outputs.lhs_latent_heat_terms[q];

              residuals[q]
                = std::abs((density * c_P + latent_heat_LHS) * (dField_dt + u_grad_field) - k_Delta_field - gamma);
            }
          else
            {
              const double field = ((scratch.old_field_values)[q] + (scratch.old_old_field_values)[q]) / 2;

              const double dreaction_term_dt =
                (this->get_old_timestep() == 0 || this->get_melt_handler().is_porosity(*scratch.advection_field))
                ?
                0.0
                :
                (scratch.material_model_outputs.reaction_terms[q][scratch.advection_field->compositional_variable]
                 / this->get_old_timestep());


              const double melt_transport_RHS = (this->get_melt_handler().is_porosity(*scratch.advection_field) ?
                                                 compute_melting_RHS (this,
                                                                      scratch,
                                                                      q,
                                                                      divergence_u)
                                                 :
                                                 0.0);

              const double melt_transport_LHS =
                (this->get_melt_handler().is_porosity(*scratch.advection_field)
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
                = std::abs(dField_dt + u_grad_field + melt_transport_LHS * field
                           - melt_transport_RHS - dreaction_term_dt);
            }
        }
      return residuals;
    }



    template <int dim>
    void
    MeltPressureRHSCompatibilityModification<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();

      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int p_f_component_index = introspection.variable("fluid pressure").first_component_index;
      const unsigned int p_c_component_index = introspection.variable("compaction pressure").first_component_index;

      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
          {
            const unsigned int component_index_i = fe.system_to_component_index(i).first;

            if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
              {
                scratch.phi_p[i_stokes] = scratch.finite_element_values[ex_p_f].value (i, q);
                data.local_pressure_shape_function_integrals(i_stokes) += scratch.phi_p[i_stokes] * scratch.finite_element_values.JxW(q);
                ++i_stokes;
              }
            ++i;
          }
    }

    template <int dim>
    void
    MeltBoundaryTraction<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      // see if any of the faces are traction boundaries for which
      // we need to assemble force terms for the right hand side
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int p_f_component_index = introspection.variable("fluid pressure").first_component_index;
      const unsigned int p_c_component_index = introspection.variable("compaction pressure").first_component_index;

      const typename DoFHandler<dim>::active_cell_iterator cell (&this->get_triangulation(),
                                                                 scratch.face_finite_element_values.get_cell()->level(),
                                                                 scratch.face_finite_element_values.get_cell()->index(),
                                                                 &this->get_dof_handler());

      unsigned int face_no = numbers::invalid_unsigned_int;
      for (face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        if (scratch.face_finite_element_values.get_face_index() == cell->face_index(face_no))
          break;
      Assert(face_no != numbers::invalid_unsigned_int,ExcInternalError());

      const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);

      if (this->get_boundary_traction()
          .find (face->boundary_id())
          !=
          this->get_boundary_traction().end())
        {
          scratch.face_finite_element_values.reinit (cell, face_no);

          for (unsigned int q=0; q<scratch.face_finite_element_values.n_quadrature_points; ++q)
            {
              const Tensor<1,dim> traction
                = this->get_boundary_traction().find(
                    face->boundary_id()
                  )->second
                  ->boundary_traction (face->boundary_id(),
                                       scratch.face_finite_element_values.quadrature_point(q),
                                       scratch.face_finite_element_values.normal_vector(q));

              for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
                {
                  const unsigned int component_index_i = fe.system_to_component_index(i).first;

                  if (is_velocity_or_pressures(introspection,p_c_component_index,p_f_component_index,component_index_i))
                    {
                      data.local_rhs(i_stokes) += scratch.face_finite_element_values[introspection.extractors.velocities].value(i,q) *
                                                  traction *
                                                  scratch.face_finite_element_values.JxW(q);
                      ++i_stokes;
                    }
                  ++i;
                }
            }
        }
    }
  }


  template <int dim>
  void
  MeltHandler<dim>::
  compute_melt_variables(LinearAlgebra::BlockSparseMatrix &system_matrix,
                         LinearAlgebra::BlockVector &solution,
                         LinearAlgebra::BlockVector &system_rhs)
  {
    if (!this->include_melt_transport())
      return;

    Assert (this->include_melt_transport(), ExcMessage ("'Include melt transport' has to be on to "
                                                        "compute melt variables"));

    LinearAlgebra::BlockVector distributed_solution (this->introspection().index_sets.system_partitioning,
                                                     this->get_mpi_communicator());

    const unsigned int por_idx = this->introspection().compositional_index_for_name("porosity");

    {
      // compute fluid_velocity
      // u_f =  u_s - K_D (nabla p_f - rho_f g) / phi  or = 0
      // by solving a mass matrix problem

      const unsigned int block_idx = this->introspection().variable("fluid velocity").block_index;
      system_matrix.block(block_idx, block_idx) = 0;
      system_rhs.block(block_idx) = 0;

      const QGauss<dim> quadrature(this->get_parameters().stokes_velocity_degree+1);
      const FiniteElement<dim> &fe = this->get_fe();

      FEValues<dim> fe_values (this->get_mapping(),
                               fe,
                               quadrature,
                               update_quadrature_points | update_values | update_gradients | update_JxW_values);

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         n_q_points    = fe_values.n_quadrature_points;

      const FEVariable<dim> &u_f_variable = this->introspection().variable("fluid velocity");
      const unsigned int fluid_velocity_dofs_per_cell = fe.element_multiplicity(u_f_variable.base_index) *
                                                        fe.base_element(u_f_variable.base_index).dofs_per_cell;

      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
      std::vector<types::global_dof_index> cell_u_f_dof_indices (fluid_velocity_dofs_per_cell);

      Vector<double> cell_vector (fluid_velocity_dofs_per_cell);
      FullMatrix<double> cell_matrix (fluid_velocity_dofs_per_cell, fluid_velocity_dofs_per_cell);

      std::vector<double> porosity_values(quadrature.size());
      std::vector<Tensor<1,dim> > grad_p_f_values(quadrature.size());
      std::vector<Tensor<1,dim> > u_s_values(quadrature.size());

      MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());

      create_material_model_outputs(out);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            cell_vector = 0;
            cell_matrix = 0;
            fe_values.reinit (cell);

            cell->get_dof_indices (cell_dof_indices);

            // Extract the dofs of this cell that are actually fluid velocity dofs
            for (unsigned int i=0, i_u_f=0; i_u_f<fluid_velocity_dofs_per_cell; /*increment at end of loop*/)
              {
                if (u_f_variable.component_mask[fe.system_to_component_index(i).first])
                  {
                    cell_u_f_dof_indices[i_u_f] = cell_dof_indices[i];
                    ++i_u_f;
                  }
                ++i;
              }

            fe_values[this->introspection().extractors.compositional_fields[por_idx]].get_function_values (
              solution, porosity_values);
            fe_values[this->introspection().extractors.velocities].get_function_values (
              solution, u_s_values);

            fe_values[this->introspection().variable("fluid pressure").extractor_scalar()].get_function_gradients (
              solution, grad_p_f_values);


            in.reinit(fe_values,
                      cell,
                      this->introspection(),
                      this->get_solution());

            this->get_material_model().evaluate(in, out);

            const bool is_melt_cell = this->get_melt_handler().is_melt_cell(in.current_cell);

            MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
            Assert(melt_outputs != nullptr, ExcMessage("Need MeltOutputs from the material model for computing the melt variables."));
            const FEValuesExtractors::Vector fluid_velocity_extractor = this->introspection().variable("fluid velocity").extractor_vector();

            std::vector<Tensor<1,dim> > phi_u_f(fluid_velocity_dofs_per_cell);

            double K_D_over_phi = 1.0;

            if (melt_parameters.average_melt_velocity)
              {
                // average the K_D and the porosity cell-wise (geometric mean)
                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    const double K_D = melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q];
                    const double K_D_over_phi_q_point = (porosity_values[q] > 0 ? K_D / porosity_values[q] : 0.0);
                    K_D_over_phi *= std::pow(K_D_over_phi_q_point, 1.0/n_q_points);
                  }
              }

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                for (unsigned int i = 0, i_u_f = 0; i_u_f < fluid_velocity_dofs_per_cell; /*increment at end of loop*/)
                  {
                    if (u_f_variable.component_mask[fe.system_to_component_index(i).first])
                      {
                        phi_u_f[i_u_f] = fe_values[fluid_velocity_extractor].value(i,q);
                        ++i_u_f;
                      }
                    ++i;
                  }

                const double JxW = fe_values.JxW(q);

                for (unsigned int i=0; i<fluid_velocity_dofs_per_cell; ++i)
                  for (unsigned int j=0; j<fluid_velocity_dofs_per_cell; ++j)
                    cell_matrix(i,j) += phi_u_f[j] * phi_u_f[i] * JxW;

                if (!melt_parameters.average_melt_velocity)
                  {
                    // use K_D without cutoff to compute u_f
                    const double phi = std::max(0.0, porosity_values[q]);
                    const double K_D = melt_outputs->permeabilities[q] / melt_outputs->fluid_viscosities[q];

                    K_D_over_phi = (phi > 0 ? K_D / phi : 0.0);
                  }

                // u_f =  u_s - K_D (nabla p_f - rho_f g) / phi  or = 0
                if (is_melt_cell)
                  {
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(in.position[q]);


                    for (unsigned int i=0; i<fluid_velocity_dofs_per_cell; ++i)
                      cell_vector(i) += (u_s_values[q] - K_D_over_phi * (grad_p_f_values[q] - melt_outputs->fluid_densities[q]*gravity))
                                        * phi_u_f[i]
                                        * JxW;
                  }
                else
                  for (unsigned int i=0; i<fluid_velocity_dofs_per_cell; ++i)
                    cell_vector(i) += (u_s_values[q])
                                      * phi_u_f[i]
                                      * JxW;

              }

            this->get_current_constraints().distribute_local_to_global (cell_matrix, cell_vector,
                                                                        cell_u_f_dof_indices, system_matrix,
                                                                        system_rhs, false);
          }

      system_rhs.compress (VectorOperation::add);
      system_matrix.compress (VectorOperation::add);

      LinearAlgebra::PreconditionAMG preconditioner;
      LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
#ifdef ASPECT_USE_PETSC
      Amg_data.symmetric_operator = false;
#else
      // Amg_data.constant_modes = constant_modes;
      Amg_data.elliptic = true;
      Amg_data.higher_order_elements = false;
      Amg_data.smoother_sweeps = 2;
      Amg_data.aggregation_threshold = 0.02;
#endif
      preconditioner.initialize(system_matrix.block(block_idx, block_idx));

      SolverControl solver_control(system_rhs.block(block_idx).size(),
                                   1e-8*system_rhs.block(block_idx).l2_norm());
      SolverCG<LinearAlgebra::Vector> cg(solver_control);

      this->get_pcout() << "   Solving fluid velocity system... " << std::flush;

      cg.solve (system_matrix.block(block_idx, block_idx),
                distributed_solution.block(block_idx),
                system_rhs.block(block_idx),
                preconditioner);

      this->get_pcout() << solver_control.last_step() <<" iterations."<< std::endl;

      this->get_current_constraints().distribute (distributed_solution);
      solution.block(block_idx) = distributed_solution.block(block_idx);
    }

    // compute solid pressure as
    // p_s = (p_c+ (1-phi)p_f) / (1-phi)
    // p_s = p_c / (1-phi) + p_f
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
                               update_quadrature_points | update_gradients | update_values);

      MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
      create_material_model_outputs(out);


      std::vector<types::global_dof_index> local_dof_indices (this->get_fe().dofs_per_cell);
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
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
            in.reinit(fe_values,
                      cell,
                      this->introspection(),
                      solution,
                      false);

            this->get_material_model().evaluate(in, out);

            const double p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(
                                       this->get_material_model()).p_c_scale(in,
                                                                             out,
                                                                             this->get_melt_handler(),
                                                                             true);

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
                if (p_c_scale > 0 && (1.0-phi) > std::numeric_limits<double>::min())
                  p = (p_c_scale*p_c_values[j] - (phi-1.0) * p_f_values[j]) / (1.0-phi);

                distributed_solution(local_dof_indices[pressure_idx]) = p;
              }
          }

      distributed_solution.block(block_p).compress(VectorOperation::insert);
      solution.block(block_p) = distributed_solution.block(block_p);
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


  namespace internal
  {
    // This is a scratch object for the setting the compaction pressure constraints
    // in cells without melt (where we do not solve the melt transport equations,
    // so we set the compaction pressure to zero).
    template <int dim>
    struct PcConstraintsAssembleData
    {
      // Standard constructor
      PcConstraintsAssembleData (const FiniteElement<dim> &finite_element,
                                 const Quadrature<dim>    &quadrature,
                                 const Mapping<dim>       &mapping,
                                 const UpdateFlags         update_flags,
                                 const unsigned int        n_compositional_fields,
                                 const unsigned int        stokes_dofs_per_cell)
        :
        finite_element_values (mapping, finite_element, quadrature,
                               update_flags),
        local_dof_indices (finite_element.dofs_per_cell),
        dof_component_indices(stokes_dofs_per_cell),
        material_model_inputs(quadrature.size(), n_compositional_fields),
        material_model_outputs(quadrature.size(), n_compositional_fields)
      {}


      // Copy constructor
      PcConstraintsAssembleData (const PcConstraintsAssembleData &scratch)
        :
        finite_element_values (scratch.finite_element_values.get_mapping(),
                               scratch.finite_element_values.get_fe(),
                               scratch.finite_element_values.get_quadrature(),
                               scratch.finite_element_values.get_update_flags()),
        local_dof_indices (scratch.local_dof_indices),
        dof_component_indices( scratch.dof_component_indices),
        material_model_inputs(scratch.material_model_inputs),
        material_model_outputs(scratch.material_model_outputs)
      {}

      virtual ~PcConstraintsAssembleData ()
      {}

      FEValues<dim>                            finite_element_values;

      std::vector<types::global_dof_index>     local_dof_indices;
      std::vector<unsigned int>                dof_component_indices;

      MaterialModel::MaterialModelInputs<dim>  material_model_inputs;
      MaterialModel::MaterialModelOutputs<dim> material_model_outputs;

    };

    template <int dim>
    struct PcConstraintsCopyData
    {
      explicit PcConstraintsCopyData (const unsigned int stokes_dofs_per_cell)
      {
        nonzero_dof_indices.reserve(stokes_dofs_per_cell);
      }

      std::vector<types::global_dof_index>   nonzero_dof_indices;
    };


    // Assembler for setting the compaction pressure constraints in cells without melt.
    // Note that the assembler does not actually assemble anything, it just uses the
    // structure of an assembler to loop over all compaction pressure dofs and to insert
    // the ones that are in melt cells into an IndexSet (nonzero_dof_indices) that is
    // later needed for setting the constraints.
    template <int dim>
    class PcNonZeroDofsAssembler : public SimulatorAccess<dim>
    {
      public:
        PcNonZeroDofsAssembler(const std::vector<bool> &is_melt_cell,
                               IndexSet &nonzero_entries)
          : is_melt_cell (is_melt_cell),
            nonzero_entries (nonzero_entries)
        {}


        void
        local_save_nonzero_pc_dofs (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    PcConstraintsAssembleData<dim> &scratch,
                                    PcConstraintsCopyData<dim> &data)
        {
          data.nonzero_dof_indices.clear();
          cell->get_dof_indices (scratch.local_dof_indices);

          const unsigned int pc_component_index = this->introspection().variable("compaction pressure").first_component_index;
          const unsigned int pc_base_index = this->introspection().variable("compaction pressure").base_index;

          Assert(cell.state() == IteratorState::valid,
                 ExcMessage("You are trying to access the property of a cell (if it is "
                            "a melt cell or not), but the cell iterator provided does not "
                            "point to a cell."));
          const bool is_melt_cell = this->get_melt_handler().is_melt_cell(cell);

          if (is_melt_cell)
            {
              data.nonzero_dof_indices.resize(this->get_fe().base_element(pc_base_index).dofs_per_cell);
              for (unsigned int j=0; j<this->get_fe().base_element(pc_base_index).dofs_per_cell; ++j)
                {
                  const unsigned int pressure_idx
                    = this->get_fe().component_to_system_index(pc_component_index,
                                                               /*dof index within component=*/ j);

                  data.nonzero_dof_indices[j] = scratch.local_dof_indices[pressure_idx];
                }
            }
        }

        void
        copy_local_to_global (const PcConstraintsCopyData<dim> &data)
        {
          nonzero_entries.add_indices(data.nonzero_dof_indices.begin(),
                                      data.nonzero_dof_indices.end());
        }

      private:
        const std::vector<bool> &is_melt_cell;
        IndexSet &nonzero_entries;

    };

  }

  template <int dim>
  void
  MeltHandler<dim>::
  add_current_constraints(AffineConstraints<double> &constraints)
  {
    IndexSet nonzero_pc_dofs(this->introspection().index_sets.system_relevant_set.size());

#if DEAL_II_VERSION_GTE(9,3,0)
    const QTrapezoid<dim> quadrature_formula;
#else
    const QTrapez<dim> quadrature_formula;
#endif
    const UpdateFlags cell_update_flags = update_quadrature_points | update_values | update_gradients;
    const FiniteElement<dim> &fe = this->get_fe();

    is_melt_cell_vector.resize(this->get_dof_handler().get_triangulation().n_active_cells());

    {
      // find the "melt cells" by looking at p_c_scale

      const unsigned int n_compositional_fields = this->introspection().n_compositional_fields;
      FEValues<dim> finite_element_values(this->get_mapping(), fe, quadrature_formula, cell_update_flags);

      MaterialModel::MaterialModelInputs<dim> material_model_inputs(quadrature_formula.size(), n_compositional_fields);
      MaterialModel::MaterialModelOutputs<dim> material_model_outputs(quadrature_formula.size(), n_compositional_fields);

      MeltHandler<dim>::create_material_model_outputs(material_model_outputs);
      Assert(Plugins::plugin_type_matches<const MaterialModel::MeltInterface<dim>>(this->get_material_model()),
             ExcMessage("Your material model does not derive from MaterialModel::MeltInterface, which is required."));

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            finite_element_values.reinit (cell);

            material_model_inputs.reinit(finite_element_values,
                                         cell,
                                         this->introspection(),
                                         this->get_current_linearization_point(),
                                         false);

            this->get_material_model().evaluate(material_model_inputs,
                                                material_model_outputs);

            const double p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(
                                       this->get_material_model()).p_c_scale(material_model_inputs,
                                                                             material_model_outputs,
                                                                             this->get_melt_handler(),
                                                                             false /*=consider_is_melt_cell*/);
            const bool is_melt_cell = (p_c_scale > 0.0);
            is_melt_cell_vector[cell->active_cell_index()] = is_melt_cell;
          }
    }


    internal::PcNonZeroDofsAssembler<dim> assembler(is_melt_cell_vector, nonzero_pc_dofs);
    assembler.initialize_simulator(this->get_simulator());

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    const unsigned int stokes_dofs_per_cell = dim * fe.base_element(this->introspection().base_elements.velocities).dofs_per_cell
                                              + fe.base_element(this->introspection().base_elements.pressure).dofs_per_cell
                                              + fe.base_element(this->introspection().base_elements.pressure).dofs_per_cell;

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::PcConstraintsAssembleData<dim> &scratch,
                      internal::PcConstraintsCopyData<dim> &data)
    {
      assembler.local_save_nonzero_pc_dofs(cell, scratch, data);
    };

    auto copier = [&](const internal::PcConstraintsCopyData<dim> &data)
    {
      assembler.copy_local_to_global(data);
    };

    // Here we call the assembler in order to save all compaction pressure dofs
    // that are in melt cells (and are nonzero) in the nonzero_pc_dofs index set.
    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     this->get_dof_handler().begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     this->get_dof_handler().end()),
         worker,
         copier,
         internal::PcConstraintsAssembleData<dim> (
           fe,
           quadrature_formula,
           this->get_mapping(),
           cell_update_flags,
           this->introspection().n_compositional_fields,
           stokes_dofs_per_cell),
         internal::PcConstraintsCopyData<dim> (stokes_dofs_per_cell));

    // For the constraints, first pick all relevant p_c dofs:
    IndexSet for_constraints = this->introspection().index_sets.system_relevant_set
                               & Utilities::extract_locally_active_dofs_with_component<dim>(this->get_dof_handler(),
                                   this->introspection().variable("compaction pressure").component_mask);

    // now subtract the ones that are nonzero as computed above:
    for_constraints.subtract_set(nonzero_pc_dofs);

    // and constrain the remaining dofs (that are not in melt cells).
    constraints.add_lines(for_constraints);
  }


  template <int dim>
  bool
  MeltHandler<dim>::
  is_melt_cell(const typename DoFHandler<dim>::active_cell_iterator &cell) const
  {
    return is_melt_cell_vector[cell->active_cell_index()];
  }


  template <int dim>
  double
  MeltHandler<dim>::
  limited_darcy_coefficient(const double K_D,
                            const bool is_melt_cell) const
  {
    const double ref_K_D = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(this->get_material_model()).reference_darcy_coefficient();
    return is_melt_cell ? std::max(K_D, melt_parameters.melt_scaling_factor_threshold*ref_K_D) : 0;
  }


  template <int dim>
  const BoundaryFluidPressure::Interface<dim> &
  MeltHandler<dim>::
  get_boundary_fluid_pressure () const
  {
    return *boundary_fluid_pressure.get();
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
                       std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree-1),
                       1,
                       0)); // same block as p_c even without a direct solver!

    variables.insert(variables.begin()+2,
                     VariableDeclaration<dim>(
                       "compaction pressure",
                       melt_parameters.use_discontinuous_p_c
                       ?
                       std::shared_ptr<FiniteElement<dim>>(
                         std::make_shared<FE_DGP<dim>>(parameters.stokes_velocity_degree-1))
                       :
                       std::shared_ptr<FiniteElement<dim>>(
                         std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree-1)),
                       1,
                       1));

    variables.insert(variables.begin()+3,
                     VariableDeclaration<dim>("fluid velocity",
                                              std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree),
                                              dim,
                                              1));

  }



  template <int dim>
  void
  MeltHandler<dim>::
  set_assemblers (Assemblers::Manager<dim> &assemblers) const
  {
    assemblers.stokes_preconditioner.push_back(std_cxx14::make_unique<Assemblers::MeltStokesPreconditioner<dim> > ());
    assemblers.stokes_system.push_back(std_cxx14::make_unique<Assemblers::MeltStokesSystem<dim> > ());

    AssertThrow((this->get_parameters().formulation_mass_conservation ==
                 Parameters<dim>::Formulation::MassConservation::isentropic_compression) ||
                (this->get_parameters().formulation_mass_conservation ==
                 Parameters<dim>::Formulation::MassConservation::incompressible),
                ExcMessage("The melt implementation currently only supports the isentropic compression "
                           "approximation or the incompressible formulation of the mass conservation equation."));

    // add the boundary integral for melt migration
    assemblers.stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data = true;
    assemblers.stokes_system_assembler_on_boundary_face_properties.needed_update_flags = (update_values  | update_quadrature_points |
        update_normal_vectors | update_gradients |
        update_JxW_values);


    assemblers.stokes_system_on_boundary_face.push_back(
      std_cxx14::make_unique<aspect::Assemblers::MeltStokesSystemBoundary<dim> >());

    // add the terms for traction boundary conditions
    if (!this->get_boundary_traction().empty())
      {
        assemblers.stokes_system_on_boundary_face.push_back(
          std_cxx14::make_unique<Assemblers::MeltBoundaryTraction<dim> > ());
      }

    // add the terms necessary to normalize the pressure
    if (this->pressure_rhs_needs_compatibility_modification())
      {
        assemblers.stokes_system.push_back(
          std_cxx14::make_unique<Assemblers::MeltPressureRHSCompatibilityModification<dim> > ());
      }

    assemblers.advection_system.push_back(
      std_cxx14::make_unique<Assemblers::MeltAdvectionSystem<dim> > ());

    if (this->get_parameters().fixed_heat_flux_boundary_indicators.size() != 0)
      {
        assemblers.advection_system_on_boundary_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::AdvectionSystemBoundaryHeatFlux<dim> >());
        assemblers.advection_system_assembler_on_face_properties[0].need_face_material_model_data = true;
        assemblers.advection_system_assembler_on_face_properties[0].need_face_finite_element_evaluation = true;
      }
  }


  namespace Melt
  {
    template <int dim>
    void
    Parameters<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Melt settings");
      {
        prm.declare_entry ("Melt scaling factor threshold", "1e-7",
                           Patterns::Double (),
                           "The factor by how much the Darcy coefficient K\\_D in a cell can be smaller than "
                           "the reference Darcy coefficient for this cell still to be considered a melt cell "
                           "(for which the melt transport equations are solved). For smaller Darcy coefficients, "
                           "the Stokes equations (without melt) are solved instead. Only used if "
                           "``Include melt transport'' is true. ");
        prm.declare_entry ("Heat advection by melt", "false",
                           Patterns::Bool (),
                           "Whether to use a porosity weighted average of the melt and solid velocity "
                           "to advect heat in the temperature equation or not. If this is set to true, "
                           "additional terms are assembled on the left-hand side of the temperature "
                           "advection equation. Only used if Include melt transport is true. "
                           "If this is set to false, only the solid velocity is used (as in models "
                           "without melt migration).");
        prm.declare_entry ("Use discontinuous compaction pressure", "true",
                           Patterns::Bool (),
                           "Whether to use a discontinuous element for the compaction pressure or not. "
                           "From our preliminary tests, continuous elements seem to work better in models "
                           "where the porosity is > 0 everywhere in the domain, and discontinuous elements "
                           "work better in models where in parts of the domain the porosity = 0.");
        prm.declare_entry ("Average melt velocity", "true",
                           Patterns::Bool (),
                           "Whether to cell-wise average the material properties that are used to "
                           "compute the melt velocity or not. The melt velocity is computed as the "
                           "sum of the solid velocity and the phase separation flux "
                           "$ - K_D / \\phi (\\nabla p_f - \\rho_f \\mathbf g)$. "
                           "If this parameter is set to true, $K_D$ and $\\phi$ will be averaged "
                           "cell-wise in the computation of the phase separation flux. "
                           "This is useful because in some models the melt velocity can have spikes "
                           "close to the interface between regions of melt and no melt, as both $K_D$ "
                           "and $\\phi$ go to zero for vanishing melt fraction. As the melt velocity is "
                           "used for computing the time step size, and in models that use heat "
                           "transport by melt or shear heating of melt, setting this parameter to true "
                           "can speed up the model and make it mode stable. In computations where "
                           "accuracy and convergence behavior of the melt velocity is important "
                           "(like in benchmark cases with an analytical solution), this parameter "
                           "should probably be set to 'false'.");
      }
      prm.leave_subsection();

      BoundaryFluidPressure::declare_parameters<dim> (prm);
    }

    template <int dim>
    void
    Parameters<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Melt settings");
      {
        melt_scaling_factor_threshold = prm.get_double("Melt scaling factor threshold");
        heat_advection_by_melt = prm.get_bool("Heat advection by melt");
        use_discontinuous_p_c = prm.get_bool("Use discontinuous compaction pressure");
        average_melt_velocity = prm.get_bool("Average melt velocity");
      }
      prm.leave_subsection();
    }
  }

  template <int dim>
  MeltHandler<dim>::MeltHandler (ParameterHandler &prm)
    :
    boundary_fluid_pressure(BoundaryFluidPressure::create_boundary_fluid_pressure<dim>(prm))
  {
    CitationInfo::add("melt");
    melt_parameters.parse_parameters(prm);
    boundary_fluid_pressure->parse_parameters(prm);
  }


  template <int dim>
  void
  MeltHandler<dim>::initialize () const
  {
    // The additional terms in the temperature systems have not been ported
    // to the DG formulation:
    AssertThrow(!this->get_parameters().use_discontinuous_temperature_discretization &&
                !this->get_parameters().use_discontinuous_composition_discretization,
                ExcMessage("Using discontinuous elements for temperature "
                           "or composition in models with melt transport is currently not implemented.") );
    if (melt_parameters.use_discontinuous_p_c)
      AssertThrow(!this->model_has_prescribed_stokes_solution(),
                  ExcMessage("You can not use a discontinuous p_c in a model "
                             "with a prescribed Stokes solution."));
    // We can not have a DG p_f.
    AssertThrow(!this->get_parameters().use_locally_conservative_discretization,
                ExcMessage ("Discontinuous elements for the fluid pressure "
                            "are not supported in models with melt transport."));
  }


  template <int dim>
  void
  MeltHandler<dim>::initialize_simulator (const Simulator<dim> &simulator_object)
  {
    this->SimulatorAccess<dim>::initialize_simulator(simulator_object);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_fluid_pressure.get()))
      sim->initialize_simulator (simulator_object);
    boundary_fluid_pressure->initialize ();
  }


  template <int dim>
  void
  MeltHandler<dim>::
  create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output)
  {
    if (output.template get_additional_output<MaterialModel::MeltOutputs<dim> >() != nullptr)
      return;

    const unsigned int n_points = output.viscosities.size();
    const unsigned int n_comp = output.reaction_terms[0].size();
    output.additional_outputs.push_back(
      std_cxx14::make_unique<MaterialModel::MeltOutputs<dim>> (n_points, n_comp));
  }


  template <int dim>
  void
  MeltHandler<dim>::
  apply_free_surface_stabilization_with_melt (const double free_surface_theta,
                                              const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                              internal::Assembly::CopyData::StokesSystem<dim>      &data) const
  {
    const std::set<types::boundary_id> free_surface_boundary_indicators = this->get_mesh_deformation_handler().get_free_surface_boundary_indicators();
    if (free_surface_boundary_indicators.empty())
      return;

    const unsigned int n_face_q_points = scratch.face_finite_element_values.n_quadrature_points;
    const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();

    // only apply on free surface faces
    if (cell->at_boundary() && cell->is_locally_owned())
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        if (cell->face(face_no)->at_boundary())
          {
            const types::boundary_id boundary_indicator
              = cell->face(face_no)->boundary_id();

            if (free_surface_boundary_indicators.find(boundary_indicator) == free_surface_boundary_indicators.end())
              continue;

            scratch.face_finite_element_values.reinit(cell, face_no);

            this->compute_material_model_input_values (this->get_solution(),
                                                       scratch.face_finite_element_values,
                                                       cell,
                                                       true,
                                                       scratch.face_material_model_inputs);

            this->get_material_model().evaluate(scratch.face_material_model_inputs, scratch.face_material_model_outputs);

            const unsigned int p_f_component_index = this->introspection().variable("fluid pressure").first_component_index;
            const unsigned int p_c_component_index = this->introspection().variable("compaction pressure").first_component_index;

            for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
              {
                for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
                  {
                    const unsigned int component_index_i = this->get_fe().system_to_component_index(i).first;

                    if (Assemblers::is_velocity_or_pressures(this->introspection(),p_c_component_index,p_f_component_index,component_index_i))
                      {
                        scratch.phi_u[i_stokes] = scratch.face_finite_element_values[this->introspection().extractors.velocities].value(i, q_point);
                        ++i_stokes;
                      }
                    ++i;
                  }

                const Tensor<1,dim>
                gravity = this->get_gravity_model().gravity_vector(scratch.face_finite_element_values.quadrature_point(q_point));
                const double g_norm = gravity.norm();

                // construct the relevant vectors
                const Tensor<1,dim> n_hat = scratch.face_finite_element_values.normal_vector(q_point);
                const Tensor<1,dim> g_hat = (g_norm == 0.0 ? Tensor<1,dim>() : gravity/g_norm);

                const double pressure_perturbation = scratch.face_material_model_outputs.densities[q_point] *
                                                     this->get_timestep() * free_surface_theta * g_norm;

                // see Kaus et al 2010 for details of the stabilization term
                for (unsigned int i=0; i< stokes_dofs_per_cell; ++i)
                  for (unsigned int j=0; j< stokes_dofs_per_cell; ++j)
                    {
                      // The fictive stabilization stress is (phi_u[i].g)*(phi_u[j].n)
                      const double stress_value = - pressure_perturbation
                                                  * (scratch.phi_u[i] * g_hat) * (scratch.phi_u[j] * n_hat)
                                                  * scratch.face_finite_element_values.JxW(q_point);

                      data.local_matrix(i,j) += stress_value;
                    }
              }
          }


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
  namespace MaterialModel \
  { \
    template \
    class \
    MeltInputs<dim>; \
  } \
  \
  namespace Melt \
  { \
    template struct Parameters<dim>; \
  } \
  namespace Assemblers \
  { \
    template class MeltInterface<dim>; \
    template class MeltStokesPreconditioner<dim>; \
    template class MeltStokesSystem<dim>; \
    template class MeltStokesSystemBoundary<dim>; \
    template class MeltAdvectionSystem<dim>; \
    template class MeltPressureRHSCompatibilityModification<dim>; \
    template class MeltBoundaryTraction<dim>; \
  } \
   
  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE

}

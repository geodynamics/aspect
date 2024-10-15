/*
  Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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

#include <aspect/simulator/assemblers/interface.h>

#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const FiniteElement<dim> &finite_element,
                              const Quadrature<dim>    &quadrature,
                              const Mapping<dim>       &mapping,
                              const UpdateFlags         update_flags,
                              const unsigned int        n_compositional_fields,
                              const unsigned int        stokes_dofs_per_cell,
                              const bool                add_compaction_pressure,
                              const bool                rebuild_matrix,
                              const bool                use_bfbt)
          :
          ScratchBase<dim>(),

          finite_element_values (mapping, finite_element, quadrature,
                                 update_flags),
          local_dof_indices (finite_element.dofs_per_cell),
          dof_component_indices(stokes_dofs_per_cell),
          grads_phi_u (stokes_dofs_per_cell, numbers::signaling_nan<SymmetricTensor<2,dim>>()),
          div_phi_u (stokes_dofs_per_cell, numbers::signaling_nan<double>()),
          phi_p (stokes_dofs_per_cell, numbers::signaling_nan<double>()),
          phi_u (stokes_dofs_per_cell,numbers::signaling_nan<Tensor<1,dim>>()),
          phi_p_c (add_compaction_pressure ? stokes_dofs_per_cell : 0, numbers::signaling_nan<double>()),
          grad_phi_p ((add_compaction_pressure || use_bfbt) ? stokes_dofs_per_cell : 0, numbers::signaling_nan<Tensor<1,dim>>()),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size(), n_compositional_fields),
          rebuild_stokes_matrix(rebuild_matrix)
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &scratch)
          :
          ScratchBase<dim>(scratch),

          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),

          local_dof_indices (scratch.local_dof_indices),
          dof_component_indices( scratch.dof_component_indices),
          grads_phi_u (scratch.grads_phi_u),
          div_phi_u (scratch.div_phi_u),
          phi_p (scratch.phi_p),
          phi_u (scratch.phi_u),
          phi_p_c (scratch.phi_p_c),
          grad_phi_p(scratch.grad_phi_p),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs),
          rebuild_stokes_matrix(scratch.rebuild_stokes_matrix)
        {}


        template <int dim>
        StokesPreconditioner<dim>::
        ~StokesPreconditioner ()
          = default;


        template <int dim>
        void
        StokesPreconditioner<dim>::
        reinit (const typename DoFHandler<dim>::active_cell_iterator &cell_ref)
        {
          this->cell = cell_ref;
          this->face_number = numbers::invalid_unsigned_int;
          finite_element_values.reinit (cell_ref);
        }


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &finite_element,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &quadrature,
                      const Quadrature<dim-1>  &face_quadrature,
                      const UpdateFlags         update_flags,
                      const UpdateFlags         face_update_flags,
                      const unsigned int        n_compositional_fields,
                      const unsigned int        stokes_dofs_per_cell,
                      const bool                add_compaction_pressure,
                      const bool                use_reference_density_profile,
                      const bool                rebuild_stokes_matrix,
                      const bool                rebuild_newton_stokes_matrix,
                      const bool                use_bfbt)
          :
          StokesPreconditioner<dim> (finite_element, quadrature,
                                     mapping,
                                     update_flags,
                                     n_compositional_fields,
                                     stokes_dofs_per_cell,
                                     add_compaction_pressure,
                                     rebuild_stokes_matrix,
                                     use_bfbt),

          face_finite_element_values (mapping,
                                      finite_element,
                                      face_quadrature,
                                      face_update_flags),

          phi_u (stokes_dofs_per_cell, numbers::signaling_nan<Tensor<1,dim>>()),
          velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          velocity_divergence(quadrature.size(), numbers::signaling_nan<double>()),
          temperature_gradients (quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_material_model_inputs(face_quadrature.size(), n_compositional_fields),
          face_material_model_outputs(face_quadrature.size(), n_compositional_fields),
          reference_densities(use_reference_density_profile ? quadrature.size() : 0, numbers::signaling_nan<double>()),
          reference_densities_depth_derivative(use_reference_density_profile ? quadrature.size() : 0, numbers::signaling_nan<double>()),
          rebuild_newton_stokes_matrix(rebuild_newton_stokes_matrix)
        {}



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &scratch)
          :
          StokesPreconditioner<dim> (scratch),

          face_finite_element_values (scratch.face_finite_element_values.get_mapping(),
                                      scratch.face_finite_element_values.get_fe(),
                                      scratch.face_finite_element_values.get_quadrature(),
                                      scratch.face_finite_element_values.get_update_flags()),

          phi_u (scratch.phi_u),
          velocity_values (scratch.velocity_values),
          velocity_divergence (scratch.velocity_divergence),
          temperature_gradients (scratch.temperature_gradients),
          face_material_model_inputs(scratch.face_material_model_inputs),
          face_material_model_outputs(scratch.face_material_model_outputs),
          reference_densities(scratch.reference_densities),
          reference_densities_depth_derivative(scratch.reference_densities_depth_derivative),
          rebuild_newton_stokes_matrix(scratch.rebuild_newton_stokes_matrix)
        {}

        template <int dim>
        void
        StokesSystem<dim>::
        reinit (const typename DoFHandler<dim>::active_cell_iterator &cell_ref,
                const unsigned face_number_ref)
        {
          this->cell = cell_ref;
          this->face_number = face_number_ref;
          face_finite_element_values.reinit(cell_ref, face_number_ref);
        }



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element,
                         const FiniteElement<dim> &advection_element,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &quadrature,
                         const Quadrature<dim-1>  &face_quadrature,
                         const UpdateFlags         update_flags,
                         const UpdateFlags         face_update_flags,
                         const unsigned int        n_compositional_fields,
                         const typename Simulator<dim>::AdvectionField &field)
          :
          ScratchBase<dim>(),

          finite_element_values (mapping,
                                 finite_element, quadrature,
                                 update_flags),
          face_finite_element_values (face_quadrature.size() > 0
                                      ?
                                      std::make_unique<FEFaceValues<dim>> (mapping,
                                                                            finite_element, face_quadrature,
                                                                            face_update_flags)
                                      :
                                      nullptr),
          neighbor_face_finite_element_values (face_quadrature.size() > 0
                                               ?
                                               std::make_unique<FEFaceValues<dim>> (mapping,
                                                                                     finite_element, face_quadrature,
                                                                                     face_update_flags)
                                               :
                                               nullptr),
          subface_finite_element_values (face_quadrature.size() > 0
                                         ?
                                         std::make_unique<FESubfaceValues<dim>> (mapping,
                                                                                  finite_element, face_quadrature,
                                                                                  face_update_flags)
                                         :
                                         nullptr),
          local_dof_indices (finite_element.dofs_per_cell),

          phi_field (advection_element.dofs_per_cell, numbers::signaling_nan<double>()),
          grad_phi_field (advection_element.dofs_per_cell, numbers::signaling_nan<Tensor<1,dim>>()),
          laplacian_phi_field (advection_element.dofs_per_cell, numbers::signaling_nan<double>()),
          face_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                          numbers::signaling_nan<double>()),
          face_grad_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                               numbers::signaling_nan<Tensor<1,dim>>()),
          neighbor_face_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                                   numbers::signaling_nan<double>()),
          neighbor_face_grad_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                                        numbers::signaling_nan<Tensor<1,dim>>()),
          old_velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_old_velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_pressure (quadrature.size(), numbers::signaling_nan<double>()),
          old_old_pressure (quadrature.size(), numbers::signaling_nan<double>()),
          old_pressure_gradients (quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_old_pressure_gradients (quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_strain_rates (quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim>>()),
          old_old_strain_rates (quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim>>()),
          old_temperature_values (quadrature.size(), numbers::signaling_nan<double>()),
          old_old_temperature_values(quadrature.size(), numbers::signaling_nan<double>()),
          old_field_values (quadrature.size(), numbers::signaling_nan<double>()),
          old_old_field_values(quadrature.size(), numbers::signaling_nan<double>()),
          old_field_grads(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_old_field_grads(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_field_laplacians(quadrature.size(), numbers::signaling_nan<double>()),
          old_old_field_laplacians(quadrature.size(), numbers::signaling_nan<double>()),
          old_composition_values(n_compositional_fields,
                                 std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
          old_old_composition_values(n_compositional_fields,
                                     std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
          current_temperature_values(quadrature.size(), numbers::signaling_nan<double>()),
          current_velocity_values(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_current_velocity_values(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          mesh_velocity_values(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_mesh_velocity_values(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          current_strain_rates(quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim>>()),
          current_composition_values(n_compositional_fields,
                                     std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
          current_velocity_divergences(quadrature.size(), numbers::signaling_nan<double>()),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size(), n_compositional_fields),
          face_material_model_inputs(face_quadrature.size(), n_compositional_fields),
          face_material_model_outputs(face_quadrature.size(), n_compositional_fields),
          neighbor_face_material_model_inputs(face_quadrature.size(), n_compositional_fields),
          neighbor_face_material_model_outputs(face_quadrature.size(), n_compositional_fields),
          heating_model_outputs(quadrature.size(), n_compositional_fields),
          face_heating_model_outputs(face_quadrature.size(), n_compositional_fields),
          neighbor_face_heating_model_outputs(face_quadrature.size(), n_compositional_fields),
          advection_field(&field),
          artificial_viscosity(numbers::signaling_nan<double>())
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const AdvectionSystem &scratch)
          :
          ScratchBase<dim>(scratch),

          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),
          face_finite_element_values (scratch.face_finite_element_values.get()
                                      ?
                                      std::make_unique<FEFaceValues<dim>> (scratch.face_finite_element_values->get_mapping(),
                                                                            scratch.face_finite_element_values->get_fe(),
                                                                            scratch.face_finite_element_values->get_quadrature(),
                                                                            scratch.face_finite_element_values->get_update_flags())
                                      :
                                      nullptr),
          neighbor_face_finite_element_values (scratch.neighbor_face_finite_element_values.get()
                                               ?
                                               std::make_unique<FEFaceValues<dim>> (scratch.neighbor_face_finite_element_values->get_mapping(),
                                                                                     scratch.neighbor_face_finite_element_values->get_fe(),
                                                                                     scratch.neighbor_face_finite_element_values->get_quadrature(),
                                                                                     scratch.neighbor_face_finite_element_values->get_update_flags())
                                               :
                                               nullptr),
          subface_finite_element_values (scratch.subface_finite_element_values.get()
                                         ?
                                         std::make_unique<FESubfaceValues<dim>> (scratch.subface_finite_element_values->get_mapping(),
                                                                                  scratch.subface_finite_element_values->get_fe(),
                                                                                  scratch.subface_finite_element_values->get_quadrature(),
                                                                                  scratch.subface_finite_element_values->get_update_flags())
                                         :
                                         nullptr),
          local_dof_indices (scratch.finite_element_values.get_fe().dofs_per_cell),

          phi_field (scratch.phi_field),
          grad_phi_field (scratch.grad_phi_field),
          laplacian_phi_field (scratch.laplacian_phi_field),
          face_phi_field (scratch.face_phi_field),
          face_grad_phi_field (scratch.face_grad_phi_field),
          neighbor_face_phi_field (scratch.neighbor_face_phi_field),
          neighbor_face_grad_phi_field (scratch.neighbor_face_grad_phi_field),
          old_velocity_values (scratch.old_velocity_values),
          old_old_velocity_values (scratch.old_old_velocity_values),
          old_pressure (scratch.old_pressure),
          old_old_pressure (scratch.old_old_pressure),
          old_pressure_gradients (scratch.old_pressure_gradients),
          old_old_pressure_gradients (scratch.old_old_pressure_gradients),
          old_strain_rates (scratch.old_strain_rates),
          old_old_strain_rates (scratch.old_old_strain_rates),
          old_temperature_values (scratch.old_temperature_values),
          old_old_temperature_values (scratch.old_old_temperature_values),
          old_field_values (scratch.old_field_values),
          old_old_field_values(scratch.old_old_field_values),
          old_field_grads (scratch.old_field_grads),
          old_old_field_grads (scratch.old_old_field_grads),
          old_field_laplacians (scratch.old_field_laplacians),
          old_old_field_laplacians (scratch.old_old_field_laplacians),
          old_composition_values(scratch.old_composition_values),
          old_old_composition_values(scratch.old_old_composition_values),
          current_temperature_values(scratch.current_temperature_values),
          current_velocity_values(scratch.current_velocity_values),
          face_current_velocity_values(scratch.face_current_velocity_values),
          mesh_velocity_values(scratch.mesh_velocity_values),
          face_mesh_velocity_values(scratch.face_mesh_velocity_values),
          current_strain_rates(scratch.current_strain_rates),
          current_composition_values(scratch.current_composition_values),
          current_velocity_divergences(scratch.current_velocity_divergences),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs),
          face_material_model_inputs(scratch.face_material_model_inputs),
          face_material_model_outputs(scratch.face_material_model_outputs),
          neighbor_face_material_model_inputs(scratch.neighbor_face_material_model_inputs),
          neighbor_face_material_model_outputs(scratch.neighbor_face_material_model_outputs),
          heating_model_outputs(scratch.heating_model_outputs),
          face_heating_model_outputs(scratch.face_heating_model_outputs),
          neighbor_face_heating_model_outputs(scratch.neighbor_face_heating_model_outputs),
          advection_field(scratch.advection_field),
          artificial_viscosity(scratch.artificial_viscosity)
        {}


        template <int dim>
        void
        AdvectionSystem<dim>::
        reinit (const typename DoFHandler<dim>::active_cell_iterator &cell_ref)
        {
          this->cell = cell_ref;
          this->face_number = numbers::invalid_unsigned_int;
          finite_element_values.reinit (cell_ref);
        }

      }


      namespace CopyData
      {

        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const unsigned int stokes_dofs_per_cell)
          :
          local_matrix (stokes_dofs_per_cell,
                        stokes_dofs_per_cell),
          local_inverse_lumped_mass_matrix (stokes_dofs_per_cell),
          local_dof_indices (stokes_dofs_per_cell)
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &data)
          :
          local_matrix (data.local_matrix),
          local_inverse_lumped_mass_matrix (data.local_inverse_lumped_mass_matrix),
          local_dof_indices (data.local_dof_indices)
        {}


        template <int dim>
        void StokesPreconditioner<dim>::
        extract_stokes_dof_indices(const std::vector<types::global_dof_index> &all_dof_indices,
                                   const Introspection<dim>                   &introspection,
                                   const FiniteElement<dim>           &finite_element)
        {
          const unsigned int dofs_per_cell = finite_element.dofs_per_cell;

          for (unsigned int i=0, i_stokes=0; i<dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(finite_element.system_to_component_index(i).first))
                {
                  this->local_dof_indices[i_stokes] = all_dof_indices[i];
                  ++i_stokes;
                }
              ++i;
            }
        }



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const unsigned int        stokes_dofs_per_cell,
                      const bool                do_pressure_rhs_compatibility_modification)
          :
          StokesPreconditioner<dim> (stokes_dofs_per_cell),
          local_rhs (stokes_dofs_per_cell),
          local_pressure_shape_function_integrals (do_pressure_rhs_compatibility_modification ?
                                                   stokes_dofs_per_cell
                                                   :
                                                   0)
        {}



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &data)
          :
          StokesPreconditioner<dim> (data),
          local_rhs (data.local_rhs),
          local_pressure_shape_function_integrals (data.local_pressure_shape_function_integrals.size())
        {}




        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element,
                         const bool                field_is_discontinuous)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_matrices_int_ext ((field_is_discontinuous
                                   ?
                                   Assemblers::n_interface_matrices(finite_element.reference_cell())
                                   :
                                   0),
                                  FullMatrix<double>(finite_element.dofs_per_cell,
                                                     finite_element.dofs_per_cell)),
          local_matrices_ext_int ((field_is_discontinuous
                                   ?
                                   Assemblers::n_interface_matrices(finite_element.reference_cell())
                                   :
                                   0),
                                  FullMatrix<double>(finite_element.dofs_per_cell,
                                                     finite_element.dofs_per_cell)),
          local_matrices_ext_ext ((field_is_discontinuous
                                   ?
                                   Assemblers::n_interface_matrices(finite_element.reference_cell())
                                   :
                                   0),
                                  FullMatrix<double>(finite_element.dofs_per_cell,
                                                     finite_element.dofs_per_cell)),
          local_rhs (finite_element.dofs_per_cell),

          assembled_matrices ((field_is_discontinuous
                               ?
                               Assemblers::n_interface_matrices(finite_element.reference_cell())
                               :
                               0), false),

          local_dof_indices (finite_element.dofs_per_cell),
          neighbor_dof_indices ((field_is_discontinuous
                                 ?
                                 Assemblers::n_interface_matrices(finite_element.reference_cell())
                                 :
                                 0),
                                std::vector<types::global_dof_index>(finite_element.dofs_per_cell))
        {}

      }
    }
  }



  namespace Assemblers
  {
    unsigned int
    n_interface_matrices (const ReferenceCell &reference_cell)
    {
      // The current implementation assumes that all faces are
      // the same; so no wedges or pyramids please.
      Assert ((reference_cell == ReferenceCells::Triangle)
              ||
              (reference_cell == ReferenceCells::Quadrilateral)
              ||
              (reference_cell == ReferenceCells::Tetrahedron)
              ||
              (reference_cell == ReferenceCells::Hexahedron),
              ExcNotImplemented());
      return (reference_cell.n_faces() *
              reference_cell.face_reference_cell(0).n_isotropic_children());
    }



    unsigned int
    nth_interface_matrix (const ReferenceCell &reference_cell,
                          const unsigned int face)
    {
      AssertIndexRange (face, reference_cell.n_faces());
      return (face *
              reference_cell.face_reference_cell(0).n_isotropic_children());
    }



    unsigned int
    nth_interface_matrix (const ReferenceCell &reference_cell,
                          const unsigned int face,
                          const unsigned int sub_face)
    {
      AssertIndexRange (face, reference_cell.n_faces());
      AssertIndexRange (sub_face,
                        reference_cell.face_reference_cell(0).n_isotropic_children());
      return (face *
              reference_cell.face_reference_cell(0).n_isotropic_children()
              + sub_face);
    }



    template <int dim>
    void
    Interface<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &) const
    {}



    template <int dim>
    std::vector<double>
    Interface<dim>::compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &) const
    {
      AssertThrow(false,
                  ExcMessage("If this is an assembler that has to compute a residual"
                             "then this function should be implemented. "
                             "If this is an assembler that does not have to compute a residual, it should not be called."));

      return {};
    }



    template <int dim>
    AdvectionStabilizationInterface<dim>::~AdvectionStabilizationInterface()
      = default;



    template <int dim>
    std::vector<double>
    AdvectionStabilizationInterface<dim>::advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      std::vector<double> prefactors(scratch.material_model_inputs.n_evaluation_points(), 1.0);

      if (scratch.advection_field->is_temperature())
        for (unsigned int i=0; i<prefactors.size(); ++i)
          prefactors[i] = scratch.material_model_outputs.densities[i] * scratch.material_model_outputs.specific_heat[i];

      return prefactors;
    }



    template <int dim>
    std::vector<double>
    AdvectionStabilizationInterface<dim>::diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const
    {
      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::AdvectionSystem<dim>&> (scratch_base);

      if (scratch.advection_field->is_temperature())
        return scratch.material_model_outputs.thermal_conductivities;

      return std::vector<double> (scratch.material_model_inputs.n_evaluation_points(), 0.0);
    }



    template <int dim>
    void Manager<dim>::reset ()
    {
      stokes_preconditioner.clear();
      stokes_system.clear();
      stokes_system_on_boundary_face.clear();
      advection_system.clear();
      advection_system_on_boundary_face.clear();
      advection_system_on_interior_face.clear();
    }

    template <int dim>
    Manager<dim>::Properties::Properties ()
      :
      need_face_material_model_data (false),
      need_face_finite_element_evaluation(false),
      need_viscosity(false),
      needed_update_flags ()
    {}
  }
} // namespace aspect

// explicit instantiation of the functions we implement in this file
namespace aspect
{

#define INSTANTIATE(dim) \
  namespace internal { \
    namespace Assembly { \
      namespace Scratch { \
        template struct StokesPreconditioner<dim>; \
        template struct StokesSystem<dim>; \
        template struct AdvectionSystem<dim>; \
      } \
      namespace CopyData { \
        template struct StokesPreconditioner<dim>; \
        template struct StokesSystem<dim>; \
        template struct AdvectionSystem<dim>; \
      } \
    } \
  } \
  namespace Assemblers { \
    template class Interface<dim>; \
    template class AdvectionStabilizationInterface<dim>; \
    template class Manager<dim>; \
  }
  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}

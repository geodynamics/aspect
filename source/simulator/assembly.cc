/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <aspect/free_surface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>

#include <limits>


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
                              const bool                add_compaction_pressure)
          :
          finite_element_values (mapping, finite_element, quadrature,
                                 update_flags),
          local_dof_indices (finite_element.dofs_per_cell),
          dof_component_indices(stokes_dofs_per_cell),
          grads_phi_u (stokes_dofs_per_cell, numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          phi_p (stokes_dofs_per_cell, numbers::signaling_nan<double>()),
          phi_p_c (add_compaction_pressure ? stokes_dofs_per_cell : 0, numbers::signaling_nan<double>()),
          grad_phi_p (add_compaction_pressure ? stokes_dofs_per_cell : 0, numbers::signaling_nan<Tensor<1,dim> >()),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size(), n_compositional_fields)
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),
          local_dof_indices (scratch.local_dof_indices),
          dof_component_indices( scratch.dof_component_indices),
          grads_phi_u (scratch.grads_phi_u),
          phi_p (scratch.phi_p),
          phi_p_c (scratch.phi_p_c),
          grad_phi_p(scratch.grad_phi_p),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs)
        {}


        template <int dim>
        StokesPreconditioner<dim>::
        ~StokesPreconditioner ()
        {}




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
                      const bool                use_reference_density_profile)
          :
          StokesPreconditioner<dim> (finite_element, quadrature,
                                     mapping,
                                     update_flags,
                                     n_compositional_fields,
                                     stokes_dofs_per_cell,
                                     add_compaction_pressure),

          face_finite_element_values (mapping,
                                      finite_element,
                                      face_quadrature,
                                      face_update_flags),

          phi_u (stokes_dofs_per_cell, numbers::signaling_nan<Tensor<1,dim> >()),
          grads_phi_u (stokes_dofs_per_cell, numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          div_phi_u (stokes_dofs_per_cell, numbers::signaling_nan<double>()),
          velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          face_material_model_inputs(face_quadrature.size(), n_compositional_fields),
          face_material_model_outputs(face_quadrature.size(), n_compositional_fields),
          reference_densities(use_reference_density_profile ? quadrature.size() : 0, numbers::signaling_nan<double>()),
          reference_densities_depth_derivative(use_reference_density_profile ? quadrature.size() : 0, numbers::signaling_nan<double>())
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
          grads_phi_u (scratch.grads_phi_u),
          div_phi_u (scratch.div_phi_u),
          velocity_values (scratch.velocity_values),
          face_material_model_inputs(scratch.face_material_model_inputs),
          face_material_model_outputs(scratch.face_material_model_outputs),
          reference_densities(scratch.reference_densities),
          reference_densities_depth_derivative(scratch.reference_densities_depth_derivative)
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element,
                         const FiniteElement<dim> &advection_element,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &quadrature,
                         const Quadrature<dim-1>  &face_quadrature,
                         const UpdateFlags         update_flags,
                         const UpdateFlags         face_update_flags,
                         const unsigned int        n_compositional_fields)
          :
          finite_element_values (mapping,
                                 finite_element, quadrature,
                                 update_flags),
          face_finite_element_values ((face_quadrature.size() > 0
                                       ?
                                       new FEFaceValues<dim> (mapping,
                                                              finite_element, face_quadrature,
                                                              face_update_flags)
                                       :
                                       NULL)),
          neighbor_face_finite_element_values ((face_quadrature.size() > 0
                                                ?
                                                new FEFaceValues<dim> (mapping,
                                                                       finite_element, face_quadrature,
                                                                       face_update_flags)
                                                :
                                                NULL)),
          subface_finite_element_values ((face_quadrature.size() > 0
                                          ?
                                          new FESubfaceValues<dim> (mapping,
                                                                    finite_element, face_quadrature,
                                                                    face_update_flags)
                                          :
                                          NULL)),
          local_dof_indices (finite_element.dofs_per_cell),

          phi_field (advection_element.dofs_per_cell, numbers::signaling_nan<double>()),
          grad_phi_field (advection_element.dofs_per_cell, numbers::signaling_nan<Tensor<1,dim> >()),
          face_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                          numbers::signaling_nan<double>()),
          face_grad_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                               numbers::signaling_nan<Tensor<1,dim> >()),
          neighbor_face_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                                   numbers::signaling_nan<double>()),
          neighbor_face_grad_phi_field ((face_quadrature.size() > 0 ? advection_element.dofs_per_cell : 0),
                                        numbers::signaling_nan<Tensor<1,dim> >()),
          old_velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          old_old_velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          old_pressure (quadrature.size(), numbers::signaling_nan<double>()),
          old_old_pressure (quadrature.size(), numbers::signaling_nan<double>()),
          old_pressure_gradients (quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          old_old_pressure_gradients (quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          old_strain_rates (quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          old_old_strain_rates (quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          old_temperature_values (quadrature.size(), numbers::signaling_nan<double>()),
          old_old_temperature_values(quadrature.size(), numbers::signaling_nan<double>()),
          old_field_values (quadrature.size(), numbers::signaling_nan<double>()),
          old_old_field_values(quadrature.size(), numbers::signaling_nan<double>()),
          old_field_grads(quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          old_old_field_grads(quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          old_field_laplacians(quadrature.size(), numbers::signaling_nan<double>()),
          old_old_field_laplacians(quadrature.size(), numbers::signaling_nan<double>()),
          old_composition_values(n_compositional_fields,
                                 std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
          old_old_composition_values(n_compositional_fields,
                                     std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
          current_temperature_values(quadrature.size(), numbers::signaling_nan<double>()),
          current_velocity_values(quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          face_current_velocity_values(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          mesh_velocity_values(quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          face_mesh_velocity_values(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          current_strain_rates(quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim> >()),
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
          neighbor_face_heating_model_outputs(face_quadrature.size(), n_compositional_fields)
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const AdvectionSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),
          face_finite_element_values (scratch.face_finite_element_values),
          neighbor_face_finite_element_values (scratch.neighbor_face_finite_element_values),
          subface_finite_element_values (scratch.subface_finite_element_values),
          local_dof_indices (scratch.finite_element_values.get_fe().dofs_per_cell),

          phi_field (scratch.phi_field),
          grad_phi_field (scratch.grad_phi_field),
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
          neighbor_face_heating_model_outputs(scratch.neighbor_face_heating_model_outputs)
        {}

      }


      namespace CopyData
      {

        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const unsigned int stokes_dofs_per_cell)
          :
          local_matrix (stokes_dofs_per_cell,
                        stokes_dofs_per_cell),
          local_dof_indices (stokes_dofs_per_cell)
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &data)
          :
          local_matrix (data.local_matrix),
          local_dof_indices (data.local_dof_indices)
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        ~StokesPreconditioner ()
        {}



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
                                   GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                   :
                                   0),
                                  FullMatrix<double>(finite_element.dofs_per_cell,
                                                     finite_element.dofs_per_cell)),
          local_matrices_ext_int ((field_is_discontinuous
                                   ?
                                   GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                   :
                                   0),
                                  FullMatrix<double>(finite_element.dofs_per_cell,
                                                     finite_element.dofs_per_cell)),
          local_matrices_ext_ext ((field_is_discontinuous
                                   ?
                                   GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                   :
                                   0),
                                  FullMatrix<double>(finite_element.dofs_per_cell,
                                                     finite_element.dofs_per_cell)),
          local_rhs (finite_element.dofs_per_cell),

          assembled_matrices ((field_is_discontinuous
                               ?
                               GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                               :
                               0), false),

          local_dof_indices (finite_element.dofs_per_cell),
          neighbor_dof_indices ((field_is_discontinuous
                                 ?
                                 GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                 :
                                 0),
                                std::vector<types::global_dof_index>(finite_element.dofs_per_cell))
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const AdvectionSystem &data)
          :
          local_matrix (data.local_matrix),
          local_matrices_int_ext (data.local_matrices_int_ext),
          local_matrices_ext_int (data.local_matrices_ext_int),
          local_matrices_ext_ext (data.local_matrices_ext_ext),
          local_rhs (data.local_rhs),

          assembled_matrices (data.assembled_matrices),

          local_dof_indices (data.local_dof_indices),
          neighbor_dof_indices (data.neighbor_dof_indices)
        {}

      }



      template <int dim>
      AssemblerLists<dim>::Properties::Properties ()
        :
        need_face_material_model_data (false),
        need_face_finite_element_evaluation(false),
        need_viscosity(false),
        needed_update_flags ()
      {}

    }
  }




  template <int dim>
  void
  Simulator<dim>::
  compute_material_model_input_values (const LinearAlgebra::BlockVector                            &input_solution,
                                       const FEValuesBase<dim>                                     &input_finite_element_values,
                                       const typename DoFHandler<dim>::active_cell_iterator        &cell,
                                       const bool                                                   compute_strainrate,
                                       MaterialModel::MaterialModelInputs<dim>                     &material_model_inputs) const
  {
    const unsigned int n_q_points = material_model_inputs.temperature.size();

    material_model_inputs.position = input_finite_element_values.get_quadrature_points();

    input_finite_element_values[introspection.extractors.temperature].get_function_values (input_solution,
        material_model_inputs.temperature);
    input_finite_element_values[introspection.extractors.pressure].get_function_values(input_solution,
        material_model_inputs.pressure);
    input_finite_element_values[introspection.extractors.velocities].get_function_values(input_solution,
        material_model_inputs.velocity);
    input_finite_element_values[introspection.extractors.pressure].get_function_gradients (input_solution,
        material_model_inputs.pressure_gradient);

    // only the viscosity in the material can depend on the strain_rate
    // if this is not needed, we can save some time here. By setting the
    // length of the strain_rate vector to 0, we signal to evaluate()
    // that we do not need to access the viscosity.
    if (compute_strainrate)
      input_finite_element_values[introspection.extractors.velocities].get_function_symmetric_gradients(input_solution,
          material_model_inputs.strain_rate);
    else
      material_model_inputs.strain_rate.resize(0);

    // the values of the compositional fields are stored as blockvectors for each field
    // we have to extract them in this structure
    std::vector<std::vector<double> > composition_values (introspection.n_compositional_fields,
                                                          std::vector<double> (n_q_points));

    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      input_finite_element_values[introspection.extractors.compositional_fields[c]].get_function_values(input_solution,
          composition_values[c]);

    // then we copy these values to exchange the inner and outer vector, because for the material
    // model we need a vector with values of all the compositional fields for every quadrature point
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
        material_model_inputs.composition[q][c] = composition_values[c][q];

    material_model_inputs.cell = &cell;
  }



  namespace Assemblers
  {
    /**
     * A namespace for the definition of functions that implement various
     * other terms that need to occasionally or always be assembled.
     */
    namespace OtherTerms
    {
      template <int dim>
      void
      pressure_rhs_compatibility_modification (const SimulatorAccess<dim>                      &simulator_access,
                                               internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                               internal::Assembly::CopyData::StokesSystem<dim> &data)
      {
        const Introspection<dim> &introspection = simulator_access.introspection();
        const FiniteElement<dim> &fe = simulator_access.get_fe();

        const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  data.local_pressure_shape_function_integrals(i_stokes) += scratch.phi_p[i_stokes] * scratch.finite_element_values.JxW(q);
                  ++i_stokes;
                }
              ++i;
            }
      }



      template <int dim>
      void
      boundary_traction (const SimulatorAccess<dim>                           &simulator_access,
                         const typename DoFHandler<dim>::active_cell_iterator &cell,
                         const unsigned int                                    face_no,
                         internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                         internal::Assembly::CopyData::StokesSystem<dim>      &data)
      {
        const Introspection<dim> &introspection = simulator_access.introspection();
        const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();

        // see if any of the faces are traction boundaries for which
        // we need to assemble force terms for the right hand side
        const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();

        if (simulator_access.get_boundary_traction()
            .find (cell->face(face_no)->boundary_id())
            !=
            simulator_access.get_boundary_traction().end())
          {
            scratch.face_finite_element_values.reinit (cell, face_no);

            for (unsigned int q=0; q<scratch.face_finite_element_values.n_quadrature_points; ++q)
              {
                const Tensor<1,dim> traction
                  = simulator_access.get_boundary_traction().find(
                      cell->face(face_no)->boundary_id()
                    )->second
                    ->boundary_traction (cell->face(face_no)->boundary_id(),
                                         scratch.face_finite_element_values.quadrature_point(q),
                                         scratch.face_finite_element_values.normal_vector(q));

                for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
                  {
                    if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
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
  }

  template <int dim>
  void
  Simulator<dim>::
  create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
  {
    typedef typename std::vector<std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> > >
    assembler_vector_t;

    for (typename assembler_vector_t::const_iterator it = assembler_objects.begin();
         it != assembler_objects.end();
         ++it)
      {
        (*it)->create_additional_material_model_outputs(outputs);
      }

  }



  template <int dim>
  void
  Simulator<dim>::
  set_assemblers ()
  {
    // create an object for the complete equations assembly; add its
    // member functions to the signals and add the object the list
    // of assembler objects
    aspect::Assemblers::StokesAssembler<dim> *stokes_assembler
      = new aspect::Assemblers::StokesAssembler<dim>();
    aspect::Assemblers::AdvectionAssembler<dim> *adv_assembler
      = new aspect::Assemblers::AdvectionAssembler<dim>();

    assemblers->advection_system_assembler_properties.resize(1+introspection.n_compositional_fields);
    assemblers->advection_system_assembler_on_face_properties.resize(1+introspection.n_compositional_fields);

    aspect::Assemblers::MeltEquations<dim> *melt_equation_assembler = NULL;
    if (parameters.include_melt_transport)
      melt_equation_assembler = new aspect::Assemblers::MeltEquations<dim>();

    if (parameters.include_melt_transport)
      assemblers->local_assemble_stokes_preconditioner
      .connect (std_cxx11::bind(&aspect::Assemblers::MeltEquations<dim>::local_assemble_stokes_preconditioner_melt,
                                std_cxx11::cref (*melt_equation_assembler),
                                std_cxx11::_1, std_cxx11::_2, std_cxx11::_3));
    else
      assemblers->local_assemble_stokes_preconditioner
      .connect (std_cxx11::bind(&aspect::Assemblers::StokesAssembler<dim>::preconditioner,
                                std_cxx11::cref (*stokes_assembler),
                                std_cxx11::_1, std_cxx11::_2, std_cxx11::_3));

    if (parameters.include_melt_transport)
      assemblers->local_assemble_stokes_system
      .connect (std_cxx11::bind(&aspect::Assemblers::MeltEquations<dim>::local_assemble_stokes_system_melt,
                                std_cxx11::cref (*melt_equation_assembler),
                                std_cxx11::_1,
                                std_cxx11::_2,
                                std_cxx11::_3,
                                std_cxx11::_4,
                                std_cxx11::_5));
    else
      {
        assemblers->local_assemble_stokes_system
        .connect (std_cxx11::bind(&aspect::Assemblers::StokesAssembler<dim>::incompressible_terms,
                                  std_cxx11::cref (*stokes_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3,
                                  std_cxx11::_4,
                                  std_cxx11::_5));

        if (material_model->is_compressible())
          assemblers->local_assemble_stokes_system
          .connect (std_cxx11::bind(&aspect::Assemblers::StokesAssembler<dim>::compressible_strain_rate_viscosity_term,
                                    std_cxx11::cref (*stokes_assembler),
                                    // discard cell,
                                    std_cxx11::_2,
                                    std_cxx11::_3,
                                    std_cxx11::_4,
                                    std_cxx11::_5));

        if (parameters.formulation_mass_conservation ==
            Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile)
          assemblers->local_assemble_stokes_system
          .connect (std_cxx11::bind(&aspect::Assemblers::StokesAssembler<dim>::implicit_reference_density_compressibility_term,
                                    std_cxx11::cref (*stokes_assembler),
                                    // discard cell,
                                    std_cxx11::_2,
                                    std_cxx11::_3,
                                    std_cxx11::_4,
                                    std_cxx11::_5,
                                    std_cxx11::cref (this->parameters)));
        else if (parameters.formulation_mass_conservation ==
                 Parameters<dim>::Formulation::MassConservation::reference_density_profile)
          {
            assemblers->local_assemble_stokes_system
            .connect (std_cxx11::bind(&aspect::Assemblers::StokesAssembler<dim>::reference_density_compressibility_term,
                                      std_cxx11::cref (*stokes_assembler),
                                      // discard cell,
                                      std_cxx11::_2,
                                      std_cxx11::_3,
                                      std_cxx11::_4,
                                      std_cxx11::_5,
                                      std_cxx11::cref (this->parameters)));
          }
        else if (parameters.formulation_mass_conservation ==
                 Parameters<dim>::Formulation::MassConservation::incompressible)
          {
            // do nothing, because we assembled div u =0 above already
          }
        else
          assemblers->local_assemble_stokes_system
          .connect (std_cxx11::bind(&aspect::Assemblers::StokesAssembler<dim>::isothermal_compression_term,
                                    std_cxx11::cref (*stokes_assembler),
                                    // discard cell,
                                    std_cxx11::_2,
                                    std_cxx11::_3,
                                    std_cxx11::_4,
                                    std_cxx11::_5,
                                    std_cxx11::cref (this->parameters)));

      }


    assembler_objects.push_back (std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >
                                 (stokes_assembler));

    assembler_objects.push_back (std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >
                                 (adv_assembler));


    if (parameters.include_melt_transport)
      assembler_objects.push_back (std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >
                                   (melt_equation_assembler));

    // add the boundary integral for melt migration
    if (parameters.include_melt_transport)
      {
        assemblers->stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data = true;
        assemblers->stokes_system_assembler_on_boundary_face_properties.needed_update_flags = (update_values  | update_quadrature_points |
            update_normal_vectors | update_gradients |
            update_JxW_values);
        assemblers->local_assemble_stokes_system_on_boundary_face
        .connect (std_cxx11::bind(&aspect::Assemblers::MeltEquations<dim>::local_assemble_stokes_system_melt_boundary,
                                  std_cxx11::cref (*melt_equation_assembler),
                                  std_cxx11::_1,
                                  std_cxx11::_2,
                                  std_cxx11::_3,
                                  // discard rebuild_stokes_matrix,
                                  std_cxx11::_5,
                                  std_cxx11::_6));
      }

    // add the terms for traction boundary conditions
    assemblers->local_assemble_stokes_system_on_boundary_face
    .connect (std_cxx11::bind(&aspect::Assemblers::OtherTerms::boundary_traction<dim>,
                              SimulatorAccess<dim>(*this),
                              std_cxx11::_1,
                              std_cxx11::_2,
                              // discard pressure_scaling,
                              // discard rebuild_stokes_matrix,
                              std_cxx11::_5,
                              std_cxx11::_6));

    // add the terms necessary to normalize the pressure
    if (do_pressure_rhs_compatibility_modification)
      {
        if (parameters.include_melt_transport)
          assemblers->local_assemble_stokes_system
          .connect (std_cxx11::bind(&aspect::Assemblers::OtherTerms::pressure_rhs_compatibility_modification_melt<dim>,
                                    SimulatorAccess<dim>(*this),
                                    // discard cell,
                                    // discard pressure_scaling,
                                    // discard rebuild_stokes_matrix,
                                    std_cxx11::_4,
                                    std_cxx11::_5));
        else
          assemblers->local_assemble_stokes_system
          .connect (std_cxx11::bind(&aspect::Assemblers::OtherTerms::pressure_rhs_compatibility_modification<dim>,
                                    SimulatorAccess<dim>(*this),
                                    // discard cell,
                                    // discard pressure_scaling,
                                    // discard rebuild_stokes_matrix,
                                    std_cxx11::_4,
                                    std_cxx11::_5));
      }

    if (parameters.include_melt_transport)
      {
        assemblers->local_assemble_advection_system
        .connect (std_cxx11::bind(&aspect::Assemblers::MeltEquations<dim>::local_assemble_advection_system_melt,
                                  std_cxx11::cref (*melt_equation_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3,
                                  std_cxx11::_4,
                                  std_cxx11::_5));

        assemblers->compute_advection_system_residual
        .connect (std_cxx11::bind(&aspect::Assemblers::MeltEquations<dim>::compute_advection_system_residual_melt,
                                  std_cxx11::cref (*melt_equation_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3));
      }
    else
      {
        assemblers->local_assemble_advection_system
        .connect (std_cxx11::bind(&aspect::Assemblers::AdvectionAssembler<dim>::local_assemble_advection_system,
                                  std_cxx11::cref (*adv_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3,
                                  std_cxx11::_4,
                                  std_cxx11::_5));

        assemblers->compute_advection_system_residual
        .connect (std_cxx11::bind(&aspect::Assemblers::AdvectionAssembler<dim>::compute_advection_system_residual,
                                  std_cxx11::cref (*adv_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3));

      }

    if (parameters.use_discontinuous_temperature_discretization ||
        parameters.use_discontinuous_composition_discretization)
      {
        assemblers->local_assemble_advection_system_on_interior_face
        .connect(std_cxx11::bind(&aspect::Assemblers::AdvectionAssembler<dim>::local_assemble_discontinuous_advection_interior_face_terms,
                                 std_cxx11::cref (*adv_assembler),
                                 std_cxx11::_1,
                                 std_cxx11::_2,
                                 std_cxx11::_3,
                                 std_cxx11::_4,
                                 std_cxx11::_5));

        assemblers->local_assemble_advection_system_on_boundary_face
        .connect(std_cxx11::bind(&aspect::Assemblers::AdvectionAssembler<dim>::local_assemble_discontinuous_advection_boundary_face_terms,
                                 std_cxx11::cref (*adv_assembler),
                                 std_cxx11::_1,
                                 std_cxx11::_2,
                                 std_cxx11::_3,
                                 std_cxx11::_4,
                                 std_cxx11::_5));

        if (parameters.use_discontinuous_temperature_discretization)
          {
            assemblers->advection_system_assembler_on_face_properties[0].need_face_material_model_data = true;
            assemblers->advection_system_assembler_on_face_properties[0].need_face_finite_element_evaluation = true;
          }

        if (parameters.use_discontinuous_composition_discretization)
          {
            for (unsigned int i = 1; i<=introspection.n_compositional_fields; ++i)
              {
                assemblers->advection_system_assembler_on_face_properties[i].need_face_material_model_data = true;
                assemblers->advection_system_assembler_on_face_properties[i].need_face_finite_element_evaluation = true;
              }
          }
      }

    // allow other assemblers to add themselves or modify the existing ones by firing the signal
    this->signals.set_assemblers(*this, *assemblers, assembler_objects);

    // ensure that all assembler objects have access to the SimulatorAccess
    // base class
    for (unsigned int i=0; i<assembler_objects.size(); ++i)
      if (SimulatorAccess<dim> *p = dynamic_cast<SimulatorAccess<dim>*>(assembler_objects[i].get()))
        p->initialize_simulator(*this);
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                        internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    // First get all dof indices of the current cell, then extract those
    // that correspond to the Stokes system we are interested in.
    // Note that assemblers below can modify this list of dofs, if they in fact
    // assemble a different system than the standard Stokes system (e.g. in
    // models with melt transport).

    cell->get_dof_indices (scratch.local_dof_indices);

    const unsigned int dofs_per_cell = finite_element.dofs_per_cell;

    for (unsigned int i=0, i_stokes=0; i<dofs_per_cell; ++i)
      if (introspection.is_stokes_component(finite_element.system_to_component_index(i).first))
        {
          data.local_dof_indices[i_stokes] = scratch.local_dof_indices[i];
          ++i_stokes;
        }

    // Prepare the data structures for assembly
    scratch.finite_element_values.reinit (cell);

    data.local_matrix = 0;

    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         cell,
                                         true,
                                         scratch.material_model_inputs);
    create_additional_material_model_outputs(scratch.material_model_outputs);

    material_model->evaluate(scratch.material_model_inputs,
                             scratch.material_model_outputs);
    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.finite_element_values.get_quadrature(),
                                               scratch.finite_element_values.get_mapping(),
                                               scratch.material_model_outputs);

    // trigger the invocation of the various functions that actually do
    // all of the assembling
    assemblers->local_assemble_stokes_preconditioner(pressure_scaling, scratch, data);
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_dof_indices,
                                                    system_preconditioner_matrix);
  }



  template <int dim>
  void
  Simulator<dim>::assemble_stokes_preconditioner ()
  {
    system_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    // determine which update flags to use for the cell integrals
    const UpdateFlags cell_update_flags
      = ((update_JxW_values |
          update_values |
          update_gradients |
          update_quadrature_points)
         |
         assemblers->stokes_preconditioner_assembler_properties.needed_update_flags);

    unsigned int stokes_dofs_per_cell = dim * finite_element.base_element(introspection.base_elements.velocities).dofs_per_cell
                                        + finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell;

    if (parameters.include_melt_transport)
      stokes_dofs_per_cell += finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx11::bind (&Simulator<dim>::
                          local_assemble_stokes_preconditioner,
                          this,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&Simulator<dim>::
                          copy_local_to_global_stokes_preconditioner,
                          this,
                          std_cxx11::_1),
         internal::Assembly::Scratch::
         StokesPreconditioner<dim> (finite_element, quadrature_formula,
                                    *mapping,
                                    cell_update_flags,
                                    introspection.n_compositional_fields,
                                    stokes_dofs_per_cell,
                                    parameters.include_melt_transport),
         internal::Assembly::CopyData::
         StokesPreconditioner<dim> (stokes_dofs_per_cell));

    system_preconditioner_matrix.compress(VectorOperation::add);
  }



  template <int dim>
  void
  Simulator<dim>::build_stokes_preconditioner ()
  {
    if (rebuild_stokes_preconditioner == false)
      return;

    if (parameters.use_direct_stokes_solver)
      return;

    computing_timer.enter_section ("   Build Stokes preconditioner");
    pcout << "   Rebuilding Stokes preconditioner..." << std::flush;

    // first assemble the raw matrices necessary for the preconditioner
    assemble_stokes_preconditioner ();

    // then extract the other information necessary to build the
    // AMG preconditioners for the A and M blocks
    std::vector<std::vector<bool> > constant_modes;
    DoFTools::extract_constant_modes (dof_handler,
                                      introspection.component_masks.velocities,
                                      constant_modes);

    Mp_preconditioner.reset (new LinearAlgebra::PreconditionILU());
    Amg_preconditioner.reset (new LinearAlgebra::PreconditionAMG());

    LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
#ifdef ASPECT_USE_PETSC
    Amg_data.symmetric_operator = false;
#else
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = true;

    // set the AMG parameters in a way that minimizes the run
    // time. compared to some of the deal.II tutorial programs, we
    // found that it pays off to set the aggregration threshold to
    // zero, especially for ill-conditioned problems with large
    // variations in the viscosity
    //
    // for extensive benchmarking of various settings of these
    // parameters and others, see
    // https://github.com/geodynamics/aspect/pull/234
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.001;
#endif

    /*  The stabilization term for the free surface (Kaus et. al., 2010)
     *  makes changes to the system matrix which are of the same form as
     *  boundary stresses.  If these stresses are not also added to the
     *  system_preconditioner_matrix, then  if fails to be very good as a
     *  preconditioner.  Instead, we just pass the system_matrix to the
     *  AMG precondition initialization so that it builds the preconditioner
     *  directly from that. However, we still need the mass matrix for the
     *  pressure block which is assembled in the preconditioner matrix.
     *  So rather than build a different preconditioner matrix which only
     *  does the mass matrix, we just reuse the same system_preconditioner_matrix
     *  for the Mp_preconditioner block.  Maybe a bit messy*/
    Mp_preconditioner->initialize (system_preconditioner_matrix.block(1,1));
    if (parameters.free_surface_enabled)
      Amg_preconditioner->initialize (system_matrix.block(0,0),
                                      Amg_data);
    else
      Amg_preconditioner->initialize (system_preconditioner_matrix.block(0,0),
                                      Amg_data);

    rebuild_stokes_preconditioner = false;

    pcout << std::endl;
    computing_timer.exit_section();
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                internal::Assembly::Scratch::StokesSystem<dim> &scratch,
                                internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    // First get all dof indices of the current cell, then extract those
    // that correspond to the Stokes system we are interested in.
    // Note that assemblers below can modify this list of dofs, if they in fact
    // assemble a different system than the standard Stokes system (e.g. in
    // models with melt transport).

    cell->get_dof_indices (scratch.local_dof_indices);

    const unsigned int dofs_per_cell = finite_element.dofs_per_cell;

    for (unsigned int i=0, i_stokes=0; i<dofs_per_cell; /*increment at end of loop*/)
      {
        if (introspection.is_stokes_component(finite_element.system_to_component_index(i).first))
          {
            data.local_dof_indices[i_stokes] = scratch.local_dof_indices[i];
            ++i_stokes;
          }
        ++i;
      }


    // Prepare the data structures for assembly
    scratch.finite_element_values.reinit (cell);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    if (do_pressure_rhs_compatibility_modification)
      data.local_pressure_shape_function_integrals = 0;

    // initialize the material model data on the cell
    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         cell,
                                         rebuild_stokes_matrix,
                                         scratch.material_model_inputs);
    create_additional_material_model_outputs(scratch.material_model_outputs);

    material_model->evaluate(scratch.material_model_inputs,
                             scratch.material_model_outputs);
    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.finite_element_values.get_quadrature(),
                                               scratch.finite_element_values.get_mapping(),
                                               scratch.material_model_outputs);

    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.velocity_values);

    const bool use_reference_density_profile = (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::reference_density_profile)
                                               || (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile);
    if (use_reference_density_profile)
      {
        for (unsigned int q=0; q<scratch.finite_element_values.n_quadrature_points; ++q)
          {
            scratch.reference_densities[q] = adiabatic_conditions->density(scratch.material_model_inputs.position[q]);
            scratch.reference_densities_depth_derivative[q] = adiabatic_conditions->density_derivative(scratch.material_model_inputs.position[q]);
          }
      }

    // trigger the invocation of the various functions that actually do
    // all of the assembling
    assemblers->local_assemble_stokes_system(cell, pressure_scaling, rebuild_stokes_matrix,
                                             scratch, data);

    // then also work on possible face terms. if necessary, initialize
    // the material model data on faces
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      if (cell->at_boundary(face_no))
        {
          scratch.face_finite_element_values.reinit (cell, face_no);

          if (assemblers->stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data)
            {
              const bool need_viscosity = rebuild_stokes_matrix |
                                          assemblers->stokes_system_assembler_on_boundary_face_properties.need_viscosity;

              compute_material_model_input_values (current_linearization_point,
                                                   scratch.face_finite_element_values,
                                                   cell,
                                                   need_viscosity,
                                                   scratch.face_material_model_inputs);
              create_additional_material_model_outputs(scratch.face_material_model_outputs);

              material_model->evaluate(scratch.face_material_model_inputs,
                                       scratch.face_material_model_outputs);

              // TODO: Currently we do not supply reference density values to Stokes face assemblers.
              // This seems acceptable for now, since the only face assemblers are the ones for the melt
              // assembly where one would not want to use the reference density anyway. In case the reference
              // density is ever needed, those assemblers can also query the adiabatic conditions for the
              // reference density.

              // TODO: the following doesn't currently compile because the get_quadrature() call returns
              //  a dim-1 dimensional quadrature
              // MaterialModel::MaterialAveraging::average (parameters.material_averaging,
              //                                            cell,
              //                                            compressibility * density *
              //                                            scratch.face_finite_element_values.get_mapping(),
              //                                            scratch.face_material_model_outputs);
            }

          assemblers->local_assemble_stokes_system_on_boundary_face(cell, face_no,
                                                                    pressure_scaling, rebuild_stokes_matrix,
                                                                    scratch, data);
        }
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    if (rebuild_stokes_matrix == true)
      current_constraints.distribute_local_to_global (data.local_matrix,
                                                      data.local_rhs,
                                                      data.local_dof_indices,
                                                      system_matrix,
                                                      system_rhs);
    else
      current_constraints.distribute_local_to_global (data.local_rhs,
                                                      data.local_dof_indices,
                                                      system_rhs);

    if (do_pressure_rhs_compatibility_modification)
      current_constraints.distribute_local_to_global (data.local_pressure_shape_function_integrals,
                                                      data.local_dof_indices,
                                                      pressure_shape_function_integrals);
  }



  template <int dim>
  void Simulator<dim>::assemble_stokes_system ()
  {
    computing_timer.enter_section ("   Assemble Stokes system");

    if (rebuild_stokes_matrix == true)
      system_matrix = 0;

    system_rhs = 0;
    if (do_pressure_rhs_compatibility_modification)
      pressure_shape_function_integrals = 0;

    const QGauss<dim>   quadrature_formula(parameters.stokes_velocity_degree+1);
    const QGauss<dim-1> face_quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    // determine which updates flags we need on cells and faces
    const UpdateFlags cell_update_flags
      = (update_values    |
         update_gradients |
         update_quadrature_points  |
         update_JxW_values)
        |
        assemblers->stokes_system_assembler_properties.needed_update_flags;
    const UpdateFlags face_update_flags
      = (
          // see if we need to assemble traction boundary conditions.
          // only if so do we actually need to have an FEFaceValues object
          parameters.prescribed_traction_boundary_indicators.size() > 0
          ?
          update_values |
          update_quadrature_points |
          update_normal_vectors |
          update_JxW_values
          :
          update_default)
        |
        (assemblers->stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data
         ?
         // if we need a material model input on the faces, we need to
         // also be able to compute the strain rate
         update_gradients
         :
         update_default)
        |
        assemblers->stokes_system_assembler_on_boundary_face_properties.needed_update_flags;

    unsigned int stokes_dofs_per_cell = dim * finite_element.base_element(introspection.base_elements.velocities).dofs_per_cell
                                        + finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell;

    if (parameters.include_melt_transport)
      stokes_dofs_per_cell += finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell;

    const bool use_reference_density_profile = (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::reference_density_profile)
                                               || (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile);

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx11::bind (&Simulator<dim>::
                          local_assemble_stokes_system,
                          this,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&Simulator<dim>::
                          copy_local_to_global_stokes_system,
                          this,
                          std_cxx11::_1),
         internal::Assembly::Scratch::
         StokesSystem<dim> (finite_element, *mapping, quadrature_formula,
                            face_quadrature_formula,
                            cell_update_flags,
                            face_update_flags,
                            introspection.n_compositional_fields,
                            stokes_dofs_per_cell,
                            parameters.include_melt_transport,
                            use_reference_density_profile),
         internal::Assembly::CopyData::
         StokesSystem<dim> (stokes_dofs_per_cell,
                            do_pressure_rhs_compatibility_modification));

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    // if the model is compressible then we need to adjust the right hand
    // side of the equation to make it compatible with the matrix on the
    // left
    if (do_pressure_rhs_compatibility_modification)
      {
        pressure_shape_function_integrals.compress(VectorOperation::add);
        make_pressure_rhs_compatible(system_rhs);
      }


    // record that we have just rebuilt the matrix
    rebuild_stokes_matrix = false;

    computing_timer.exit_section();
  }

  template <int dim>
  void
  Simulator<dim>::build_advection_preconditioner(const AdvectionField &advection_field,
                                                 std_cxx11::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner)
  {
    switch (advection_field.field_type)
      {
        case AdvectionField::temperature_field:
        {
          computing_timer.enter_section ("   Build temperature preconditioner");
          break;
        }

        case AdvectionField::compositional_field:
        {
          computing_timer.enter_section ("   Build composition preconditioner");
          break;
        }

        default:
          Assert (false, ExcNotImplemented());
      }

    const unsigned int block_idx = advection_field.block_index(introspection);
    preconditioner.reset (new LinearAlgebra::PreconditionILU());
    preconditioner->initialize (system_matrix.block(block_idx, block_idx));
    computing_timer.exit_section();
  }


  template <int dim>
  void Simulator<dim>::
  local_assemble_advection_system (const AdvectionField     &advection_field,
                                   const Vector<double>           &viscosity_per_cell,
                                   const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                   internal::Assembly::CopyData::AdvectionSystem<dim> &data)
  {
    // also have the number of dofs that correspond just to the element for
    // the system we are currently trying to assemble
    const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

    Assert (advection_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
    Assert (scratch.grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
    Assert (scratch.phi_field.size() == advection_dofs_per_cell, ExcInternalError());

    const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

    const unsigned int solution_component = advection_field.component_index(introspection);

    scratch.finite_element_values.reinit (cell);

    // get all dof indices on the current cell, then extract those
    // that correspond to the solution_field we are interested in
    cell->get_dof_indices (scratch.local_dof_indices);
    for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
      {
        if (finite_element.system_to_component_index(i).first == solution_component)
          {
            data.local_dof_indices[i_advection] = scratch.local_dof_indices[i];
            ++i_advection;
          }
        ++i;
      }

    data.local_matrix = 0;
    data.local_rhs = 0;

    scratch.finite_element_values[solution_field].get_function_values (old_solution,
                                                                       scratch.old_field_values);
    scratch.finite_element_values[solution_field].get_function_values (old_old_solution,
                                                                       scratch.old_old_field_values);


    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.current_velocity_values);

    if (parameters.include_melt_transport)
      scratch.finite_element_values[introspection.extractors.velocities].get_function_divergences(current_linearization_point,
          scratch.current_velocity_divergences);

    // get the mesh velocity, as we need to subtract it off of the advection systems
    if (parameters.free_surface_enabled)
      scratch.finite_element_values[introspection.extractors.velocities].get_function_values(free_surface->mesh_velocity,
          scratch.mesh_velocity_values);

    // compute material properties and heating terms
    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         cell,
                                         true,
                                         scratch.material_model_inputs);
    create_additional_material_model_outputs(scratch.material_model_outputs);

    material_model->evaluate(scratch.material_model_inputs,
                             scratch.material_model_outputs);
    if (parameters.formulation_temperature_equation ==
        Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
      {
        const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            scratch.material_model_outputs.densities[q] = adiabatic_conditions->density(scratch.material_model_inputs.position[q]);
          }
      }

    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.finite_element_values.get_quadrature(),
                                               scratch.finite_element_values.get_mapping(),
                                               scratch.material_model_outputs);

    heating_model_manager.evaluate(scratch.material_model_inputs,
                                   scratch.material_model_outputs,
                                   scratch.heating_model_outputs);

    // TODO: Compute artificial viscosity once per timestep instead of each time
    // temperature system is assembled (as this might happen more than once per
    // timestep for iterative solvers)
    const double nu = viscosity_per_cell[cell->active_cell_index()];
    Assert (nu >= 0, ExcMessage ("The artificial viscosity needs to be a non-negative quantity."));

    // trigger the invocation of the various functions that actually do
    // all of the assembling
    assemblers->local_assemble_advection_system(cell, advection_field, nu, scratch, data);

    // then also work on possible face terms. if necessary, initialize
    // the material model data on faces
    const bool has_boundary_face_assemblers = !assemblers->local_assemble_advection_system_on_boundary_face.empty()
                                              && assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;
    const bool has_interior_face_assemblers = !assemblers->local_assemble_advection_system_on_interior_face.empty()
                                              && assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;

    // skip the remainder if no work needs to be done on faces
    if (!has_boundary_face_assemblers && !has_interior_face_assemblers)
      return;

    if (has_interior_face_assemblers)
      {
        // for interior face contributions loop over all possible
        // subfaces of the cell, and reset their matrices.
        for (unsigned int f = 0; f < GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell; ++f)
          {
            data.local_matrices_ext_int[f] = 0;
            data.local_matrices_int_ext[f] = 0;
            data.local_matrices_ext_ext[f] = 0;
            data.assembled_matrices[f] = false;
          }
      }

    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        const typename DoFHandler<dim>::face_iterator face = cell->face (face_no);

        if ((has_boundary_face_assemblers && face->at_boundary()) ||
            (has_interior_face_assemblers && !face->at_boundary()))
          {
            (*scratch.face_finite_element_values).reinit (cell, face_no);

            (*scratch.face_finite_element_values)[introspection.extractors.velocities].get_function_values(current_linearization_point,
                scratch.face_current_velocity_values);

            //get the mesh velocity, as we need to subtract it off of the advection systems
            if (parameters.free_surface_enabled)
              (*scratch.face_finite_element_values)[introspection.extractors.velocities].get_function_values(free_surface->mesh_velocity,
                  scratch.face_mesh_velocity_values);

            if (assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_material_model_data)
              {
                compute_material_model_input_values (current_linearization_point,
                                                     *scratch.face_finite_element_values,
                                                     cell,
                                                     true,
                                                     scratch.face_material_model_inputs);

                create_additional_material_model_outputs(scratch.face_material_model_outputs);

                material_model->evaluate(scratch.face_material_model_inputs,
                                         scratch.face_material_model_outputs);

                heating_model_manager.evaluate(scratch.face_material_model_inputs,
                                               scratch.face_material_model_outputs,
                                               scratch.face_heating_model_outputs);
                //TODO: the following doesn't currently compile because the get_quadrature() call returns
                //  a dim-1 dimensional quadrature
                // MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                //                                            cell,
                //                                            scratch.face_finite_element_values.get_quadrature(),
                //                                            scratch.face_finite_element_values.get_mapping(),
                //                                            scratch.face_material_model_outputs);
              }

            if (face->at_boundary())
              assemblers->local_assemble_advection_system_on_boundary_face(cell, face_no,
                                                                           advection_field,
                                                                           scratch, data);
            else
              assemblers->local_assemble_advection_system_on_interior_face(cell, face_no,
                                                                           advection_field,
                                                                           scratch, data);
          }
      }
  }

  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_advection_system (const AdvectionField &advection_field,
                                         const internal::Assembly::CopyData::AdvectionSystem<dim> &data)
  {
    // copy entries into the global matrix. note that these local contributions
    // only correspond to the advection dofs, as assembled above
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);

    /* In the following, we copy DG contributions element by element. This
     * is allowed since there are no constraints imposed on discontinuous fields.
     */
    if (!assemblers->local_assemble_advection_system_on_interior_face.empty() &&
        assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation)
      {
        for (unsigned int f=0; f<GeometryInfo<dim>::max_children_per_face
             * GeometryInfo<dim>::faces_per_cell; ++f)
          {
            if (data.assembled_matrices[f])
              {
                for (unsigned int i=0; i<data.local_dof_indices.size(); ++i)
                  for (unsigned int j=0; j<data.neighbor_dof_indices[f].size(); ++j)
                    {
                      system_matrix.add (data.local_dof_indices[i],
                                         data.neighbor_dof_indices[f][j],
                                         data.local_matrices_int_ext[f](i,j));
                      system_matrix.add (data.neighbor_dof_indices[f][j],
                                         data.local_dof_indices[i],
                                         data.local_matrices_ext_int[f](j,i));
                    }

                for (unsigned int i=0; i<data.neighbor_dof_indices[f].size(); ++i)
                  for (unsigned int j=0; j<data.neighbor_dof_indices[f].size(); ++j)
                    system_matrix.add (data.neighbor_dof_indices[f][i],
                                       data.neighbor_dof_indices[f][j],
                                       data.local_matrices_ext_ext[f](i,j));
              }
          }
      }
  }



  template <int dim>
  void Simulator<dim>::assemble_advection_system (const AdvectionField &advection_field)
  {
    if (advection_field.is_temperature())
      computing_timer.enter_section ("   Assemble temperature system");
    else
      computing_timer.enter_section ("   Assemble composition system");

    const unsigned int block_idx = advection_field.block_index(introspection);
    system_matrix.block(block_idx, block_idx) = 0;
    system_rhs = 0;

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    Vector<double> viscosity_per_cell;
    viscosity_per_cell.reinit(triangulation.n_active_cells());
    get_artificial_viscosity(viscosity_per_cell, advection_field);

    // We have to assemble the term u.grad phi_i * phi_j, which is
    // of total polynomial degree
    //   stokes_deg + 2*temp_deg -1
    // (or similar for comp_deg). This suggests using a Gauss
    // quadrature formula of order
    //   temp_deg + stokes_deg/2
    // rounded up (note that x/2 rounded up
    // equals (x+1)/2 using integer division.)
    //
    // (Note: All compositional fields have the same base element and therefore
    // the same composition_degree. Thus, we do not need to find out the degree
    // of the current field, but use the global instead)
    const unsigned int advection_quadrature_degree = advection_field.polynomial_degree(introspection)
                                                     +
                                                     (parameters.stokes_velocity_degree+1)/2;

    const bool allocate_face_quadrature = (!assemblers->local_assemble_advection_system_on_boundary_face.empty() ||
                                           !assemblers->local_assemble_advection_system_on_interior_face.empty()) &&
                                          assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;
    const bool allocate_neighbor_contributions = !assemblers->local_assemble_advection_system_on_interior_face.empty() &&
                                                 assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;;

    const UpdateFlags update_flags = update_values |
                                     update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values;

    const UpdateFlags face_update_flags = (allocate_face_quadrature ?
                                           update_values |
                                           update_gradients |
                                           update_quadrature_points |
                                           update_normal_vectors |
                                           update_JxW_values
                                           :
                                           update_default);

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx11::bind (&Simulator<dim>::
                          local_assemble_advection_system,
                          this,
                          advection_field,
                          std_cxx11::cref(viscosity_per_cell),
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&Simulator<dim>::
                          copy_local_to_global_advection_system,
                          this,
                          std_cxx11::cref(advection_field),
                          std_cxx11::_1),
         internal::Assembly::Scratch::
         AdvectionSystem<dim> (finite_element,
                               finite_element.base_element(advection_field.base_element(introspection)),
                               *mapping,
                               QGauss<dim>(advection_quadrature_degree),
                               /* Only generate a valid face quadrature if necessary.
                                * Otherwise, generate invalid face quadrature rule.
                                */
                               (allocate_face_quadrature ?
                                QGauss<dim-1>(advection_quadrature_degree) :
                                Quadrature<dim-1> ()),
                               update_flags,
                               face_update_flags,
                               introspection.n_compositional_fields),
         internal::Assembly::CopyData::
         AdvectionSystem<dim> (finite_element.base_element(advection_field.base_element(introspection)),
                               allocate_neighbor_contributions));

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    computing_timer.exit_section();
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  namespace internal { \
    namespace Assembly { \
      namespace Scratch { \
        template struct AdvectionSystem<dim>; \
        template struct StokesSystem<dim>; \
      } \
      template struct AssemblerLists<dim>::Properties; \
    } \
  } \
  template void Simulator<dim>::set_assemblers (); \
  template void Simulator<dim>::local_assemble_stokes_preconditioner ( \
                                                                       const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                       internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch, \
                                                                       internal::Assembly::CopyData::StokesPreconditioner<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_stokes_preconditioner ( \
                                                                             const internal::Assembly::CopyData::StokesPreconditioner<dim> &data); \
  template void Simulator<dim>::assemble_stokes_preconditioner (); \
  template void Simulator<dim>::build_stokes_preconditioner (); \
  template void Simulator<dim>::local_assemble_stokes_system ( \
                                                               const DoFHandler<dim>::active_cell_iterator &cell, \
                                                               internal::Assembly::Scratch::StokesSystem<dim>  &scratch, \
                                                               internal::Assembly::CopyData::StokesSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_stokes_system ( \
                                                                     const internal::Assembly::CopyData::StokesSystem<dim> &data); \
  template void Simulator<dim>::assemble_stokes_system (); \
  template void Simulator<dim>::build_advection_preconditioner (const AdvectionField &, \
                                                                std_cxx11::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner); \
  template void Simulator<dim>::local_assemble_advection_system ( \
                                                                  const AdvectionField          &advection_field, \
                                                                  const Vector<double>           &viscosity_per_cell, \
                                                                  const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                  internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch, \
                                                                  internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_advection_system ( \
                                                                        const AdvectionField          &advection_field, \
                                                                        const internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::assemble_advection_system (const AdvectionField     &advection_field); \
  template void Simulator<dim>::compute_material_model_input_values ( \
                                                                      const LinearAlgebra::BlockVector                      &input_solution, \
                                                                      const FEValuesBase<dim,dim>                           &input_finite_element_values, \
                                                                      const DoFHandler<dim>::active_cell_iterator  &cell, \
                                                                      const bool                                             compute_strainrate, \
                                                                      MaterialModel::MaterialModelInputs<dim>               &material_model_inputs) const; \
  template void Simulator<dim>::create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const; \
   


  ASPECT_INSTANTIATE(INSTANTIATE)
}

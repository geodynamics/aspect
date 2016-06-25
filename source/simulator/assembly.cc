/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
                              const bool                add_compaction_pressure)
          :
          finite_element_values (mapping, finite_element, quadrature,
                                 update_flags),
          grads_phi_u (finite_element.dofs_per_cell, numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          phi_p (finite_element.dofs_per_cell, numbers::signaling_nan<double>()),
          phi_p_c (add_compaction_pressure ? finite_element.dofs_per_cell : 0, numbers::signaling_nan<double>()),
          grad_phi_p (add_compaction_pressure ? finite_element.dofs_per_cell : 0, numbers::signaling_nan<Tensor<1,dim> >()),
          temperature_values (quadrature.size(), numbers::signaling_nan<double>()),
          pressure_values (quadrature.size(), numbers::signaling_nan<double>()),
          strain_rates (quadrature.size(), numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          composition_values(n_compositional_fields,
                             std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
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
          grads_phi_u (scratch.grads_phi_u),
          phi_p (scratch.phi_p),
          phi_p_c (scratch.phi_p_c),
          grad_phi_p(scratch.grad_phi_p),
          temperature_values (scratch.temperature_values),
          pressure_values (scratch.pressure_values),
          strain_rates (scratch.strain_rates),
          composition_values(scratch.composition_values),
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
                      const bool                add_compaction_pressure)
          :
          StokesPreconditioner<dim> (finite_element, quadrature,
                                     mapping,
                                     update_flags,
                                     n_compositional_fields,
                                     add_compaction_pressure),

          face_finite_element_values (mapping,
                                      finite_element,
                                      face_quadrature,
                                      face_update_flags),

          phi_u (finite_element.dofs_per_cell, numbers::signaling_nan<Tensor<1,dim> >()),
          grads_phi_u (finite_element.dofs_per_cell, numbers::signaling_nan<SymmetricTensor<2,dim> >()),
          div_phi_u (finite_element.dofs_per_cell, numbers::signaling_nan<double>()),
          velocity_values (quadrature.size(), numbers::signaling_nan<Tensor<1,dim> >()),
          face_material_model_inputs(face_quadrature.size(), n_compositional_fields),
          face_material_model_outputs(face_quadrature.size(), n_compositional_fields)
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
          face_material_model_outputs(scratch.face_material_model_outputs)
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element,
                         const FiniteElement<dim> &advection_element,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &quadrature,
                         const Quadrature<dim-1>  &face_quadrature,
                         const unsigned int        n_compositional_fields)
          :
          finite_element_values (mapping,
                                 finite_element, quadrature,
                                 update_values    |
                                 update_gradients |
                                 update_hessians  |
                                 update_quadrature_points |
                                 update_JxW_values),
          face_finite_element_values ((face_quadrature.size() > 0
                                       ?
                                       new FEFaceValues<dim> (mapping,
                                                              finite_element, face_quadrature,
                                                              update_values    |
                                                              update_gradients |
                                                              update_quadrature_points |
                                                              update_normal_vectors |
                                                              update_JxW_values)
                                       :
                                       NULL)),
          neighbor_face_finite_element_values ((face_quadrature.size() > 0
                                                ?
                                                new FEFaceValues<dim> (mapping,
                                                                       finite_element, face_quadrature,
                                                                       update_values    |
                                                                       update_gradients |
                                                                       update_quadrature_points |
                                                                       update_normal_vectors |
                                                                       update_JxW_values)
                                                :
                                                NULL)),
          subface_finite_element_values ((face_quadrature.size() > 0
                                          ?
                                          new FESubfaceValues<dim> (mapping,
                                                                    finite_element, face_quadrature,
                                                                    update_values    |
                                                                    update_gradients |
                                                                    update_quadrature_points |
                                                                    update_normal_vectors |
                                                                    update_JxW_values)
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
        StokesPreconditioner (const FiniteElement<dim> &finite_element)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_dof_indices (finite_element.dofs_per_cell)
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
        StokesSystem (const FiniteElement<dim> &finite_element,
                      const bool                do_pressure_rhs_compatibility_modification)
          :
          StokesPreconditioner<dim> (finite_element),
          local_rhs (finite_element.dofs_per_cell),
          local_pressure_shape_function_integrals (do_pressure_rhs_compatibility_modification ?
                                                   finite_element.dofs_per_cell
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
        needed_update_flags ()
      {}

    }
  }


  /**
   * Compute the variation in the entropy needed in the definition of the
   * artificial viscosity used to stabilize the composition/temperature equation.
   */
  template <int dim>
  double
  Simulator<dim>::get_entropy_variation (const double average_field,
                                         const AdvectionField &advection_field) const
  {
    // only do this if we really need entropy
    // variation. otherwise return something that's obviously
    // nonsensical
    if (parameters.stabilization_alpha != 2)
      return std::numeric_limits<double>::quiet_NaN();

    // record maximal entropy on Gauss quadrature
    // points
    const QGauss<dim> quadrature_formula (parameters.temperature_degree+1);
    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    FEValues<dim> fe_values (finite_element, quadrature_formula,
                             update_values | update_JxW_values);
    std::vector<double> old_field_values(n_q_points);
    std::vector<double> old_old_field_values(n_q_points);

    double min_entropy = std::numeric_limits<double>::max(),
           max_entropy = -std::numeric_limits<double>::max(),
           area = 0,
           entropy_integrated = 0;

    // loop over all locally owned cells and evaluate the entropy
    // at all quadrature points. keep a running tally of the
    // integral over the entropy as well as the area and the
    // maximal and minimal entropy
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[field].get_function_values (old_solution,
                                                old_field_values);
          fe_values[field].get_function_values (old_old_solution,
                                                old_old_field_values);
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double T = (old_field_values[q] +
                                old_old_field_values[q]) / 2;
              const double entropy = ((T-average_field) *
                                      (T-average_field));

              min_entropy = std::min (min_entropy, entropy);
              max_entropy = std::max (max_entropy, entropy);

              area += fe_values.JxW(q);
              entropy_integrated += fe_values.JxW(q) * entropy;
            }
        }

    // do MPI data exchange: we need to sum over
    // the two integrals (area,
    // entropy_integrated), and get the extrema
    // for maximum and minimum. combine
    // MPI_Allreduce for two values since that is
    // an expensive operation
    const double local_for_sum[2] = { entropy_integrated, area },
                                    local_for_max[2] = { -min_entropy, max_entropy };
    double global_for_sum[2], global_for_max[2];

    dealii::Utilities::MPI::sum (local_for_sum, mpi_communicator, global_for_sum);
    dealii::Utilities::MPI::max (local_for_max, mpi_communicator, global_for_max);

    const double average_entropy = global_for_sum[0] / global_for_sum[1];

    // return the maximal deviation of the entropy everywhere from the
    // average value
    return std::max(global_for_max[1] - average_entropy,
                    average_entropy - (-global_for_max[0]));
  }


  template <int dim>
  double
  Simulator<dim>::
  compute_viscosity (internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                     const double                        global_u_infty,
                     const double                        global_field_variation,
                     const double                        average_field,
                     const double                        global_entropy_variation,
                     const double                        cell_diameter,
                     const AdvectionField               &advection_field) const
  {
    // discontinuous Galerkin doesn't require an artificial viscosity
    if (advection_field.is_discontinuous(introspection))
      return 0.;

    std::vector<double> residual = assemblers->compute_advection_system_residual(*scratch.material_model_inputs.cell,
                                                                                 advection_field,
                                                                                 scratch);

    double max_residual = 0;
    double max_velocity = 0;
    double max_density = (advection_field.is_temperature()) ? 0.0 : 1.0;
    double max_specific_heat = (advection_field.is_temperature()) ? 0.0 : 1.0;
    double max_conductivity = 0;

    for (unsigned int q=0; q < scratch.finite_element_values.n_quadrature_points; ++q)
      {
        const Tensor<1,dim> velocity = (scratch.old_velocity_values[q] +
                                        scratch.old_old_velocity_values[q]) / 2;

        if (parameters.stabilization_alpha == 2)
          {
            const double field = (scratch.old_field_values[q] + scratch.old_old_field_values[q]) / 2;
            residual[q] *= std::abs(field - average_field);
          }

        max_residual = std::max (residual[q],     max_residual);
        max_velocity = std::max (std::sqrt (velocity*velocity), max_velocity);

        if (advection_field.is_temperature())
          {
            max_density = std::max       (scratch.material_model_outputs.densities[q],              max_density);
            max_specific_heat = std::max (scratch.material_model_outputs.specific_heat[q],          max_specific_heat);
            max_conductivity = std::max  (scratch.material_model_outputs.thermal_conductivities[q], max_conductivity);
          }
      }

    // If the velocity is 0 we have to assume a sensible velocity to calculate
    // an artificial diffusion. We choose similar to nondimensional
    // formulations: v ~ thermal_diffusivity / length_scale, which cancels
    // the density and specific heat from the entropy formulation. It seems
    // surprising at first that only the conductivity remains, but remember
    // that this actually *is* an additional artificial diffusion.
    if (std::abs(global_u_infty) < 1e-50)
      return parameters.stabilization_beta *
             max_conductivity / geometry_model->length_scale() *
             cell_diameter;

    const double max_viscosity = parameters.stabilization_beta *
                                 max_density *
                                 max_specific_heat *
                                 max_velocity * cell_diameter;

    if (timestep_number <= 1
        || std::abs(global_entropy_variation) < 1e-50
        || std::abs(global_field_variation) < 1e-50)
      // we don't have sensible timesteps during the first two iterations
      // and we can not divide by the entropy_variation if it is zero
      return max_viscosity;
    else
      {
        Assert (old_time_step > 0, ExcInternalError());

        double entropy_viscosity;
        if (parameters.stabilization_alpha == 2)
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * cell_diameter *
                               max_residual /
                               global_entropy_variation);
        else
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * global_Omega_diameter *
                               max_velocity * max_residual /
                               (global_u_infty * global_field_variation));


        return std::min (max_viscosity, entropy_viscosity);
      }
  }

  template <int dim>
  template <typename T>
  void
  Simulator<dim>::
  get_artificial_viscosity (Vector<T> &viscosity_per_cell,
                            const AdvectionField &advection_field) const
  {
    Assert(viscosity_per_cell.size()==triangulation.n_active_cells(), ExcInternalError());

    if (advection_field.field_type == AdvectionField::compositional_field)
      Assert(parameters.n_compositional_fields > advection_field.compositional_variable, ExcInternalError());

    viscosity_per_cell = 0.0;

    //discontinuous Galerkin doesn't require an artificial viscosity
    if (advection_field.is_discontinuous(introspection))
      return;

    const std::pair<double,double>
    global_field_range = get_extrapolated_advection_field_range (advection_field);
    double global_entropy_variation = get_entropy_variation ((global_field_range.first +
                                                              global_field_range.second) / 2,
                                                             advection_field);
    double global_max_velocity = get_maximal_velocity(old_solution);


    internal::Assembly::Scratch::
    AdvectionSystem<dim> scratch (finite_element,
                                  finite_element.base_element(advection_field.base_element(introspection)),
                                  *mapping,
                                  QGauss<dim>((advection_field.is_temperature()
                                               ?
                                               parameters.temperature_degree
                                               :
                                               parameters.composition_degree)
                                              +
                                              (parameters.stokes_velocity_degree+1)/2),
                                  /* Because we can only get here in the continuous case, which never requires
                                   * face integrals, we supply an invalid face_quadrature to the scratch object
                                   * to reduce the initialization cost.
                                   */
                                  Quadrature<dim-1> (),
                                  parameters.n_compositional_fields);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    for (; cell<dof_handler.end(); ++cell)
      {
        if (!cell->is_locally_owned()
            || (parameters.use_artificial_viscosity_smoothing  == true  &&  cell->is_artificial()))
          {
            viscosity_per_cell[cell->active_cell_index()]=-1;
            continue;
          }

        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

        // also have the number of dofs that correspond just to the element for
        // the system we are currently trying to assemble
        const unsigned int advection_dofs_per_cell = scratch.phi_field.size();
        (void)advection_dofs_per_cell;
        Assert (advection_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
        Assert (scratch.grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
        Assert (scratch.phi_field.size() == advection_dofs_per_cell, ExcInternalError());

        const FEValuesExtractors::Scalar solution_field
          = (advection_field.is_temperature()
             ?
             introspection.extractors.temperature
             :
             introspection.extractors.compositional_fields[advection_field.compositional_variable]);

        scratch.finite_element_values.reinit (cell);

        // get all dof indices on the current cell, then extract those
        // that correspond to the solution_field we are interested in
        cell->get_dof_indices (scratch.local_dof_indices);

        // initialize all of the scratch fields for further down
        scratch.finite_element_values[introspection.extractors.temperature].get_function_values (old_solution,
            scratch.old_temperature_values);
        scratch.finite_element_values[introspection.extractors.temperature].get_function_values (old_old_solution,
            scratch.old_old_temperature_values);

        scratch.finite_element_values[introspection.extractors.velocities].get_function_symmetric_gradients (old_solution,
            scratch.old_strain_rates);
        scratch.finite_element_values[introspection.extractors.velocities].get_function_symmetric_gradients (old_old_solution,
            scratch.old_old_strain_rates);

        scratch.finite_element_values[introspection.extractors.pressure].get_function_values (old_solution,
            scratch.old_pressure);
        scratch.finite_element_values[introspection.extractors.pressure].get_function_values (old_old_solution,
            scratch.old_old_pressure);

        for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
          {
            scratch.finite_element_values[introspection.extractors.compositional_fields[c]].get_function_values(old_solution,
                scratch.old_composition_values[c]);
            scratch.finite_element_values[introspection.extractors.compositional_fields[c]].get_function_values(old_old_solution,
                scratch.old_old_composition_values[c]);
          }

        scratch.finite_element_values[introspection.extractors.velocities].get_function_values (old_solution,
            scratch.old_velocity_values);
        scratch.finite_element_values[introspection.extractors.velocities].get_function_values (old_old_solution,
            scratch.old_old_velocity_values);
        scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
            scratch.current_velocity_values);

        scratch.finite_element_values[introspection.extractors.pressure].get_function_gradients (old_solution,
            scratch.old_pressure_gradients);
        scratch.finite_element_values[introspection.extractors.pressure].get_function_gradients (old_old_solution,
            scratch.old_old_pressure_gradients);


        scratch.old_field_values = (advection_field.is_temperature()
                                    ?
                                    scratch.old_temperature_values
                                    :
                                    scratch.old_composition_values[advection_field.compositional_variable]);
        scratch.old_old_field_values = (advection_field.is_temperature()
                                        ?
                                        scratch.old_old_temperature_values
                                        :
                                        scratch.old_old_composition_values[advection_field.compositional_variable]);

        scratch.finite_element_values[solution_field].get_function_gradients (old_solution,
                                                                              scratch.old_field_grads);
        scratch.finite_element_values[solution_field].get_function_gradients (old_old_solution,
                                                                              scratch.old_old_field_grads);

        scratch.finite_element_values[solution_field].get_function_laplacians (old_solution,
                                                                               scratch.old_field_laplacians);
        scratch.finite_element_values[solution_field].get_function_laplacians (old_old_solution,
                                                                               scratch.old_old_field_laplacians);

        if (parameters.include_melt_transport && melt_handler->is_porosity(advection_field))
          {
            scratch.finite_element_values[introspection.extractors.velocities].get_function_divergences (current_linearization_point,
                scratch.current_velocity_divergences);
          }

        /**
         * Explicit material model inputs and outputs.
         */
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            scratch.material_model_inputs.temperature[q] = (scratch.old_temperature_values[q] + scratch.old_old_temperature_values[q]) / 2;
            scratch.material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
            scratch.material_model_inputs.pressure[q] = (scratch.old_pressure[q] + scratch.old_old_pressure[q]) / 2;
            scratch.material_model_inputs.velocity[q] = (scratch.old_velocity_values[q] + scratch.old_old_velocity_values[q]) / 2;
            scratch.material_model_inputs.pressure_gradient[q] = (scratch.old_pressure_gradients[q] + scratch.old_old_pressure_gradients[q]) / 2;

            for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
              scratch.material_model_inputs.composition[q][c] = (scratch.old_composition_values[c][q] + scratch.old_old_composition_values[c][q]) / 2;
            scratch.material_model_inputs.strain_rate[q] = (scratch.old_strain_rates[q] + scratch.old_old_strain_rates[q]) / 2;
          }
        scratch.material_model_inputs.cell = &cell;
        create_additional_material_model_outputs(scratch.material_model_outputs);

        material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);
        MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                                   cell,
                                                   scratch.finite_element_values.get_quadrature(),
                                                   scratch.finite_element_values.get_mapping(),
                                                   scratch.material_model_outputs);

        viscosity_per_cell[cell->active_cell_index()] = compute_viscosity(scratch,
                                                                          global_max_velocity,
                                                                          global_field_range.second - global_field_range.first,
                                                                          0.5 * (global_field_range.second + global_field_range.first),
                                                                          global_entropy_variation,
                                                                          cell->diameter(),
                                                                          advection_field);
      }

    // if set to true, the maximum of the artificial viscosity in the cell as well
    // as the neighbors of the cell is computed and used instead
    if (parameters.use_artificial_viscosity_smoothing  == true)
      {
        Vector<T> viscosity_per_cell_temp;
        viscosity_per_cell_temp.reinit(triangulation.n_active_cells());

        viscosity_per_cell_temp = viscosity_per_cell;
        typename DoFHandler<dim>::active_cell_iterator
        cell,
        end_cell = dof_handler.end();
        for (cell = dof_handler.begin_active(); cell!=end_cell; ++cell)
          {
            if (cell->is_locally_owned())
              for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
                if (cell->at_boundary(face_no) == false)
                  {
                    if (cell->neighbor(face_no)->active())
                      viscosity_per_cell[cell->active_cell_index()] = std::max(viscosity_per_cell[cell->active_cell_index()],
                                                                               viscosity_per_cell_temp[cell->neighbor(face_no)->active_cell_index()]);
                    else
                      for (unsigned int l=0; l<cell->neighbor(face_no)->n_children(); l++)
                        if (cell->neighbor(face_no)->child(l)->active())
                          viscosity_per_cell[cell->active_cell_index()] = std::max(viscosity_per_cell[cell->active_cell_index()],
                                                                                   viscosity_per_cell_temp[cell->neighbor(face_no)->child(l)->active_cell_index()]);
                  }
          }
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
    std::vector<std::vector<double> > composition_values (parameters.n_compositional_fields,
                                                          std::vector<double> (n_q_points));

    for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
      input_finite_element_values[introspection.extractors.compositional_fields[c]].get_function_values(input_solution,
          composition_values[c]);

    // then we copy these values to exchange the inner and outer vector, because for the material
    // model we need a vector with values of all the compositional fields for every quadrature point
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
        material_model_inputs.composition[q][c] = composition_values[c][q];

    material_model_inputs.cell = &cell;
  }



  namespace Assemblers
  {

    /**
     * A class for the definition of functions that implement the
     * linear system terms for the *complete* compressible or
     * incompressible equations.
     */
    template <int dim>
    class CompleteEquations : public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
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


        void
        local_assemble_stokes_system_compressible (const double                                     pressure_scaling,
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
                  scratch.phi_p[k] = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
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

              const SymmetricTensor<4,dim> &stress_strain_director =
                scratch.material_model_outputs.stress_strain_directors[q];
              const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

              const Tensor<1,dim>
              gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

              const double compressibility
                = scratch.material_model_outputs.compressibilities[q];
              const double density = scratch.material_model_outputs.densities[q];

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
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * scratch.finite_element_values.JxW(q);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                data.local_rhs(i) += (
                                       (density * gravity * scratch.phi_u[i])
                                       +
                                       // add the term that results from the compressibility. compared
                                       // to the manual, this term seems to have the wrong sign, but this
                                       // is because we negate the entire equation to make sure we get
                                       // -div(u) as the adjoint operator of grad(p) (see above where
                                       // we assemble the matrix)
                                       (pressure_scaling *
                                        compressibility * density *
                                        (scratch.velocity_values[q] * gravity) *
                                        scratch.phi_p[i])
                                     )
                                     * scratch.finite_element_values.JxW(q);
            }

        }


        void
        local_assemble_stokes_system_incompressible (const double                                     pressure_scaling,
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
                  scratch.phi_p[k] = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
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

              const SymmetricTensor<4,dim> &stress_strain_director =
                scratch.material_model_outputs.stress_strain_directors[q];
              const bool use_tensor = (stress_strain_director !=  dealii::identity_tensor<dim> ());

              const Tensor<1,dim>
              gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

              const double density = scratch.material_model_outputs.densities[q];

              if (rebuild_stokes_matrix)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    data.local_matrix(i,j) += ( (use_tensor ?
                                                 eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                                 :
                                                 eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                                                - (pressure_scaling *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // finally the term -div(u). note the negative sign to make this
                                                // operator adjoint to the grad(p) term
                                                - (pressure_scaling *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * scratch.finite_element_values.JxW(q);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                data.local_rhs(i) += (density * gravity * scratch.phi_u[i])
                                     * scratch.finite_element_values.JxW(q);
            }
        }

        void
        local_assemble_advection_system (const typename Simulator<dim>::AdvectionField &advection_field,
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

          const unsigned int solution_component = advection_field.component_index(introspection);

          const FEValuesExtractors::Scalar solution_field
            = (advection_field.is_temperature()
               ?
               introspection.extractors.temperature
               :
               introspection.extractors.compositional_fields[advection_field.compositional_variable]
              );

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

              const double density_c_P              =
                ((advection_field.is_temperature())
                 ?
                 scratch.material_model_outputs.densities[q] *
                 scratch.material_model_outputs.specific_heat[q]
                 :
                 1.0);

              Assert (density_c_P >= 0,
                      ExcMessage ("The product of density and c_P needs to be a "
                                  "non-negative quantity."));

              const double conductivity =
                ((advection_field.is_temperature())
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
                ((advection_field.is_temperature())
                 ?
                 0.0
                 :
                 scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]);

              const double field_term_for_rhs
                = (use_bdf2_scheme ?
                   (scratch.old_field_values[q] *
                    (1 + time_step/old_time_step)
                    -
                    scratch.old_old_field_values[q] *
                    (time_step * time_step) /
                    (old_time_step * (time_step + old_time_step)))
                   :
                   scratch.old_field_values[q])
                  *
                  (density_c_P + latent_heat_LHS);

              Tensor<1,dim> current_u = scratch.current_velocity_values[q];
              //Subtract off the mesh velocity for ALE corrections if necessary
              if (this->get_parameters().free_surface_enabled)
                current_u -= scratch.mesh_velocity_values[q];

              const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                        (time_step + old_time_step)) : 1.0;

              // do the actual assembly. note that we only need to loop over the advection
              // shape functions because these are the only contributions we compute here
              for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                {
                  data.local_rhs(i)
                  += (field_term_for_rhs * scratch.phi_field[i]
                      + time_step *
                      scratch.phi_field[i]
                      * gamma
                      + scratch.phi_field[i]
                      * reaction_term)
                     *
                     scratch.finite_element_values.JxW(q);

                  for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                    {
                      data.local_matrix(i,j)
                      += (
                           (this->get_timestep() * (conductivity + artificial_viscosity)
                            * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                           + ((time_step * (scratch.phi_field[i] * (current_u * scratch.grad_phi_field[j])))
                              + (factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                           (density_c_P + latent_heat_LHS)
                         )
                         * scratch.finite_element_values.JxW(q);
                    }
                }
            }
        }


        std::vector<double>
        compute_advection_system_residual(const typename Simulator<dim>::AdvectionField     &advection_field,
                                          internal::Assembly::Scratch::AdvectionSystem<dim> &scratch) const
        {
          const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
          std::vector<double> residuals(n_q_points);

          this->get_heating_model_manager().evaluate(scratch.material_model_inputs,
                                                     scratch.material_model_outputs,
                                                     scratch.heating_model_outputs);

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

              const double density       = ((advection_field.is_temperature()) ? scratch.material_model_outputs.densities[q] : 1.0);
              const double conductivity  = ((advection_field.is_temperature()) ? scratch.material_model_outputs.thermal_conductivities[q] : 0.0);
              const double c_P           = ((advection_field.is_temperature()) ? scratch.material_model_outputs.specific_heat[q] : 1.0);
              const double k_Delta_field = conductivity
                                           * (scratch.old_field_laplacians[q] +
                                              scratch.old_old_field_laplacians[q]) / 2;

              const double gamma =
                ((advection_field.is_temperature())
                 ?
                 scratch.heating_model_outputs.heating_source_terms[q]
                 :
                 0.0);

              const double latent_heat_LHS =
                ((advection_field.is_temperature())
                 ?
                 scratch.heating_model_outputs.lhs_latent_heat_terms[q]
                 :
                 0.0);

              const double dreaction_term_dt =
                (advection_field.is_temperature() || this->get_old_timestep() == 0)
                ?
                0.0
                :
                (scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]
                 / this->get_old_timestep());

              residuals[q]
                = std::abs((density * c_P + latent_heat_LHS) * (dField_dt + u_grad_field) - k_Delta_field - gamma
                           - dreaction_term_dt);
            }
          return residuals;
        }


        void
        local_assemble_discontinuous_advection_boundary_face_terms(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                   const unsigned int face_no,
                                                                   const typename Simulator<dim>::AdvectionField &advection_field,
                                                                   internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                                                   internal::Assembly::CopyData::AdvectionSystem<dim> &data) const
        {
          const Parameters<dim> &parameters = this->get_parameters();
          const Introspection<dim> &introspection = this->introspection();
          const unsigned int n_q_points    = scratch.face_finite_element_values->n_quadrature_points;

          const double time_step = this->get_timestep();

          if (!advection_field.is_discontinuous(introspection))
            return;

          // also have the number of dofs that correspond just to the element for
          // the system we are currently trying to assemble
          const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

          Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
          Assert (scratch.face_grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
          Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());

          const unsigned int solution_component = advection_field.component_index(introspection);

          const FEValuesExtractors::Scalar solution_field
            = (advection_field.is_temperature()
               ?
               introspection.extractors.temperature
               :
               introspection.extractors.compositional_fields[advection_field.compositional_variable]
              );

          typename DoFHandler<dim>::face_iterator face = cell->face (face_no);

          if (((parameters.fixed_temperature_boundary_indicators.find(
                  face->boundary_id()
                )
                != parameters.fixed_temperature_boundary_indicators.end())
               && (advection_field.is_temperature()))
              ||
              (( parameters.fixed_composition_boundary_indicators.find(
                   face->boundary_id()
                 )
                 != parameters.fixed_composition_boundary_indicators.end())
               && (!advection_field.is_temperature())))
            {
              /*
               * We are in the case of a Dirichlet temperature or composition boundary.
               * In the temperature case, impose the Dirichlet value weakly using a matrix term
               * and RHS term. In the composition case, Dirichlet conditions can only be imposed
               * on inflow boundaries, and we only have the flow-dependent terms, so we only
               * assemble the corresponding flow-dependent and matrix and RHS terms
               * if we are on an inflow boundary.
               */

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  // precompute the values of shape functions and their gradients.
                  // We only need to look up values of shape functions if they
                  // belong to 'our' component. They are zero otherwise anyway.
                  // Note that we later only look at the values that we do set here.
                  for (unsigned int k=0; k<advection_dofs_per_cell; ++k)
                    {
                      scratch.face_grad_phi_field[k] = (*scratch.face_finite_element_values)[solution_field].gradient (scratch.face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                      scratch.face_phi_field[k]      = (*scratch.face_finite_element_values)[solution_field].value (scratch.face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                    }
                  const double density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.densities[q] *
                     scratch.face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  Assert (density_c_P >= 0,
                          ExcMessage ("The product of density and c_P needs to be a "
                                      "non-negative quantity."));

                  const double conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  Assert (density_c_P + latent_heat_LHS >= 0,
                          ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                      "to the left hand side needs to be a non-negative quantity."));

                  const double penalty = (advection_field.is_temperature()
                                          ?
                                          parameters.discontinuous_penalty
                                          * parameters.temperature_degree
                                          * parameters.temperature_degree
                                          / face->measure()
                                          * conductivity
                                          / (density_c_P + latent_heat_LHS)
                                          :
                                          0.0);

                  const double dirichlet_value = (advection_field.is_temperature()
                                                  ?
                                                  this->get_boundary_temperature().boundary_temperature(
                                                    cell->face(face_no)->boundary_id(),
                                                    scratch.face_finite_element_values->quadrature_point(q))
                                                  :
                                                  this->get_boundary_composition().boundary_composition(
                                                    cell->face(face_no)->boundary_id(),
                                                    scratch.face_finite_element_values->quadrature_point(q),
                                                    advection_field.compositional_variable));

                  Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                  //Subtract off the mesh velocity for ALE corrections if necessary
                  if (parameters.free_surface_enabled)
                    current_u -= scratch.face_mesh_velocity_values[q];

                  /**
                   * The discontinuous Galerkin method uses 2 types of jumps over edges:
                   * undirected and directed jumps. Undirected jumps are dependent only
                   * on the order of the numbering of cells. Directed jumps are dependent
                   * on the direction of the flow. Thus the flow-dependent terms below are
                   * only calculated if the edge is an inflow edge.
                   */
                  const bool inflow = ((current_u * scratch.face_finite_element_values->normal_vector(q)) < 0.);

                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    {
                      data.local_rhs(i)
                      += (- time_step *  conductivity
                          * scratch.face_grad_phi_field[i]
                          * scratch.face_finite_element_values->normal_vector(q)
                          * dirichlet_value

                          + time_step
                          * (density_c_P + latent_heat_LHS)
                          * penalty
                          * scratch.face_phi_field[i]
                          * dirichlet_value

                          + (inflow
                             ?
                             - (density_c_P + latent_heat_LHS)
                             * time_step
                             * (current_u
                                * scratch.face_finite_element_values->normal_vector(q))
                             * dirichlet_value
                             * scratch.face_phi_field[i]
                             :
                             0.)
                         )
                         *
                         scratch.face_finite_element_values->JxW(q);

                      // local_matrix terms
                      for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j)
                          += (- time_step *  conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - time_step *  conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + time_step
                              * (density_c_P + latent_heat_LHS)
                              * penalty
                              * scratch.face_phi_field[i]
                              * scratch.face_phi_field[j]

                              + (inflow
                                 ?
                                 - (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);
                        }
                    }
                }
            }
          else
            {
              //Neumann temperature term - no non-zero contribution as only homogeneous Neumann boundary conditions are implemented elsewhere for temperature
            }
        }

        void
        local_assemble_discontinuous_advection_interior_face_terms(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                   const unsigned int face_no,
                                                                   const typename Simulator<dim>::AdvectionField &advection_field,
                                                                   internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                                                   internal::Assembly::CopyData::AdvectionSystem<dim> &data) const
        {
          const Parameters<dim> &parameters = this->get_parameters();
          const Introspection<dim> &introspection = this->introspection();
          const unsigned int n_q_points    = scratch.face_finite_element_values->n_quadrature_points;

          const double time_step = this->get_timestep();

          if (!advection_field.is_discontinuous(introspection))
            return;

          // also have the number of dofs that correspond just to the element for
          // the system we are currently trying to assemble
          const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

          Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
          Assert (scratch.face_grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
          Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());

          const unsigned int solution_component = advection_field.component_index(introspection);

          const FEValuesExtractors::Scalar solution_field
            = (advection_field.is_temperature()
               ?
               introspection.extractors.temperature
               :
               introspection.extractors.compositional_fields[advection_field.compositional_variable]
              );

          typename DoFHandler<dim>::face_iterator face = cell->face (face_no);

          //interior face - no contribution on RHS
          Assert (cell->neighbor(face_no).state() == IteratorState::valid,
                  ExcInternalError());
          const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no); //note: NOT active_cell_iterator, so this includes cells that are refined.

          if (!(face->has_children()))
            {
              if (!cell->neighbor_is_coarser(face_no) &&
                  (((neighbor->is_locally_owned()) && (cell->index() < neighbor->index()))
                   ||
                   ((!neighbor->is_locally_owned()) && (cell->subdomain_id() < neighbor->subdomain_id()))))
                {
                  Assert (cell->is_locally_owned(), ExcInternalError());
                  //cell and neighbor are equal-sized, and cell has been chosen to assemble this face, so calculate from cell

                  //how does the neighbor talk about this cell?
                  const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

                  //set up neighbor values
                  scratch.neighbor_face_finite_element_values->reinit (neighbor, neighbor2);

                  this->compute_material_model_input_values (this->get_current_linearization_point(),
                                                             *scratch.neighbor_face_finite_element_values,
                                                             neighbor,
                                                             true,
                                                             scratch.neighbor_face_material_model_inputs);
                  this->get_material_model().evaluate(scratch.neighbor_face_material_model_inputs,
                                                      scratch.neighbor_face_material_model_outputs);

                  this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                             scratch.neighbor_face_material_model_outputs,
                                                             scratch.neighbor_face_heating_model_outputs);

                  std::vector<types::global_dof_index> neighbor_dof_indices (scratch.face_finite_element_values->get_fe().dofs_per_cell);
                  // get all dof indices on the neighbor, then extract those
                  // that correspond to the solution_field we are interested in
                  neighbor->get_dof_indices (neighbor_dof_indices);
                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face][i] = neighbor_dof_indices[scratch.face_finite_element_values->get_fe().component_to_system_index(solution_component, i)];
                  data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face] = true;

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      // precompute the values of shape functions and their gradients.
                      // We only need to look up values of shape functions if they
                      // belong to 'our' component. They are zero otherwise anyway.
                      // Note that we later only look at the values that we do set here.
                      for (unsigned int k=0; k<advection_dofs_per_cell; ++k)
                        {
                          scratch.face_grad_phi_field[k]          = (*scratch.face_finite_element_values)[solution_field].gradient (scratch.face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                          scratch.face_phi_field[k]               = (*scratch.face_finite_element_values)[solution_field].value (scratch.face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                          scratch.neighbor_face_grad_phi_field[k] = (*scratch.neighbor_face_finite_element_values)[solution_field].gradient (scratch.neighbor_face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                          scratch.neighbor_face_phi_field[k]      = (*scratch.neighbor_face_finite_element_values)[solution_field].value (scratch.neighbor_face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                        }

                      const double density_c_P              =
                        ((advection_field.is_temperature())
                         ?
                         scratch.face_material_model_outputs.densities[q] *
                         scratch.face_material_model_outputs.specific_heat[q]
                         :
                         1.0);

                      Assert (density_c_P >= 0,
                              ExcMessage ("The product of density and c_P needs to be a "
                                          "non-negative quantity."));

                      const double conductivity =
                        ((advection_field.is_temperature())
                         ?
                         scratch.face_material_model_outputs.thermal_conductivities[q]
                         :
                         0.0);
                      const double latent_heat_LHS =
                        ((advection_field.is_temperature())
                         ?
                         scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                         :
                         0.0);
                      Assert (density_c_P + latent_heat_LHS >= 0,
                              ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                          "to the left hand side needs to be a non-negative quantity."));

                      const double penalty = (advection_field.is_temperature()
                                              ?
                                              parameters.discontinuous_penalty
                                              * parameters.temperature_degree
                                              * parameters.temperature_degree
                                              / face->measure()
                                              * conductivity
                                              / (density_c_P + latent_heat_LHS)
                                              :
                                              0.0);

                      Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                      //Subtract off the mesh velocity for ALE corrections if necessary
                      if (parameters.free_surface_enabled)
                        current_u -= scratch.face_mesh_velocity_values[q];

                      const double neighbor_density_c_P              =
                        ((advection_field.is_temperature())
                         ?
                         scratch.neighbor_face_material_model_outputs.densities[q] *
                         scratch.neighbor_face_material_model_outputs.specific_heat[q]
                         :
                         1.0);

                      Assert (neighbor_density_c_P >= 0,
                              ExcMessage ("The product of density and c_P on the neighbor needs to be a "
                                          "non-negative quantity."));

                      const double neighbor_conductivity =
                        ((advection_field.is_temperature())
                         ?
                         scratch.neighbor_face_material_model_outputs.thermal_conductivities[q]
                         :
                         0.0);
                      const double neighbor_latent_heat_LHS =
                        ((advection_field.is_temperature())
                         ?
                         scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q]
                         :
                         0.0);
                      Assert (neighbor_density_c_P + neighbor_latent_heat_LHS >= 0,
                              ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                          "to the left hand side on the neighbor needs to be a non-negative quantity."));

                      const double neighbor_penalty = (advection_field.is_temperature()
                                                       ?
                                                       parameters.discontinuous_penalty
                                                       * parameters.temperature_degree
                                                       * parameters.temperature_degree
                                                       / neighbor->face(neighbor2)->measure()
                                                       * neighbor_conductivity
                                                       / (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                                       :
                                                       0.0);

                      const double max_penalty = std::max(penalty, neighbor_penalty);

                      const double max_density_c_P_and_latent_heat =
                        std::max(density_c_P + latent_heat_LHS,
                                 neighbor_density_c_P + neighbor_latent_heat_LHS);

                      Assert (numbers::is_finite(max_density_c_P_and_latent_heat),
                              ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a finite quantity."));
                      Assert (max_density_c_P_and_latent_heat >= 0,
                              ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a "
                                          "non-negative quantity."));

                      /**
                       * The discontinuous Galerkin method uses 2 types of jumps over edges:
                       * undirected and directed jumps. Undirected jumps are dependent only
                       * on the order of the numbering of cells. Directed jumps are dependent
                       * on the direction of the flow. Thus the flow-dependent terms below are
                       * only calculated if the edge is an inflow edge.
                       */
                      const bool inflow = ((current_u * scratch.face_finite_element_values->normal_vector(q)) < 0.);

                      for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                        {
                          for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                            {
                              data.local_matrix(i,j)
                              += (- 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[i]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[j]

                                  - 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[j]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[i]

                                  + time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.face_phi_field[i]
                                  * scratch.face_phi_field[j]

                                  - (inflow
                                     ?
                                     (density_c_P + latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.face_finite_element_values->normal_vector(q))
                                     * scratch.face_phi_field[i]
                                     * scratch.face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.face_finite_element_values->JxW(q);

                              data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face](i,j)
                              += (- 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[j]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[i]

                                  + 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[i]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[j]

                                  - time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.neighbor_face_phi_field[j]
                                  * scratch.face_phi_field[i]

                                  + (inflow
                                     ?
                                     (density_c_P + latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.face_finite_element_values->normal_vector(q))
                                     * scratch.face_phi_field[i]
                                     * scratch.neighbor_face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.face_finite_element_values->JxW(q);

                              data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face](i,j)
                              += (+ 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[j]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[i]

                                  - 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[i]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[j]

                                  - time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.face_phi_field[j]
                                  * scratch.neighbor_face_phi_field[i]

                                  - (!inflow
                                     ?
                                     (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.face_finite_element_values->normal_vector(q))
                                     * scratch.neighbor_face_phi_field[i]
                                     * scratch.face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.face_finite_element_values->JxW(q);

                              data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face](i,j)
                              += (+ 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[i]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[j]

                                  + 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[j]
                                  * scratch.face_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[i]

                                  + time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.neighbor_face_phi_field[i]
                                  * scratch.neighbor_face_phi_field[j]

                                  + (!inflow
                                     ?
                                     (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.face_finite_element_values->normal_vector(q))
                                     * scratch.neighbor_face_phi_field[i]
                                     * scratch.neighbor_face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.face_finite_element_values->JxW(q);
                            }
                        }
                    }
                }
              else
                {
                  /* neighbor is taking responsibility for assembly of this face, because
                   * either (1) neighbor is coarser, or
                   *        (2) neighbor is equally-sized and
                   *           (a) neighbor is on a different subdomain, with lower subdmain_id(), or
                   *           (b) neighbor is on the same subdomain and has lower index().
                  */
                }
            }
          else //face->has_children(), so always assemble from here.
            {
              //how does the neighbor talk about this cell?
              const unsigned int neighbor2 = cell->neighbor_face_no(face_no);

              //loop over subfaces
              for (unsigned int subface_no=0; subface_no<face->number_of_children(); ++subface_no)
                {
                  const typename DoFHandler<dim>::active_cell_iterator neighbor_child
                    = cell->neighbor_child_on_subface (face_no, subface_no);

                  //set up subface values
                  scratch.subface_finite_element_values->reinit (cell, face_no, subface_no);

                  //subface->face
                  (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_current_linearization_point(),
                      scratch.face_current_velocity_values);

                  //get the mesh velocity, as we need to subtract it off of the advection systems
                  if (parameters.free_surface_enabled)
                    (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_mesh_velocity(),
                        scratch.face_mesh_velocity_values);

                  //get the mesh velocity, as we need to subtract it off of the advection systems
                  if (parameters.free_surface_enabled)
                    (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_mesh_velocity(),
                        scratch.face_mesh_velocity_values);

                  this->compute_material_model_input_values (this->get_current_linearization_point(),
                                                             *scratch.subface_finite_element_values,
                                                             cell,
                                                             true,
                                                             scratch.face_material_model_inputs);
                  this->get_material_model().evaluate(scratch.face_material_model_inputs,
                                                      scratch.face_material_model_outputs);

                  this->get_heating_model_manager().evaluate(scratch.face_material_model_inputs,
                                                             scratch.face_material_model_outputs,
                                                             scratch.face_heating_model_outputs);

                  //set up neighbor values
                  scratch.neighbor_face_finite_element_values->reinit (neighbor_child, neighbor2);

                  this->compute_material_model_input_values (this->get_current_linearization_point(),
                                                             *scratch.neighbor_face_finite_element_values,
                                                             neighbor_child,
                                                             true,
                                                             scratch.neighbor_face_material_model_inputs);
                  this->get_material_model().evaluate(scratch.neighbor_face_material_model_inputs,
                                                      scratch.neighbor_face_material_model_outputs);

                  this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                             scratch.neighbor_face_material_model_outputs,
                                                             scratch.neighbor_face_heating_model_outputs);

                  std::vector<types::global_dof_index> neighbor_dof_indices (scratch.face_finite_element_values->get_fe().dofs_per_cell);
                  // get all dof indices on the neighbor, then extract those
                  // that correspond to the solution_field we are interested in
                  neighbor_child->get_dof_indices (neighbor_dof_indices);
                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no][i] = neighbor_dof_indices[scratch.subface_finite_element_values->get_fe().component_to_system_index(solution_component, i)];
                  data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no] = true;

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      // precompute the values of shape functions and their gradients.
                      // We only need to look up values of shape functions if they
                      // belong to 'our' component. They are zero otherwise anyway.
                      // Note that we later only look at the values that we do set here.
                      for (unsigned int k=0; k<advection_dofs_per_cell; ++k)
                        {
                          scratch.face_grad_phi_field[k]          = (*scratch.subface_finite_element_values)[solution_field].gradient (scratch.subface_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                          scratch.face_phi_field[k]               = (*scratch.subface_finite_element_values)[solution_field].value (scratch.subface_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                          scratch.neighbor_face_grad_phi_field[k] = (*scratch.neighbor_face_finite_element_values)[solution_field].gradient (scratch.neighbor_face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                          scratch.neighbor_face_phi_field[k]      = (*scratch.neighbor_face_finite_element_values)[solution_field].value (scratch.neighbor_face_finite_element_values->get_fe().component_to_system_index(solution_component, k), q);
                        }

                      const double density_c_P              =
                        ((advection_field.is_temperature())
                         ?
                         scratch.face_material_model_outputs.densities[q] *
                         scratch.face_material_model_outputs.specific_heat[q]
                         :
                         1.0);

                      Assert (density_c_P >= 0,
                              ExcMessage ("The product of density and c_P needs to be a "
                                          "non-negative quantity."));

                      const double conductivity =
                        ((advection_field.is_temperature())
                         ?
                         scratch.face_material_model_outputs.thermal_conductivities[q]
                         :
                         0.0);
                      const double latent_heat_LHS =
                        ((advection_field.is_temperature())
                         ?
                         scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                         :
                         0.0);
                      Assert (density_c_P + latent_heat_LHS >= 0,
                              ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                          "to the left hand side needs to be a non-negative quantity."));

                      const double penalty = (advection_field.is_temperature()
                                              ?
                                              parameters.discontinuous_penalty
                                              * parameters.temperature_degree
                                              * parameters.temperature_degree
                                              / face->measure()
                                              * conductivity
                                              / (density_c_P + latent_heat_LHS)
                                              :
                                              0.0);

                      Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                      //Subtract off the mesh velocity for ALE corrections if necessary
                      if (parameters.free_surface_enabled)
                        current_u -= scratch.face_mesh_velocity_values[q];

                      const double neighbor_density_c_P              =
                        ((advection_field.is_temperature())
                         ?
                         scratch.neighbor_face_material_model_outputs.densities[q] *
                         scratch.neighbor_face_material_model_outputs.specific_heat[q]
                         :
                         1.0);

                      Assert (neighbor_density_c_P >= 0,
                              ExcMessage ("The product of density and c_P on the neighbor needs to be a "
                                          "non-negative quantity."));

                      const double neighbor_conductivity =
                        ((advection_field.is_temperature())
                         ?
                         scratch.neighbor_face_material_model_outputs.thermal_conductivities[q]
                         :
                         0.0);
                      const double neighbor_latent_heat_LHS =
                        ((advection_field.is_temperature())
                         ?
                         scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q]
                         :
                         0.0);
                      Assert (neighbor_density_c_P + neighbor_latent_heat_LHS >= 0,
                              ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                          "to the left hand side on the neighbor needs to be a non-negative quantity."));

                      const double neighbor_penalty = (advection_field.is_temperature()
                                                       ?
                                                       parameters.discontinuous_penalty
                                                       * parameters.temperature_degree
                                                       * parameters.temperature_degree
                                                       / neighbor_child->face(neighbor2)->measure()
                                                       * neighbor_conductivity
                                                       / (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                                       :
                                                       0.0);

                      const double max_penalty = std::max(penalty, neighbor_penalty);

                      const double max_density_c_P_and_latent_heat =
                        std::max(density_c_P + latent_heat_LHS,
                                 neighbor_density_c_P + neighbor_latent_heat_LHS);

                      Assert (numbers::is_finite(max_density_c_P_and_latent_heat),
                              ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a finite quantity."));
                      Assert (max_density_c_P_and_latent_heat >= 0,
                              ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a "
                                          "non-negative quantity."));

                      /**
                       * The discontinuous Galerkin method uses 2 types of jumps over edges:
                       * undirected and directed jumps. Undirected jumps are dependent only
                       * on the order of the numbering of cells. Directed jumps are dependent
                       * on the direction of the flow. Thus the flow-dependent terms below are
                       * only calculated if the edge is an inflow edge.
                       */
                      const bool inflow = ((current_u * scratch.subface_finite_element_values->normal_vector(q)) < 0.);

                      for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                        {
                          for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                            {
                              data.local_matrix(i,j)
                              += (- 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[i]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[j]

                                  - 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[j]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[i]

                                  + time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.face_phi_field[i]
                                  * scratch.face_phi_field[j]

                                  - (inflow
                                     ?
                                     (density_c_P + latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.subface_finite_element_values->normal_vector(q))
                                     * scratch.face_phi_field[i]
                                     * scratch.face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.subface_finite_element_values->JxW(q);

                              data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                              += (- 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[j]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[i]

                                  + 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[i]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[j]

                                  - time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.neighbor_face_phi_field[j]
                                  * scratch.face_phi_field[i]

                                  + (inflow
                                     ?
                                     (density_c_P + latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.subface_finite_element_values->normal_vector(q))
                                     * scratch.face_phi_field[i]
                                     * scratch.neighbor_face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.subface_finite_element_values->JxW(q);

                              data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                              += (+ 0.5 * time_step * conductivity
                                  * scratch.face_grad_phi_field[j]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[i]

                                  - 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[i]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.face_phi_field[j]

                                  - time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.face_phi_field[j]
                                  * scratch.neighbor_face_phi_field[i]

                                  - (!inflow
                                     ?
                                     (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.subface_finite_element_values->normal_vector(q))
                                     * scratch.neighbor_face_phi_field[i]
                                     * scratch.face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.subface_finite_element_values->JxW(q);

                              data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                              += (+ 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[i]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[j]

                                  + 0.5 * time_step * neighbor_conductivity
                                  * scratch.neighbor_face_grad_phi_field[j]
                                  * scratch.subface_finite_element_values->normal_vector(q)
                                  * scratch.neighbor_face_phi_field[i]

                                  + time_step
                                  * max_density_c_P_and_latent_heat
                                  * max_penalty
                                  * scratch.neighbor_face_phi_field[i]
                                  * scratch.neighbor_face_phi_field[j]

                                  + (!inflow
                                     ?
                                     (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                     * time_step
                                     * (current_u
                                        * scratch.subface_finite_element_values->normal_vector(q))
                                     * scratch.neighbor_face_phi_field[i]
                                     * scratch.neighbor_face_phi_field[j]
                                     :
                                     0.)
                                 )
                                 * scratch.subface_finite_element_values->JxW(q);
                            }
                        }
                    }
                }
            }
        }
    };


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

        const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              scratch.phi_p[i] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
              data.local_pressure_shape_function_integrals(i) += scratch.phi_p[i] * scratch.finite_element_values.JxW(q);
            }
      }



      template <int dim>
      void
      traction_boundary_conditions (const SimulatorAccess<dim>                           &simulator_access,
                                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    const unsigned int                                    face_no,
                                    internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                    internal::Assembly::CopyData::StokesSystem<dim>      &data)
      {
        const Introspection<dim> &introspection = simulator_access.introspection();


        // see if any of the faces are traction boundaries for which
        // we need to assemble force terms for the right hand side
        const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
        if (simulator_access.get_traction_boundary_conditions()
            .find (cell->face(face_no)->boundary_id())
            !=
            simulator_access.get_traction_boundary_conditions().end())
          {
            scratch.face_finite_element_values.reinit (cell, face_no);

            for (unsigned int q=0; q<scratch.face_finite_element_values.n_quadrature_points; ++q)
              {
                const Tensor<1,dim> traction
                  = simulator_access.get_traction_boundary_conditions().find(
                      cell->face(face_no)->boundary_id()
                    )->second
                    ->boundary_traction (cell->face(face_no)->boundary_id(),
                                         scratch.face_finite_element_values.quadrature_point(q),
                                         scratch.face_finite_element_values.normal_vector(q));
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  data.local_rhs(i) += scratch.face_finite_element_values[introspection.extractors.velocities].value(i,q) *
                                       traction *
                                       scratch.face_finite_element_values.JxW(q);
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
    aspect::Assemblers::CompleteEquations<dim> *complete_equation_assembler
      = new aspect::Assemblers::CompleteEquations<dim>();

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
      .connect (std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::local_assemble_stokes_preconditioner,
                                std_cxx11::cref (*complete_equation_assembler),
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
    else if (material_model->is_compressible())
      assemblers->local_assemble_stokes_system
      .connect (std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::local_assemble_stokes_system_compressible,
                                std_cxx11::cref (*complete_equation_assembler),
                                // discard cell,
                                std_cxx11::_2,
                                std_cxx11::_3,
                                std_cxx11::_4,
                                std_cxx11::_5));
    else
      assemblers->local_assemble_stokes_system
      .connect (std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::local_assemble_stokes_system_incompressible,
                                std_cxx11::cref (*complete_equation_assembler),
                                // discard cell,
                                std_cxx11::_2,
                                std_cxx11::_3,
                                std_cxx11::_4,
                                std_cxx11::_5));

    assembler_objects.push_back (std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >
                                 (complete_equation_assembler));

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
    .connect (std_cxx11::bind(&aspect::Assemblers::OtherTerms::traction_boundary_conditions<dim>,
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
        .connect (std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::local_assemble_advection_system,
                                  std_cxx11::cref (*complete_equation_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3,
                                  std_cxx11::_4,
                                  std_cxx11::_5));

        assemblers->compute_advection_system_residual
        .connect (std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::compute_advection_system_residual,
                                  std_cxx11::cref (*complete_equation_assembler),
                                  // discard cell,
                                  std_cxx11::_2,
                                  std_cxx11::_3));
      }

    if (parameters.use_discontinuous_temperature_discretization ||
        parameters.use_discontinuous_composition_discretization)
      {
        assemblers->local_assemble_advection_system_on_interior_face
        .connect(std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::local_assemble_discontinuous_advection_interior_face_terms,
                                 std_cxx11::cref (*complete_equation_assembler),
                                 std_cxx11::_1,
                                 std_cxx11::_2,
                                 std_cxx11::_3,
                                 std_cxx11::_4,
                                 std_cxx11::_5));

        assemblers->local_assemble_advection_system_on_boundary_face
        .connect(std_cxx11::bind(&aspect::Assemblers::CompleteEquations<dim>::local_assemble_discontinuous_advection_boundary_face_terms,
                                 std_cxx11::cref (*complete_equation_assembler),
                                 std_cxx11::_1,
                                 std_cxx11::_2,
                                 std_cxx11::_3,
                                 std_cxx11::_4,
                                 std_cxx11::_5));

        assemblers->advection_system_assembler_on_face_properties.need_face_material_model_data = true;
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

    cell->get_dof_indices (data.local_dof_indices);
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
                                    parameters.n_compositional_fields,
                                    parameters.include_melt_transport),
         internal::Assembly::CopyData::
         StokesPreconditioner<dim> (finite_element));

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
              compute_material_model_input_values (current_linearization_point,
                                                   scratch.face_finite_element_values,
                                                   cell,
                                                   rebuild_stokes_matrix,
                                                   scratch.face_material_model_inputs);
              create_additional_material_model_outputs(scratch.face_material_model_outputs);

              material_model->evaluate(scratch.face_material_model_inputs,
                                       scratch.face_material_model_outputs);
//TODO: the following doesn't currently compile because the get_quadrature() call returns
//  a dim-1 dimensional quadrature
              // MaterialModel::MaterialAveraging::average (parameters.material_averaging,
              //                                            cell,
              //                                            scratch.face_finite_element_values.get_quadrature(),
              //                                            scratch.face_finite_element_values.get_mapping(),
              //                                            scratch.face_material_model_outputs);
            }

          assemblers->local_assemble_stokes_system_on_boundary_face(cell, face_no,
                                                                    pressure_scaling, rebuild_stokes_matrix,
                                                                    scratch, data);
        }

    cell->get_dof_indices (data.local_dof_indices);
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
          UpdateFlags(0))
        |
        (assemblers->stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data
         ?
         // if we need a material model input on the faces, we need to
         // also be able to compute the strain rate
         update_gradients
         :
         UpdateFlags(0))
        |
        assemblers->stokes_system_assembler_on_boundary_face_properties.needed_update_flags;

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
                            parameters.n_compositional_fields,
                            parameters.include_melt_transport),
         internal::Assembly::CopyData::
         StokesSystem<dim> (finite_element,
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

    const FEValuesExtractors::Scalar solution_field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    const unsigned int solution_component = advection_field.component_index(introspection);

    scratch.finite_element_values.reinit (cell);

    // get all dof indices on the current cell, then extract those
    // that correspond to the solution_field we are interested in
    cell->get_dof_indices (scratch.local_dof_indices);
    for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
      data.local_dof_indices[i] = scratch.local_dof_indices[scratch.finite_element_values.get_fe().component_to_system_index(solution_component, i)];

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
    const bool has_boundary_face_assemblers = !assemblers->local_assemble_advection_system_on_boundary_face.empty();
    const bool has_interior_face_assemblers = !assemblers->local_assemble_advection_system_on_interior_face.empty();

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

            if (assemblers->advection_system_assembler_on_face_properties.need_face_material_model_data)
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
  copy_local_to_global_advection_system (const AdvectionField &/*advection_field*/,
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
    if (!assemblers->local_assemble_advection_system_on_interior_face.empty())
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
    const unsigned int advection_quadrature_degree = (advection_field.is_temperature()
                                                      ?
                                                      parameters.temperature_degree
                                                      :
                                                      parameters.composition_degree)
                                                     +
                                                     (parameters.stokes_velocity_degree+1)/2;

    const bool allocate_face_quadrature = !assemblers->local_assemble_advection_system_on_boundary_face.empty() ||
                                          !assemblers->local_assemble_advection_system_on_interior_face.empty();
    const bool allocate_neighbor_contributions = !assemblers->local_assemble_advection_system_on_interior_face.empty();

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
                               parameters.n_compositional_fields),
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
  template void Simulator<dim>::get_artificial_viscosity (Vector<double> &viscosity_per_cell,  \
                                                          const AdvectionField &advection_field) const; \
  template void Simulator<dim>::get_artificial_viscosity (Vector<float> &viscosity_per_cell,  \
                                                          const AdvectionField &advection_field) const; \
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

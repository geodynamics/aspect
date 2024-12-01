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
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/melt.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
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
      return numbers::signaling_nan<double>();

    // record maximal entropy on Gauss quadrature points
    const Quadrature<dim> &quadrature_formula =
      (advection_field.is_temperature() ?
       introspection.quadratures.temperature :
       introspection.quadratures.compositional_fields[advection_field.compositional_variable]);
    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar field = advection_field.scalar_extractor(introspection);

    FEValues<dim> fe_values (finite_element, quadrature_formula,
                             update_values | update_JxW_values);
    std::vector<double> old_field_values(n_q_points);
    std::vector<double> old_old_field_values(n_q_points);

    double min_entropy = std::numeric_limits<double>::max(),
           max_entropy = std::numeric_limits<double>::lowest(),
           area = 0,
           entropy_integrated = 0;

    // loop over all locally owned cells and evaluate the entropy
    // at all quadrature points. keep a running tally of the
    // integral over the entropy as well as the area and the
    // maximal and minimal entropy
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[field].get_function_values (old_solution,
                                                old_field_values);
          fe_values[field].get_function_values (old_old_solution,
                                                old_old_field_values);
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double field_value = (old_field_values[q] +
                                          old_old_field_values[q]) / 2;
              const double entropy = ((field_value-average_field) *
                                      (field_value-average_field));

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

    double max_residual = 0.0;
    double max_velocity = 0.0;
    double max_advection_prefactor = 0.0;
    double max_conductivity = 0.0;

    std::vector<double> residual (scratch.finite_element_values.n_quadrature_points,0.0);

    for (unsigned int i=0; i<assemblers->advection_system[advection_field.field_index()].size(); ++i)
      {
        const std::vector<double> new_residual = assemblers->advection_system[advection_field.field_index()][i]->compute_residual(scratch);
        for (unsigned int j=0; j<residual.size(); ++j)
          residual[j] += new_residual[j];

        auto *stabilization_assembler =
          dynamic_cast<Assemblers::AdvectionStabilizationInterface<dim>*> ((assemblers->advection_system[advection_field.field_index()][i]).get());

        Assert (stabilization_assembler != nullptr,
                ExcMessage("Entropy viscosity can only be computed for assemblers that "
                           "are derived from the advection stabilization assembler interface."));

        // Ensure no other assembler has set max_advection_prefactor or max_conductivity before,
        // otherwise we dont know which one to use.
        Assert (max_advection_prefactor == 0.0 && max_conductivity == 0.0,
                ExcMessage("More than one assembler has provided scaling factors for the entropy "
                           "viscosity stabilization, which is not supported. Make sure only one active advection "
                           "assembler is derived from the class AdvectionStabilizationInterface."));

        const std::vector<double> advection_prefactors = stabilization_assembler->advection_prefactors(scratch);
        const std::vector<double> conductivities = stabilization_assembler->diffusion_prefactors(scratch);

        max_advection_prefactor = *std::max_element(advection_prefactors.begin(),advection_prefactors.end());
        max_conductivity = *std::max_element(conductivities.begin(),conductivities.end());
      }

    std::vector<Tensor<1,dim>> old_fluid_velocity_values(scratch.finite_element_values.n_quadrature_points);
    std::vector<Tensor<1,dim>> old_old_fluid_velocity_values(scratch.finite_element_values.n_quadrature_points);
    if (parameters.include_melt_transport)
      {
        const FEValuesExtractors::Vector ex_u_f = introspection.variable("fluid velocity").extractor_vector();
        scratch.finite_element_values[ex_u_f].get_function_values (old_solution,old_fluid_velocity_values);
        scratch.finite_element_values[ex_u_f].get_function_values (old_old_solution,old_old_fluid_velocity_values);
      }

    for (unsigned int q=0; q < scratch.finite_element_values.n_quadrature_points; ++q)
      {
        const Tensor<1,dim> velocity = (scratch.old_velocity_values[q] +
                                        scratch.old_old_velocity_values[q]) / 2;
        double velocity_norm = velocity.norm();

        if (parameters.include_melt_transport)
          {
            const Tensor<1,dim> fluid_velocity = (old_fluid_velocity_values[q] +
                                                  old_old_fluid_velocity_values[q]) / 2;
            velocity_norm = std::max (fluid_velocity.norm(), velocity_norm);
          }

        const double strain_rate = ((scratch.old_strain_rates[q]
                                     + scratch.old_old_strain_rates[q]) / 2).norm();

        if (parameters.stabilization_alpha == 2)
          {
            const double field = (scratch.old_field_values[q] + scratch.old_old_field_values[q]) / 2;
            residual[q] *= std::abs(field - average_field);
          }

        max_residual = std::max (residual[q],     max_residual);
        max_velocity = std::max (velocity_norm
                                 + parameters.stabilization_gamma * strain_rate * cell_diameter,
                                 max_velocity);
      }

    // If the velocity is 0 we have to assume a sensible velocity to calculate
    // an artificial diffusion. We choose similar to nondimensional
    // formulations: v ~ thermal_diffusivity / length_scale, which cancels
    // the density and specific heat from the entropy formulation. It seems
    // surprising at first that only the conductivity remains, but remember
    // that this actually *is* an additional artificial diffusion.
    if (std::abs(global_u_infty) < 1e-50)
      return parameters.stabilization_beta[advection_field.field_index()] *
             max_conductivity / geometry_model->length_scale() *
             cell_diameter;

    const double maximum_viscosity = parameters.stabilization_beta[advection_field.field_index()] *
                                     max_advection_prefactor *
                                     max_velocity * cell_diameter;

    if (timestep_number <= 1
        || std::abs(global_entropy_variation) < 1e-50
        || std::abs(global_field_variation) < 1e-50)
      // we don't have sensible time-steps during the first two iterations
      // and we can not divide by the entropy_variation if it is zero
      return maximum_viscosity;
    else
      {
        Assert (old_time_step > 0, ExcInternalError());

        double entropy_viscosity;
        if (parameters.stabilization_alpha == 2)
          entropy_viscosity = (parameters.stabilization_c_R[advection_field.field_index()] *
                               cell_diameter * cell_diameter *
                               max_residual /
                               global_entropy_variation);
        else
          entropy_viscosity = (parameters.stabilization_c_R[advection_field.field_index()] *
                               cell_diameter * global_Omega_diameter *
                               max_velocity * max_residual /
                               (global_u_infty * global_field_variation));


        return std::min (maximum_viscosity, entropy_viscosity);
      }
  }



  template <int dim>
  template <typename T>
  void
  Simulator<dim>::
  get_artificial_viscosity (Vector<T> &viscosity_per_cell,
                            const AdvectionField &advection_field,
                            const bool skip_interior_cells) const
  {
    Assert(viscosity_per_cell.size()==triangulation.n_active_cells(), ExcInternalError());

    if (advection_field.field_type == AdvectionField::compositional_field)
      Assert(introspection.n_compositional_fields > advection_field.compositional_variable, ExcInternalError());

    viscosity_per_cell = 0.0;

    // discontinuous Galerkin doesn't require an artificial viscosity
    if (advection_field.is_discontinuous(introspection))
      return;

    bool skip_EV_dirichlet_boundary_cells = false;

    if (advection_field.is_temperature())
      skip_EV_dirichlet_boundary_cells = true;
    else
      {
        const std::string field_name = introspection.name_for_compositional_index(advection_field.compositional_variable);
        if (std::find(parameters.compositional_fields_with_disabled_boundary_entropy_viscosity.begin(),
                      parameters.compositional_fields_with_disabled_boundary_entropy_viscosity.end(),
                      field_name)
            !=
            parameters.compositional_fields_with_disabled_boundary_entropy_viscosity.end())
          skip_EV_dirichlet_boundary_cells = true;
      }

    const std::pair<double,double>
    global_field_range = get_extrapolated_advection_field_range (advection_field);
    const double global_entropy_variation = get_entropy_variation ((global_field_range.first +
                                                                    global_field_range.second) / 2,
                                                                   advection_field);
    const double global_max_velocity = get_maximal_velocity(old_solution);

    UpdateFlags update_flags = update_values |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values;

    if (advection_field.is_temperature())
      update_flags = update_flags | update_hessians;

    for (unsigned int i = 0; i < assemblers->advection_system_assembler_properties.size(); ++i)
      update_flags = update_flags | assemblers->advection_system_assembler_properties[i].needed_update_flags;

    // We need the face integrals to determine if a Dirichlet boundary
    // is conduction dominated in which case we disable stabilization
    // to get the most accurate heat flux, or still advection
    // dominated (e.g. a boundary with normal flow) in which case we
    // still need stabilization.
    const UpdateFlags face_update_flags = update_values |
                                          update_quadrature_points |
                                          update_normal_vectors |
                                          update_JxW_values;

    internal::Assembly::Scratch::
    AdvectionSystem<dim> scratch (finite_element,
                                  finite_element.base_element(advection_field.base_element(introspection)),
                                  *mapping,
                                  QGauss<dim>(advection_field.polynomial_degree(introspection)
                                              +
                                              (parameters.stokes_velocity_degree+1)/2),
                                  QTrapezoid<dim-1> (),
                                  update_flags,
                                  face_update_flags,
                                  introspection.n_compositional_fields,
                                  advection_field);

    std::vector<Tensor<1,dim>> face_old_velocity_values (scratch.face_finite_element_values->n_quadrature_points);
    std::vector<Tensor<1,dim>> face_old_old_velocity_values (scratch.face_finite_element_values->n_quadrature_points);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        // Skip cells for which we can not/do not need to compute the
        // stabilization. We need to compute the artificial viscosity
        // on all locally owned cells, but if we want to
        // smooth/average it over a neighborhood of a locally owned
        // cell, then we also need it on ghost cells; we could get it
        // there through parallel communication, but the easier way is
        // to simply compute it there as well
        if (cell->is_artificial()
            ||
            (cell->is_ghost() &&
             parameters.use_artificial_viscosity_smoothing == false))
          {
            viscosity_per_cell[cell->active_cell_index()] = numbers::signaling_nan<T>();
            continue;
          }
        // Also skip all interior cells if we are asked to do so. Do not
        // skip neighbor cells of boundary cells if smoothing is on, because
        // the smoothing uses both the boundary cell and its neighbor.
        else if (skip_interior_cells && !cell->at_boundary())
          {
            bool neighbor_at_boundary = false;
            for (const unsigned int face_no : cell->face_indices())
              if (cell->neighbor(face_no)->at_boundary() == true)
                neighbor_at_boundary = true;

            if (parameters.use_artificial_viscosity_smoothing == false ||
                neighbor_at_boundary == false)
              {
                viscosity_per_cell[cell->active_cell_index()] = numbers::signaling_nan<T>();
                continue;
              }
          }

        // For fields that have physical diffusion (e.g. temperature),
        // we can disable artificial viscosity stabilization at
        // Dirichlet boundaries, because the boundary is conduction
        // dominated anyway. Moreover, the residual we would compute
        // would be erroneously large, because it does not take into
        // account the boundary constraints. This would lead to
        // unnecessary large diffusion in the cells that matter most
        // for the overall energy balance of the system. However, we
        // sometimes have Dirichlet temperature boundary conditions
        // with prescribed non-tangential velocities, in these cases
        // we need the stabilization, because the boundary cells can
        // be advection dominated. Hence, only disable artificial
        // viscosity if flow through the boundary is slow, or
        // tangential. Also disable artificial viscosity for compositions
        // for which this is requested (e.g. compositions with
        // physical diffusion).
        if (parameters.advection_stabilization_method
            == Parameters<dim>::AdvectionStabilizationMethod::entropy_viscosity
            && skip_EV_dirichlet_boundary_cells == true)
          {
            const std::set<types::boundary_id> &dirichlet_boundaries =
              (advection_field.is_temperature() == true)
              ?
              boundary_temperature_manager.get_fixed_temperature_boundary_indicators()
              :
              boundary_composition_manager.get_fixed_composition_boundary_indicators();
            const std::set<types::boundary_id> &tangential_velocity_boundaries =
              boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
            const std::set<types::boundary_id> &zero_velocity_boundaries =
              boundary_velocity_manager.get_zero_boundary_velocity_indicators();

            bool cell_at_conduction_dominated_dirichlet_boundary = false;
            for (const unsigned int face_no : cell->face_indices())
              if (cell->at_boundary(face_no) == true &&
                  dirichlet_boundaries.find(cell->face(face_no)->boundary_id()) != dirichlet_boundaries.end())
                {
                  // If the velocity is tangential or zero we can always disable stabilization, except if there is another
                  // face at a different boundary. Therefore continue with the next face rather than break the loop.
                  if ((tangential_velocity_boundaries.find(cell->face(face_no)->boundary_id())
                       != tangential_velocity_boundaries.end())
                      ||
                      (zero_velocity_boundaries.find(cell->face(face_no)->boundary_id())
                       != zero_velocity_boundaries.end()))
                    {
                      cell_at_conduction_dominated_dirichlet_boundary = true;
                      continue;   // test next face
                    }

                  scratch.face_finite_element_values->reinit (cell, face_no);
                  (*scratch.face_finite_element_values)[introspection.extractors.velocities].get_function_values(old_solution,
                      face_old_velocity_values);
                  (*scratch.face_finite_element_values)[introspection.extractors.velocities].get_function_values(old_old_solution,
                      face_old_old_velocity_values);

                  // ... check if the face is a boundary with normal flow by integrating the normal velocities
                  // (flux through the boundary) as: int u*n ds = Sum_q u(x_q)*n(x_q) JxW(x_q)...
                  double normal_flow = 0.0;
                  double flow = 0.0;
                  double area = 0.0;
                  for (unsigned int q=0; q<scratch.face_finite_element_values->n_quadrature_points; ++q)
                    {
                      normal_flow += ((face_old_velocity_values[q]+face_old_old_velocity_values[q])/2.0 *
                                      scratch.face_finite_element_values->normal_vector(q)) *
                                     scratch.face_finite_element_values->JxW(q);
                      flow += ((face_old_velocity_values[q]+face_old_old_velocity_values[q])/2.0).norm() *
                              scratch.face_finite_element_values->JxW(q);
                      area += scratch.face_finite_element_values->JxW(q);
                    }

                  // Disable stabilization for boundaries with slow flow, or tangential flow.
                  // Break the loop in case a face is at multiple boundaries, some with flow, some without.
                  // In those cases we can not disable stabilization.
                  if ((std::abs(flow/area) * time_step
                       < std::sqrt(std::numeric_limits<double>::epsilon()) * cell->diameter())
                      ||
                      (std::abs(normal_flow)
                       < std::sqrt(std::numeric_limits<double>::epsilon()) * std::abs(flow)))
                    {
                      cell_at_conduction_dominated_dirichlet_boundary = true;
                    }
                  else
                    {
                      cell_at_conduction_dominated_dirichlet_boundary = false;
                      break; // no need to check any other face
                    }
                }

            if (cell_at_conduction_dominated_dirichlet_boundary)
              {
                // If we set the viscosity to zero, we don't need any further computation on this cell
                viscosity_per_cell[cell->active_cell_index()] = 0.0;
                continue;   // next cell
              }
          }

        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

        // also have the number of dofs that correspond just to the element for
        // the system we are currently trying to assemble
        const unsigned int advection_dofs_per_cell = scratch.phi_field.size();
        (void)advection_dofs_per_cell;
        Assert (advection_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
        Assert (scratch.grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
        Assert (scratch.phi_field.size() == advection_dofs_per_cell, ExcInternalError());

        const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

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

        for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
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

        if (update_flags & update_hessians)
          {
            scratch.finite_element_values[solution_field].get_function_laplacians (old_solution,
                                                                                   scratch.old_field_laplacians);
            scratch.finite_element_values[solution_field].get_function_laplacians (old_old_solution,
                                                                                   scratch.old_old_field_laplacians);
          }

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

            for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
              scratch.material_model_inputs.composition[q][c] = (scratch.old_composition_values[c][q] + scratch.old_old_composition_values[c][q]) / 2;
            scratch.material_model_inputs.strain_rate[q] = (scratch.old_strain_rates[q] + scratch.old_old_strain_rates[q]) / 2;
          }
        scratch.material_model_inputs.current_cell = cell;

        for (unsigned int i=0; i<assemblers->advection_system[advection_field.field_index()].size(); ++i)
          assemblers->advection_system[advection_field.field_index()][i]->create_additional_material_model_outputs(scratch.material_model_outputs);
        heating_model_manager.create_additional_material_model_inputs_and_outputs(scratch.material_model_inputs,
                                                                                  scratch.material_model_outputs);

        material_model->fill_additional_material_model_inputs(scratch.material_model_inputs,
                                                              solution,
                                                              scratch.finite_element_values,
                                                              introspection);
        material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);
        heating_model_manager.evaluate(scratch.material_model_inputs,scratch.material_model_outputs,scratch.heating_model_outputs);

        if (parameters.formulation_temperature_equation
            == Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
          {
            // Overwrite the density by the reference density coming from the
            // adiabatic conditions as required by the formulation
            for (unsigned int q=0; q<n_q_points; ++q)
              scratch.material_model_outputs.densities[q] = adiabatic_conditions->density(scratch.material_model_inputs.position[q]);
          }
        else if (parameters.formulation_temperature_equation
                 == Parameters<dim>::Formulation::TemperatureEquation::real_density)
          {
            // use real density
          }
        else
          AssertThrow(false, ExcNotImplemented());

        MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                                   cell,
                                                   scratch.finite_element_values.get_quadrature(),
                                                   scratch.finite_element_values.get_mapping(),
                                                   scratch.material_model_inputs.requested_properties,
                                                   scratch.material_model_outputs);

        if (parameters.advection_stabilization_method == Parameters<dim>::AdvectionStabilizationMethod::entropy_viscosity)
          {
            viscosity_per_cell[cell->active_cell_index()] = compute_viscosity(scratch,
                                                                              global_max_velocity,
                                                                              global_field_range.second - global_field_range.first,
                                                                              0.5 * (global_field_range.second + global_field_range.first),
                                                                              global_entropy_variation,
                                                                              cell->diameter(),
                                                                              advection_field);
          }
        else if (parameters.advection_stabilization_method == Parameters<dim>::AdvectionStabilizationMethod::supg)
          {
            double norm_of_advection_term = 0.0;
            double max_conductivity_on_cell = 0.0;

            {
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  if (advection_field.is_temperature())
                    {
                      norm_of_advection_term =
                        std::max(scratch.current_velocity_values[q].norm()*
                                 (scratch.material_model_outputs.densities[q] *
                                  scratch.material_model_outputs.specific_heat[q] +
                                  scratch.heating_model_outputs.lhs_latent_heat_terms[q]),
                                 norm_of_advection_term);

                      max_conductivity_on_cell =
                        std::max(scratch.material_model_outputs.thermal_conductivities[q],max_conductivity_on_cell);
                    }
                  else
                    {
                      norm_of_advection_term =
                        std::max(scratch.current_velocity_values[q].norm(),norm_of_advection_term);

                      max_conductivity_on_cell = 0.0;
                    }
                }
            }

            const double fe_order = advection_field.polynomial_degree(introspection);
            const double h = cell->diameter();
            const double eps = max_conductivity_on_cell;

            // SUPG parameter design from "On Discontinuity-Capturing Methods
            // for Convection-Diffusion Equations" by Volker John and Petr
            // Knobloch. Also see deal.II step-63:
            // delta_k = h / (2 \|u\| k) * (coth(Pe) - 1/Pe)
            // Pe = \| u \| h/(2 p eps)
            const double peclet_times_eps = norm_of_advection_term * h / (2.0 * fe_order);

            // Instead of Pe < 1, we check Pe*eps < eps as eps can be ==0:
            if (peclet_times_eps==0.0 || peclet_times_eps < eps)
              {
                // Diffusion dominant case, no stabilization needed:
                viscosity_per_cell[cell->active_cell_index()] = 0.0;
              }
            else
              {
                // To avoid a division by zero, increase eps slightly. The actual value is not
                // important, as long as the result is still a valid number. Note that this
                // is only important if \|u\| and eps are zero.
                const double peclet = peclet_times_eps / (eps + 1e-100);
                const double coth_of_peclet = (1.0 + std::exp(-2.0*peclet)) / (1.0 - std::exp(-2.0*peclet));
                const double delta = h/(2.0*norm_of_advection_term*fe_order) * (coth_of_peclet - 1.0/peclet);
                viscosity_per_cell[cell->active_cell_index()] = delta;
              }
            Assert (viscosity_per_cell[cell->active_cell_index()] >= 0, ExcMessage ("tau for SUPG needs to be a nonnegative constant."));
          }
        else
          AssertThrow(false, ExcNotImplemented());
      }

    // if set to true, the maximum of the artificial viscosity in the cell as well
    // as the neighbors of the cell is computed and used instead
    if (parameters.use_artificial_viscosity_smoothing  == true)
      {
        Vector<T> viscosity_per_cell_temp;
        viscosity_per_cell_temp.reinit(triangulation.n_active_cells());

        viscosity_per_cell_temp = viscosity_per_cell;
        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              if (skip_interior_cells && !cell->at_boundary())
                continue;

              for (const unsigned int face_no : cell->face_indices())
                if (cell->at_boundary(face_no) == false)
                  {
                    if (cell->neighbor(face_no)->is_active())
                      viscosity_per_cell[cell->active_cell_index()] = std::max(viscosity_per_cell[cell->active_cell_index()],
                                                                               viscosity_per_cell_temp[cell->neighbor(face_no)->active_cell_index()]);
                    else
                      for (unsigned int l=0; l<cell->neighbor(face_no)->n_children(); ++l)
                        if (cell->neighbor(face_no)->child(l)->is_active())
                          viscosity_per_cell[cell->active_cell_index()] = std::max(viscosity_per_cell[cell->active_cell_index()],
                                                                                   viscosity_per_cell_temp[cell->neighbor(face_no)->child(l)->active_cell_index()]);
                  }
            }
      }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::get_artificial_viscosity (Vector<double> &viscosity_per_cell,  \
                                                          const AdvectionField &advection_field, \
                                                          const bool skip_interior_cells) const; \
  template void Simulator<dim>::get_artificial_viscosity (Vector<float> &viscosity_per_cell,  \
                                                          const AdvectionField &advection_field, \
                                                          const bool skip_interior_cells) const;


  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}

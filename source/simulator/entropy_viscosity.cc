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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/simulator.h>
#include <aspect/assembly.h>
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

    // record maximal entropy on Gauss quadrature
    // points
    const QGauss<dim> quadrature_formula (parameters.temperature_degree+1);
    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar field = advection_field.scalar_extractor(introspection);

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

    std::vector<double> residual = assemblers->compute_advection_system_residual(scratch.material_model_inputs.current_cell,
                                                                                 advection_field,
                                                                                 scratch);

    double max_residual = 0;
    double max_velocity = 0;
    double max_density = (advection_field.is_temperature()) ? 0.0 : 1.0;
    double max_specific_heat = (advection_field.is_temperature()) ? 0.0 : 1.0;
    double max_conductivity = 0;

    std::vector<Tensor<1,dim> > old_fluid_velocity_values(scratch.finite_element_values.n_quadrature_points);
    std::vector<Tensor<1,dim> > old_old_fluid_velocity_values(scratch.finite_element_values.n_quadrature_points);
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
      // we don't have sensible time-steps during the first two iterations
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
      Assert(introspection.n_compositional_fields > advection_field.compositional_variable, ExcInternalError());

    viscosity_per_cell = 0.0;

    // discontinuous Galerkin doesn't require an artificial viscosity
    if (advection_field.is_discontinuous(introspection))
      return;

    const std::pair<double,double>
    global_field_range = get_extrapolated_advection_field_range (advection_field);
    double global_entropy_variation = get_entropy_variation ((global_field_range.first +
                                                              global_field_range.second) / 2,
                                                             advection_field);
    double global_max_velocity = get_maximal_velocity(old_solution);

    const UpdateFlags update_flags = update_values |
                                     update_gradients |
                                     (advection_field.is_temperature() ? update_hessians : update_default) |
                                     update_quadrature_points |
                                     update_JxW_values;

    /* Because we can only get here in the continuous case, which never requires
     * face integrals, we supply no face_update_flags and an invalid
     * face_quadrature to the scratch object to reduce the initialization cost.
     */
    const UpdateFlags face_update_flags = update_default;

    internal::Assembly::Scratch::
    AdvectionSystem<dim> scratch (finite_element,
                                  finite_element.base_element(advection_field.base_element(introspection)),
                                  *mapping,
                                  QGauss<dim>(advection_field.polynomial_degree(introspection)
                                              +
                                              (parameters.stokes_velocity_degree+1)/2),
                                  Quadrature<dim-1> (),
                                  update_flags,
                                  face_update_flags,
                                  introspection.n_compositional_fields);

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

        if (advection_field.is_temperature())
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
        create_additional_material_model_outputs(scratch.material_model_outputs);

        material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);

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
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::get_artificial_viscosity (Vector<double> &viscosity_per_cell,  \
                                                          const AdvectionField &advection_field) const; \
  template void Simulator<dim>::get_artificial_viscosity (Vector<float> &viscosity_per_cell,  \
                                                          const AdvectionField &advection_field) const; \
   

  ASPECT_INSTANTIATE(INSTANTIATE)
}

/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
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
        struct StokesPreconditioner
        {
          StokesPreconditioner (const FiniteElement<dim> &finite_element,
                                const Quadrature<dim>    &quadrature,
                                const Mapping<dim>       &mapping,
                                const UpdateFlags         update_flags,
                                const unsigned int        n_compositional_fields);
          StokesPreconditioner (const StokesPreconditioner &data);

          FEValues<dim>               finite_element_values;

          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  phi_p;

          std::vector<double>                  temperature_values;
          std::vector<double>                  pressure_values;
          std::vector<SymmetricTensor<2,dim> > strain_rates;
          std::vector<std::vector<double> >     composition_values;

          typename MaterialModel::Interface<dim>::MaterialModelInputs material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs material_model_outputs;
        };



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const FiniteElement<dim> &finite_element,
                              const Quadrature<dim>    &quadrature,
                              const Mapping<dim>       &mapping,
                              const UpdateFlags         update_flags,
                              const unsigned int        n_compositional_fields)
          :
          finite_element_values (mapping, finite_element, quadrature,
                                 update_flags),
          grads_phi_u (finite_element.dofs_per_cell),
          phi_p (finite_element.dofs_per_cell),
          temperature_values (quadrature.size()),
          pressure_values (quadrature.size()),
          strain_rates (quadrature.size()),
          composition_values(n_compositional_fields,
                             std::vector<double>(quadrature.size())),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size())
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
          temperature_values (scratch.temperature_values),
          pressure_values (scratch.pressure_values),
          strain_rates (scratch.strain_rates),
          composition_values(scratch.composition_values),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs)
        {}



        // We derive the StokesSystem scratch array from the
        // StokesPreconditioner array. We do this because all the objects that
        // are necessary for the assembly of the preconditioner are also
        // needed for the actual matrix system and right hand side, plus some
        // extra data that we need for the time stepping right hand side.
        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &finite_element,
                        const Mapping<dim>       &mapping,
                        const Quadrature<dim>    &quadrature,
                        const UpdateFlags         update_flags,
                        const unsigned int        n_compositional_fields);

          StokesSystem (const StokesSystem<dim> &data);

          std::vector<Tensor<1,dim> >          phi_u;
          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  div_phi_u;
          std::vector<Tensor<1,dim> >          velocity_values;
          std::vector<std::vector<double> >     composition_values;

          typename MaterialModel::Interface<dim>::MaterialModelInputs material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs material_model_outputs;
        };



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &finite_element,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        n_compositional_fields)
          :
          StokesPreconditioner<dim> (finite_element, quadrature,
                                     mapping,
                                     update_flags, n_compositional_fields),
          phi_u (finite_element.dofs_per_cell),
          grads_phi_u (finite_element.dofs_per_cell),
          div_phi_u (finite_element.dofs_per_cell),
          velocity_values (quadrature.size()),
          composition_values(n_compositional_fields,
                             std::vector<double>(quadrature.size())),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size())
        {}



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &scratch)
          :
          StokesPreconditioner<dim> (scratch),
          phi_u (scratch.phi_u),
          grads_phi_u (scratch.grads_phi_u),
          div_phi_u (scratch.div_phi_u),
          velocity_values (scratch.velocity_values),
          composition_values(scratch.composition_values),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs)
        {}

        template <int dim>
        struct AdvectionSystem
        {
          AdvectionSystem (const FiniteElement<dim> &finite_element,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &quadrature,
                           const unsigned int        n_compositional_fields);
          AdvectionSystem (const AdvectionSystem &data);

          FEValues<dim>               finite_element_values;

          std::vector<double>         phi_field;
          std::vector<Tensor<1,dim> > grad_phi_field;

          std::vector<Tensor<1,dim> > old_velocity_values;
          std::vector<Tensor<1,dim> > old_old_velocity_values;

          std::vector<double>         old_pressure;
          std::vector<double>         old_old_pressure;

          std::vector<SymmetricTensor<2,dim> > old_strain_rates;
          std::vector<SymmetricTensor<2,dim> > old_old_strain_rates;

          std::vector<double>         old_temperature_values;
          std::vector<double>         old_old_temperature_values;

          std::vector<double>        *old_field_values;
          std::vector<double>        *old_old_field_values;
          std::vector<Tensor<1,dim> > old_field_grads;
          std::vector<Tensor<1,dim> > old_old_field_grads;
          std::vector<double>         old_field_laplacians;
          std::vector<double>         old_old_field_laplacians;

          std::vector<std::vector<double> > old_composition_values;
          std::vector<std::vector<double> > old_old_composition_values;

          std::vector<double>         current_temperature_values;
          std::vector<Tensor<1,dim> > current_velocity_values;
          std::vector<SymmetricTensor<2,dim> > current_strain_rates;
          std::vector<double>         current_pressure_values;
          std::vector<std::vector<double> > current_composition_values;

          typename MaterialModel::Interface<dim>::MaterialModelInputs material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs material_model_outputs;

          typename MaterialModel::Interface<dim>::MaterialModelInputs explicit_material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs explicit_material_model_outputs;
        };



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &quadrature,
                         const unsigned int        n_compositional_fields)
          :
          finite_element_values (mapping,
                                 finite_element, quadrature,
                                 update_values    |
                                 update_gradients |
                                 update_hessians  |
                                 update_quadrature_points |
                                 update_JxW_values),

          phi_field (finite_element.dofs_per_cell),
          grad_phi_field (finite_element.dofs_per_cell),
          old_velocity_values (quadrature.size()),
          old_old_velocity_values (quadrature.size()),
          old_pressure (quadrature.size()),
          old_old_pressure (quadrature.size()),
          old_strain_rates (quadrature.size()),
          old_old_strain_rates (quadrature.size()),
          old_temperature_values (quadrature.size()),
          old_old_temperature_values(quadrature.size()),
          old_field_grads(quadrature.size()),
          old_old_field_grads(quadrature.size()),
          old_field_laplacians(quadrature.size()),
          old_old_field_laplacians(quadrature.size()),
          old_composition_values(n_compositional_fields,
                                 std::vector<double>(quadrature.size())),
          old_old_composition_values(n_compositional_fields,
                                     std::vector<double>(quadrature.size())),
          current_temperature_values(quadrature.size()),
          current_velocity_values(quadrature.size()),
          current_strain_rates(quadrature.size()),
          current_pressure_values(quadrature.size()),
          current_composition_values(n_compositional_fields,
                                     std::vector<double>(quadrature.size())),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size()),
          explicit_material_model_inputs(quadrature.size(), n_compositional_fields),
          explicit_material_model_outputs(quadrature.size())
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const AdvectionSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),

          phi_field (scratch.phi_field),
          grad_phi_field (scratch.grad_phi_field),
          old_velocity_values (scratch.old_velocity_values),
          old_old_velocity_values (scratch.old_old_velocity_values),
          old_pressure (scratch.old_pressure),
          old_old_pressure (scratch.old_old_pressure),
          old_strain_rates (scratch.old_strain_rates),
          old_old_strain_rates (scratch.old_old_strain_rates),
          old_temperature_values (scratch.old_temperature_values),
          old_old_temperature_values (scratch.old_old_temperature_values),
          old_field_grads (scratch.old_field_grads),
          old_old_field_grads (scratch.old_old_field_grads),
          old_field_laplacians (scratch.old_field_laplacians),
          old_old_field_laplacians (scratch.old_old_field_laplacians),
          old_composition_values(scratch.old_composition_values),
          old_old_composition_values(scratch.old_old_composition_values),
          current_temperature_values(scratch.current_temperature_values),
          current_velocity_values(scratch.current_velocity_values),
          current_strain_rates(scratch.current_strain_rates),
          current_pressure_values(scratch.current_pressure_values),
          current_composition_values(scratch.current_composition_values),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs),
          explicit_material_model_inputs(scratch.explicit_material_model_inputs),
          explicit_material_model_outputs(scratch.explicit_material_model_outputs)
        {}

      }


      // The CopyData arrays are similar to the
      // Scratch arrays. They provide a
      // constructor, a copy operation, and
      // some arrays for local matrix, local
      // vectors and the relation between local
      // and global degrees of freedom (a.k.a.
      // <code>local_dof_indices</code>).
      namespace CopyData
      {
        template <int dim>
        struct StokesPreconditioner
        {
          StokesPreconditioner (const FiniteElement<dim> &finite_element);
          StokesPreconditioner (const StokesPreconditioner &data);

          FullMatrix<double>          local_matrix;
          std::vector<unsigned int>   local_dof_indices;
        };



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
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &finite_element);
          StokesSystem (const StokesSystem<dim> &data);

          Vector<double> local_rhs;
          Vector<double> local_pressure_shape_function_integrals;
        };



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &finite_element)
          :
          StokesPreconditioner<dim> (finite_element),
          local_rhs (finite_element.dofs_per_cell),
          local_pressure_shape_function_integrals (finite_element.dofs_per_cell)
        {}



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &data)
          :
          StokesPreconditioner<dim> (data),
          local_rhs (data.local_rhs),
          local_pressure_shape_function_integrals (data.local_pressure_shape_function_integrals)
        {}



        template <int dim>
        struct AdvectionSystem
        {
          AdvectionSystem (const FiniteElement<dim> &finite_element);
          AdvectionSystem (const AdvectionSystem &data);

          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;
          std::vector<unsigned int>   local_dof_indices;
        };



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_rhs (finite_element.dofs_per_cell),
          local_dof_indices (finite_element.dofs_per_cell)
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const AdvectionSystem &data)
          :
          local_matrix (data.local_matrix),
          local_rhs (data.local_rhs),
          local_dof_indices (data.local_dof_indices)
        {}

      }
    }



    template <class T>
    inline
    T
    bdf2_extrapolate (const bool use_bdf_scheme,
                      const double old_time_step,
                      const double time_step,
                      const T &old_data,
                      const T &new_data)
    {
      return (use_bdf_scheme) ?
             (new_data * (1 + time_step/old_time_step)
              - old_data * time_step/old_time_step)
             : new_data;
    }
  }


  /**
   * Compute the variation in the entropy needed in the definition of the
   * artificial viscosity used to stabilize the composition/temperature equation.
   */
  template <int dim>
  double
  Simulator<dim>::get_entropy_variation (const double average_field,
                                         const TemperatureOrComposition &temperature_or_composition) const
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
      = (temperature_or_composition.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[temperature_or_composition.compositional_variable]
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

    Utilities::MPI::sum (local_for_sum, mpi_communicator, global_for_sum);
    Utilities::MPI::max (local_for_max, mpi_communicator, global_for_max);

    const double average_entropy = global_for_sum[0] / global_for_sum[1];

    // return the maximal deviation of the entropy everywhere from the
    // average value
    return std::max(global_for_max[1] - average_entropy,
                    average_entropy - (-global_for_max[0]));
  }

  template <int dim>
  void
  Simulator<dim>::
  compute_advection_system_residual(internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                    const double                        average_field,
                                    const TemperatureOrComposition     &temperature_or_composition,
                                    double                             &max_residual,
                                    double                             &max_velocity,
                                    double                             &max_density,
                                    double                             &max_specific_heat) const
  {
    const unsigned int n_q_points = scratch.old_field_values->size();

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (scratch.old_velocity_values[q] +
                                 scratch.old_old_velocity_values[q]) / 2;

        const SymmetricTensor<2,dim> strain_rate = scratch.material_model_inputs.strain_rate[q];

        const double dField_dt = ((*scratch.old_field_values)[q] - (*scratch.old_old_field_values)[q])
                                 / old_time_step;
        const double u_grad_field = u * (scratch.old_field_grads[q] +
                                         scratch.old_old_field_grads[q]) / 2;

        const double density              = ((temperature_or_composition.is_temperature())
                                             ? scratch.explicit_material_model_outputs.densities[q] : 1.0);
        const double conductivity = ((temperature_or_composition.is_temperature()) ? scratch.explicit_material_model_outputs.thermal_conductivities[q] : 0.0);
        const double c_P                  = ((temperature_or_composition.is_temperature()) ? scratch.explicit_material_model_outputs.specific_heat[q] : 1.0);
        const double k_Delta_field = conductivity
                                     * (scratch.old_field_laplacians[q] +
                                        scratch.old_old_field_laplacians[q]) / 2;

        const double field = ((*scratch.old_field_values)[q] + (*scratch.old_old_field_values)[q]) / 2;

        const double gamma
          = compute_heating_term(scratch,scratch.explicit_material_model_inputs,scratch.explicit_material_model_outputs,temperature_or_composition,q);
        double residual
          = std::abs(density * c_P * (dField_dt + u_grad_field) - k_Delta_field - gamma);

        if (parameters.stabilization_alpha == 2)
          residual *= std::abs(field - average_field);

        max_residual = std::max      (residual,        max_residual);
        max_velocity = std::max      (std::sqrt (u*u), max_velocity);
        max_density  = std::max      (density,         max_density);
        max_specific_heat = std::max (c_P,   max_specific_heat);
      }
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
                     const TemperatureOrComposition     &temperature_or_composition) const
  {
    if (global_u_infty == 0)
      return 5e-3 * cell_diameter;

    double max_residual = 0;
    double max_velocity = 0;
    double max_density = 0;
    double max_specific_heat = 0;

    compute_advection_system_residual(scratch,
                                      average_field,
                                      temperature_or_composition,
                                      max_residual,
                                      max_velocity,
                                      max_density,
                                      max_specific_heat);

    const double max_viscosity = parameters.stabilization_beta *
                                 max_density *
                                 max_specific_heat *
                                 max_velocity * cell_diameter;
    if (timestep_number <= 1)
      // we don't have sensible timesteps during the first two iterations
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
  void
  Simulator<dim>::
  compute_material_model_input_values (const TrilinosWrappers::MPI::BlockVector                    &input_solution,
                                       const FEValues<dim>                                         &input_finite_element_values,
                                       const bool                                                   compute_strainrate,
                                       typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs) const
  {
    unsigned int n_q_points = material_model_inputs.temperature.size();
    for (unsigned int q=0; q<n_q_points; ++q)
      material_model_inputs.position[q] = input_finite_element_values.quadrature_point(q);

    input_finite_element_values[introspection.extractors.temperature].get_function_values (input_solution,
        material_model_inputs.temperature);
    input_finite_element_values[introspection.extractors.pressure].get_function_values(input_solution,
        material_model_inputs.pressure);
    if (compute_strainrate)
      input_finite_element_values[introspection.extractors.velocities].get_function_symmetric_gradients(input_solution,
          material_model_inputs.strain_rate);

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
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                        internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    const unsigned int   dofs_per_cell   = finite_element.dofs_per_cell;
    const unsigned int   n_q_points      = scratch.finite_element_values.n_quadrature_points;

    scratch.finite_element_values.reinit (cell);

    data.local_matrix = 0;

    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         true,
                                         scratch.material_model_inputs);

    material_model->compute_parameters(scratch.material_model_inputs,scratch.material_model_outputs);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grads_phi_u[k] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(k,q);
            scratch.phi_p[k]       = scratch.finite_element_values[introspection.extractors.pressure].value (k, q);
          }

        const double eta = scratch.material_model_outputs.viscosities[q];

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            if (finite_element.system_to_component_index(i).first
                ==
                finite_element.system_to_component_index(j).first)
              data.local_matrix(i,j) += (eta *
                                         (scratch.grads_phi_u[i] *
                                          scratch.grads_phi_u[j])
                                         +
                                         (1./eta) *
                                         pressure_scaling *
                                         pressure_scaling *
                                         (scratch.phi_p[i] * scratch.phi_p[j]))
                                        * scratch.finite_element_values.JxW(q);
      }

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

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_stokes_preconditioner,
                          this,
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_stokes_preconditioner,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         StokesPreconditioner<dim> (finite_element, quadrature_formula,
                                    mapping,
                                    update_JxW_values |
                                    update_values |
                                    update_gradients |
                                    update_quadrature_points,
                                    parameters.n_compositional_fields),
         internal::Assembly::CopyData::
         StokesPreconditioner<dim> (finite_element));

    system_preconditioner_matrix.compress();
  }



  template <int dim>
  void
  Simulator<dim>::build_stokes_preconditioner ()
  {
    if (rebuild_stokes_preconditioner == false)
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
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = true;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;

    Mp_preconditioner->initialize (system_preconditioner_matrix.block(1,1));
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
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    scratch.finite_element_values.reinit (cell);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    if (material_model->is_compressible())
      data.local_pressure_shape_function_integrals = 0;

    // we only need the strain rates for the viscosity,
    // which we only need when rebuilding the matrix
    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         rebuild_stokes_matrix,
                                         scratch.material_model_inputs);

    material_model->compute_parameters(scratch.material_model_inputs,scratch.material_model_outputs);

    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.velocity_values);

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

        const double eta = (rebuild_stokes_matrix
                            ?
                            scratch.material_model_outputs.viscosities[q]
                            :
                            std::numeric_limits<double>::quiet_NaN());

        const Tensor<1,dim>
        gravity = gravity_model->gravity_vector (scratch.finite_element_values.quadrature_point(q));

        const double compressibility
          = (scratch.material_model_outputs.is_compressible
             ?
             scratch.material_model_outputs.compressibilities[q]
             :
             std::numeric_limits<double>::quiet_NaN() );
        const double density = scratch.material_model_outputs.densities[q];

        if (rebuild_stokes_matrix)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              data.local_matrix(i,j) += ( eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j])
                                          - (scratch.material_model_outputs.is_compressible
                                             ?
                                             eta * 2.0/3.0 * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                             :
                                             0)
                                          - (pressure_scaling *
                                             scratch.div_phi_u[i] * scratch.phi_p[j])
                                          - (pressure_scaling *
                                             scratch.phi_p[i] * scratch.div_phi_u[j]))
                                        * scratch.finite_element_values.JxW(q);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs(i) += (
                                 (density * gravity * scratch.phi_u[i])
                                 + (scratch.material_model_outputs.is_compressible
                                    ?
                                    (pressure_scaling *
                                     compressibility * density *
                                     (scratch.velocity_values[q] * gravity) *
                                     scratch.phi_p[i])
                                    :
                                    0)
                               )
                               * scratch.finite_element_values.JxW(q);
        if (scratch.material_model_outputs.is_compressible)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            data.local_pressure_shape_function_integrals(i) += scratch.phi_p[i] * scratch.finite_element_values.JxW(q);
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

    if (material_model->is_compressible())
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
    if (material_model->is_compressible())
      pressure_shape_function_integrals = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_stokes_system,
                          this,
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_stokes_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         StokesSystem<dim> (finite_element, mapping, quadrature_formula,
                            (update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             (rebuild_stokes_matrix == true
                              ?
                              update_gradients
                              :
                              UpdateFlags(0))),
                            parameters.n_compositional_fields),
         internal::Assembly::CopyData::
         StokesSystem<dim> (finite_element));

    system_matrix.compress();
    system_rhs.compress(VectorOperation::add);

    if (material_model->is_compressible())
      pressure_shape_function_integrals.compress(VectorOperation::add);


    // if we got here and have rebuilt the matrix, make sure that
    // the flag for the preconditioner had previously also been
    // set
    if (rebuild_stokes_matrix == true)
      Assert (rebuild_stokes_preconditioner == true,
              ExcInternalError());


    // record that we have just rebuilt the matrix
    rebuild_stokes_matrix = false;

    computing_timer.exit_section();
  }

  template <int dim>
  void
  Simulator<dim>::build_advection_preconditioner(const TemperatureOrComposition &temperature_or_composition,
                                                 std_cxx1x::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner)
  {
    switch (temperature_or_composition.field_type)
      {
        case TemperatureOrComposition::temperature_field:
        {
          computing_timer.enter_section ("   Build temperature preconditioner");

          preconditioner.reset (new TrilinosWrappers::PreconditionILU());
          preconditioner->initialize (system_matrix.block(2,2));

          computing_timer.exit_section();

          break;
        }

        case TemperatureOrComposition::compositional_field:
        {
          computing_timer.enter_section ("   Build composition preconditioner");

          const unsigned int block_number
            = 3+temperature_or_composition.compositional_variable;
          preconditioner.reset (new TrilinosWrappers::PreconditionILU());
          preconditioner->initialize (system_matrix.block(block_number,
                                                          block_number));

          computing_timer.exit_section();

          break;
        }

        default:
          Assert (false, ExcNotImplemented());
      }
  }


  template <int dim>
  double
  Simulator<dim>::compute_heating_term(const internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                       typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
                                       typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
                                       const TemperatureOrComposition     &temperature_or_composition,
                                       const unsigned int q) const
  {

    if (temperature_or_composition.field_type == TemperatureOrComposition::compositional_field)
      return 0.0;

    const double current_T = material_model_inputs.temperature[q];
    const SymmetricTensor<2,dim> current_strain_rate = material_model_inputs.strain_rate[q];
    const Tensor<1,dim> current_u = scratch.current_velocity_values[q];

    const double alpha                = material_model_outputs.thermal_expansion_coefficients[q];
    const double density              = material_model_outputs.densities[q];
    const double viscosity            = material_model_outputs.viscosities[q];
    const bool is_compressible        = material_model_outputs.is_compressible;
    const double compressibility      = (is_compressible
                                         ?
                                         material_model_outputs.compressibilities[q]
                                         :
                                         std::numeric_limits<double>::quiet_NaN() );

    const Tensor<1,dim>
    gravity = gravity_model->gravity_vector (scratch.finite_element_values.quadrature_point(q));

    const double gamma
      = (parameters.radiogenic_heating_rate * density
         +
         // add the term 2*eta*(eps - 1/3*(tr eps)1):(eps - 1/3*(tr eps)1)
         //
         // we can multiply this out to obtain
         //   2*eta*(eps:eps - 1/3*(tr eps)^2)
         // and can then use that in the compressible case we have
         //   tr eps = div u
         //          = -1/rho u . grad rho
         // and by the usual approximation we make,
         //   tr eps = -1/rho drho/dp u . grad p
         //          = -1/rho drho/dp rho (u . g)
         //          = - drho/dp (u . g)
         //          = - compressibility rho (u . g)
         // to yield the final form of the term:
         //   2*eta [eps:eps - 1/3 (compressibility * rho * (u.g))^2]
         (parameters.include_shear_heating
          ?
          2 * viscosity *
          current_strain_rate * current_strain_rate
          -
          (is_compressible
           ?
           2./3.*viscosity*std::pow(compressibility * density * (current_u * gravity),
                                    2)
           :
           0)
            :
            0)
           +
           // finally add the term from adiabatic compression heating
           //   - drho/dT T (u . g)
           // where we use the definition of
           //   alpha = - 1/rho drho/dT
           (parameters.include_adiabatic_heating
            ?
            alpha * density * (current_u*gravity) * current_T
            :
            0)
          );

    return gamma;
  }


  template <int dim>
  void Simulator<dim>::
  local_assemble_advection_system (const TemperatureOrComposition     &temperature_or_composition,
                                   const std::pair<double,double> global_field_range,
                                   const double                   global_max_velocity,
                                   const double                   global_entropy_variation,
                                   const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                   internal::Assembly::CopyData::AdvectionSystem<dim> &data)
{
    const bool use_bdf2_scheme = (timestep_number > 1);

    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    const FEValuesExtractors::Scalar solution_field
      = (temperature_or_composition.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[temperature_or_composition.compositional_variable]
        );

    scratch.finite_element_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_matrix = 0;
    data.local_rhs = 0;

    if (temperature_or_composition.is_temperature())
      {
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
      }
    else
      {
        scratch.finite_element_values[introspection.extractors.compositional_fields[temperature_or_composition.compositional_variable]].get_function_values(old_solution,
            scratch.old_composition_values[temperature_or_composition.compositional_variable]);
        scratch.finite_element_values[introspection.extractors.compositional_fields[temperature_or_composition.compositional_variable]].get_function_values(old_old_solution,
            scratch.old_old_composition_values[temperature_or_composition.compositional_variable]);
      }

    scratch.finite_element_values[introspection.extractors.velocities].get_function_values (old_solution,
        scratch.old_velocity_values);
    scratch.finite_element_values[introspection.extractors.velocities].get_function_values (old_old_solution,
        scratch.old_old_velocity_values);
    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.current_velocity_values);


    scratch.old_field_values = ((temperature_or_composition.is_temperature()) ? &scratch.old_temperature_values : &scratch.old_composition_values[temperature_or_composition.compositional_variable]);
    scratch.old_old_field_values = ((temperature_or_composition.is_temperature()) ? &scratch.old_old_temperature_values : &scratch.old_old_composition_values[temperature_or_composition.compositional_variable]);

    scratch.finite_element_values[solution_field].get_function_gradients (old_solution,
                                                                          scratch.old_field_grads);
    scratch.finite_element_values[solution_field].get_function_gradients (old_old_solution,
                                                                          scratch.old_old_field_grads);

    scratch.finite_element_values[solution_field].get_function_laplacians (old_solution,
                                                                           scratch.old_field_laplacians);
    scratch.finite_element_values[solution_field].get_function_laplacians (old_old_solution,
                                                                           scratch.old_old_field_laplacians);

    if (temperature_or_composition.is_temperature())
      {
        compute_material_model_input_values (current_linearization_point,
                                             scratch.finite_element_values,
                                             true,
                                             scratch.material_model_inputs);
        material_model->compute_parameters(scratch.material_model_inputs,scratch.material_model_outputs);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            scratch.explicit_material_model_inputs.temperature[q] = (scratch.old_temperature_values[q] + scratch.old_old_temperature_values[q]) / 2;
            scratch.explicit_material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
            scratch.explicit_material_model_inputs.pressure[q] = (scratch.old_pressure[q] + scratch.old_old_pressure[q]) / 2;
            for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
              scratch.explicit_material_model_inputs.composition[q][c] = (scratch.old_composition_values[c][q] + scratch.old_old_composition_values[c][q]) / 2;
            scratch.explicit_material_model_inputs.strain_rate[q] = (scratch.old_strain_rates[q] + scratch.old_old_strain_rates[q]) / 2;
          }
        material_model->compute_parameters(scratch.explicit_material_model_inputs,scratch.explicit_material_model_outputs);
      }

    // TODO: Compute artificial viscosity once per timestep instead of each time
    // temperature system is assembled (as this might happen more than once per
    // timestep for iterative solvers)
    const double nu
      = compute_viscosity (scratch,
                           global_max_velocity,
                           global_field_range.second - global_field_range.first,
                           0.5 * (global_field_range.second + global_field_range.first),
                           global_entropy_variation,
                           cell->diameter(),
                           temperature_or_composition);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_field[k] = scratch.finite_element_values[solution_field].gradient (k,q);
            scratch.phi_field[k]      = scratch.finite_element_values[solution_field].value (k, q);
          }

        const double density_c_P              =
          ((temperature_or_composition.is_temperature())
           ?
           scratch.material_model_outputs.densities[q] *
           scratch.material_model_outputs.specific_heat[q]
           :
           1.0);
        const double conductivity =
          ((temperature_or_composition.is_temperature())
           ?
           scratch.material_model_outputs.thermal_conductivities[q]
           :
           0.0);
        const double gamma = compute_heating_term(scratch,
                                                  scratch.material_model_inputs,
                                                  scratch.material_model_outputs,
                                                  temperature_or_composition,
                                                  q);

        const double field_term_for_rhs
          = (use_bdf2_scheme ?
             ((*scratch.old_field_values)[q] *
              (1 + time_step/old_time_step)
              -
              (*scratch.old_old_field_values)[q] *
              (time_step * time_step) /
              (old_time_step * (time_step + old_time_step)))
             :
             (*scratch.old_field_values)[q])
            *
            density_c_P;


        const Tensor<1,dim> current_u = scratch.current_velocity_values[q];

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            data.local_rhs(i) += (field_term_for_rhs * scratch.phi_field[i]
                                  + time_step *
                                  gamma
                                  * scratch.phi_field[i])
                                 *
                                 scratch.finite_element_values.JxW(q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                          (time_step + old_time_step)) : 1.0;
                data.local_matrix(i,j)
                += (
                     (time_step * (conductivity + nu)
                      * scratch.grad_phi_field[i] * scratch.grad_phi_field[j])
                     + ((time_step * (current_u * scratch.grad_phi_field[j] * scratch.phi_field[i]))
                        + (factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                     density_c_P
                   )
                   * scratch.finite_element_values.JxW(q);
              }
          }
      }
  }

  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_advection_system (const internal::Assembly::CopyData::AdvectionSystem<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);
  }



  template <int dim>
  void Simulator<dim>::assemble_advection_system (const TemperatureOrComposition &temperature_or_composition)
  {
    if (temperature_or_composition.is_temperature())
      {
        computing_timer.enter_section ("   Assemble temperature system");
        system_matrix.block (2,2) = 0;
      }
    else
      {
        computing_timer.enter_section ("   Assemble composition system");
        system_matrix.block(3+temperature_or_composition.compositional_variable,
                            3+temperature_or_composition.compositional_variable) = 0;
      }
    system_rhs = 0;

    const std::pair<double,double>
    global_field_range = get_extrapolated_temperature_or_composition_range (temperature_or_composition);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_advection_system,
                          this,
                          temperature_or_composition,
                          global_field_range,
                          get_maximal_velocity(old_solution),
                          // use the mid-value of the advected field instead of the
                          // integral mean. results are not very
                          // sensitive to this and this is far simpler
                          get_entropy_variation ((global_field_range.first +
                                                  global_field_range.second) / 2,
                                                 temperature_or_composition),
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_advection_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         AdvectionSystem<dim> (finite_element, mapping, QGauss<dim>(parameters.composition_degree+2),
                               parameters.n_compositional_fields),
         internal::Assembly::CopyData::
         AdvectionSystem<dim> (finite_element));

    system_matrix.compress();
    system_rhs.compress(VectorOperation::add);

    computing_timer.exit_section();
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
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
  template void Simulator<dim>::build_advection_preconditioner (const TemperatureOrComposition &, \
                                                                std_cxx1x::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner); \
  template void Simulator<dim>::local_assemble_advection_system ( \
                                                                  const TemperatureOrComposition     &temperature_or_composition, \
                                                                  const std::pair<double,double> global_field_range, \
                                                                  const double                   global_max_velocity, \
                                                                  const double                   global_entropy_variation, \
                                                                  const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                  internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch, \
                                                                  internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_advection_system ( \
                                                                        const internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::assemble_advection_system (const TemperatureOrComposition     &temperature_or_composition); \
   


  ASPECT_INSTANTIATE(INSTANTIATE)
}

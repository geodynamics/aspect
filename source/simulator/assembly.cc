/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

          virtual ~StokesPreconditioner ();

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
          velocity_values (quadrature.size())
        {}



        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &scratch)
          :
          StokesPreconditioner<dim> (scratch),
          phi_u (scratch.phi_u),
          grads_phi_u (scratch.grads_phi_u),
          div_phi_u (scratch.div_phi_u),
          velocity_values (scratch.velocity_values)
        {}



        template <int dim>
        struct AdvectionSystem
        {
          AdvectionSystem (const FiniteElement<dim> &finite_element,
                           const FiniteElement<dim> &advection_element,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &quadrature,
                           const unsigned int        n_compositional_fields);
          AdvectionSystem (const AdvectionSystem &data);

          FEValues<dim>               finite_element_values;

          std::vector<types::global_dof_index>   local_dof_indices;

          /**
           * Variables describing the values and gradients of the
           * shape functions at the quadrature points, as they are
           * used in the advection assembly function. note that the sizes
           * of these arrays are equal to the number of shape functions
           * corresponding to the advected field (and not of the overall
           * field!), and that they are also correspondingly indexed.
           */
          std::vector<double>         phi_field;
          std::vector<Tensor<1,dim> > grad_phi_field;

          std::vector<Tensor<1,dim> > old_velocity_values;
          std::vector<Tensor<1,dim> > old_old_velocity_values;

          std::vector<double>         old_pressure;
          std::vector<double>         old_old_pressure;
          std::vector<Tensor<1,dim> > old_pressure_gradients;
          std::vector<Tensor<1,dim> > old_old_pressure_gradients;

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
          std::vector<Tensor<1,dim> > mesh_velocity_values;

          std::vector<SymmetricTensor<2,dim> > current_strain_rates;
          std::vector<double>         current_pressure_values;
          std::vector<Tensor<1,dim> > current_pressure_gradients;
          std::vector<std::vector<double> > current_composition_values;

          typename MaterialModel::Interface<dim>::MaterialModelInputs material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs material_model_outputs;

          typename MaterialModel::Interface<dim>::MaterialModelInputs explicit_material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs explicit_material_model_outputs;
        };



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const FiniteElement<dim> &finite_element,
                         const FiniteElement<dim> &advection_element,
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

          local_dof_indices (finite_element.dofs_per_cell),

          phi_field (advection_element.dofs_per_cell),
          grad_phi_field (advection_element.dofs_per_cell),
          old_velocity_values (quadrature.size()),
          old_old_velocity_values (quadrature.size()),
          old_pressure (quadrature.size()),
          old_old_pressure (quadrature.size()),
          old_pressure_gradients (quadrature.size()),
          old_old_pressure_gradients (quadrature.size()),
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
          mesh_velocity_values(quadrature.size()),
          current_strain_rates(quadrature.size()),
          current_pressure_values(quadrature.size()),
          current_pressure_gradients(quadrature.size()),
          current_composition_values(n_compositional_fields,
                                     std::vector<double>(quadrature.size())),
          material_model_inputs(quadrature.size(), n_compositional_fields),
          material_model_outputs(quadrature.size(), n_compositional_fields),
          explicit_material_model_inputs(quadrature.size(), n_compositional_fields),
          explicit_material_model_outputs(quadrature.size(), n_compositional_fields)
        {}



        template <int dim>
        AdvectionSystem<dim>::
        AdvectionSystem (const AdvectionSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),

          local_dof_indices (scratch.finite_element_values.get_fe().dofs_per_cell),

          phi_field (scratch.phi_field),
          grad_phi_field (scratch.grad_phi_field),
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
          old_field_grads (scratch.old_field_grads),
          old_old_field_grads (scratch.old_old_field_grads),
          old_field_laplacians (scratch.old_field_laplacians),
          old_old_field_laplacians (scratch.old_old_field_laplacians),
          old_composition_values(scratch.old_composition_values),
          old_old_composition_values(scratch.old_old_composition_values),
          current_temperature_values(scratch.current_temperature_values),
          current_velocity_values(scratch.current_velocity_values),
          mesh_velocity_values(scratch.mesh_velocity_values),
          current_strain_rates(scratch.current_strain_rates),
          current_pressure_values(scratch.current_pressure_values),
          current_pressure_gradients(scratch.current_pressure_gradients),
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

          virtual ~StokesPreconditioner ();

          FullMatrix<double>          local_matrix;
          std::vector<types::global_dof_index>   local_dof_indices;
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
        StokesPreconditioner<dim>::
        ~StokesPreconditioner ()
        {}



        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &finite_element,
                        const bool                do_pressure_rhs_compatibility_modification);
          StokesSystem (const StokesSystem<dim> &data);

          Vector<double> local_rhs;
          Vector<double> local_pressure_shape_function_integrals;
        };



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
        struct AdvectionSystem
        {
          /**
           * Constructor.
           * @param finite_element The element that describes the field for which we
           *    are trying to assemble a linear system. <b>Not</b> the global finite
           *    element.
           */
          AdvectionSystem (const FiniteElement<dim> &finite_element);
          AdvectionSystem (const AdvectionSystem &data);

          /**
           * Local contributions to the global matrix and right hand side
           * that correspond only to the variables listed in local_dof_indices
           */
          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;

          /**
           * Indices of those degrees of freedom that actually correspond
           * to the temperature or compositional field. since this structure
           * is used to represent just contributions to the advection
           * systems, there will be no contributions to other parts of the
           * system and consequently, we do not need to list here indices
           * that correspond to velocity or pressure degrees (or, in fact
           * any other variable outside the block we are currently considering)
           */
          std::vector<types::global_dof_index>   local_dof_indices;
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
                                    const AdvectionField               &advection_field,
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

        const double dField_dt = ((*scratch.old_field_values)[q] - (*scratch.old_old_field_values)[q])
                                 / old_time_step;
        const double u_grad_field = u * (scratch.old_field_grads[q] +
                                         scratch.old_old_field_grads[q]) / 2;

        const double density              = ((advection_field.is_temperature())
                                             ? scratch.explicit_material_model_outputs.densities[q] : 1.0);
        const double conductivity = ((advection_field.is_temperature()) ? scratch.explicit_material_model_outputs.thermal_conductivities[q] : 0.0);
        const double c_P                  = ((advection_field.is_temperature()) ? scratch.explicit_material_model_outputs.specific_heat[q] : 1.0);
        const double k_Delta_field = conductivity
                                     * (scratch.old_field_laplacians[q] +
                                        scratch.old_old_field_laplacians[q]) / 2;

        const double field = ((*scratch.old_field_values)[q] + (*scratch.old_old_field_values)[q]) / 2;

        const double gamma
          = compute_heating_term(scratch,
                                 scratch.explicit_material_model_inputs,
                                 scratch.explicit_material_model_outputs,
                                 advection_field,
                                 q);
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
                     const AdvectionField     &advection_field) const
  {
    if (std::abs(global_u_infty) < 1e-50
        || std::abs(global_entropy_variation) < 1e-50
        || std::abs(global_field_variation) < 1e-50)
      return 5e-3 * cell_diameter;

    double max_residual = 0;
    double max_velocity = 0;
    double max_density = 0;
    double max_specific_heat = 0;

    compute_advection_system_residual(scratch,
                                      average_field,
                                      advection_field,
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
  get_artificial_viscosity (Vector<float> &viscosity_per_cell) const
  {
    Assert(viscosity_per_cell.size()==triangulation.n_active_cells(), ExcInternalError());

    viscosity_per_cell = 0.0;

    // this function computes the artificial viscosity for the temperature
    // equation only. create an object that signifies this.
    const AdvectionField torc = AdvectionField::temperature_field;
    const std::pair<double,double>
    global_field_range = get_extrapolated_advection_field_range (torc);
    double global_entropy_variation = get_entropy_variation ((global_field_range.first +
                                                              global_field_range.second) / 2,
                                                             torc);
    double global_max_velocity = get_maximal_velocity(old_solution);


    internal::Assembly::Scratch::
    AdvectionSystem<dim> scratch (finite_element,
                                  finite_element.base_element(introspection.block_indices.temperature),
                                  mapping,
                                  QGauss<dim>(parameters.temperature_degree
                                              +
                                              (parameters.stokes_velocity_degree+1)/2),
                                  parameters.n_compositional_fields);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    for (unsigned int cellidx=0; cellidx<triangulation.n_active_cells(); ++cellidx, ++cell)
      {
        if (!cell->is_locally_owned())
          {
            viscosity_per_cell[cellidx]=-1;
            continue;
          }

        const bool use_bdf2_scheme = (timestep_number > 1);

        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

        // also have the number of dofs that correspond just to the element for
        // the system we are currently trying to assemble
        const unsigned int advection_dofs_per_cell = scratch.phi_field.size();

        Assert (advection_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
        Assert (scratch.grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
        Assert (scratch.phi_field.size() == advection_dofs_per_cell, ExcInternalError());

        const unsigned int solution_component
          = introspection.component_indices.temperature;

        const FEValuesExtractors::Scalar solution_field
          = introspection.extractors.temperature;

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


        scratch.old_field_values = &scratch.old_temperature_values;
        scratch.old_old_field_values = &scratch.old_old_temperature_values;

        scratch.finite_element_values[solution_field].get_function_gradients (old_solution,
                                                                              scratch.old_field_grads);
        scratch.finite_element_values[solution_field].get_function_gradients (old_old_solution,
                                                                              scratch.old_old_field_grads);

        scratch.finite_element_values[solution_field].get_function_laplacians (old_solution,
                                                                               scratch.old_field_laplacians);
        scratch.finite_element_values[solution_field].get_function_laplacians (old_old_solution,
                                                                               scratch.old_old_field_laplacians);

        compute_material_model_input_values (current_linearization_point,
                                             scratch.finite_element_values,
                                             true,
                                             scratch.material_model_inputs);
        material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            scratch.explicit_material_model_inputs.temperature[q] = (scratch.old_temperature_values[q] + scratch.old_old_temperature_values[q]) / 2;
            scratch.explicit_material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
            scratch.explicit_material_model_inputs.pressure[q] = (scratch.old_pressure[q] + scratch.old_old_pressure[q]) / 2;
            for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
              scratch.explicit_material_model_inputs.composition[q][c] = (scratch.old_composition_values[c][q] + scratch.old_old_composition_values[c][q]) / 2;
            scratch.explicit_material_model_inputs.strain_rate[q] = (scratch.old_strain_rates[q] + scratch.old_old_strain_rates[q]) / 2;
          }
        material_model->evaluate(scratch.explicit_material_model_inputs,scratch.explicit_material_model_outputs);

        viscosity_per_cell[cellidx] = compute_viscosity(scratch,
                                                        global_max_velocity,
                                                        global_field_range.second - global_field_range.first,
                                                        0.5 * (global_field_range.second + global_field_range.first),
                                                        global_entropy_variation,
                                                        cell->diameter(),
                                                        torc);
      }


  }



  template <int dim>
  void
  Simulator<dim>::
  compute_material_model_input_values (const LinearAlgebra::BlockVector                    &input_solution,
                                       const FEValues<dim>                                         &input_finite_element_values,
                                       const bool                                                   compute_strainrate,
                                       typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs) const
  {
    const unsigned int n_q_points = material_model_inputs.temperature.size();

    material_model_inputs.position = input_finite_element_values.get_quadrature_points();

    input_finite_element_values[introspection.extractors.temperature].get_function_values (input_solution,
        material_model_inputs.temperature);
    input_finite_element_values[introspection.extractors.pressure].get_function_values(input_solution,
        material_model_inputs.pressure);

    // only the viscosity in the material can depend on the strain_rate
    // if this is not needed, we can same some time here. By setting the
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

    material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);

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
    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
    const bool is_compressible = material_model->is_compressible();

    scratch.finite_element_values.reinit (cell);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    if (do_pressure_rhs_compatibility_modification)
      data.local_pressure_shape_function_integrals = 0;

    // we only need the strain rates for the viscosity,
    // which we only need when rebuilding the matrix
    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         rebuild_stokes_matrix,
                                         scratch.material_model_inputs);

    material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);

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
          = (is_compressible
             ?
             scratch.material_model_outputs.compressibilities[q]
             :
             std::numeric_limits<double>::quiet_NaN() );
        const double density = scratch.material_model_outputs.densities[q];

        if (rebuild_stokes_matrix)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              data.local_matrix(i,j) += ( eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j])
                                          - (is_compressible
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
                                 + (is_compressible
                                    ?
                                    (pressure_scaling *
                                     compressibility * density *
                                     (scratch.velocity_values[q] * gravity) *
                                     scratch.phi_p[i])
                                    :
                                    0)
                               )
                               * scratch.finite_element_values.JxW(q);
        if (do_pressure_rhs_compatibility_modification)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            data.local_pressure_shape_function_integrals(i) += scratch.phi_p[i] * scratch.finite_element_values.JxW(q);
      }

    //Add stabilization terms if necessary.
    if (parameters.free_surface_enabled)
      free_surface->apply_stabilization(cell, data.local_matrix);

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
         StokesSystem<dim> (finite_element,
                            do_pressure_rhs_compatibility_modification));

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    if (do_pressure_rhs_compatibility_modification)
      pressure_shape_function_integrals.compress(VectorOperation::add);


    // record that we have just rebuilt the matrix
    rebuild_stokes_matrix = false;

    computing_timer.exit_section();
  }

  template <int dim>
  void
  Simulator<dim>::build_advection_preconditioner(const AdvectionField &advection_field,
                                                 std_cxx1x::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner)
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
  double
  Simulator<dim>::compute_heating_term(const internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                       typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
                                       typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
                                       const AdvectionField     &advection_field,
                                       const unsigned int q) const
  {

    if (advection_field.field_type == AdvectionField::compositional_field)
      return 0.0;

    const double current_T = material_model_inputs.temperature[q];
    const SymmetricTensor<2,dim> current_strain_rate = material_model_inputs.strain_rate[q];
    const Tensor<1,dim> current_u = scratch.current_velocity_values[q];
    const Tensor<1,dim> current_grad_p = scratch.current_pressure_gradients[q];

    const double alpha                = material_model_outputs.thermal_expansion_coefficients[q];
    const double density              = material_model_outputs.densities[q];
    const double viscosity            = material_model_outputs.viscosities[q];
    const bool is_compressible        = material_model->is_compressible();
    const double specific_radiogenic_heating_rate = heating_model->specific_heating_rate(material_model_inputs.temperature[q],
                                                    material_model_inputs.pressure[q],
                                                    material_model_inputs.composition[q],
                                                    material_model_inputs.position[q]);
    const double compressibility      = (is_compressible
                                         ?
                                         material_model_outputs.compressibilities[q]
                                         :
                                         std::numeric_limits<double>::quiet_NaN() );
    const double entropy_gradient     = material_model_outputs.entropy_derivative_pressure[q];

    const Tensor<1,dim>
    gravity = gravity_model->gravity_vector (scratch.finite_element_values.quadrature_point(q));

    const double gamma
      = (specific_radiogenic_heating_rate * density
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
         // add the term from adiabatic compression heating
         //   + alpha T (u . nabla p)
         // where we use the definition of
         //   alpha = - 1/rho drho/dT
         // Note: this term is often simplified using the relationship
         //   rho g = -nabla p
         // to yield
         //   - alpha rho T (u . g)
         // However, we do not use this simplification here, see the
         // comment in the manual in the section on the governing
         // equations
         (parameters.include_adiabatic_heating
          ?
          (current_u * current_grad_p) * alpha * current_T
          :
          0)
         +
         // finally add the right-hand side term from latent heating
         //   DeltaS dLambda/dpi T rho (v . grad p)
         // DeltaS:      change of entropy across phase transition
         // dLambda/dpi: derivative of the phase function
         // pi:          excess pressure (argument of the phase function)
         // formulation modified after Christensen & Yuen, 1985
         (parameters.include_latent_heat
          ?
          current_T * density * entropy_gradient * (current_u * current_grad_p)
          :
          0)
        );

    return gamma;
  }


  template <int dim>
  void Simulator<dim>::
  local_assemble_advection_system (const AdvectionField     &advection_field,
                                   const std::pair<double,double> global_field_range,
                                   const double                   global_max_velocity,
                                   const double                   global_entropy_variation,
                                   const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                   internal::Assembly::CopyData::AdvectionSystem<dim> &data)
  {
    const bool use_bdf2_scheme = (timestep_number > 1);

    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    // also have the number of dofs that correspond just to the element for
    // the system we are currently trying to assemble
    const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

    Assert (advection_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
    Assert (scratch.grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
    Assert (scratch.phi_field.size() == advection_dofs_per_cell, ExcInternalError());

    const unsigned int solution_component
      = (advection_field.is_temperature()
         ?
         introspection.component_indices.temperature
         :
         introspection.component_indices.compositional_fields[advection_field.compositional_variable]
        );

    const FEValuesExtractors::Scalar solution_field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    scratch.finite_element_values.reinit (cell);

    // get all dof indices on the current cell, then extract those
    // that correspond to the solution_field we are interested in
    cell->get_dof_indices (scratch.local_dof_indices);
    for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
      data.local_dof_indices[i] = scratch.local_dof_indices[scratch.finite_element_values.get_fe().component_to_system_index(solution_component, i)];

    data.local_matrix = 0;
    data.local_rhs = 0;

    if (advection_field.is_temperature())
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

        scratch.finite_element_values[introspection.extractors.pressure].get_function_gradients (old_solution,
            scratch.old_pressure_gradients);
        scratch.finite_element_values[introspection.extractors.pressure].get_function_gradients (old_old_solution,
            scratch.old_old_pressure_gradients);
        scratch.finite_element_values[introspection.extractors.pressure].get_function_gradients (current_linearization_point,
            scratch.current_pressure_gradients);

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
        scratch.finite_element_values[introspection.extractors.compositional_fields[advection_field.compositional_variable]].get_function_values(old_solution,
            scratch.old_composition_values[advection_field.compositional_variable]);
        scratch.finite_element_values[introspection.extractors.compositional_fields[advection_field.compositional_variable]].get_function_values(old_old_solution,
            scratch.old_old_composition_values[advection_field.compositional_variable]);
      }

    scratch.finite_element_values[introspection.extractors.velocities].get_function_values (old_solution,
        scratch.old_velocity_values);
    scratch.finite_element_values[introspection.extractors.velocities].get_function_values (old_old_solution,
        scratch.old_old_velocity_values);
    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.current_velocity_values);

    //get the mesh velocity, as we need to subtract it off of the advection systems
    if (parameters.free_surface_enabled)
      scratch.finite_element_values[introspection.extractors.velocities].get_function_values(free_surface->mesh_velocity,
          scratch.mesh_velocity_values);


    scratch.old_field_values = ((advection_field.is_temperature()) ? &scratch.old_temperature_values : &scratch.old_composition_values[advection_field.compositional_variable]);
    scratch.old_old_field_values = ((advection_field.is_temperature()) ? &scratch.old_old_temperature_values : &scratch.old_old_composition_values[advection_field.compositional_variable]);

    scratch.finite_element_values[solution_field].get_function_gradients (old_solution,
                                                                          scratch.old_field_grads);
    scratch.finite_element_values[solution_field].get_function_gradients (old_old_solution,
                                                                          scratch.old_old_field_grads);

    scratch.finite_element_values[solution_field].get_function_laplacians (old_solution,
                                                                           scratch.old_field_laplacians);
    scratch.finite_element_values[solution_field].get_function_laplacians (old_old_solution,
                                                                           scratch.old_old_field_laplacians);

    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         true,
                                         scratch.material_model_inputs);
    material_model->evaluate(scratch.material_model_inputs,scratch.material_model_outputs);

    if (advection_field.is_temperature())
      {
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            scratch.explicit_material_model_inputs.temperature[q] = (scratch.old_temperature_values[q] + scratch.old_old_temperature_values[q]) / 2;
            scratch.explicit_material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
            scratch.explicit_material_model_inputs.pressure[q] = (scratch.old_pressure[q] + scratch.old_old_pressure[q]) / 2;
            for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
              scratch.explicit_material_model_inputs.composition[q][c] = (scratch.old_composition_values[c][q] + scratch.old_old_composition_values[c][q]) / 2;
            scratch.explicit_material_model_inputs.strain_rate[q] = (scratch.old_strain_rates[q] + scratch.old_old_strain_rates[q]) / 2;
          }
        material_model->evaluate(scratch.explicit_material_model_inputs,scratch.explicit_material_model_outputs);
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
                           advection_field);
    Assert (nu >= 0, ExcMessage ("The artificial viscosity needs to be a non-negative quantity."));

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
        Assert (density_c_P >= 0, ExcMessage ("The product of density and c_P needs to be a non-negative quantity."));

        const double conductivity =
          ((advection_field.is_temperature())
           ?
           scratch.material_model_outputs.thermal_conductivities[q]
           :
           0.0);
        const double latent_heat_LHS =
          ((parameters.include_latent_heat && advection_field.is_temperature())
           ?
           - scratch.material_model_outputs.densities[q] *
           scratch.material_model_inputs.temperature[q] *
           scratch.material_model_outputs.entropy_derivative_temperature[q]
           :
           0.0);
        Assert (density_c_P + latent_heat_LHS >= 0,
                ExcMessage ("The sum of density times c_P and the latent heat contribution "
                            "to the left hand side needs to be a non-negative quantity."));
        const double gamma = compute_heating_term(scratch,
                                                  scratch.material_model_inputs,
                                                  scratch.material_model_outputs,
                                                  advection_field,
                                                  q);
        const double reaction_term =
          ((advection_field.is_temperature())
           ?
           0.0
           :
           scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]);

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
            (density_c_P + latent_heat_LHS);

        Tensor<1,dim> current_u = scratch.current_velocity_values[q];
        //Subtract off the mesh velocity for ALE corrections if necessary
        if (parameters.free_surface_enabled)
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
                gamma
                * scratch.phi_field[i]
                + reaction_term
                * scratch.phi_field[i])
               *
               scratch.finite_element_values.JxW(q);

            for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
              {
                data.local_matrix(i,j)
                += (
                     (time_step * (conductivity + nu)
                      * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                     + ((time_step * (current_u * scratch.grad_phi_field[j] * scratch.phi_field[i]))
                        + (factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                     (density_c_P + latent_heat_LHS)
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
    // copy entries into the global matrix. note that these local contributions
    // only correspond to the advection dofs, as assembled above
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);
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

    const std::pair<double,double>
    global_field_range = get_extrapolated_advection_field_range (advection_field);

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
                          advection_field,
                          global_field_range,
                          get_maximal_velocity(old_solution),
                          // use the mid-value of the advected field instead of the
                          // integral mean. results are not very
                          // sensitive to this and this is far simpler
                          get_entropy_variation ((global_field_range.first +
                                                  global_field_range.second) / 2,
                                                 advection_field),
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_advection_system,
                          this,
                          std_cxx1x::_1),

         // we have to assemble the term u.grad phi_i * phi_j, which is
         // of total polynomial degree
         //   stokes_deg + 2*temp_deg -1
         // (or similar for comp_deg). this suggests using a Gauss
         // quadrature formula of order
         //   temp_deg + stokes_deg/2
         // rounded up. do so. (note that x/2 rounded up
         // equals (x+1)/2 using integer division.)
         //
         // (note: we need to get at the advection element in
         // use for the scratch and copy objects below. the
         // base element for the compositional fields exists
         // only once, with multiplicity, so only query
         // introspection.block_indices.compositional_fields[0]
         // instead of subscripting with the correct compositional
         // field index.)
         internal::Assembly::Scratch::
         AdvectionSystem<dim> (finite_element,
                               finite_element.base_element(advection_field.base_element(introspection)),
                               mapping,
                               QGauss<dim>((advection_field.is_temperature()
                                            ?
                                            parameters.temperature_degree
                                            :
                                            parameters.composition_degree)
                                           +
                                           (parameters.stokes_velocity_degree+1)/2),
                               parameters.n_compositional_fields),
         internal::Assembly::CopyData::
         AdvectionSystem<dim> (finite_element.base_element(advection_field.base_element(introspection))));

    system_matrix.compress(VectorOperation::add);
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
  template void Simulator<dim>::get_artificial_viscosity (Vector<float> &viscosity_per_cell) const; \
  template void Simulator<dim>::build_advection_preconditioner (const AdvectionField &, \
                                                                std_cxx1x::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner); \
  template void Simulator<dim>::local_assemble_advection_system ( \
                                                                  const AdvectionField          &advection_field, \
                                                                  const std::pair<double,double> global_field_range, \
                                                                  const double                   global_max_velocity, \
                                                                  const double                   global_entropy_variation, \
                                                                  const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                  internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch, \
                                                                  internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_advection_system ( \
                                                                        const internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::assemble_advection_system (const AdvectionField     &advection_field); \
   


  ASPECT_INSTANTIATE(INSTANTIATE)
}

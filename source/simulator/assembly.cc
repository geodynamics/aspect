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
          std::vector<std::vector<double>>     composition_values;

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



        // Observe that we derive the
        // StokesSystem scratch array from the
        // StokesPreconditioner array. We do this
        // because all the objects that are
        // necessary for the assembly of the
        // preconditioner are also needed for the
        // actual matrix system and right hand
        // side, plus some extra data. This makes
        // the program more compact. Note also
        // that the assembly of the Stokes system
        // and the temperature right hand side
        // further down requires data from
        // temperature and velocity,
        // respectively, so we actually need two
        // FEValues objects for those two cases.
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
          std::vector<std::vector<double>>     composition_values;

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
        struct TemperatureSystem
        {
          TemperatureSystem (const FiniteElement<dim> &finite_element,
                             const Mapping<dim>       &mapping,
                             const Quadrature<dim>    &quadrature,
                             const unsigned int        n_compositional_fields);
          TemperatureSystem (const TemperatureSystem &data);

          FEValues<dim>               finite_element_values;

          std::vector<double>         phi_T;
          std::vector<Tensor<1,dim> > grad_phi_T;

          std::vector<Tensor<1,dim> > old_velocity_values;
          std::vector<Tensor<1,dim> > old_old_velocity_values;

          std::vector<double>         old_pressure;
          std::vector<double>         old_old_pressure;

          std::vector<SymmetricTensor<2,dim> > old_strain_rates;
          std::vector<SymmetricTensor<2,dim> > old_old_strain_rates;

          std::vector<double>         old_temperature_values;
          std::vector<double>         old_old_temperature_values;
          std::vector<Tensor<1,dim> > old_temperature_grads;
          std::vector<Tensor<1,dim> > old_old_temperature_grads;
          std::vector<double>         old_temperature_laplacians;
          std::vector<double>         old_old_temperature_laplacians;

          std::vector<std::vector<double>> old_composition_values;
          std::vector<std::vector<double>> old_old_composition_values;

          std::vector<double>         current_temperature_values;
          std::vector<Tensor<1,dim> > current_velocity_values;
          std::vector<SymmetricTensor<2,dim> > current_strain_rates;
          std::vector<double>         current_pressure_values;
          std::vector<std::vector<double>> current_composition_values;


          typename MaterialModel::Interface<dim>::MaterialModelInputs material_model_inputs;
          typename MaterialModel::Interface<dim>::MaterialModelOutputs material_model_outputs;
        };



        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const FiniteElement<dim> &finite_element,
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

          phi_T (finite_element.dofs_per_cell),
          grad_phi_T (finite_element.dofs_per_cell),
          old_velocity_values (quadrature.size()),
          old_old_velocity_values (quadrature.size()),
          old_pressure (quadrature.size()),
          old_old_pressure (quadrature.size()),
          old_strain_rates (quadrature.size()),
          old_old_strain_rates (quadrature.size()),
          old_temperature_values (quadrature.size()),
          old_old_temperature_values(quadrature.size()),
          old_temperature_grads(quadrature.size()),
          old_old_temperature_grads(quadrature.size()),
          old_temperature_laplacians(quadrature.size()),
          old_old_temperature_laplacians(quadrature.size()),
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
          material_model_outputs(quadrature.size())
        {}



        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const TemperatureSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),

          phi_T (scratch.phi_T),
          grad_phi_T (scratch.grad_phi_T),
          old_velocity_values (scratch.old_velocity_values),
          old_old_velocity_values (scratch.old_old_velocity_values),
          old_pressure (scratch.old_pressure),
          old_old_pressure (scratch.old_old_pressure),
          old_strain_rates (scratch.old_strain_rates),
          old_old_strain_rates (scratch.old_old_strain_rates),
          old_temperature_values (scratch.old_temperature_values),
          old_old_temperature_values (scratch.old_old_temperature_values),
          old_temperature_grads (scratch.old_temperature_grads),
          old_old_temperature_grads (scratch.old_old_temperature_grads),
          old_temperature_laplacians (scratch.old_temperature_laplacians),
          old_old_temperature_laplacians (scratch.old_old_temperature_laplacians),
          old_composition_values(scratch.old_composition_values),
          old_old_composition_values(scratch.old_old_composition_values),
          current_temperature_values(scratch.current_temperature_values),
          current_velocity_values(scratch.current_velocity_values),
          current_strain_rates(scratch.current_strain_rates),
          current_pressure_values(scratch.current_pressure_values),
          current_composition_values(scratch.current_composition_values),
          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs)
        {}


        template <int dim>
        struct CompositionSystem
        {
          CompositionSystem (const FiniteElement<dim> &finite_element,
                             const Mapping<dim>       &mapping,
                             const Quadrature<dim>    &quadrature,
                             const unsigned int        n_compositional_fields);
          CompositionSystem (const CompositionSystem &data);

          FEValues<dim>               finite_element_values;

          std::vector<double>         phi_C;
          std::vector<Tensor<1,dim> > grad_phi_C;

          std::vector<Tensor<1,dim> > old_velocity_values;
          std::vector<Tensor<1,dim> > old_old_velocity_values;

          std::vector<double>         old_composition_values;
          std::vector<double>         old_old_composition_values;
          std::vector<Tensor<1,dim> > old_composition_grads;
          std::vector<Tensor<1,dim> > old_old_composition_grads;
          std::vector<double>         old_composition_laplacians;
          std::vector<double>         old_old_composition_laplacians;

          std::vector<Tensor<1,dim> > current_velocity_values;
        };



        template <int dim>
        CompositionSystem<dim>::
        CompositionSystem (const FiniteElement<dim> &finite_element,
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

          phi_C (finite_element.dofs_per_cell),
          grad_phi_C (finite_element.dofs_per_cell),
          old_velocity_values (quadrature.size()),
          old_old_velocity_values (quadrature.size()),
          old_composition_values(quadrature.size()),
          old_old_composition_values(quadrature.size()),
          old_composition_grads(quadrature.size()),
          old_old_composition_grads(quadrature.size()),
          old_composition_laplacians(quadrature.size()),
          old_old_composition_laplacians(quadrature.size()),
          current_velocity_values(quadrature.size())
        {}



        template <int dim>
        CompositionSystem<dim>::
        CompositionSystem (const CompositionSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),

          phi_C (scratch.phi_C),
          grad_phi_C (scratch.grad_phi_C),
          old_velocity_values (scratch.old_velocity_values),
          old_old_velocity_values (scratch.old_old_velocity_values),
          old_composition_values (scratch.old_composition_values),
          old_old_composition_values (scratch.old_old_composition_values),
          old_composition_grads (scratch.old_composition_grads),
          old_old_composition_grads (scratch.old_old_composition_grads),
          old_composition_laplacians (scratch.old_composition_laplacians),
          old_old_composition_laplacians (scratch.old_old_composition_laplacians),
          current_velocity_values(scratch.current_velocity_values)
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
        struct TemperatureSystem
        {
          TemperatureSystem (const FiniteElement<dim> &finite_element);
          TemperatureSystem (const TemperatureSystem &data);

          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;
          std::vector<unsigned int>   local_dof_indices;
        };



        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const FiniteElement<dim> &finite_element)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_rhs (finite_element.dofs_per_cell),
          local_dof_indices (finite_element.dofs_per_cell)
        {}



        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const TemperatureSystem &data)
          :
          local_matrix (data.local_matrix),
          local_rhs (data.local_rhs),
          local_dof_indices (data.local_dof_indices)
        {}



        template <int dim>
        struct CompositionSystem
        {
          CompositionSystem (const FiniteElement<dim> &finite_element);
          CompositionSystem (const CompositionSystem &data);

          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;
          std::vector<unsigned int>   local_dof_indices;
        };



        template <int dim>
        CompositionSystem<dim>::
        CompositionSystem (const FiniteElement<dim> &finite_element)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_rhs (finite_element.dofs_per_cell),
          local_dof_indices (finite_element.dofs_per_cell)
        {}



        template <int dim>
        CompositionSystem<dim>::
        CompositionSystem (const CompositionSystem &data)
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
  Simulator<dim>::get_entropy_variation (const double average_field, const unsigned int index) const
  {
    // only do this if we really need entropy
    // variation. otherwise return something that's obviously
    // nonsensical
    if (parameters.stabilization_alpha != 2)
      return std::numeric_limits<double>::quiet_NaN();

    // we should only calculate the entropy of a temperature or
    //  composition field
    AssertIndexRange(index,parameters.n_compositional_fields+1);

    // record maximal entropy on Gauss quadrature
    // points
    const QGauss<dim> quadrature_formula (parameters.temperature_degree+1);
    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar field (dim+1+index);

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
  compute_temperature_system_residual(const std::vector<double>          &old_temperature,
                                      const std::vector<double>          &old_old_temperature,
                                      const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                                      const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                                      const std::vector<double>          &old_temperature_laplacians,
                                      const std::vector<double>          &old_old_temperature_laplacians,
                                      const std::vector<Tensor<1,dim> >  &old_velocity_values,
                                      const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                                      const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                                      const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                                      const std::vector<double>          &old_pressure,
                                      const std::vector<double>          &old_old_pressure,
                                      const std::vector<std::vector<double>> &old_composition,
                                      const std::vector<std::vector<double>> &old_old_composition,
                                      const double                        average_temperature,
                                      const std::vector<Point<dim> >     &evaluation_points,
                                      double                             &max_residual,
                                      double                             &max_velocity) const
    {

    const unsigned int n_q_points = old_temperature.size();

    typename MaterialModel::Interface<dim>::MaterialModelInputs material_model_inputs (n_q_points, parameters.n_compositional_fields);
    typename MaterialModel::Interface<dim>::MaterialModelOutputs material_model_outputs (n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q) {
      material_model_inputs.temperature[q] = (old_temperature[q] + old_old_temperature[q]) / 2;
      material_model_inputs.position[q] = evaluation_points[q];
      material_model_inputs.pressure[q] = (old_pressure[q] + old_old_pressure[q]) / 2;
      for (unsigned int i=0; i<parameters.n_compositional_fields; ++i)
        material_model_inputs.composition[i][q] = (old_composition[i][q] + old_old_composition[i][q]) / 2;
      material_model_inputs.strain_rate[q] = (old_strain_rates[q] + old_old_strain_rates[q]) / 2;
    }

    material_model->compute_parameters(material_model_inputs,material_model_outputs);

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (old_velocity_values[q] +
                                 old_old_velocity_values[q]) / 2;

        const SymmetricTensor<2,dim> strain_rate = material_model_inputs.strain_rate[q];

        const double dT_dt = (old_temperature[q] - old_old_temperature[q])
                             / old_time_step;
        const double u_grad_T = u * (old_temperature_grads[q] +
                                     old_old_temperature_grads[q]) / 2;

        const double alpha                = material_model_outputs.thermal_expansion_coefficients[q];
        const double density              = material_model_outputs.densities[q];
        const double thermal_conductivity = material_model_outputs.thermal_conductivities[q];
        const double c_P                  = material_model_outputs.specific_heat[q];
        const double k_Delta_T = thermal_conductivity
                                 * (old_temperature_laplacians[q] +
                                    old_old_temperature_laplacians[q]) / 2;

        // verify correctness of the heating term
        const double viscosity =  material_model_outputs.viscosities[q];
        const bool is_compressible = material_model_outputs.is_compressible;
        const double compressibility
          = (is_compressible
             ?
             material_model_outputs.compressibilities[q]
             :
             std::numeric_limits<double>::quiet_NaN() );
        const Tensor<1,dim> gravity = gravity_model->gravity_vector (evaluation_points[q] );
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
              strain_rate * strain_rate
              -
              (is_compressible
               ?
               2./3.*viscosity*std::pow(compressibility * density * (u * gravity),
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
                alpha * density * (u*gravity) * material_model_inputs.temperature[q]
                :
                0)
              );
        double residual
          = std::abs(density * c_P * (dT_dt + u_grad_T) - k_Delta_T - gamma);
        if (parameters.stabilization_alpha == 2)
          residual *= std::abs(material_model_inputs.temperature[q] - average_temperature);

        max_residual = std::max (residual,        max_residual);
        max_velocity = std::max (std::sqrt (u*u), max_velocity);
      }
    }

  template <int dim>
  void
  Simulator<dim>::
  compute_composition_system_residual(const std::vector<double>          &old_composition,
                                      const std::vector<double>          &old_old_composition,
                                      const std::vector<Tensor<1,dim> >  &old_composition_grads,
                                      const std::vector<Tensor<1,dim> >  &old_old_composition_grads,
                                      const std::vector<double>          &old_composition_laplacians,
                                      const std::vector<double>          &old_old_composition_laplacians,
                                      const std::vector<Tensor<1,dim> >  &old_velocity_values,
                                      const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                                      const double                        average_composition,
                                      const unsigned int                  composition_index,
                                      double                             &max_residual,
                                      double                             &max_velocity) const
    {

    const unsigned int n_q_points = old_composition.size();

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (old_velocity_values[q] +
                                 old_old_velocity_values[q]) / 2;
        const double C = (old_composition[q] + old_old_composition[q]) / 2;

        const double dC_dt = (old_composition[q] - old_old_composition[q])
                             / old_time_step;
        const double u_grad_C = u * (old_composition_grads[q] +
                                     old_old_composition_grads[q]) / 2;

        const double kappa = parameters.chemical_diffusivities[composition_index];

        const double kappa_Delta_C = kappa
                                 * (old_composition_laplacians[q] +
                                    old_old_composition_laplacians[q]) / 2;

        double residual
          = std::abs(dC_dt + u_grad_C - kappa_Delta_C);
        if (parameters.stabilization_alpha == 2)
          residual *= std::abs(C - average_composition);

        max_residual = std::max (residual,        max_residual);
        max_velocity = std::max (std::sqrt (u*u), max_velocity);
      }
    }



  template <int dim>
  double
  Simulator<dim>::
  compute_viscosity (const std::vector<double>          &old_field,
                     const std::vector<double>          &old_old_field,
                     const std::vector<Tensor<1,dim> >  &old_field_grads,
                     const std::vector<Tensor<1,dim> >  &old_old_field_grads,
                     const std::vector<double>          &old_field_laplacians,
                     const std::vector<double>          &old_old_field_laplacians,
                     const std::vector<Tensor<1,dim> >  &old_velocity_values,
                     const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                     const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                     const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                     const std::vector<double>          &old_pressure,
                     const std::vector<double>          &old_old_pressure,
                     const std::vector<std::vector<double>> &old_composition,
                     const std::vector<std::vector<double>> &old_old_composition,
                     const double                        global_u_infty,
                     const double                        global_T_variation,
                     const double                        average_field,
                     const double                        global_entropy_variation,
                     const std::vector<Point<dim> >     &evaluation_points,
                     const double                        cell_diameter,
                     const unsigned int                  index) const
  {
    if (global_u_infty == 0)
      return 5e-3 * cell_diameter;

    double max_residual = 0;
    double max_velocity = 0;

    AssertIndexRange(index,parameters.n_compositional_fields+1);

    if(index == 0) {
      //make sure that all arguments we need for computing the residual are passed
      Assert (old_strain_rates.size() > 0 && old_old_strain_rates.size() > 0
              && old_pressure.size() > 0 && old_old_pressure.size() > 0,
              ExcMessage ("Not enough parameters to calculate artificial viscosity "
                          "for the temperature equation."));

      compute_temperature_system_residual(old_field,
                                          old_old_field,
                                          old_field_grads,
                                          old_old_field_grads,
                                          old_field_laplacians,
                                          old_old_field_laplacians,
                                          old_velocity_values,
                                          old_old_velocity_values,
                                          old_strain_rates,
                                          old_old_strain_rates,
                                          old_pressure,
                                          old_old_pressure,
                                          old_composition,
                                          old_old_composition,
                                          average_field,
                                          evaluation_points,
                                          max_residual,
                                          max_velocity);
    }
    else
      compute_composition_system_residual(old_field,
                                          old_old_field,
                                          old_field_grads,
                                          old_old_field_grads,
                                          old_field_laplacians,
                                          old_old_field_laplacians,
                                          old_velocity_values,
                                          old_old_velocity_values,
                                          average_field,
                                          index-1,   //index of compositional field
                                          max_residual,
                                          max_velocity);

    const double max_viscosity = (parameters.stabilization_beta *
                                  max_velocity * cell_diameter);
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
                               (global_u_infty * global_T_variation));

        return std::min (max_viscosity, entropy_viscosity);
      }
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

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);
    std::vector<FEValuesExtractors::Scalar> compositional_fields;

    for (unsigned int q=0;q<parameters.n_compositional_fields;++q)
      {
      const FEValuesExtractors::Scalar temp(dim+2+q);
      compositional_fields.push_back(temp);
      }

    scratch.finite_element_values.reinit (cell);

    scratch.finite_element_values[temperature].get_function_values (current_linearization_point,
                                                                    scratch.temperature_values);
    scratch.finite_element_values[pressure].get_function_values(current_linearization_point,
                                                                scratch.pressure_values);
    scratch.finite_element_values[velocities].get_function_symmetric_gradients(current_linearization_point,
                                                                               scratch.strain_rates);
    for(unsigned int q=0;q<parameters.n_compositional_fields;++q)
      scratch.finite_element_values[compositional_fields[q]].get_function_values(current_linearization_point,
                                                                            scratch.composition_values[q]);

    data.local_matrix = 0;

    scratch.material_model_inputs.temperature = scratch.temperature_values;
    for (unsigned int q=0; q<n_q_points; ++q)
      scratch.material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
    scratch.material_model_inputs.pressure = scratch.pressure_values;
    for (unsigned int i=0; i<parameters.n_compositional_fields; ++i)
      scratch.material_model_inputs.composition[i] = scratch.composition_values[i];
    scratch.material_model_inputs.strain_rate = scratch.strain_rates;

    material_model->compute_parameters(scratch.material_model_inputs,scratch.material_model_outputs);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        const double current_temperature = scratch.temperature_values[q];
        const double current_pressure = scratch.pressure_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grads_phi_u[k] = scratch.finite_element_values[velocities].symmetric_gradient(k,q);
            scratch.phi_p[k]       = scratch.finite_element_values[pressure].value (k, q);
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
    std::vector<bool>  velocity_components (dim+2+parameters.n_compositional_fields,true);
    for (unsigned int i=dim; i<dim+2+parameters.n_compositional_fields; ++i)
      velocity_components[i] = false;
    DoFTools::extract_constant_modes (dof_handler, velocity_components,
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
    const bool use_bdf2_scheme       = (timestep_number > 1);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);
    std::vector<FEValuesExtractors::Scalar> compositional_fields;

    for (unsigned int q=0;q<parameters.n_compositional_fields;++q)
      {
      const FEValuesExtractors::Scalar temp(dim+2+q);
      compositional_fields.push_back(temp);
      }

    scratch.finite_element_values.reinit (cell);
    //scratch.finite_element_values[temperature].get_function_values (old_solution,
    //                                              scratch.old_temperature_values);
    // Assuming we already have the temperature for the current time step:
    scratch.finite_element_values[temperature].get_function_values (current_linearization_point,
                                                                    scratch.temperature_values);
    scratch.finite_element_values[pressure].get_function_values(current_linearization_point,
                                                                scratch.pressure_values);
    scratch.finite_element_values[velocities].get_function_values(current_linearization_point,
                                                                  scratch.velocity_values);
    for(unsigned int q=0;q<parameters.n_compositional_fields;++q)
      scratch.finite_element_values[compositional_fields[q]].get_function_values(current_linearization_point,
                                                                            scratch.composition_values[q]);

    // we only need the strain rates for the viscosity,
    // which we only need when rebuilding the matrix
    if (rebuild_stokes_matrix)
      scratch.finite_element_values[velocities]
      .get_function_symmetric_gradients(current_linearization_point,
                                        scratch.strain_rates);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    if (material_model->is_compressible())
      data.local_pressure_shape_function_integrals = 0;

    scratch.material_model_inputs.temperature = scratch.temperature_values;
    for (unsigned int q=0; q<n_q_points; ++q)
      scratch.material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
    scratch.material_model_inputs.pressure = scratch.pressure_values;
    for (unsigned int i=0; i<parameters.n_compositional_fields; ++i)
      scratch.material_model_inputs.composition[i] = scratch.composition_values[i];
    scratch.material_model_inputs.strain_rate = scratch.strain_rates;

    material_model->compute_parameters(scratch.material_model_inputs,scratch.material_model_outputs);


    for (unsigned int q=0; q<n_q_points; ++q)
      {
        const double current_temperature = scratch.temperature_values[q];
        const double current_pressure    = scratch.pressure_values[q];
        const Tensor<1,dim> current_u    = scratch.velocity_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_u[k] = scratch.finite_element_values[velocities].value (k,q);
            scratch.phi_p[k] = scratch.finite_element_values[pressure].value (k, q);
            if (rebuild_stokes_matrix)
              {
                scratch.grads_phi_u[k] = scratch.finite_element_values[velocities].symmetric_gradient(k,q);
                scratch.div_phi_u[k]   = scratch.finite_element_values[velocities].divergence (k, q);
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
    system_rhs.compress(Add);

    if (material_model->is_compressible())
      pressure_shape_function_integrals.compress(Add);

    rebuild_stokes_matrix = false;

    computing_timer.exit_section();
  }



  template <int dim>
  void
  Simulator<dim>::build_temperature_preconditioner ()
  {
    computing_timer.enter_section ("   Build temperature preconditioner");
    {
      T_preconditioner.reset (new TrilinosWrappers::PreconditionILU());
      T_preconditioner->initialize (system_matrix.block(2,2));
    }
    computing_timer.exit_section();
  }



  template <int dim>
  void Simulator<dim>::
  local_assemble_temperature_system (const std::pair<double,double> global_T_range,
                                     const double                   global_max_velocity,
                                     const double                   global_entropy_variation,
                                     const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     internal::Assembly::Scratch::TemperatureSystem<dim> &scratch,
                                     internal::Assembly::CopyData::TemperatureSystem<dim> &data)
  {
    const bool use_bdf2_scheme = (timestep_number > 1);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);
    std::vector<FEValuesExtractors::Scalar> compositional_fields;

    for (unsigned int q=0;q<parameters.n_compositional_fields;++q)
      {
      const FEValuesExtractors::Scalar temp(dim+2+q);
      compositional_fields.push_back(temp);
      }

    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    scratch.finite_element_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_matrix = 0;
    data.local_rhs = 0;

    scratch.finite_element_values[temperature].get_function_values (old_solution,
                                                                    scratch.old_temperature_values);
    scratch.finite_element_values[temperature].get_function_values (old_old_solution,
                                                                    scratch.old_old_temperature_values);

    scratch.finite_element_values[temperature].get_function_gradients (old_solution,
                                                                       scratch.old_temperature_grads);
    scratch.finite_element_values[temperature].get_function_gradients (old_old_solution,
                                                                       scratch.old_old_temperature_grads);

    scratch.finite_element_values[temperature].get_function_laplacians (old_solution,
                                                                        scratch.old_temperature_laplacians);
    scratch.finite_element_values[temperature].get_function_laplacians (old_old_solution,
                                                                        scratch.old_old_temperature_laplacians);

    scratch.finite_element_values[velocities].get_function_values (old_solution,
                                                                   scratch.old_velocity_values);
    scratch.finite_element_values[velocities].get_function_values (old_old_solution,
                                                                   scratch.old_old_velocity_values);
    scratch.finite_element_values[velocities].get_function_symmetric_gradients (old_solution,
                                                                                scratch.old_strain_rates);
    scratch.finite_element_values[velocities].get_function_symmetric_gradients (old_old_solution,
                                                                                scratch.old_old_strain_rates);

    scratch.finite_element_values[pressure].get_function_values (old_solution,
                                                                 scratch.old_pressure);
    scratch.finite_element_values[pressure].get_function_values (old_old_solution,
                                                                 scratch.old_old_pressure);

    scratch.finite_element_values[temperature].get_function_values(current_linearization_point,
                                                                   scratch.current_temperature_values);
    scratch.finite_element_values[velocities].get_function_values(current_linearization_point,
                                                                  scratch.current_velocity_values);
    scratch.finite_element_values[velocities].get_function_symmetric_gradients(current_linearization_point,
                                                                               scratch.current_strain_rates);
    scratch.finite_element_values[pressure].get_function_values(current_linearization_point,
                                                                scratch.current_pressure_values);

    for(unsigned int q=0;q<parameters.n_compositional_fields;++q) {
      scratch.finite_element_values[compositional_fields[q]].get_function_values(current_linearization_point,
                                                                            scratch.current_composition_values[q]);
      scratch.finite_element_values[compositional_fields[q]].get_function_values(old_solution,
                                                                                  scratch.old_composition_values[q]);
      scratch.finite_element_values[compositional_fields[q]].get_function_values(old_old_solution,
                                                                                  scratch.old_old_composition_values[q]);
    }


    // TODO: Compute artificial viscosity once per timestep instead of each time
    // temperature system is assembled (as this might happen more than once per
    // timestep for iterative solvers)
    const double nu
      = compute_viscosity (scratch.old_temperature_values,
                           scratch.old_old_temperature_values,
                           scratch.old_temperature_grads,
                           scratch.old_old_temperature_grads,
                           scratch.old_temperature_laplacians,
                           scratch.old_old_temperature_laplacians,
                           scratch.old_velocity_values,
                           scratch.old_old_velocity_values,
                           scratch.old_strain_rates,
                           scratch.old_old_strain_rates,
                           scratch.old_pressure,
                           scratch.old_old_pressure,
                           scratch.old_composition_values,
                           scratch.old_old_composition_values,
                           global_max_velocity,
                           global_T_range.second - global_T_range.first,
                           0.5 * (global_T_range.second + global_T_range.first),
                           global_entropy_variation,
                           scratch.finite_element_values.get_quadrature_points(),
                           cell->diameter(),
                           0);  //index for temperature field

    scratch.material_model_inputs.temperature = scratch.current_temperature_values;
    for (unsigned int q=0; q<n_q_points; ++q)
      scratch.material_model_inputs.position[q] = scratch.finite_element_values.quadrature_point(q);
    scratch.material_model_inputs.pressure = scratch.current_pressure_values;
    for (unsigned int i=0; i<parameters.n_compositional_fields; ++i)
      scratch.material_model_inputs.composition[i] = scratch.current_composition_values[i];
    scratch.material_model_inputs.strain_rate = scratch.current_strain_rates;

    material_model->compute_parameters(scratch.material_model_inputs,scratch.material_model_outputs);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_T[k] = scratch.finite_element_values[temperature].gradient (k,q);
            scratch.phi_T[k]      = scratch.finite_element_values[temperature].value (k, q);
          }

        const double T_term_for_rhs
          = (use_bdf2_scheme ?
             (scratch.old_temperature_values[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_temperature_values[q] *
              (time_step * time_step) /
              (old_time_step * (time_step + old_time_step)))
             :
             scratch.old_temperature_values[q]);

        const double current_T = scratch.current_temperature_values[q];
        const SymmetricTensor<2,dim> current_strain_rate = scratch.current_strain_rates[q];
        const Tensor<1,dim> current_u = scratch.current_velocity_values[q];
        const double current_p = scratch.current_pressure_values[q];

        const double alpha                = scratch.material_model_outputs.thermal_expansion_coefficients[q];
        const double density              = scratch.material_model_outputs.densities[q];
        const double thermal_conductivity = scratch.material_model_outputs.thermal_conductivities[q];
        const double c_P                  = scratch.material_model_outputs.specific_heat[q];
        const double viscosity            = scratch.material_model_outputs.viscosities[q];
        const bool is_compressible        = scratch.material_model_outputs.is_compressible;
        const double compressibility      = (is_compressible
                                            ?
                                            scratch.material_model_outputs.compressibilities[q]
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

        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            data.local_rhs(i) += (T_term_for_rhs * density * c_P * scratch.phi_T[i]
                                  +
                                  time_step *
                                  gamma * scratch.phi_T[i])
                                 *
                                 scratch.finite_element_values.JxW(q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                          (time_step + old_time_step)) : 1.0;
                data.local_matrix(i,j)
                += (
                     (time_step * (thermal_conductivity +
                                   nu * density * c_P) * scratch.grad_phi_T[i] * scratch.grad_phi_T[j])
                     + ((time_step * (current_u * scratch.grad_phi_T[j] * scratch.phi_T[i]))
                        + (factor * scratch.phi_T[i] * scratch.phi_T[j])) * density * c_P
                   )
                   * scratch.finite_element_values.JxW(q);
              }
          }
      }
  }

  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_temperature_system (const internal::Assembly::CopyData::TemperatureSystem<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs
                                                   );

  }


  template <int dim>
  void Simulator<dim>::assemble_temperature_system ()
  {
    computing_timer.enter_section ("   Assemble temperature system");

    // Reset only temperature block (might reuse Stokes block)
    system_matrix.block(2,2) = 0;
    system_rhs = 0;

    const std::pair<double,double>
    global_T_range = get_extrapolated_temperature_or_composition_range (0);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_temperature_system,
                          this,
                          global_T_range,
                          get_maximal_velocity(old_solution),
                          // use the mid temperature instead of the
                          // integral mean. results are not very
                          // sensitive to this and this is far simpler
                          get_entropy_variation ((global_T_range.first +
                                                  global_T_range.second) / 2, 0),
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_temperature_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         TemperatureSystem<dim> (finite_element, mapping, QGauss<dim>(parameters.temperature_degree+2),
                                 parameters.n_compositional_fields),
         internal::Assembly::CopyData::
         TemperatureSystem<dim> (finite_element));

    system_matrix.compress();
    system_rhs.compress(Add);

    computing_timer.exit_section();
  }


  template <int dim>
  void
  Simulator<dim>::build_composition_preconditioner (unsigned int composition_index)
  {
    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(composition_index,parameters.n_compositional_fields);

    computing_timer.enter_section ("   Build composition preconditioner");
    C_preconditioner.reset (new TrilinosWrappers::PreconditionILU());
    C_preconditioner->initialize (system_matrix.block(3+composition_index,3+composition_index));

    computing_timer.exit_section();
  }


  template <int dim>
  void Simulator<dim>::
  local_assemble_composition_system (const unsigned int             composition_index,
                                     const std::pair<double,double> global_C_range,
                                     const double                   global_max_velocity,
                                     const double                   global_entropy_variation,
                                     const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     internal::Assembly::Scratch::CompositionSystem<dim> &scratch,
                                     internal::Assembly::CopyData::CompositionSystem<dim> &data)
  {
    const bool use_bdf2_scheme = (timestep_number > 1);

    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(composition_index,parameters.n_compositional_fields);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);
    const FEValuesExtractors::Scalar composition (dim+2+composition_index);

    const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

    scratch.finite_element_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_matrix = 0;
    data.local_rhs = 0;

    scratch.finite_element_values[composition].get_function_values (old_solution,
                                                                    scratch.old_composition_values);
    scratch.finite_element_values[composition].get_function_values (old_old_solution,
                                                                    scratch.old_old_composition_values);

    scratch.finite_element_values[composition].get_function_gradients (old_solution,
                                                                       scratch.old_composition_grads);
    scratch.finite_element_values[composition].get_function_gradients (old_old_solution,
                                                                       scratch.old_old_composition_grads);

    scratch.finite_element_values[composition].get_function_laplacians (old_solution,
                                                                        scratch.old_composition_laplacians);
    scratch.finite_element_values[composition].get_function_laplacians (old_old_solution,
                                                                        scratch.old_old_composition_laplacians);

    scratch.finite_element_values[velocities].get_function_values (old_solution,
                                                                   scratch.old_velocity_values);
    scratch.finite_element_values[velocities].get_function_values (old_old_solution,
                                                                   scratch.old_old_velocity_values);

    scratch.finite_element_values[velocities].get_function_values(current_linearization_point,
                                                                  scratch.current_velocity_values);

    // TODO: Compute artificial viscosity once per timestep instead of each time
    // temperature system is assembled (as this might happen more than once per
    // timestep for iterative solvers)
    const double nu
      = compute_viscosity (scratch.old_composition_values,
                           scratch.old_old_composition_values,
                           scratch.old_composition_grads,
                           scratch.old_old_composition_grads,
                           scratch.old_composition_laplacians,
                           scratch.old_old_composition_laplacians,
                           scratch.old_velocity_values,
                           scratch.old_old_velocity_values,
                           std::vector<SymmetricTensor<2,dim> >(), //we do not need the strain rate for calculating the artificial viscosity
                           std::vector<SymmetricTensor<2,dim> >(), //strain rate
                           std::vector<double>(),                  //we also do not need the pressure
                           std::vector<double>(),                  //pressure
                           std::vector<std::vector<double>> (),    //we also do not need the other compositional fields
                           std::vector<std::vector<double>> (),    //other compositional fields
                           global_max_velocity,
                           global_C_range.second - global_C_range.first,
                           0.5 * (global_C_range.second + global_C_range.first),
                           global_entropy_variation,
                           std::vector<Point<dim, double> >(),     //we also do not need the position for calculating the artificial viscosity
                           cell->diameter(),
                           composition_index+1);                   //index for the compositional field (0 is temperature)

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_C[k] = scratch.finite_element_values[composition].gradient (k,q);
            scratch.phi_C[k]      = scratch.finite_element_values[composition].value (k, q);
          }

        const double C_term_for_rhs
          = (use_bdf2_scheme ?
             (scratch.old_composition_values[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_composition_values[q] *
              (time_step * time_step) /
              (old_time_step * (time_step + old_time_step)))
             :
             scratch.old_composition_values[q]);

        const Tensor<1,dim> current_u = scratch.current_velocity_values[q];

        const double kappa = parameters.chemical_diffusivities[composition_index];

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            data.local_rhs(i) += C_term_for_rhs * scratch.phi_C[i]
                                 *
                                 scratch.finite_element_values.JxW(q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                          (time_step + old_time_step)) : 1.0;
                data.local_matrix(i,j)
                += (
                     time_step * (kappa + nu) * scratch.grad_phi_C[i] * scratch.grad_phi_C[j]
                     + time_step * current_u * scratch.grad_phi_C[j] * scratch.phi_C[i]
                     + factor * scratch.phi_C[i] * scratch.phi_C[j]
                   )
                   * scratch.finite_element_values.JxW(q);
              }
          }
      }
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_composition_system (const internal::Assembly::CopyData::CompositionSystem<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs
                                                   );

  }


  template <int dim>
  void Simulator<dim>::assemble_composition_system (unsigned int composition_index)
  {
    computing_timer.enter_section ("   Assemble composition system");

    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(composition_index,parameters.n_compositional_fields);

    // Reset only composition block (might reuse Stokes block)
    system_matrix.block(3+composition_index,3+composition_index) = 0;
    system_rhs = 0;

    const std::pair<double,double>
    global_C_range = get_extrapolated_temperature_or_composition_range (1+composition_index);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_composition_system,
                          this,
                          composition_index,
                          global_C_range,
                          get_maximal_velocity(old_solution),
                          // use the mid-value of the composition instead of the
                          // integral mean. results are not very
                          // sensitive to this and this is far simpler
                          get_entropy_variation ((global_C_range.first +
                                                  global_C_range.second) / 2, 1+composition_index),
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_composition_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         CompositionSystem<dim> (finite_element, mapping, QGauss<dim>(parameters.composition_degree+2),
                                 parameters.n_compositional_fields),
         internal::Assembly::CopyData::
         CompositionSystem<dim> (finite_element));

    system_matrix.compress();
    system_rhs.compress(Add);

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
  template void Simulator<dim>::copy_local_to_global_temperature_system ( \
                                                                          const internal::Assembly::CopyData::TemperatureSystem<dim> &data); \
  template void Simulator<dim>::build_temperature_preconditioner (); \
  template void Simulator<dim>::assemble_temperature_system (); \
  template void Simulator<dim>::copy_local_to_global_composition_system ( \
                                                                          const internal::Assembly::CopyData::CompositionSystem<dim> &data); \
  template void Simulator<dim>::build_composition_preconditioner (unsigned int composition_index); \
  template void Simulator<dim>::assemble_composition_system (unsigned int composition_index); \
  template void Simulator<dim>::local_assemble_temperature_system ( \
                                                                    const std::pair<double,double> global_T_range, \
                                                                    const double                   global_max_velocity, \
                                                                    const double                   global_entropy_variation, \
                                                                    const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                    internal::Assembly::Scratch::TemperatureSystem<dim>  &scratch, \
                                                                    internal::Assembly::CopyData::TemperatureSystem<dim> &data); \
  template void Simulator<dim>::local_assemble_composition_system ( \
                                                                    const unsigned int             composition_index, \
                                                                    const std::pair<double,double> global_C_range, \
                                                                    const double                   global_max_velocity, \
                                                                    const double                   global_entropy_variation, \
                                                                    const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                    internal::Assembly::Scratch::CompositionSystem<dim>  &scratch, \
                                                                    internal::Assembly::CopyData::CompositionSystem<dim> &data);

  ASPECT_INSTANTIATE(INSTANTIATE)
}

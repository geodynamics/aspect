/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University,
     Timo Heister, University of Goettingen, 2008-2011 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011, 2012 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

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
          StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
                                const Quadrature<dim>    &stokes_quadrature,
                                const Mapping<dim>       &mapping,
                                const UpdateFlags         stokes_update_flags,
                                const FiniteElement<dim> &temperature_fe,
                                const UpdateFlags         temperature_update_flags);
          StokesPreconditioner (const StokesPreconditioner &data);

          FEValues<dim>               stokes_fe_values;
          FEValues<dim>               temperature_fe_values;

          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  phi_p;

          std::vector<double>                  old_temperature_values;
          std::vector<double>                  old_pressure_values;
        };

        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
                              const Quadrature<dim>    &stokes_quadrature,
                              const Mapping<dim>       &mapping,
                              const UpdateFlags         stokes_update_flags,
                              const FiniteElement<dim> &temperature_fe,
                              const UpdateFlags         temperature_update_flags)
          :
          stokes_fe_values (mapping, stokes_fe, stokes_quadrature,
                            stokes_update_flags),
          temperature_fe_values (mapping,
                                 temperature_fe, stokes_quadrature,
                                 temperature_update_flags),
          grads_phi_u (stokes_fe.dofs_per_cell),
          phi_p (stokes_fe.dofs_per_cell),
          old_temperature_values (stokes_quadrature.size()),
          old_pressure_values (stokes_quadrature.size())
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &scratch)
          :
          stokes_fe_values (scratch.stokes_fe_values.get_mapping(),
                            scratch.stokes_fe_values.get_fe(),
                            scratch.stokes_fe_values.get_quadrature(),
                            scratch.stokes_fe_values.get_update_flags()),
          temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                                 scratch.temperature_fe_values.get_fe(),
                                 scratch.temperature_fe_values.get_quadrature(),
                                 scratch.temperature_fe_values.get_update_flags()),
          grads_phi_u (scratch.grads_phi_u),
          phi_p (scratch.phi_p),
          old_temperature_values (scratch.old_temperature_values),
          old_pressure_values (scratch.old_pressure_values)
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
          StokesSystem (const FiniteElement<dim> &stokes_fe,
                        const Mapping<dim>       &mapping,
                        const Quadrature<dim>    &stokes_quadrature,
                        const UpdateFlags         stokes_update_flags,
                        const FiniteElement<dim> &temperature_fe,
                        const UpdateFlags         temperature_update_flags);

          StokesSystem (const StokesSystem<dim> &data);

          std::vector<Tensor<1,dim> >          phi_u;
          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  div_phi_u;
          std::vector<Tensor<1,dim> > old_velocity_values;
        };


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &stokes_fe,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &stokes_quadrature,
                      const UpdateFlags         stokes_update_flags,
                      const FiniteElement<dim> &temperature_fe,
                      const UpdateFlags         temperature_update_flags)
          :
          StokesPreconditioner<dim> (stokes_fe, stokes_quadrature,
                                     mapping,
                                     stokes_update_flags,
                                     temperature_fe,
                                     temperature_update_flags),
          phi_u (stokes_fe.dofs_per_cell),
          grads_phi_u (stokes_fe.dofs_per_cell),
          div_phi_u (stokes_fe.dofs_per_cell),
          old_velocity_values (stokes_quadrature.size())
        {}


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &scratch)
          :
          StokesPreconditioner<dim> (scratch),
          phi_u (scratch.phi_u),
          grads_phi_u (scratch.grads_phi_u),
          div_phi_u (scratch.div_phi_u),
          old_velocity_values (scratch.old_velocity_values)
        {}



        template <int dim>
        struct TemperatureSystem
        {
          TemperatureSystem (const FiniteElement<dim> &temperature_fe,
                             const FiniteElement<dim> &stokes_fe,
                             const Mapping<dim>       &mapping,
                             const Quadrature<dim>    &quadrature);
          TemperatureSystem (const TemperatureSystem &data);

          FEValues<dim>               temperature_fe_values;
          FEValues<dim>               stokes_fe_values;

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
        };

        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const FiniteElement<dim> &temperature_fe,
                           const FiniteElement<dim> &stokes_fe,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &quadrature)
          :
          temperature_fe_values (mapping,
                                 temperature_fe, quadrature,
                                 update_values    |
                                 update_gradients |
                                 update_hessians  |
                                 update_quadrature_points |
                                 update_JxW_values),
          stokes_fe_values(mapping, stokes_fe, quadrature,
                           update_values    | update_gradients |
                           update_JxW_values),
          phi_T (temperature_fe.dofs_per_cell),
          grad_phi_T (temperature_fe.dofs_per_cell),
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
          old_old_temperature_laplacians(quadrature.size())
        {}


        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const TemperatureSystem &scratch)
          :
          temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                                 scratch.temperature_fe_values.get_fe(),
                                 scratch.temperature_fe_values.get_quadrature(),
                                 scratch.temperature_fe_values.get_update_flags()),
          stokes_fe_values (scratch.stokes_fe_values.get_mapping(),
                            scratch.stokes_fe_values.get_fe(),
                            scratch.stokes_fe_values.get_quadrature(),
                            scratch.stokes_fe_values.get_update_flags()),
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
          old_old_temperature_laplacians (scratch.old_old_temperature_laplacians)
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
          StokesPreconditioner (const FiniteElement<dim> &stokes_fe);
          StokesPreconditioner (const StokesPreconditioner &data);

          FullMatrix<double>          local_matrix;
          std::vector<unsigned int>   local_dof_indices;
        };

        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const FiniteElement<dim> &stokes_fe)
          :
          local_matrix (stokes_fe.dofs_per_cell,
                        stokes_fe.dofs_per_cell),
          local_dof_indices (stokes_fe.dofs_per_cell)
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
          StokesSystem (const FiniteElement<dim> &stokes_fe);
          StokesSystem (const StokesSystem<dim> &data);

          Vector<double> local_rhs;
          Vector<double> local_pressure_shape_function_integrals;
        };


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &stokes_fe)
          :
          StokesPreconditioner<dim> (stokes_fe),
          local_rhs (stokes_fe.dofs_per_cell),
          local_pressure_shape_function_integrals (stokes_fe.dofs_per_cell)
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
          TemperatureSystem (const FiniteElement<dim> &temperature_fe);
          TemperatureSystem (const TemperatureSystem &data);

          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;
          std::vector<unsigned int>   local_dof_indices;
        };

        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const FiniteElement<dim> &temperature_fe)
          :
          local_matrix (temperature_fe.dofs_per_cell,
                        temperature_fe.dofs_per_cell),
          local_rhs (temperature_fe.dofs_per_cell),
          local_dof_indices (temperature_fe.dofs_per_cell)
        {}


        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const TemperatureSystem &data)
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
   * artificial viscosity used to stabilize the temperature equation.
   */
  template <int dim>
  double
  Simulator<dim>::get_entropy_variation (const double average_temperature) const
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

    FEValues<dim> fe_values (temperature_fe, quadrature_formula,
                             update_values | update_JxW_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    double min_entropy = std::numeric_limits<double>::max(),
           max_entropy = -std::numeric_limits<double>::max(),
           area = 0,
           entropy_integrated = 0;

    // loop over all locally owned cells and evaluate the entrop
    // at all quadrature points. keep a running tally of the
    // integral over the entropy as well as the area and the
    // maximal and minimal entropy
    typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values.get_function_values (old_temperature_solution,
                                         old_temperature_values);
          fe_values.get_function_values (old_old_temperature_solution,
                                         old_old_temperature_values);
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double T = (old_temperature_values[q] +
                                old_old_temperature_values[q]) / 2;
              const double entropy = ((T-average_temperature) *
                                      (T-average_temperature));

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

    Utilities::MPI::sum (local_for_sum, MPI_COMM_WORLD, global_for_sum);
    Utilities::MPI::max (local_for_max, MPI_COMM_WORLD, global_for_max);

    const double average_entropy = global_for_sum[0] / global_for_sum[1];

    // return the maximal deviation of the entropy everywhere from the
    // average value
    return std::max(global_for_max[1] - average_entropy,
                    average_entropy - (-global_for_max[0]));
  }



  template <int dim>
  double
  Simulator<dim>::
  compute_viscosity (const std::vector<double>          &old_temperature,
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
                     const double                        global_u_infty,
                     const double                        global_T_variation,
                     const double                        average_temperature,
                     const double                        global_entropy_variation,
                     const std::vector<Point<dim> >     &evaluation_points,
                     const double                        cell_diameter) const
  {
    if (global_u_infty == 0)
      return 5e-3 * cell_diameter;

    const unsigned int n_q_points = old_temperature.size();

    double max_residual = 0;
    double max_velocity = 0;

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (old_velocity_values[q] +
                                 old_old_velocity_values[q]) / 2;

        const SymmetricTensor<2,dim> strain_rate = (old_strain_rates[q] +
                                                    old_old_strain_rates[q]) / 2;

        const double T = (old_temperature[q] + old_old_temperature[q]) / 2;
        const double p = (old_pressure[q] + old_old_pressure[q]) / 2;
        const double dT_dt = (old_temperature[q] - old_old_temperature[q])
                             / old_time_step;
        const double u_grad_T = u * (old_temperature_grads[q] +
                                     old_old_temperature_grads[q]) / 2;

        const double density              = material_model->density(T, p, evaluation_points[q]);
        const double thermal_conductivity = material_model->thermal_conductivity(T, p, evaluation_points[q]);
        const double c_P                  = material_model->specific_heat(T, p, evaluation_points[q]);
        const double k_Delta_T = thermal_conductivity
                                 * (old_temperature_laplacians[q] +
                                    old_old_temperature_laplacians[q]) / 2;

        // verify correctness of the heating term
        const double gamma
          = (parameters.radiogenic_heating_rate * density
             +
             (parameters.include_shear_heating
              ?
              2 * material_model->viscosity(T, p, evaluation_points[q]) *
              strain_rate * strain_rate
              :
              0));

        double residual
          = std::abs(density * c_P * (dT_dt + u_grad_T) - k_Delta_T - gamma);
        if (parameters.stabilization_alpha == 2)
          residual *= std::abs(T - average_temperature);

        max_residual = std::max (residual,        max_residual);
        max_velocity = std::max (std::sqrt (u*u), max_velocity);
      }

    const double max_viscosity = (parameters.stabilization_beta *
                                  max_velocity * cell_diameter);
    if (timestep_number <= 1)
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
    const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
    const unsigned int   n_q_points      = scratch.stokes_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    scratch.stokes_fe_values.reinit (cell);

    typename DoFHandler<dim>::active_cell_iterator
    temperature_cell (&triangulation,
                      cell->level(),
                      cell->index(),
                      &temperature_dof_handler);
    scratch.temperature_fe_values.reinit (temperature_cell);
    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.stokes_fe_values[pressure].get_function_values(old_stokes_solution,
                                                           scratch.old_pressure_values);

    data.local_matrix = 0;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
//TODO: make the temperature be something useful here
        const double old_temperature = scratch.old_temperature_values[q];
        const double old_pressure = scratch.old_pressure_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grads_phi_u[k] = scratch.stokes_fe_values[velocities].symmetric_gradient(k,q);
            scratch.phi_p[k]       = scratch.stokes_fe_values[pressure].value (k, q);
          }

        double eta = material_model->viscosity(old_temperature,
                                               old_pressure,
                                               scratch.stokes_fe_values.quadrature_point(q) );

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            if (stokes_fe.system_to_component_index(i).first
                ==
                stokes_fe.system_to_component_index(j).first)
              data.local_matrix(i,j) += (eta *
                                         (scratch.grads_phi_u[i] *
                                          scratch.grads_phi_u[j])
                                         +
                                         (1./eta) *
                                         pressure_scaling *
                                         pressure_scaling *
                                         (scratch.phi_p[i] * scratch.phi_p[j]))
                                        * scratch.stokes_fe_values.JxW(q);
      }

    cell->get_dof_indices (data.local_dof_indices);
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    current_stokes_constraints.distribute_local_to_global (data.local_matrix,
                                                           data.local_dof_indices,
                                                           stokes_preconditioner_matrix);
  }



  template <int dim>
  void
  Simulator<dim>::assemble_stokes_preconditioner ()
  {
    stokes_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.end()),
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
         StokesPreconditioner<dim> (stokes_fe, quadrature_formula,
                                    mapping,
                                    update_JxW_values |
                                    update_values |
                                    update_gradients |
                                    update_quadrature_points,
                                    temperature_fe,
                                    update_values),
         internal::Assembly::CopyData::
         StokesPreconditioner<dim> (stokes_fe));

    stokes_preconditioner_matrix.compress();
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
    std::vector<bool>  velocity_components (dim+1,true);
    velocity_components[dim] = false;
    DoFTools::extract_constant_modes (stokes_dof_handler, velocity_components,
                                      constant_modes);

    Mp_preconditioner.reset (new TrilinosWrappers::PreconditionILU());
    Amg_preconditioner.reset (new TrilinosWrappers::PreconditionAMG());

    TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = true;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;

    Mp_preconditioner->initialize (stokes_preconditioner_matrix.block(1,1));
    Amg_preconditioner->initialize (stokes_preconditioner_matrix.block(0,0),
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
    const unsigned int dofs_per_cell = scratch.stokes_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.stokes_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    scratch.stokes_fe_values.reinit (cell);

    typename DoFHandler<dim>::active_cell_iterator
    temperature_cell (&triangulation,
                      cell->level(),
                      cell->index(),
                      &temperature_dof_handler);
    scratch.temperature_fe_values.reinit (temperature_cell);
    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.stokes_fe_values[pressure].get_function_values(old_stokes_solution,
                                                           scratch.old_pressure_values);
    scratch.stokes_fe_values[velocities].get_function_values(old_stokes_solution,
                                                             scratch.old_velocity_values);

    // cache whether the model is compressible or not
    const bool is_compressible = material_model->is_compressible ();

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    if (is_compressible)
      data.local_pressure_shape_function_integrals = 0;


    for (unsigned int q=0; q<n_q_points; ++q)
      {
        const double old_temperature = scratch.old_temperature_values[q];
        const double old_pressure = scratch.old_pressure_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_u[k] = scratch.stokes_fe_values[velocities].value (k,q);
            scratch.phi_p[k]       = scratch.stokes_fe_values[pressure].value (k, q);
            if (rebuild_stokes_matrix)
              {
                scratch.grads_phi_u[k] = scratch.stokes_fe_values[velocities].symmetric_gradient(k,q);
                scratch.div_phi_u[k]   = scratch.stokes_fe_values[velocities].divergence (k, q);
              }
          }

        const double eta = material_model->viscosity(old_temperature,
                                                     old_pressure,
                                                     scratch.stokes_fe_values.quadrature_point(q));

        const Tensor<1,dim>
        gravity = gravity_model->gravity_vector (scratch.stokes_fe_values.quadrature_point(q));

        const double compressibility
          = (is_compressible
             ?
             material_model->compressibility(old_temperature,
                                             old_pressure,
                                             scratch.stokes_fe_values
                                             .quadrature_point(q))
             :
             1./0.);
        const double density = material_model->density(old_temperature,
                                                       old_pressure,
                                                       scratch.stokes_fe_values
                                                       .quadrature_point(q));

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
                                        * scratch.stokes_fe_values.JxW(q);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs(i) += (  // TODO: extrapolation von old_velocity
                                 (density * gravity * scratch.phi_u[i])
                                 + (is_compressible
                                    ?
                                    (pressure_scaling *
                                     compressibility * density *
                                     (scratch.old_velocity_values[q] * gravity) *
                                     scratch.phi_p[i])
                                    :
                                    0)
                               )
                               * scratch.stokes_fe_values.JxW(q);

        if (is_compressible)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            data.local_pressure_shape_function_integrals(i) += scratch.phi_p[i] * scratch.stokes_fe_values.JxW(q);
      }

    cell->get_dof_indices (data.local_dof_indices);
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    if (rebuild_stokes_matrix == true)
      current_stokes_constraints.distribute_local_to_global (data.local_matrix,
                                                             data.local_rhs,
                                                             data.local_dof_indices,
                                                             stokes_matrix,
                                                             stokes_rhs);
    else
      current_stokes_constraints.distribute_local_to_global (data.local_rhs,
                                                             data.local_dof_indices,
                                                             stokes_rhs);

    if (material_model->is_compressible())
      current_stokes_constraints.distribute_local_to_global (data.local_pressure_shape_function_integrals,
                                                             data.local_dof_indices,
                                                             pressure_shape_function_integrals);
  }



  template <int dim>
  void Simulator<dim>::assemble_stokes_system ()
  {
    computing_timer.enter_section ("   Assemble Stokes system");

    if (rebuild_stokes_matrix == true)
      stokes_matrix = 0;

    stokes_rhs = 0;
    if (material_model->is_compressible())
      pressure_shape_function_integrals = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.end()),
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
         StokesSystem<dim> (stokes_fe, mapping, quadrature_formula,
                            (update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             (rebuild_stokes_matrix == true
                              ?
                              update_gradients
                              :
                              UpdateFlags(0))),
                            temperature_fe,
                            update_values),
         internal::Assembly::CopyData::
         StokesSystem<dim> (stokes_fe));

    stokes_matrix.compress();
    stokes_rhs.compress(Add);

    if (material_model->is_compressible())
      pressure_shape_function_integrals.compress(Add);

    rebuild_stokes_matrix = false;

    pcout << std::endl;
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
    const bool use_bdf2_scheme = (timestep_number != 0);



    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    const unsigned int dofs_per_cell = scratch.temperature_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.temperature_fe_values.n_quadrature_points;

    scratch.temperature_fe_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_matrix = 0;
    data.local_rhs = 0;

    typename DoFHandler<dim>::active_cell_iterator
    stokes_cell (&triangulation,
                 cell->level(),
                 cell->index(),
                 &stokes_dof_handler);
    scratch.stokes_fe_values.reinit (stokes_cell);

    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.temperature_fe_values.get_function_values (old_old_temperature_solution,
                                                       scratch.old_old_temperature_values);

    scratch.temperature_fe_values.get_function_gradients (old_temperature_solution,
                                                          scratch.old_temperature_grads);
    scratch.temperature_fe_values.get_function_gradients (old_old_temperature_solution,
                                                          scratch.old_old_temperature_grads);

    scratch.temperature_fe_values.get_function_laplacians (old_temperature_solution,
                                                           scratch.old_temperature_laplacians);
    scratch.temperature_fe_values.get_function_laplacians (old_old_temperature_solution,
                                                           scratch.old_old_temperature_laplacians);

    scratch.stokes_fe_values[velocities].get_function_values (stokes_solution,
                                                              scratch.old_velocity_values);
    scratch.stokes_fe_values[velocities].get_function_values (old_stokes_solution,
                                                              scratch.old_old_velocity_values);
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (stokes_solution,
                                                                           scratch.old_strain_rates);
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (old_stokes_solution,
                                                                           scratch.old_old_strain_rates);

    scratch.stokes_fe_values[pressure].get_function_values (stokes_solution,
                                                            scratch.old_pressure);
    scratch.stokes_fe_values[pressure].get_function_values (old_stokes_solution,
                                                            scratch.old_old_pressure);

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
                           global_max_velocity,
                           global_T_range.second - global_T_range.first,
                           0.5 * (global_T_range.second + global_T_range.first),
                           global_entropy_variation,
                           scratch.temperature_fe_values.get_quadrature_points(),
                           cell->diameter());

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_T[k] = scratch.temperature_fe_values.shape_grad (k,q);
            scratch.phi_T[k]      = scratch.temperature_fe_values.shape_value (k, q);
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

        const double ext_T
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_temperature_values[q], scratch.old_temperature_values[q]);

        const Tensor<1,dim> extrapolated_u
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_velocity_values[q], scratch.old_velocity_values[q]);

        const SymmetricTensor<2,dim> extrapolated_strain_rate
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_strain_rates[q], scratch.old_strain_rates[q]);

        const double ext_pressure
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_pressure[q], scratch.old_pressure[q]);

        const double density              = material_model->density(ext_T,
                                                                    ext_pressure,
                                                                    scratch.temperature_fe_values.quadrature_point(q));
        const double thermal_conductivity = material_model->thermal_conductivity(ext_T,
                                                                                 ext_pressure,
                                                                                 scratch.temperature_fe_values.quadrature_point(q));
        const double c_P                  = material_model->specific_heat (ext_T,
                                                                           ext_pressure,
                                                                           scratch.temperature_fe_values.quadrature_point(q));

        //TODO: this is the wrong formula for the compressible case
        const double gamma
          = (parameters.radiogenic_heating_rate * density
             +
             (parameters.include_shear_heating
              ?
              2 * material_model->viscosity(ext_T,
                                            ext_pressure,
                                            scratch.temperature_fe_values.quadrature_point(q)) *
              extrapolated_strain_rate * extrapolated_strain_rate
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
                                 scratch.temperature_fe_values.JxW(q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                          (time_step + old_time_step)) : 1.0;
                data.local_matrix(i,j)
                += (
                     (time_step * (thermal_conductivity +
                                   nu * density * c_P) * scratch.grad_phi_T[i] * scratch.grad_phi_T[j])
                     + ((time_step * (extrapolated_u * scratch.grad_phi_T[j] * scratch.phi_T[i]))
                        + (factor * scratch.phi_T[i] * scratch.phi_T[j])) * density * c_P
                   )
                   * scratch.temperature_fe_values.JxW(q);
              }
          }
      }
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_temperature_system (const internal::Assembly::CopyData::TemperatureSystem<dim> &data)
  {
    temperature_constraints.distribute_local_to_global (data.local_matrix,
                                                        data.local_rhs,
                                                        data.local_dof_indices,
                                                        temperature_matrix,
                                                        temperature_rhs
                                                       );
  }


  template <int dim>
  void Simulator<dim>::assemble_temperature_system ()
  {
    computing_timer.enter_section ("   Assemble temperature system");

    temperature_matrix = 0;
    temperature_rhs = 0;

    const std::pair<double,double>
    global_T_range = get_extrapolated_temperature_range();

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_temperature_system,
                          this,
                          global_T_range,
                          get_maximal_velocity(),
                          // use the mid temperature instead of the
                          // integral mean. results are not very
                          // sensitive to this and this is far simpler
                          get_entropy_variation ((global_T_range.first +
                                                  global_T_range.second) / 2),
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_temperature_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         TemperatureSystem<dim> (temperature_fe, stokes_fe, mapping, QGauss<dim>(parameters.temperature_degree+2)),
         internal::Assembly::CopyData::
         TemperatureSystem<dim> (temperature_fe));

    temperature_matrix.compress();
    temperature_rhs.compress(Add);

    computing_timer.exit_section();

    computing_timer.enter_section ("   Build temperature preconditioner");
    T_preconditioner.reset (new TrilinosWrappers::PreconditionILU());
    T_preconditioner->initialize (temperature_matrix);

    computing_timer.exit_section();
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
  template void Simulator<deal_II_dimension>::local_assemble_stokes_preconditioner (
    const DoFHandler<deal_II_dimension>::active_cell_iterator &cell,
    internal::Assembly::Scratch::StokesPreconditioner<deal_II_dimension> &scratch,
    internal::Assembly::CopyData::StokesPreconditioner<deal_II_dimension> &data);

  template void Simulator<deal_II_dimension>::copy_local_to_global_stokes_preconditioner (
    const internal::Assembly::CopyData::StokesPreconditioner<deal_II_dimension> &data);

  template void Simulator<deal_II_dimension>::assemble_stokes_preconditioner ();

  template void Simulator<deal_II_dimension>::build_stokes_preconditioner ();

  template void Simulator<deal_II_dimension>::local_assemble_stokes_system (
    const DoFHandler<deal_II_dimension>::active_cell_iterator &cell,
    internal::Assembly::Scratch::StokesSystem<deal_II_dimension>  &scratch,
    internal::Assembly::CopyData::StokesSystem<deal_II_dimension> &data);

  template void Simulator<deal_II_dimension>::copy_local_to_global_stokes_system (
    const internal::Assembly::CopyData::StokesSystem<deal_II_dimension> &data);

  template void Simulator<deal_II_dimension>::assemble_stokes_system ();

  template void Simulator<deal_II_dimension>::copy_local_to_global_temperature_system (
    const internal::Assembly::CopyData::TemperatureSystem<deal_II_dimension> &data);

  template void Simulator<deal_II_dimension>::assemble_temperature_system ();

  template void Simulator<deal_II_dimension>::local_assemble_temperature_system (
    const std::pair<double,double> global_T_range,
    const double                   global_max_velocity,
    const double                   global_entropy_variation,
    const DoFHandler<deal_II_dimension>::active_cell_iterator &cell,
    internal::Assembly::Scratch::TemperatureSystem<deal_II_dimension>  &scratch,
    internal::Assembly::CopyData::TemperatureSystem<deal_II_dimension> &data);
}

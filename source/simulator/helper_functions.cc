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
#include <aspect/global.h>

#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vectors.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <string>


namespace aspect
{

  template <int dim>
  void Simulator<dim>::output_statistics()
  {
    if (!aspect::output_parallel_statistics)
      return;

    Utilities::System::MemoryStats stats;
    Utilities::System::get_memory_stats(stats);
    pcout << "VmPeak (proc0): " << stats.VmPeak/1024 << " mb" << std::endl;

    // memory consumption:
    const double mb = 1024*1024; //convert from bytes into mb
    pcout << "memory in MB:" << std::endl
          << "* tria " << triangulation.memory_consumption()/mb << std::endl
          << "  - p4est " << triangulation.memory_consumption_p4est()/mb << std::endl
          << "* DoFHandler " << dof_handler.memory_consumption()/mb <<std::endl
          << "* ConstraintMatrix " << constraints.memory_consumption()/mb << std::endl
          << "* current_constraints " << current_constraints.memory_consumption()/mb << std::endl
          << "* Matrix " << system_matrix.memory_consumption()/mb << std::endl
          << "* 5 Vectors " << 5*solution.memory_consumption()/mb << std::endl
          << "* preconditioner " << (system_preconditioner_matrix.memory_consumption()
                                     + Amg_preconditioner->memory_consumption()
                                     /*+Mp_preconditioner->memory_consumption()
                                                                      +T_preconditioner->memory_consumption()*/)/mb
          << std::endl
          << "  - matrix " << system_preconditioner_matrix.memory_consumption()/mb << std::endl
          << "  - prec vel " << Amg_preconditioner->memory_consumption()/mb << std::endl
          << "  - prec mass " << 0/*Mp_preconditioner->memory_consumption()/mb*/ << std::endl
          << "  - prec T " << 0/*T_preconditioner->memory_consumption()/mb*/ << std::endl
          << std::endl;
  }



  /**
   * Find the largest velocity throughout the domain.
   **/
  template <int dim>
  double Simulator<dim>::get_maximal_velocity (
    const LinearAlgebra::BlockVector &solution) const
  {
    // use a quadrature formula that has one point at
    // the location of each degree of freedom in the
    // velocity element
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();


    FEValues<dim> fe_values (mapping, finite_element, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_velocity = 0;

    // loop over all locally owned cells and evaluate the velocities at each
    // quadrature point (i.e. each node). keep a running tally of the largest
    // such velocity
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (solution,
                                                     velocity_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
        }

    // return the largest value over all processors
    return Utilities::MPI::max (max_local_velocity, mpi_communicator);
  }



  template <int dim>
  double Simulator<dim>::compute_time_step () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, finite_element, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_speed_over_meshsize = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (solution,
                                                     velocity_values);

          double max_local_velocity = 0;
          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
          max_local_speed_over_meshsize = std::max(max_local_speed_over_meshsize,
                                                   max_local_velocity
                                                   /
                                                   cell->minimum_vertex_distance());
        }

    const double max_global_speed_over_meshsize
      = Utilities::MPI::max (max_local_speed_over_meshsize, mpi_communicator);

    // if the velocity is zero, then it is somewhat arbitrary what time step
    // we should choose. in that case, do as if the velocity was one
    if (max_global_speed_over_meshsize == 0)
      return  (parameters.CFL_number / (parameters.temperature_degree *
                                        1));
    else
      return (parameters.CFL_number / (parameters.temperature_degree *
                                       max_global_speed_over_meshsize));
  }



  template <int dim>
  std::pair<double,double>
  Simulator<dim>::get_extrapolated_temperature_range () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.temperature_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar temperature (dim+1);

    FEValues<dim> fe_values (mapping, finite_element, quadrature_formula,
                             update_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    // This presets the minimum with a bigger
    // and the maximum with a smaller number
    // than one that is going to appear. Will
    // be overwritten in the cell loop or in
    // the communication step at the
    // latest.
    double min_local_temperature = std::numeric_limits<double>::max(),
           max_local_temperature = -std::numeric_limits<double>::max();

    if (timestep_number != 0)
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[temperature].get_function_values (old_solution,
                                                          old_temperature_values);
              fe_values[temperature].get_function_values (old_old_solution,
                                                          old_old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double extrapolated_temperature =
                    (1. + time_step/old_time_step) * old_temperature_values[q]-
                    time_step/old_time_step * old_old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    extrapolated_temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    extrapolated_temperature);
                }
            }
      }
    else
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[temperature].get_function_values (old_solution,
                                                          old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double extrapolated_temperature = old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    extrapolated_temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    extrapolated_temperature);
                }
            }
      }

    return std::make_pair(-Utilities::MPI::max (-min_local_temperature,
                                                mpi_communicator),
                          Utilities::MPI::max (max_local_temperature,
                                               mpi_communicator));
  }



  /*
   * normalize the pressure by calculating the surface integral of the pressure on the outer
   * shell and subtracting this from all pressure nodes.
   */
  template <int dim>
  void Simulator<dim>::normalize_pressure(LinearAlgebra::BlockVector &vector)
  {
    if (parameters.pressure_normalization == "no")
      return;

    double my_pressure = 0.0;
    double my_area = 0.0;
    if (parameters.pressure_normalization == "surface")
      {
        QGauss < dim - 1 > quadrature (parameters.stokes_velocity_degree + 1);

        const unsigned int n_q_points = quadrature.size();
        FEFaceValues<dim> fe_face_values (mapping, finite_element,  quadrature,
                                          update_JxW_values | update_values);
        const FEValuesExtractors::Scalar pressure (dim);

        std::vector<double> pressure_values(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                {
                  const typename DoFHandler<dim>::face_iterator face = cell->face (face_no);
                  if (face->at_boundary()
                      &&
                      (geometry_model->depth (face->center()) <
                       (face->diameter() / std::sqrt(1.*dim-1) / 3)))
                    {
                      fe_face_values.reinit (cell, face_no);
                      fe_face_values[pressure].get_function_values(vector,
                                                                   pressure_values);

                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          my_pressure += pressure_values[q]
                                         * fe_face_values.JxW (q);
                          my_area += fe_face_values.JxW (q);
                        }
                    }
                }
            }
      }
    else if (parameters.pressure_normalization=="volume")
      {
        QGauss<dim> quadrature (parameters.stokes_velocity_degree + 1);

        const unsigned int n_q_points = quadrature.size();
        FEValues<dim> fe_values (mapping, finite_element,  quadrature,
                                 update_JxW_values | update_values);
        const FEValuesExtractors::Scalar pressure (dim);

        std::vector<double> pressure_values(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[pressure].get_function_values(vector,
                                                      pressure_values);

              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  my_pressure += pressure_values[q]
                                 * fe_values.JxW (q);
                  my_area += fe_values.JxW (q);
                }
            }
      }
    else
      AssertThrow (false, ExcMessage("Invalid pressure normalization method: " +
                                     parameters.pressure_normalization));

    double adjust = 0;
    // sum up the integrals from each processor
    {
      const double my_temp[2] = {my_pressure, my_area};
      double temp[2];
      Utilities::MPI::sum (my_temp, mpi_communicator, temp);
      if (parameters.pressure_normalization == "surface")
        {

          adjust = -temp[0]/temp[1] + parameters.surface_pressure;
        }
      else if (parameters.pressure_normalization == "volume")
        {
          adjust = -temp[0];
        }
      else
        AssertThrow(false, ExcNotImplemented());
    }

    if (parameters.use_locally_conservative_discretization == false)
      vector.block(1).add(adjust);
    else
      {
        // this case is a bit more complicated: if the condition above is false
        // then we use the FE_DGP element for which the shape functions do not
        // add up to one; consequently, adding a constant to all degrees of
        // freedom does not alter the overall function by that constant, but
        // by something different
        //
        // we can work around this by using the documented property of the
        // FE_DGP element that the first shape function is constant.
        // consequently, adding the adjustment to the global function is
        // achieved by adding the adjustment to the first pressure degree
        // of freedom on each cell.
        //
        // an additional complication is that we can't modify individual
        // elements of the solution vector since that one has ghost element.
        // rather, we first need to localize it and then distribute back
//TODO: It isn't strictly necessary to copy system_rhs here. all we want
//      is to get the partitioning
        LinearAlgebra::BlockVector distributed_vector (system_rhs);
        distributed_vector = vector;

        Assert (dynamic_cast<const FE_DGP<dim>*>(&finite_element.base_element(1)) != 0,
                ExcInternalError());
        std::vector<unsigned int> local_dof_indices (finite_element.dofs_per_cell);
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = finite_element.component_to_system_index (dim, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (dof_handler.locally_owned_dofs().is_element(first_pressure_dof),
                      ExcInternalError());
              Assert (local_dof_indices[first_pressure_dof] >= vector.block(0).size(),
                      ExcInternalError());

              // then adjust its value
              distributed_vector(local_dof_indices[first_pressure_dof]) += adjust;
            }

        // now get back to the original vector
        vector = distributed_vector;
      }
  }



  /**
   * This routine adjusts the second block of the right hand side of the
   * system containing the compressibility, so that the system becomes
   * compatible. See the general documentation of this class for more
   * information.
   */
  template <int dim>
  void Simulator<dim>::make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector)
  {
    if (parameters.use_locally_conservative_discretization)
      AssertThrow(false, ExcNotImplemented());

    const double mean       = vector.block(1).mean_value();
    const double correction = -mean*vector.block(1).size()/global_volume;

    vector.block(1).add(correction, pressure_shape_function_integrals.block(1));
  }


  template <int dim>
  void Simulator<dim>::compute_running_average(std::vector<double> &values, const int ncells) const
  {

    std::vector<double> temp(values.size());

    for (unsigned int idx=0; idx<values.size(); idx++)
      {
        double sum = 0;
        for (int isum=-ncells; isum<=ncells; isum++)
          {
            sum += values[std::max(0,std::min((int) values.size()-1,(int) idx+isum))];
          }
        temp[idx] = sum/((double) (ncells*2+1));
      }
    values = temp;

  }


//TODO: unify the following functions
  template <int dim>
  void Simulator<dim>::compute_depth_average_temperature(std::vector<double> &values) const
  {
    const unsigned int num_slices = 100;
    values.resize(num_slices);
    std::vector<unsigned int> counts(num_slices);

    // this yields 10^dim quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();
    const double max_depth = geometry_model->maximal_depth();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values | update_quadrature_points);
    const FEValuesExtractors::Scalar temperature (dim+1);
    std::vector<double> temperature_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[temperature].get_function_values (this->solution,
                                                      temperature_values);
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const Point<dim> &p = fe_values.quadrature_point(q);
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              ++counts[idx];
              values[idx]+= temperature_values[q];
            }
        }

    std::vector<double> values_all(num_slices);
    std::vector<unsigned int> counts_all(num_slices);
    Utilities::MPI::sum(counts, mpi_communicator, counts_all);
    Utilities::MPI::sum(values, mpi_communicator, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(counts_all[i])+1e-20);
  }


  template <int dim>
  void Simulator<dim>::compute_depth_average_viscosity(std::vector<double> &values) const
  {
    const unsigned int num_slices = 100;
    values.resize(num_slices);
    std::vector<unsigned int> counts(num_slices);

    // this yields 100 quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();
    const double max_depth = geometry_model->maximal_depth();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values | update_gradients | update_quadrature_points);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);

    std::vector<SymmetricTensor<2,dim> > strain_rates(n_q_points);
    std::vector<double> pressure_values(n_q_points);
    std::vector<double> temperature_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // compute the integral quantities by quadrature
    unsigned int cell_index = 0;
    for (; cell!=endc; ++cell,++cell_index)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[pressure].get_function_values (this->solution,
                                                   pressure_values);
          fe_values[temperature].get_function_values (this->solution,
                                                      temperature_values);
          fe_values[velocities].get_function_symmetric_gradients (this->solution,
                                                                  strain_rates);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              const double viscosity = material_model->viscosity(temperature_values[q],
                                                                 pressure_values[q],
                                                                 strain_rates[q],
                                                                 fe_values.quadrature_point(q));
              ++counts[idx];
              values[idx] += viscosity;
            }
        }
    std::vector<double> values_all(num_slices);
    std::vector<unsigned int> counts_all(num_slices);
    Utilities::MPI::sum(counts, MPI_COMM_WORLD, counts_all);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(counts_all[i])+1e-20);
  }



  template <int dim>
  void Simulator<dim>::compute_depth_average_velocity_magnitude(std::vector<double> &values) const
  {
    const unsigned int num_slices = 100;
    values.resize(num_slices);
    std::vector<unsigned int> counts(num_slices);

    // this yields 100 quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();
    const double max_depth = geometry_model->maximal_depth();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values | update_quadrature_points | update_JxW_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (this->solution,
                                                     velocity_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const Point<dim> &p = fe_values.quadrature_point(q);
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              ++counts[idx];
              values[idx]+= ((velocity_values[q] * velocity_values[q]) *
                             fe_values.JxW(q));
            }
        }
    std::vector<double> values_all(num_slices);
    std::vector<unsigned int> counts_all(num_slices);
    Utilities::MPI::sum(counts, MPI_COMM_WORLD, counts_all);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(counts_all[i])+1e-20)*year_in_seconds;
  }



  template <int dim>
  void Simulator<dim>::compute_depth_average_sinking_velocity(std::vector<double> &values) const
  {
    const unsigned int num_slices = 100;
    values.resize(num_slices);
    std::vector<unsigned int> counts(num_slices);

    // this yields 100 quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();
    const double max_depth = geometry_model->maximal_depth();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values   |
                             update_quadrature_points |
                             update_JxW_values);

    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);


    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (this->solution,
                                                     velocity_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const Point<dim> p = fe_values.quadrature_point(q);
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              ++counts[idx];
              values[idx]+= std::fabs(std::min(-1e-16,p*velocity_values[q]/p.norm()))*year_in_seconds;
            }
        }
    std::vector<double> values_all(num_slices);
    std::vector<unsigned int> counts_all(num_slices);
    Utilities::MPI::sum(counts, MPI_COMM_WORLD, counts_all);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(counts_all[i])+1e-20);
  }



  template <int dim>
  void Simulator<dim>::compute_depth_average_Vs(std::vector<double> &values) const
  {

    std::vector<double> average_temperature;
    compute_depth_average_temperature(average_temperature);

    values.resize(average_temperature.size());
    std::vector<unsigned int> counts(average_temperature.size());

    const unsigned int num_slices = average_temperature.size();
    const double max_depth = geometry_model->maximal_depth();

    // this yields 100 quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values | update_quadrature_points);

    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<double> pressure_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // compute the integral quantities by quadrature
    unsigned int cell_index = 0;
    for (; cell!=endc; ++cell,++cell_index)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[pressure].get_function_values (this->solution,
                                                   pressure_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              const double Vs_depth_average = material_model->seismic_Vs(average_temperature[idx], pressure_values[q]);
              ++counts[idx];
              values[idx] += Vs_depth_average;
            }
        }
    std::vector<double> values_all(num_slices);
    std::vector<unsigned int> counts_all(num_slices);
    Utilities::MPI::sum(counts, MPI_COMM_WORLD, counts_all);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(counts_all[i])+1e-20);
  }

  template <int dim>
  void Simulator<dim>::compute_depth_average_Vp(std::vector<double> &values) const
  {

    std::vector<double> average_temperature;

    compute_depth_average_temperature(average_temperature);

    values.resize(average_temperature.size());
    std::vector<unsigned int> counts(average_temperature.size());

    const unsigned int num_slices = average_temperature.size();
    const double max_depth = geometry_model->maximal_depth();

    // this yields 100 quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values   |
                             update_quadrature_points );

    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<double> pressure_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // compute the integral quantities by quadrature
    unsigned int cell_index = 0;
    for (; cell!=endc; ++cell,++cell_index)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[pressure].get_function_values (this->solution,
                                                   pressure_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              const double Vp_depth_average = material_model->seismic_Vp(average_temperature[idx], pressure_values[q]);
              ++counts[idx];
              values[idx] += Vp_depth_average;
            }
        }
    std::vector<double> values_all(num_slices);
    std::vector<unsigned int> counts_all(num_slices);
    Utilities::MPI::sum(counts, MPI_COMM_WORLD, counts_all);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(counts_all[i])+1e-20);
  }



  template <int dim>
  void Simulator<dim>::compute_Vs_anomaly(Vector<float> &values) const
  {
    const int npoints = 2; // npoint in running average half-width of window
    std::vector<double> Vs_depth_average;

    compute_depth_average_Vs(Vs_depth_average);
    compute_running_average(Vs_depth_average, npoints);

    const unsigned int num_slices = Vs_depth_average.size();
    const double max_depth = geometry_model->maximal_depth();

    // evaluate a single point per cell
    const QMidpoint<dim> quadrature_formula;
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values   |
                             update_quadrature_points );

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);

    std::vector<double> pressure_values(n_q_points);
    std::vector<double> temperature_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    unsigned int cell_index = 0;
    for (; cell!=endc; ++cell,++cell_index)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[pressure].get_function_values (this->solution,
                                                   pressure_values);
          fe_values[temperature].get_function_values (this->solution,
                                                      temperature_values);

          const double Vs = material_model->seismic_Vs(temperature_values[0], pressure_values[0]);
          const double depth = geometry_model->depth(fe_values.quadrature_point(0));
          const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
          Assert(idx<num_slices, ExcInternalError());

          values(cell_index) = (Vs - Vs_depth_average[idx])/Vs_depth_average[idx]*1e2;
        }
  }



  template <int dim>
  void Simulator<dim>::compute_Vp_anomaly(Vector<float> &values) const
  {

    const int npoints = 2; // npoint in running average half-width of window
    std::vector<double> Vp_depth_average;

    compute_depth_average_Vp(Vp_depth_average);
    compute_running_average(Vp_depth_average, npoints);

    const unsigned int num_slices = Vp_depth_average.size();
    const double max_depth = geometry_model->maximal_depth();

    const QMidpoint<dim> quadrature_formula;
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values   |
                             update_quadrature_points );


    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Scalar temperature (dim+1);

    std::vector<double> pressure_values(n_q_points);
    std::vector<double> temperature_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // compute the integral quantities by quadrature
    unsigned int cell_index = 0;
    for (; cell!=endc; ++cell,++cell_index)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[pressure].get_function_values (this->solution,
                                                   pressure_values);
          fe_values[temperature].get_function_values (this->solution,
                                                      temperature_values);

          const double Vp = material_model->seismic_Vp(temperature_values[0], pressure_values[0]);
          const double depth = geometry_model->depth(fe_values.quadrature_point(0));
          const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
          Assert(idx<num_slices, ExcInternalError());
          values(cell_index) = (Vp - Vp_depth_average[idx])/Vp_depth_average[idx]*1e2;
        }
  }
}

// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::normalize_pressure(LinearAlgebra::BlockVector &vector); \
  template double Simulator<dim>::get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const; \
  template std::pair<double,double> Simulator<dim>::get_extrapolated_temperature_range () const; \
  template double Simulator<dim>::compute_time_step () const; \
  template void Simulator<dim>::make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector); \
  template void Simulator<dim>::compute_depth_average_temperature(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_viscosity(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_velocity_magnitude(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_sinking_velocity(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_Vs(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_Vp(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_Vs_anomaly(Vector<float> &values) const; \
  template void Simulator<dim>::compute_Vp_anomaly(Vector<float> &values) const; \
  template void Simulator<dim>::output_statistics();

  ASPECT_INSTANTIATE(INSTANTIATE)
}

/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University,
     Timo Heister, University of Goettingen, 2008-2011 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

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

  /**
   * Find the largest velocity throughout the domain.
   **/
  template <int dim>
  double Simulator<dim>::get_maximal_velocity () const
  {
    // use a quadrature formula that has one point at
    // the location of each degree of freedom in the
    // velocity element
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();


    FEValues<dim> fe_values (mapping, stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_velocity = 0;

    // loop over all locally owned cells and evaluate the velocities at each
    // quadrature point (i.e. each node). keep a running tally of the largest
    // such velocity
    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (stokes_solution,
                                                     velocity_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
        }

    // return the largest value over all processors
    return Utilities::MPI::max (max_local_velocity, MPI_COMM_WORLD);
  }



  /**
   * Similar function to before, but we now
   * compute the cfl number, i.e., maximal
   * velocity on a cell divided by the cell
   * diameter.
   */
  template <int dim>
  double Simulator<dim>::compute_time_step () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_speed_over_meshsize = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (stokes_solution,
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
      = Utilities::MPI::max (max_local_speed_over_meshsize, MPI_COMM_WORLD);

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

    FEValues<dim> fe_values (mapping, temperature_fe, quadrature_formula,
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
                  const double temperature =
                    (1. + time_step/old_time_step) * old_temperature_values[q]-
                    time_step/old_time_step * old_old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    temperature);
                }
            }
      }
    else
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values.get_function_values (old_temperature_solution,
                                             old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double temperature = old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    temperature);
                }
            }
      }

    return std::make_pair(-Utilities::MPI::max (-min_local_temperature,
                                                MPI_COMM_WORLD),
                          Utilities::MPI::max (max_local_temperature,
                                               MPI_COMM_WORLD));
  }



  /*
   * normalize the pressure by calculating the surface integral of the pressure on the outer
   * shell and subtracting this from all pressure nodes.
   */
  template <int dim>
  void Simulator<dim>::normalize_pressure(TrilinosWrappers::MPI::BlockVector &vector)
  {
    // TODO: somehow parameterize based on the geometry model
    // on which parts of the boundary the pressure should be
    // zero
    if (dynamic_cast<const GeometryModel::SphericalShell<dim>*>(&*geometry_model) == 0)
      return;

    double my_pressure = 0.0;
    double my_area = 0.0;
    {
      QGauss < dim - 1 > quadrature (parameters.stokes_velocity_degree + 1);

      const unsigned int n_q_points = quadrature.size();
      FEFaceValues<dim> fe_face_values (mapping, stokes_fe,  quadrature,
                                        update_JxW_values | update_values);
      const FEValuesExtractors::Scalar pressure (dim);

      std::vector<double> pressure_values(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator
      cell = stokes_dof_handler.begin_active(),
      endc = stokes_dof_handler.end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
              {
                const typename DoFHandler<dim>::face_iterator face = cell->face (face_no);
                if (face->at_boundary() && face->boundary_indicator() == 1) // outer shell boundary
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

    double surf_pressure = 0;
    // sum up the surface integrals from each processor
    {
      const double my_temp[2] = {my_pressure, my_area};
      double temp[2];
      Utilities::MPI::sum (my_temp, MPI_COMM_WORLD, temp);
      surf_pressure = temp[0]/temp[1];
    }

    const double adjust = -surf_pressure + 1e7;
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
        Assert (dynamic_cast<const FE_DGP<dim>*>(&stokes_fe.base_element(1)) != 0,
                ExcInternalError());
        std::vector<unsigned int> local_dof_indices (stokes_fe.dofs_per_cell);
        typename DoFHandler<dim>::active_cell_iterator
        cell = stokes_dof_handler.begin_active(),
        endc = stokes_dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = stokes_fe.component_to_system_index (dim, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (stokes_dof_handler.locally_owned_dofs().is_element(first_pressure_dof),
                      ExcInternalError());
              Assert (local_dof_indices[first_pressure_dof] >= vector.block(0).size(),
                      ExcInternalError());

              // then adjust its value
              vector(local_dof_indices[first_pressure_dof]) += adjust;
            }
      }
  }


  /**
  * This routine adjusts the second block of the right hand side of the
  * system containing the compressibility, so that the system becomes
  * compatible: 0=\int div u = \int g
  * the helper vector h contains h_i=(q_i,1) with the pressure functions q_i
  * and we adjust the right hand side g by h_i \int g / |\Omega|
  */
  template <int dim>
  void Simulator<dim>::make_pressure_rhs_compatible(TrilinosWrappers::MPI::BlockVector &vector)
  {
    if (parameters.use_locally_conservative_discretization)
      throw ExcNotImplemented();

    const double mean       = vector.block(1).mean_value();
    const double correction = -mean*vector.block(1).size()/global_volume;

    vector.block(1).add(correction, pressure_shape_function_integrals.block(1));
  }

}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  template void Simulator<deal_II_dimension>::normalize_pressure(TrilinosWrappers::MPI::BlockVector &vector);

  template double Simulator<deal_II_dimension>::get_maximal_velocity () const;

  template std::pair<double,double> Simulator<deal_II_dimension>::get_extrapolated_temperature_range () const;

  template double Simulator<deal_II_dimension>::compute_time_step () const;

  template void Simulator<deal_II_dimension>::make_pressure_rhs_compatible(TrilinosWrappers::MPI::BlockVector &vector);
}

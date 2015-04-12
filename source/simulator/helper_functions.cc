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
#include <aspect/global.h>


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
#include <deal.II/numerics/vector_tools.h>

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
  Simulator<dim>::AdvectionField::
  AdvectionField (const FieldType field_type,
                  const unsigned int compositional_variable)
    :
    field_type (field_type),
    compositional_variable (compositional_variable)
  {
    if (field_type == temperature_field)
      Assert (compositional_variable == numbers::invalid_unsigned_int,
              ExcMessage ("You can't specify a compositional variable if you "
                          "have in fact selected the temperature."));
  }



  template <int dim>
  typename Simulator<dim>::AdvectionField
  Simulator<dim>::AdvectionField::temperature ()
  {
    return AdvectionField(temperature_field);
  }



  template <int dim>
  typename Simulator<dim>::AdvectionField
  Simulator<dim>::AdvectionField::composition (const unsigned int compositional_variable)
  {
    return AdvectionField(compositional_field,
                          compositional_variable);
  }


  template <int dim>
  bool
  Simulator<dim>::AdvectionField::is_temperature() const
  {
    return (field_type == temperature_field);
  }


  template <int dim>
  bool
  Simulator<dim>::AdvectionField::is_porosity(const Introspection<dim> &introspection) const
  {
    if (field_type != compositional_field)
      return false;
    else
	  return (introspection.name_for_compositional_index(compositional_variable) == "porosity");
  }


  template <int dim>
  unsigned int
  Simulator<dim>::AdvectionField::block_index(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.block_indices.temperature;
    else
      return introspection.block_indices.compositional_fields[compositional_variable];
  }

  template <int dim>
  unsigned int
  Simulator<dim>::AdvectionField::component_index(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.component_indices.temperature;
    else
      return introspection.component_indices.compositional_fields[compositional_variable];
  }

  template <int dim>
  unsigned int
  Simulator<dim>::AdvectionField::base_element(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.base_elements.temperature;
    else
      return introspection.base_elements.compositional_fields;
  }


  template <int dim>
  void Simulator<dim>::output_program_stats()
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
//                                     + Amg_preconditioner->memory_consumption()
                                     /*+Mp_preconditioner->memory_consumption()
                                                                      +T_preconditioner->memory_consumption()*/)/mb
          << std::endl
          << "  - matrix " << system_preconditioner_matrix.memory_consumption()/mb << std::endl
//          << "  - prec vel " << Amg_preconditioner->memory_consumption()/mb << std::endl
          << "  - prec mass " << 0/*Mp_preconditioner->memory_consumption()/mb*/ << std::endl
          << "  - prec T " << 0/*T_preconditioner->memory_consumption()/mb*/ << std::endl
          << std::endl;
  }



  namespace
  {
    /**
     * A function that writes the statistics object into a file.
     *
     * @param stat_file_name The name of the file into which the result
     * should go
     * @param copy_of_table A copy of the table that we're to write. Since
     * this function is called in the background on a separate thread,
     * the actual table might be modified while we are about to write
     * it, so we need to work on a copy. This copy is deleted at the end
     * of this function.
     */
    void do_output_statistics (const std::string stat_file_name,
                               const TableHandler *copy_of_table)
    {
      // write into a temporary file for now so that we don't
      // interrupt anyone who might want to look at the real
      // statistics file while the program is still running
      const std::string tmp_file_name = stat_file_name + " tmp";

      std::ofstream stat_file (tmp_file_name.c_str());
      copy_of_table->write_text (stat_file,
                                 TableHandler::table_with_separate_column_description);
      stat_file.close();

      // now move the temporary file into place
      std::rename(tmp_file_name.c_str(), stat_file_name.c_str());

      // delete the copy now:
      delete copy_of_table;
    }
  }


  template <int dim>
  void Simulator<dim>::output_statistics()
  {
    // only write the statistics file from processor zero
    if (Utilities::MPI::this_mpi_process(mpi_communicator)!=0)
      return;

    if (parameters.convert_to_years == true)
      {
        statistics.set_scientific("Time (years)", true);
        statistics.set_scientific("Time step size (years)", true);
      }
    else
      {
        statistics.set_scientific("Time (seconds)", true);
        statistics.set_scientific("Time step size (seconds)", true);
      }

    // formatting the table we're about to output and writing the
    // actual file may take some time, so do it on a separate
    // thread. we pass a pointer to a copy of the statistics
    // object which the called function then has to destroy
    //
    // before we can start working on a new thread, we need to
    // make sure that the previous thread is done or they'll
    // stomp on each other's feet
    output_statistics_thread.join();
    output_statistics_thread = Threads::new_thread (&do_output_statistics,
                                                    parameters.output_directory+"statistics",
                                                    new TableHandler(statistics));
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
          fe_values[introspection.extractors.velocities].get_function_values (solution,
                                                                              velocity_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
        }

    // return the largest value over all processors
    return Utilities::MPI::max (max_local_velocity, mpi_communicator);
  }



  template <int dim>
  std::pair<double,bool> Simulator<dim>::compute_time_step () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, finite_element, quadrature_formula,
                             update_values | (parameters.use_conduction_timestep || parameters.include_melt_transport
                                              ?
                                               (parameters.include_melt_transport
                                                ?
                                                update_gradients | update_quadrature_points
                                                :
                                                update_quadrature_points)
                                              :
                                              update_default));
    std::vector<Tensor<1,dim> > velocity_values(n_q_points), fluid_pressure_gradients(n_q_points);
    std::vector<double> pressure_values(n_q_points), fluid_pressure_values(n_q_points), temperature_values(n_q_points);
    std::vector<std::vector<double> > composition_values (parameters.n_compositional_fields,std::vector<double> (n_q_points));
    std::vector<double> composition_values_at_q_point (parameters.n_compositional_fields);

    double new_time_step;
    bool convection_dominant;

    double max_local_speed_over_meshsize = 0;
    double min_local_conduction_timestep = std::numeric_limits<double>::max(), min_conduction_timestep;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[introspection.extractors.velocities].get_function_values (solution,
                                                                              velocity_values);

          double max_local_velocity = 0;
          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
          max_local_speed_over_meshsize = std::max(max_local_speed_over_meshsize,
                                                   max_local_velocity
                                                   /
                                                   cell->minimum_vertex_distance());
          if (parameters.use_conduction_timestep)
            {
              fe_values[introspection.extractors.pressure].get_function_values (solution,
                                                                                pressure_values);
              fe_values[introspection.extractors.temperature].get_function_values (solution,
                                                                                   temperature_values);
              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                fe_values[introspection.extractors.compositional_fields[c]].get_function_values (solution,
                    composition_values[c]);

              typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points, parameters.n_compositional_fields);
              typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_q_points, parameters.n_compositional_fields);

              in.strain_rate.resize(0);// we are not reading the viscosity

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
                    composition_values_at_q_point[k] = composition_values[k][q];

                  in.position[q] = fe_values.quadrature_point(q);
                  in.temperature[q] = temperature_values[q];
                  in.pressure[q] = pressure_values[q];
                  for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                    in.composition[q][c] = composition_values_at_q_point[c];
                }
              in.cell = cell;

              material_model->evaluate(in, out);

              // In the future we may want to evaluate thermal diffusivity at
              // each point in the mesh, but for now we just use the reference value
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  double k = out.thermal_conductivities[q];
                  double rho = out.densities[q];
                  double c_p = out.specific_heat[q];

                  double thermal_diffusivity = k/(rho*c_p);

                  min_local_conduction_timestep = std::min(min_local_conduction_timestep,
                                                           parameters.CFL_number*pow(cell->minimum_vertex_distance(),2)
                                                           / thermal_diffusivity);
                }
            }

          if (parameters.include_melt_transport)
            {
              fe_values[introspection.extractors.pressure].get_function_values (solution,
                                                                                pressure_values);
              fe_values[introspection.extractors.compaction_pressure].get_function_values (solution,
                                                                                           fluid_pressure_values);
              fe_values[introspection.extractors.compaction_pressure].get_function_gradients (solution,
                                                                                   fluid_pressure_gradients);
              fe_values[introspection.extractors.temperature].get_function_values (solution,
                                                                                   temperature_values);
              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                fe_values[introspection.extractors.compositional_fields[c]].get_function_values (solution,
                    composition_values[c]);

              typename MaterialModel::MeltInterface<dim>::MaterialModelInputs in(n_q_points, parameters.n_compositional_fields);
              typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs out(n_q_points, parameters.n_compositional_fields);

              const unsigned int porosity_idx = introspection.compositional_index_for_name("porosity");

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
                    composition_values_at_q_point[k] = composition_values[k][q];

                  in.position[q] = fe_values.quadrature_point(q);
                  in.temperature[q] = temperature_values[q];
                  in.pressure[q] = pressure_values[q];
                  for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                    in.composition[q][c] = composition_values_at_q_point[c];
                }

              in.cell = cell;
              const typename MaterialModel::MeltInterface<dim> * melt_mat = dynamic_cast<const MaterialModel::MeltInterface<dim>*> (&*material_model);
              AssertThrow(melt_mat != NULL, ExcMessage("Need MeltMaterial if include_melt_transport is on."));
              melt_mat->evaluate_with_melt(in, out);

              // melt velocity = v_f =  v_s - K_D (nabla p_f - rho_f g) / porosity  or = v_s, if porosity = 0
              double max_local_melt_velocity = 0;
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  Tensor<1,dim> melt_velocity;
                  double porosity = in.composition[q][porosity_idx];

                  double p_s = in.pressure[q];
                  double p_f = fluid_pressure_values[q];

                  if (porosity < parameters.melt_transport_threshold)
                    melt_velocity = velocity_values[q];
                  else
                    {
                      double K_D = out.permeabilities[q] / out.fluid_viscosities[q];
                      Tensor<1,dim> grad_p_f = fluid_pressure_gradients[q];
                      const Tensor<1,dim> gravity = gravity_model->gravity_vector(in.position[q]);

                      // v_f =  v_s - K_D (nabla p_f - rho_f g) / porosity
                      melt_velocity = velocity_values[q]  - K_D * (grad_p_f - out.fluid_densities[q] * gravity) / porosity;
                    }
                  max_local_melt_velocity = std::max (max_local_melt_velocity,
                                                      melt_velocity.norm());

                }
              max_local_speed_over_meshsize = std::max(max_local_speed_over_meshsize,
                                                       max_local_melt_velocity
                                                       /
                                                       cell->minimum_vertex_distance());
            }

        }

    const double max_global_speed_over_meshsize
      = Utilities::MPI::max (max_local_speed_over_meshsize, mpi_communicator);
    if (parameters.use_conduction_timestep)
      MPI_Allreduce (&min_local_conduction_timestep, &min_conduction_timestep, 1, MPI_DOUBLE, MPI_MIN, mpi_communicator);
    else
      min_conduction_timestep = std::numeric_limits<double>::max();

    // if the velocity is zero, then it is somewhat arbitrary what time step
    // we should choose. in that case, do as if the velocity was one
    if (max_global_speed_over_meshsize == 0 && !parameters.use_conduction_timestep)
      {
        new_time_step = (parameters.CFL_number /
                         (parameters.temperature_degree * 1));
        convection_dominant = false;
      }
    else
      {
        new_time_step = std::min(min_conduction_timestep,
                                 (parameters.CFL_number / (parameters.temperature_degree * max_global_speed_over_meshsize)));
        convection_dominant = (new_time_step < min_conduction_timestep);
      }

    return std::make_pair(new_time_step, convection_dominant);
  }


  namespace
  {
    void
    extract_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                           const unsigned int                      q,
                                           std::vector<double>                    &composition_values_at_q_point)
    {
      for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
        composition_values_at_q_point[k] = composition_values[k][q];
    }
  }


  template <int dim>
  std::pair<double,double>
  Simulator<dim>::
  get_extrapolated_advection_field_range (const AdvectionField &advection_field) const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             (advection_field.is_temperature() ?
                                              parameters.temperature_degree :
                                              parameters.composition_degree));

    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    FEValues<dim> fe_values (mapping, finite_element, quadrature_formula,
                             update_values);
    std::vector<double> old_field_values(n_q_points);
    std::vector<double> old_old_field_values(n_q_points);

    // This presets the minimum with a bigger
    // and the maximum with a smaller number
    // than one that is going to appear. Will
    // be overwritten in the cell loop or in
    // the communication step at the
    // latest.
    double min_local_field = std::numeric_limits<double>::max(),
           max_local_field = -std::numeric_limits<double>::max();

    if (timestep_number > 1)
      {
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
                  const double extrapolated_field =
                    (1. + time_step/old_time_step) * old_field_values[q]-
                    time_step/old_time_step * old_old_field_values[q];

                  min_local_field = std::min (min_local_field,
                                              extrapolated_field);
                  max_local_field = std::max (max_local_field,
                                              extrapolated_field);
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
              fe_values[field].get_function_values (old_solution,
                                                    old_field_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double extrapolated_field = old_field_values[q];

                  min_local_field = std::min (min_local_field,
                                              extrapolated_field);
                  max_local_field = std::max (max_local_field,
                                              extrapolated_field);
                }
            }
      }

    return std::make_pair(-Utilities::MPI::max (-min_local_field,
                                                mpi_communicator),
                          Utilities::MPI::max (max_local_field,
                                               mpi_communicator));
  }


  template <int dim>
  void Simulator<dim>::interpolate_onto_velocity_system(const TensorFunction<1,dim> &func,
                                                        LinearAlgebra::Vector &vec)
  {
    ConstraintMatrix hanging_constraints(introspection.index_sets.system_relevant_set);
    DoFTools::make_hanging_node_constraints(dof_handler, hanging_constraints);
    hanging_constraints.close();

    Assert(introspection.block_indices.velocities == 0, ExcNotImplemented());
    const std::vector<Point<dim> > mesh_support_points = finite_element.base_element(introspection.base_elements.velocities).get_unit_support_points();
    FEValues<dim> mesh_points (mapping, finite_element, mesh_support_points, update_quadrature_points);
    std::vector<unsigned int> cell_dof_indices (finite_element.dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          mesh_points.reinit(cell);
          cell->get_dof_indices (cell_dof_indices);
          for (unsigned int j=0; j<finite_element.base_element(introspection.base_elements.velocities).dofs_per_cell; ++j)
            for (unsigned int dir=0; dir<dim; ++dir)
              {
                unsigned int support_point_index
                  = finite_element.component_to_system_index(/*velocity component=*/ introspection.component_indices.velocities[dir],
                                                                                     /*dof index within component=*/ j);
                Assert(introspection.block_indices.velocities == 0, ExcNotImplemented());
                vec[cell_dof_indices[support_point_index]] = func.value(mesh_points.quadrature_point(j))[dir];
              }
        }

    vec.compress(VectorOperation::insert);
    hanging_constraints.distribute(vec);
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
                      fe_face_values[introspection.extractors.pressure].get_function_values(vector,
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
        const QGauss<dim> quadrature (parameters.stokes_velocity_degree + 1);

        const unsigned int n_q_points = quadrature.size();
        FEValues<dim> fe_values (mapping, finite_element,  quadrature,
                                 update_JxW_values | update_values);

        std::vector<double> pressure_values(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[introspection.extractors.pressure].get_function_values(vector,
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

    pressure_adjustment = 0.0;
    // sum up the integrals from each processor
    {
      const double my_temp[2] = {my_pressure, my_area};
      double temp[2];
      Utilities::MPI::sum (my_temp, mpi_communicator, temp);
      if (parameters.pressure_normalization == "surface")
        pressure_adjustment = -temp[0]/temp[1] + parameters.surface_pressure;
      else if (parameters.pressure_normalization == "volume")
        pressure_adjustment = -temp[0]/temp[1];
      else
        AssertThrow(false, ExcNotImplemented());
    }

    // A complication is that we can't modify individual
    // elements of the solution vector since that one has ghost element.
    // rather, we first need to localize it and then distribute back
    LinearAlgebra::BlockVector distributed_vector (introspection.index_sets.system_partitioning,
                                                   mpi_communicator);
    distributed_vector = vector;

    if (parameters.use_locally_conservative_discretization == false)
      {
        if ((introspection.block_indices.velocities != introspection.block_indices.pressure)
            &&
            (introspection.component_indices.compaction_pressure == numbers::invalid_unsigned_int))
          distributed_vector.block(introspection.block_indices.pressure).add(pressure_adjustment);
        else
          {
            // velocity and pressure are in the same block, so we have to modify the values manually
            std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
            typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
            for (; cell != endc; ++cell)
              if (cell->is_locally_owned())
                {
                  cell->get_dof_indices (local_dof_indices);
                  for (unsigned int j=0; j<finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell; ++j)
                    {
                      unsigned int support_point_index
                        = finite_element.component_to_system_index(introspection.component_indices.pressure,
                                                                   /*dof index within component=*/ j);

                      Assert (introspection.block_indices.velocities == introspection.block_indices.pressure
                              || local_dof_indices[support_point_index] >= vector.block(0).size(),
                              ExcInternalError());

                      // then adjust its value. Note that because we end up touching
                      // entries more than once, we are not simply incrementing
                      // distributed_vector but copy from the unchanged vector.
                      distributed_vector(local_dof_indices[support_point_index]) = vector(local_dof_indices[support_point_index]) + pressure_adjustment;
                    }
                }
            distributed_vector.compress(VectorOperation::insert);
          }
      }
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
        Assert (dynamic_cast<const FE_DGP<dim>*>(&finite_element.base_element(introspection.base_elements.pressure)) != 0,
                ExcInternalError());
        std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = finite_element.component_to_system_index (introspection.component_indices.pressure, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (dof_handler.locally_owned_dofs().is_element(local_dof_indices[first_pressure_dof]),
                      ExcInternalError());

              Assert (introspection.block_indices.velocities == introspection.block_indices.pressure
                      || local_dof_indices[first_pressure_dof] >= vector.block(0).size(),
                      ExcInternalError());

              // then adjust its value
              distributed_vector(local_dof_indices[first_pressure_dof]) += pressure_adjustment;
            }
        distributed_vector.compress(VectorOperation::insert);
      }

    // now get back to the original vector
    vector = distributed_vector;
  }


  /*
   * inverse to normalize_pressure.
   */
  template <int dim>
  void Simulator<dim>::denormalize_pressure (LinearAlgebra::BlockVector &vector, const LinearAlgebra::BlockVector &relevant_vector)
  {
    if (parameters.pressure_normalization == "no")
      return;

    if (parameters.use_locally_conservative_discretization == false)
      {
        if ((introspection.block_indices.velocities != introspection.block_indices.pressure)
            &&
            (introspection.component_indices.compaction_pressure == numbers::invalid_unsigned_int))
          vector.block(introspection.block_indices.pressure).add(-1.0 * pressure_adjustment);
        else
          {
            // velocity and pressure are in the same block, or we have several pressures in one block,
            // so we have to modify the values manually

            const unsigned int block_idx = introspection.block_indices.pressure;

            std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
            typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
            for (; cell != endc; ++cell)
              if (cell->is_locally_owned())
                {
                  cell->get_dof_indices (local_dof_indices);
                  for (unsigned int j=0; j<finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell; ++j)
                    {
                      const unsigned int local_dof_index
                        = finite_element.component_to_system_index(introspection.component_indices.pressure,
                                                                   /*dof index within component=*/ j);

                      // then adjust its value. Note that because we end up touching
                      // entries more than once, we are not simply incrementing
                      // distributed_vector but copy from the unchanged vector.
                      vector(local_dof_indices[local_dof_index])
                        = relevant_vector(local_dof_indices[local_dof_index]) - pressure_adjustment;
                    }
                }
            vector.compress(VectorOperation::insert);
          }
      }
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
        Assert (dynamic_cast<const FE_DGP<dim>*>(&finite_element.base_element(introspection.base_elements.pressure)) != 0,
                ExcInternalError());
        std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
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
              Assert (dof_handler.locally_owned_dofs().is_element(local_dof_indices[first_pressure_dof]),
                      ExcInternalError());
              Assert (local_dof_indices[first_pressure_dof] >= vector.block(0).size(),
                      ExcInternalError());

              // then adjust its value
              vector (local_dof_indices[first_pressure_dof]) -= pressure_adjustment;
            }

        vector.compress(VectorOperation::insert);
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

    // In the following we integrate the normal velocity over every surface
    // of the model. This integral is part of the correction term that needs
    // to be added to the pressure right hand side. To calculate the normal
    // velocity we need the positions and normals of every quadrature point on
    // the surface.

    const QGauss<dim-1> quadrature_formula (parameters.stokes_velocity_degree+1);
    FEFaceValues<dim> fe_face_values (mapping,
                                      finite_element,
                                      quadrature_formula,
                                      update_normal_vectors |
                                      update_q_points |
                                      update_JxW_values);

    double local_normal_velocity_integral = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    // for every surface face that is part of a geometry boundary with
    // prescribed velocity and that is owned by this processor,
    // integrate the normal velocity magnitude.
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          {
            if (!cell->face(f)->at_boundary())
              continue;

            fe_face_values.reinit (cell, f);

            if (parameters.prescribed_velocity_boundary_indicators.find(cell->face(f)->boundary_indicator())!=
                parameters.prescribed_velocity_boundary_indicators.end())
              {
                                // for each of the quadrature points, evaluate the
                                // normal velocity by calling the boundary conditions and add
                                // it to the integral
                                for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                                      {
                                        const Tensor<1,dim> velocity =
                                              velocity_boundary_conditions.find(cell->face(f)->boundary_indicator())->second->boundary_velocity(fe_face_values.quadrature_point(q));

                                        const double normal_velocity = velocity * fe_face_values.normal_vector(q);
                                        local_normal_velocity_integral += normal_velocity * fe_face_values.JxW(q);
                                      }
               }

            // correct u by -K_D grad p_f = -K_D f
            if (parameters.include_melt_transport)
              {
                typename MaterialModel::MeltInterface<dim>::MaterialModelInputs melt_inputs(quadrature_formula.size(), parameters.n_compositional_fields);
                typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs melt_outputs(quadrature_formula.size(), parameters.n_compositional_fields);

                compute_material_model_input_values (current_linearization_point,
                                                     fe_face_values,
                                                     rebuild_stokes_matrix,
                                                     melt_inputs);
                // TODO: fill other inputs
                melt_inputs.cell = cell;
                typename MaterialModel::MeltInterface<dim> * melt_mat = dynamic_cast<MaterialModel::MeltInterface<dim>*> (&*material_model);
                AssertThrow(melt_mat != NULL, ExcMessage("Need MeltMaterial if include_melt_transport is on."));
                melt_mat->evaluate_with_melt(melt_inputs, melt_outputs);

                std::vector<Tensor<1,dim> > grad_p_f(quadrature_formula.size());
                fluid_pressure_boundary_conditions->fluid_pressure_gradient(melt_inputs,
                    melt_outputs,
                    grad_p_f);

                for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                  {
                    const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");
                    const double porosity = std::max(melt_inputs.composition[q][porosity_index],0.000);
                    const double K_D = (porosity > parameters.melt_transport_threshold
                                      ?
                                      melt_outputs.permeabilities[q] / melt_outputs.fluid_viscosities[q]
                                      :
                                      0.0);
                    const double K_D_grad_p_f = - K_D * grad_p_f[q] * fe_face_values.normal_vector(q);
                    local_normal_velocity_integral += K_D_grad_p_f * fe_face_values.JxW(q);
                  }
                }
            }

    const double global_normal_velocity_integral =
      Utilities::MPI::sum (local_normal_velocity_integral,mpi_communicator);

    if (!parameters.include_melt_transport && introspection.block_indices.velocities != introspection.block_indices.pressure)
      {
        const double mean       = vector.block(introspection.block_indices.pressure).mean_value();
        const double correction = (global_normal_velocity_integral - mean * vector.block(introspection.block_indices.pressure).size()) / global_volume;

        vector.block(introspection.block_indices.pressure).add(correction, pressure_shape_function_integrals.block(introspection.block_indices.pressure));
      }
    else
      {
        // we need to operate only on p_f not on p_c
        double pressure_sum = 0.0;

        for (unsigned int i=0; i < introspection.index_sets.locally_owned_fluid_pressure_dofs.n_elements(); ++i)
           {
             types::global_dof_index idx = introspection.index_sets.locally_owned_fluid_pressure_dofs.nth_index_in_set(i);
             pressure_sum += vector(idx);
           }

        const double global_pressure_sum = Utilities::MPI::sum(pressure_sum, mpi_communicator);
        const double correction = (global_normal_velocity_integral - global_pressure_sum) / global_volume;

        for (unsigned int i=0; i < introspection.index_sets.locally_owned_fluid_pressure_dofs.n_elements(); ++i)
           {
             types::global_dof_index idx = introspection.index_sets.locally_owned_fluid_pressure_dofs.nth_index_in_set(i);
             vector(idx) += correction * pressure_shape_function_integrals(idx);
           }

        vector.compress(VectorOperation::insert);
      }
  }


  template <int dim>
  void Simulator<dim>::compute_melt_variables(const LinearAlgebra::BlockVector &input_solution,
                                               LinearAlgebra::BlockVector       &output_solution)
    {
      if (!parameters.include_melt_transport)
        return;

      const unsigned int por_idx = introspection.compositional_index_for_name("porosity");


      //   * fluid_velocity
{
      // u_f =  u_s - K_D (nabla p_f - rho_f g) / phi  or = 0

          const Quadrature<dim> quadrature(finite_element.base_element(introspection.base_elements.pressure).get_unit_support_points());

          typename MaterialModel::MeltInterface<dim>::MaterialModelInputs in(quadrature.size(), parameters.n_compositional_fields);
          typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs out(quadrature.size(), parameters.n_compositional_fields);

          std::vector<double> porosity_values(quadrature.size());
          std::vector<Tensor<1,dim> > grad_p_f_values(quadrature.size());
          std::vector<Tensor<1,dim> > u_s_values(quadrature.size());

          FEValues<dim> fe_values (mapping,
                                   finite_element,
                                   quadrature,
                                   update_quadrature_points | update_values);

          std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
          typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);
              cell->get_dof_indices (local_dof_indices);
              fe_values[introspection.extractors.compositional_fields[por_idx]].get_function_values (
                  input_solution, porosity_values);
              fe_values[introspection.extractors.velocities].get_function_values (
                  input_solution, u_s_values);
              fe_values[introspection.extractors.fluid_pressure].get_function_gradients (
                  input_solution, grad_p_f_values);

              compute_material_model_input_values (input_solution,
                                       fe_values,
                                       true, // TODO: use rebuild_stokes_matrix here?
                                       in);
              in.cell = cell;
              typename MaterialModel::MeltInterface<dim> * melt_mat = dynamic_cast<MaterialModel::MeltInterface<dim>*> (&*material_model);
              melt_mat->evaluate_with_melt(in, out);

              for (unsigned int j=0;j < finite_element.dofs_per_cell; ++j)
              {
                  std::pair<unsigned int, unsigned int> base_index
                   = finite_element.system_to_base_index(j).first;
                  if (base_index.first != introspection.base_elements.fluid_velocities)
                      continue;

                  const unsigned int d = base_index.second;

                  // skip entries that are not locally owned:
                  if (!dof_handler.locally_owned_dofs().is_element(local_dof_indices[j]))
                    continue;

                  const double phi = std::max(0.0, porosity_values[j]);

                  double value = 0.0;
                  // u_f =  u_s - K_D (nabla p_f - rho_f g) / phi  or = 0
                  if (phi > parameters.melt_transport_threshold)
                  {
                          const double K_D = out.permeabilities[j] / out.fluid_viscosities[j];
                          const double gravity_d = this->gravity_model->gravity_vector(in.position[j])[d];

                          // v_f =  v_s - K_D (nabla p_f - rho_f g) / phi
                          value = u_s_values[j][d] - K_D * (grad_p_f_values[j][d] - out.fluid_densities[j] * gravity_d) / phi;
                    }

                  output_solution(local_dof_indices[j]) = value;
                }
            }
        output_solution.block(introspection.block_indices.fluid_velocities).compress(VectorOperation::insert);
  }

      //* pressure

      {

      const unsigned int block_p = introspection.block_indices.pressure;

      // Think what we need to do if the pressure is not an FE_Q...
      Assert(parameters.use_locally_conservative_discretization == false, ExcNotImplemented());

      const Quadrature<dim> quadrature(finite_element.base_element(introspection.base_elements.pressure).get_unit_support_points());
      std::vector<double> porosity_values(quadrature.size());
      std::vector<double> p_c_values(quadrature.size());
      std::vector<double> p_f_values(quadrature.size());
      FEValues<dim> fe_values (mapping,
                               finite_element,
                               quadrature,
                               update_quadrature_points | update_values);

      std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            cell->get_dof_indices (local_dof_indices);
            fe_values[introspection.extractors.compositional_fields[por_idx]].get_function_values (
                current_linearization_point, porosity_values);
            fe_values[introspection.extractors.compaction_pressure].get_function_values (
                current_linearization_point, p_c_values);
            fe_values[introspection.extractors.fluid_pressure].get_function_values (
                current_linearization_point, p_f_values);

            for (unsigned int j=0; j<finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell; ++j)
              {
                const unsigned int pressure_idx
                = finite_element.component_to_system_index(introspection.component_indices.pressure,
                    /*dof index within component=*/ j);

                // skip entries that are not locally owned:
                if (!dof_handler.locally_owned_dofs().is_element(local_dof_indices[pressure_idx]))
                  continue;

                const double phi = std::max(0.0, porosity_values[j]);

                double p = p_f_values[j];
                if (phi < 1.0-parameters.melt_transport_threshold)
                    p = (p_c_values[j] - (phi-1.0) * p_f_values[j]) / (1.0-phi);

                output_solution(local_dof_indices[pressure_idx]) = p;
              }
          }
      output_solution.block(block_p).compress(VectorOperation::insert);
      }
    }

  template <int dim>
  void Simulator<dim>::convert_pressure_blocks(const LinearAlgebra::BlockVector &input_solution,
                                               const bool                       solid_to_fluid_pressure,
                                               LinearAlgebra::BlockVector       &output_solution)
  {
    if (!parameters.include_melt_transport)
      return;

    // for the direct solver we have to copy the whole block,
    // because the velocity is included as well.
    const unsigned int block_p = introspection.block_indices.pressure;
    output_solution.block(block_p) = input_solution.block(block_p);

    // Think what we need to do if the pressure is not an FE_Q...
    Assert(parameters.use_locally_conservative_discretization == false, ExcNotImplemented());

    const unsigned int por_idx = introspection.compositional_index_for_name("porosity");
    const Quadrature<dim> quadrature(finite_element.base_element(introspection.base_elements.pressure).get_unit_support_points());
    std::vector<double> porosity_values(quadrature.size());
    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature,
                             update_quadrature_points | update_values);

    std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell->get_dof_indices (local_dof_indices);
          fe_values[introspection.extractors.compositional_fields[por_idx]].get_function_values (
              current_linearization_point, porosity_values);

          for (unsigned int j=0; j<finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell; ++j)
            {
              const unsigned int pressure_idx
              = finite_element.component_to_system_index(introspection.component_indices.pressure,
                  /*dof index within component=*/ j);

              // skip entries that are not locally owned:
              if (!dof_handler.locally_owned_dofs().is_element(local_dof_indices[pressure_idx]))
                continue;

              const unsigned int p_c_idx
              = finite_element.component_to_system_index(introspection.component_indices.compaction_pressure,
                  /*dof index within component=*/ j);

              double phi = porosity_values[j];

              if (solid_to_fluid_pressure)
                // (p_s, p_f) -> (p_f, p_c)
                {
                  double p_s = input_solution(local_dof_indices[pressure_idx]);
                  double p_f = input_solution(local_dof_indices[p_c_idx]);
                  double p_c = (1-phi)*(p_s-p_f);

                  output_solution(local_dof_indices[pressure_idx]) = p_f;
                  output_solution(local_dof_indices[p_c_idx]) = p_c;
                }
              else
                // (p_f, p_c) -> (p_s, p_f)
                {
                  double p_f = input_solution(local_dof_indices[pressure_idx]);
                  double p_s;
                  if (phi >(1.0-parameters.melt_transport_threshold) || (phi <= 0.0))
                    p_s = p_f;
                  else
                  {
                    double p_c = input_solution(local_dof_indices[p_c_idx]);
                    p_s = (p_c - (phi-1.0) * p_f) / (1.0-phi);
                  }

                  output_solution(local_dof_indices[pressure_idx]) = p_s;
                  output_solution(local_dof_indices[p_c_idx]) = p_f;
                }
            }
        }
    output_solution.block(block_p).compress(VectorOperation::insert);
  }


  template <int dim>
  double
  Simulator<dim>::compute_initial_stokes_residual()
  {
    LinearAlgebra::BlockVector remap (introspection.index_sets.stokes_partitioning, mpi_communicator);
    LinearAlgebra::BlockVector residual (introspection.index_sets.stokes_partitioning, mpi_communicator);
    const unsigned int block_p = introspection.block_indices.pressure;

    // if velocity and pressure are in the same block, we have to copy the
    // pressure to the solution and RHS vector with a zero velocity
    if(introspection.block_indices.pressure == introspection.block_indices.velocities)
      {
        const IndexSet & idxset = (parameters.include_melt_transport) ?
            introspection.index_sets.locally_owned_fluid_pressure_dofs
            :
            introspection.index_sets.locally_owned_pressure_dofs;

        for (unsigned int i=0; i < idxset.n_elements(); ++i)
          {
            types::global_dof_index idx = idxset.nth_index_in_set(i);
            remap(idx)        = current_linearization_point(idx);
          }
        remap.block(0).compress(VectorOperation::insert);
      }
    else
      remap.block (block_p) = current_linearization_point.block (block_p);

    // we have to do the same conversions and rescaling we do before solving
    // the Stokes system: conversion of the pressure blocks (p_s, p_f --> p_f, p_c),
    // denormalizing the pressure and applying the pressure scaling
    residual = remap;
    if (parameters.include_melt_transport)
      convert_pressure_blocks(residual, true, remap);
    residual = remap;

    // TODO: we don't have .stokes_relevant_partitioning so I am creating a much
    // bigger vector here, oh well.
    LinearAlgebra::BlockVector ghosted (introspection.index_sets.system_partitioning,
        introspection.index_sets.system_relevant_partitioning,
        mpi_communicator);
    // TODO for Timo: can we create the ghost vector inside of denormalize_pressure
    // (only in cases where we need it)
    ghosted.block(block_p) = remap.block(block_p);
    denormalize_pressure (remap, ghosted);
    current_constraints.set_zero (remap);

    remap.block (block_p) /= pressure_scaling;
    if (do_pressure_rhs_compatibility_modification)
      make_pressure_rhs_compatible(system_rhs);

    // we calculate the velocity residual with a zero velocity,
    // computing only the part of the RHS not balanced by the static pressure
    if(introspection.block_indices.pressure == introspection.block_indices.velocities)
      {
        // we can use the whole block here because we set the velocity to zero above
        return system_matrix.block(0,0).residual (residual.block(0),
                                                  remap.block(0),
                                                  system_rhs.block(0));
      }
    else
      {
        double residual_u = system_matrix.block(0,1).residual (residual.block(0),
                                                        remap.block(1),
                                                        system_rhs.block(0));
        double residual_p = system_rhs.block(block_p).l2_norm();
        return sqrt(residual_u*residual_u+residual_p*residual_p);
      }
  }


  template <int dim>
  template<class FUNCTOR>
  void Simulator<dim>::compute_depth_average(std::vector<double> &values,
                                             FUNCTOR &fctr) const
  {
    Assert (values.size() > 0,
            ExcMessage ("To call this function, you need to request a positive "
                        "number of depth slices."));
    const unsigned int num_slices = values.size();
    std::vector<double> volume(num_slices);

    // this yields 10^dim quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();
    const double max_depth = geometry_model->maximal_depth();

    FEValues<dim> fe_values (mapping,
                             finite_element,
                             quadrature_formula,
                             update_values | update_gradients | update_quadrature_points);

    std::vector<std::vector<double> > composition_values (parameters.n_compositional_fields,std::vector<double> (n_q_points));

    std::vector<double> output_values(quadrature_formula.size());

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points,
                                                                   parameters.n_compositional_fields);
    typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_q_points,
                                                                     parameters.n_compositional_fields);

    fctr.setup(quadrature_formula.size());

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          if (fctr.need_material_properties())
            {
              fe_values[introspection.extractors.pressure].get_function_values (this->solution,
                                                                                in.pressure);
              fe_values[introspection.extractors.temperature].get_function_values (this->solution,
                                                                                   in.temperature);
              fe_values[introspection.extractors.velocities].get_function_symmetric_gradients (this->solution,
                  in.strain_rate);
              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                fe_values[introspection.extractors.compositional_fields[c]].get_function_values(this->solution,
                                                                                                composition_values[c]);

              for (unsigned int i=0; i<n_q_points; ++i)
                {
                  in.position[i] = fe_values.quadrature_point(i);
                  for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                    in.composition[i][c] = composition_values[c][i];
                }
              in.cell = cell;

              material_model->evaluate(in, out);
            }

          fctr(in, out, fe_values, this->solution, output_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double depth = geometry_model->depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              values[idx] += output_values[q] * fe_values.JxW(q);
              volume[idx] += fe_values.JxW(q);
            }
        }

    std::vector<double> values_all(num_slices);
    std::vector<double> volume_all(num_slices);
    Utilities::MPI::sum(volume, mpi_communicator, volume_all);
    Utilities::MPI::sum(values, mpi_communicator, values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(volume_all[i])+1e-20);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageField
    {
      public:
        FunctorDepthAverageField(const FEValuesExtractors::Scalar &field)
          : field_(field) {}

        bool need_material_properties() const
        {
          return false;
        }

        void setup(const unsigned int q_points)
        {
        }

        void operator()(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &out,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution, std::vector<double> &output)
        {
          fe_values[field_].get_function_values (solution, output);
        }

        const FEValuesExtractors::Scalar field_;
    };
  }

  template <int dim>
  void Simulator<dim>::compute_depth_average_field(const AdvectionField &advection_field,
                                                   std::vector<double> &values) const
  {
    const FEValuesExtractors::Scalar field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    FunctorDepthAverageField<dim> f(field);
    compute_depth_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageViscosity
    {
      public:
        bool need_material_properties() const
        {
          return true;
        }

        void setup(const unsigned int q_points)
        {}

        void operator()(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &out,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution, std::vector<double> &output)
        {
          output = out.viscosities;
        }
    };
  }

  template <int dim>
  void Simulator<dim>::compute_depth_average_viscosity(std::vector<double> &values) const
  {
    FunctorDepthAverageViscosity<dim> f;
    compute_depth_average(values, f);
  }


  namespace
  {
    template <int dim>
    class FunctorDepthAverageVelocityMagnitude
    {
      public:
        FunctorDepthAverageVelocityMagnitude(const FEValuesExtractors::Vector &field)
          : field_(field) {}

        bool need_material_properties() const
        {
          return false;
        }

        void setup(const unsigned int q_points)
        {
          velocity_values.resize(q_points);
        }

        void operator()(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &out,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution, std::vector<double> &output)
        {
          fe_values[field_].get_function_values (solution, velocity_values);
          for (unsigned int q=0; q<output.size(); ++q)
            output[q] = velocity_values[q] * velocity_values[q];
        }

        std::vector<Tensor<1,dim> > velocity_values;
        const FEValuesExtractors::Vector field_;
    };
  }


  template <int dim>
  void Simulator<dim>::compute_depth_average_velocity_magnitude(std::vector<double> &values) const
  {
    FunctorDepthAverageVelocityMagnitude<dim> f(introspection.extractors.velocities);

    compute_depth_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageSinkingVelocity
    {
      public:
        FunctorDepthAverageSinkingVelocity(const FEValuesExtractors::Vector &field, GravityModel::Interface<dim> *gravity)
          : field_(field), gravity_(gravity) {}

        bool need_material_properties() const
        {
          return false;
        }

        void setup(const unsigned int q_points)
        {
          velocity_values.resize(q_points);
        }

        void operator()(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &out,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution, std::vector<double> &output)
        {
          fe_values[field_].get_function_values (solution, velocity_values);
          for (unsigned int q=0; q<output.size(); ++q)
            {
              Tensor<1,dim> g = gravity_->gravity_vector(in.position[q]);
              output[q] = std::fabs(std::min(-1e-16,g*velocity_values[q]/g.norm()))*year_in_seconds;
            }
        }

        std::vector<Tensor<1,dim> > velocity_values;
        const FEValuesExtractors::Vector field_;
        GravityModel::Interface<dim> *gravity_;
    };
  }


  template <int dim>
  void Simulator<dim>::compute_depth_average_sinking_velocity(std::vector<double> &values) const
  {
    FunctorDepthAverageSinkingVelocity<dim> f(introspection.extractors.velocities,
                                              this->gravity_model.get());

    compute_depth_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageVsVp
    {
      public:
        FunctorDepthAverageVsVp(const MaterialModel::Interface<dim> *mm, bool vs)
          : material_model(mm), vs_(vs)
        {}

        bool need_material_properties() const
        {
          return true;
        }

        void setup(const unsigned int q_points)
        {}

        void operator()(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &out,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution, std::vector<double> &output)
        {
          if (vs_)
            for (unsigned int q=0; q<output.size(); ++q)
              output[q] = material_model->seismic_Vs(
                            in.temperature[q], in.pressure[q], in.composition[q],
                            in.position[q]);
          else
            for (unsigned int q=0; q<output.size(); ++q)
              output[q] = material_model->seismic_Vp(
                            in.temperature[q], in.pressure[q], in.composition[q],
                            in.position[q]);
        }

        const MaterialModel::Interface<dim> *material_model;
        bool vs_;
    };

  }

  template <int dim>
  void Simulator<dim>::compute_depth_average_Vs(std::vector<double> &values) const
  {
    FunctorDepthAverageVsVp<dim> f(this->material_model.get(), true /* Vs */);

    compute_depth_average(values, f);
  }

  template <int dim>
  void Simulator<dim>::compute_depth_average_Vp(std::vector<double> &values) const
  {
    FunctorDepthAverageVsVp<dim> f(this->material_model.get(), false /* Vp */);

    compute_depth_average(values, f);
  }



  template <int dim>
  bool
  Simulator<dim>::stokes_matrix_depends_on_solution() const
  {
    // currently, the only coefficient that really appears on the
    // left hand side of the Stokes equation is the viscosity. note
    // that our implementation of compressible materials makes sure
    // that the density does not appear on the lhs.
    // if melt transport is included in the simulation, we have an
    // additional equation with more coefficients on the left hand
    // side.
    return (material_model->viscosity_depends_on (MaterialModel::NonlinearDependence::any_variable)
            || parameters.include_melt_transport);
  }
}

// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Simulator<dim>::AdvectionField; \
  template void Simulator<dim>::normalize_pressure(LinearAlgebra::BlockVector &vector); \
  template void Simulator<dim>::denormalize_pressure(LinearAlgebra::BlockVector &vector, const LinearAlgebra::BlockVector &relevant_vector); \
  template double Simulator<dim>::get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const; \
  template std::pair<double,double> Simulator<dim>::get_extrapolated_advection_field_range (const AdvectionField &advection_field) const; \
  template std::pair<double,bool> Simulator<dim>::compute_time_step () const; \
  template void Simulator<dim>::make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector); \
  template void Simulator<dim>::compute_depth_average_field(const AdvectionField &advection_field, std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_viscosity(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_velocity_magnitude(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_sinking_velocity(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_Vs(std::vector<double> &values) const; \
  template void Simulator<dim>::compute_depth_average_Vp(std::vector<double> &values) const; \
  template void Simulator<dim>::output_program_stats(); \
  template void Simulator<dim>::output_statistics(); \
  template void Simulator<dim>::convert_pressure_blocks(const LinearAlgebra::BlockVector &input_solution, const bool solid_to_fluid_pressure, LinearAlgebra::BlockVector &output_solution); \
  template void Simulator<dim>::compute_melt_variables(const LinearAlgebra::BlockVector &input_solution,\
                                                 LinearAlgebra::BlockVector       &output_solution); \
  template double Simulator<dim>::compute_initial_stokes_residual(); \
  template bool Simulator<dim>::stokes_matrix_depends_on_solution() const; \
  template void Simulator<dim>::interpolate_onto_velocity_system(const TensorFunction<1,dim> &func, LinearAlgebra::Vector &vec);

  ASPECT_INSTANTIATE(INSTANTIATE)
}

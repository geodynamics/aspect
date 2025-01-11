/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/melt.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/newton.h>
#include <aspect/global.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/heating_model/interface.h>
#include <aspect/heating_model/adiabatic_heating.h>
#include <aspect/material_model/interface.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/particle/generator/interface.h>
#include <aspect/particle/integrator/interface.h>
#include <aspect/particle/interpolator/interface.h>
#include <aspect/particle/property/interface.h>
#include <aspect/postprocess/visualization.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/sundials/arkode.h>

#include <fstream>
#include <iostream>
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
  Simulator<dim>::AdvectionField::is_discontinuous(const Introspection<dim> &introspection) const
  {
    if (field_type == temperature_field)
      return introspection.use_discontinuous_temperature_discretization;
    else if (field_type == compositional_field)
      return introspection.use_discontinuous_composition_discretization[compositional_variable];

    Assert (false, ExcInternalError());
    return false;
  }

  template <int dim>
  typename Parameters<dim>::AdvectionFieldMethod::Kind
  Simulator<dim>::AdvectionField::advection_method(const Introspection<dim> &introspection) const
  {
    if (field_type == temperature_field)
      return introspection.temperature_method;
    else if (field_type == compositional_field)
      return introspection.compositional_field_methods[compositional_variable];

    Assert (false, ExcInternalError());
    return Parameters<dim>::AdvectionFieldMethod::fem_field;
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
  Simulator<dim>::AdvectionField::sparsity_pattern_block_index(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.block_indices.temperature;
    else
      return introspection.block_indices.compositional_field_sparsity_pattern[compositional_variable];
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
  Simulator<dim>::AdvectionField::field_index() const
  {
    if (this->is_temperature())
      return 0;
    else
      return compositional_variable + 1;
  }

  template <int dim>
  unsigned int
  Simulator<dim>::AdvectionField::base_element(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.base_elements.temperature;
    else
      return introspection.base_elements.compositional_fields[compositional_variable];
  }

  template <int dim>
  FEValuesExtractors::Scalar
  Simulator<dim>::AdvectionField::scalar_extractor(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.extractors.temperature;
    else
      {
        Assert(compositional_variable < introspection.n_compositional_fields,
               ExcMessage("Invalid AdvectionField."));
        return introspection.extractors.compositional_fields[compositional_variable];
      }
  }

  template <int dim>
  unsigned int
  Simulator<dim>::AdvectionField::polynomial_degree(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.polynomial_degree.temperature;
    else
      return introspection.polynomial_degree.compositional_fields[compositional_variable];
  }



  template <int dim>
  std::string
  Simulator<dim>::AdvectionField::name(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return "temperature";
    else
      return "composition " + std::to_string(compositional_variable) + " (" + introspection.name_for_compositional_index(compositional_variable) + ")";
  }



  template <int dim>
  void Simulator<dim>::write_plugin_graph (std::ostream &out) const
  {
    // write the preamble
    out << "digraph Plugins\n"
        "{\n"
        "  splines=line;\n"
        "  splines=true;\n"
        "  overlap=false;\n"
        "  edge [fontname=\"FreeSans\",\n"
        "        fontsize=\"10\",\n"
        "        labelfontname=\"FreeSans\",\n"
        "        labelfontsize=\"10\",\n"
        "        color=\"black\",\n"
        "        style=\"solid\"];\n"
        "  node [fontname=\"FreeSans\",\n"
        "        fontsize=\"10\",\n"
        "        shape=\"rectangle\",\n"
        "        height=0.2,\n"
        "        width=0.4,\n"
        "        color=\"black\",\n"
        "        fillcolor=\"white\",\n"
        "        style=\"filled\"];\n"
        "  layout=neato;\n"
        "\n";

    // then also write nodes for the Simulator and SimulatorAccess classes,
    // and an arrow from the former to the latter to indicate flow of
    // information
    out << "  Simulator [height=1.5,width=2,shape=\"octagon\",fillcolor=\"yellow\"];\n";
    out << "  SimulatorAccess [height=1.2,width=1.2,shape=\"rect\",fillcolor=\"yellow\"];\n";
    out << "  Simulator -> SimulatorAccess [len=1, weight=100];\n";

    // then go through all plugin systems and output everything we have
    AdiabaticConditions::write_plugin_graph<dim>(out);
    BoundaryComposition::Manager<dim>::write_plugin_graph(out);
    BoundaryFluidPressure::write_plugin_graph<dim>(out);
    BoundaryTemperature::Manager<dim>::write_plugin_graph(out);
    BoundaryTraction::Manager<dim>::write_plugin_graph(out);
    BoundaryVelocity::Manager<dim>::write_plugin_graph(out);
    InitialTopographyModel::write_plugin_graph<dim>(out);
    GeometryModel::write_plugin_graph<dim>(out);
    GravityModel::write_plugin_graph<dim>(out);
    HeatingModel::Manager<dim>::write_plugin_graph(out);
    InitialComposition::Manager<dim>::write_plugin_graph(out);
    InitialTemperature::Manager<dim>::write_plugin_graph(out);
    MaterialModel::write_plugin_graph<dim>(out);
    MeshRefinement::Manager<dim>::write_plugin_graph(out);
    Particle::Generator::write_plugin_graph<dim>(out);
    Particle::Integrator::write_plugin_graph<dim>(out);
    Particle::Interpolator::write_plugin_graph<dim>(out);
    Particle::Property::Manager<dim>::write_plugin_graph(out);
    Postprocess::Manager<dim>::write_plugin_graph(out);
    Postprocess::Visualization<dim>::write_plugin_graph(out);
    PrescribedStokesSolution::write_plugin_graph<dim>(out);
    TerminationCriteria::Manager<dim>::write_plugin_graph(out);

    // end the graph
    out << '}'
        << std::endl;
  }



  template <int dim>
  void Simulator<dim>::output_statistics()
  {
    // Only write the statistics file from processor zero
    if (Utilities::MPI::this_mpi_process(mpi_communicator)!=0)
      return;

    // Formatting the table we're about to output and writing the
    // actual file may take some time, so do it on a separate
    // thread. We do this using a lambda function that takes
    // a copy of the statistics object to make sure that whatever
    // we do to the 'real' statistics object at the time of
    // writing data doesn't affect what we write.
    //
    // Before we can start working on a new thread, we need to
    // make sure that the previous thread is done or they'll
    // step on each other's feet.
    if (output_statistics_thread.joinable())
      output_statistics_thread.join();

    // Write data in the background through a lambda function.
    // This happening in the background means that we have
    // to create a copy of the statistics table, since whatever is
    // running in the foreground may continue to add entries to the
    // statistics table at the same time.
    output_statistics_thread = std::thread (
                                 [statistics_copy_ptr = std::make_unique<TableHandler>(statistics),this]()
    {
      // First write everything into a string in memory
      std::ostringstream stream;
      statistics_copy_ptr->write_text (stream,
                                       TableHandler::table_with_separate_column_description);
      stream.flush();

      const std::string statistics_contents = stream.str();

      // Next find out whether we need to write everything into
      // the statistics file, or whether it is enough to just write
      // the last few bytes that were added since we wrote to that
      // file again. The way we do that is by checking whether the
      // first few bytes of the string we just created match what we
      // had previously written. One might think that they always should,
      // but the statistics object automatically sizes the column widths
      // of its output to match what is being written, and so if a later
      // entry requires more width, then even the first columns are
      // changed -- in that case, we will have to write everything,
      // not just append one line.
      const bool write_everything
        = ( // We may have never written anything. More precisely, this
            // case happens if the statistics_last_write_size is at the
            // value initialized by the Simulator::Simulator()
            // constructor, and this can happen in two situations:
            // (i) At the end of the first time step; and (ii) upon restart
            // since the variable we query here is not serialized. It is clear
            // that in both situations, we want to write the
            // entire contents of the statistics object. For the second
            // case, this is also appropriate since someone may have
            // previously restarted from a checkpoint, run a couple of
            // time steps that have added to the statistics file, but then
            // aborted the run again; a later restart from the same
            // checkpoint then requires overwriting the statistics file
            // on disk with what we have when this function is called for
            // the first time after the restart. The same situation
            // happens if the simulation kept running for some time after
            // a checkpoint, but is resumed from that checkpoint (i.e.,
            // at an earlier time step than when the statistics file was
            // written to last). In these situations, we effectively want
            // to "truncate" the file to the state stored in the checkpoint,
            // and we do that by just overwriting the entire file.
            (statistics_last_write_size == 0)
            ||
            // Or the size of the statistics file may have
            // shrunk mysteriously -- this shouldn't happen
            // but if it did we'd get into trouble with the
            // .substr() call in the next check.
            (statistics_last_write_size > statistics_contents.size())
            ||
            // Or the hash of what we wrote last time doesn't match
            // the hash of the first part of what we want to write
            (statistics_last_hash
             !=
             std::hash<std::string>()(statistics_contents.substr(0, statistics_last_write_size))) );

      const std::string stat_file_name = parameters.output_directory + "statistics";
      if (write_everything)
        {
          // Write what we have into a tmp file, then move that into
          // place
          const std::string tmp_file_name = stat_file_name + ".tmp";
          {
            std::ofstream tmp_file (tmp_file_name);
            tmp_file << statistics_contents;
          }
          std::rename(tmp_file_name.c_str(), stat_file_name.c_str());
        }
      else
        {
          // If we don't have to write everything, then the first part of what
          // we want to write matches what's already on disk. In that case,
          // we just have to append what's new.
          std::ofstream stat_file (stat_file_name, std::ios::app);
          stat_file << statistics_contents.substr(statistics_last_write_size, std::string::npos);
        }

      // Now update the size and hash of what we just wrote so that
      // we can compare against it next time we get here. Note that we do
      // not need to guard access to these variables with a mutex because
      // this is the only function that touches the variables, and
      // this function runs only once at a time (on a different
      // thread, but it's not started a second time while the previous
      // run hasn't finished).
      statistics_last_write_size = statistics_contents.size();
      statistics_last_hash       = std::hash<std::string>()(statistics_contents);
    });
  }



  template <int dim>
  double
  Simulator<dim>::
  compute_pressure_scaling_factor() const
  {
    // Determine how to treat the pressure. We have to scale it for the solver
    // to make velocities and pressures of roughly the same (numerical) size,
    // and we may have to fix up the right hand side vector before solving for
    // compressible models if there are no in-/outflow boundaries
    //
    // We do this by scaling the divergence equation by a constant factor that
    // is equal to a reference viscosity divided by a length scale.
    // We get the latter from the geometry model, and the former
    // by looping over all cells and averaging the "order of magnitude"
    // of the viscosity. The order of magnitude is the logarithm of
    // the viscosity, so
    //
    //   \eta_{ref} = exp( 1/N * (log(eta_1) + log(eta_2) + ... + log(eta_N))
    //
    // where the \eta_i are typical viscosities on the cells of the mesh.
    // For this, we just take the viscosity at the cell center.
    //
    // The formula above computes the exponential of the average of the
    // logarithms. It is easy to verify that this is equivalent to
    // computing the *geometric* mean of the viscosities, but the
    // formula above is numerically stable.

    // Do return NaN if we do not solve the Stokes equation. Technically, there is
    // no harm in computing the factor anyway, but the viscosity has to be positive
    // for this function to work, and some models may set viscosity to 0 if not solving
    // the Stokes equation. Returning NaN guarantees this value cannot be used.
    if (parameters.nonlinear_solver == NonlinearSolver::no_Advection_no_Stokes ||
        parameters.nonlinear_solver == NonlinearSolver::single_Advection_no_Stokes)
      return numbers::signaling_nan<double>();

    const QMidpoint<dim> quadrature_formula;

    FEValues<dim> fe_values (*mapping,
                             finite_element,
                             quadrature_formula,
                             update_values |
                             update_gradients |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();

    double local_integrated_viscosity_logarithm = 0.;
    double local_volume = 0.;

    MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                               introspection.n_compositional_fields);
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                 introspection.n_compositional_fields);

    // We do not need to compute anything but the viscosity
    in.requested_properties = MaterialModel::MaterialProperties::viscosity;

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          in.reinit(fe_values,
                    cell,
                    introspection,
                    solution);

          // We do not call the cell-wise average function of the
          // material model, because we average globally below
          material_model->evaluate(in, out);

          // Evaluate viscosity at the mid-point of each cell and
          // calculate the volume weighted harmonic average of all cells
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              Assert(out.viscosities[q] > 0,
                     ExcMessage ("The viscosity needs to be a "
                                 "positive quantity."));

              const double JxW = fe_values.JxW(q);
              local_integrated_viscosity_logarithm += std::log(out.viscosities[q]) * JxW;
              local_volume += JxW;
            }
        }

    // vector for packing local values before MPI summing them
    double values[2] = {local_integrated_viscosity_logarithm, local_volume};

    Utilities::MPI::sum(values, mpi_communicator, values);

    const double reference_viscosity = std::exp(values[0]/values[1]);
    const double length_scale = geometry_model->length_scale();

    // Allow the user to inspect and/or overwrite our result:
    const double result = reference_viscosity / length_scale;
    if (boost::optional<double> value = signals.modify_pressure_scaling(result, reference_viscosity, length_scale))
      return *value;
    else
      return result;
  }



  template <int dim>
  double
  Simulator<dim>::
  get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const
  {
    // use a quadrature formula that has one point at
    // the location of each degree of freedom in the
    // velocity element
    const QIterated<dim> quadrature_formula (QTrapezoid<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();


    FEValues<dim> fe_values (*mapping, finite_element, quadrature_formula, update_values);
    std::vector<Tensor<1,dim>> velocity_values(n_q_points);

    double max_local_velocity = 0;

    // loop over all locally owned cells and evaluate the velocities at each
    // quadrature point (i.e. each node). keep a running tally of the largest
    // such velocity
    for (const auto &cell : dof_handler.active_cell_iterators())
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
  bool Simulator<dim>::maybe_do_initial_refinement (const unsigned int max_refinement_level)
  {
    if (pre_refinement_step < parameters.initial_adaptive_refinement)
      {
        if (parameters.timing_output_frequency == 0)
          {
            computing_timer.print_summary ();
            pcout << "-- Total wallclock time elapsed including restarts: "
                  << std::round(wall_timer.wall_time()+total_walltime_until_last_snapshot)
                  << 's' << std::endl;
          }

        output_statistics();

        // we only want to do the postprocessing here if it is not already done in
        // the nonlinear iteration scheme, which is the case if we run postprocessors
        // on all nonlinear iterations
        if (parameters.run_postprocessors_on_initial_refinement && (!parameters.run_postprocessors_on_nonlinear_iterations))
          postprocess ();

        refine_mesh (max_refinement_level);
        ++pre_refinement_step;
        return true;
      }
    else
      {
        // invalidate the value of pre_refinement_step since it will no longer be used from here on
        pre_refinement_step = std::numeric_limits<unsigned int>::max();
        return false;
      }
  }



  template <int dim>
  void Simulator<dim>::exchange_refinement_flags ()
  {
    // Communicate refinement flags on ghost cells from the owner of the
    // cell. This is necessary to get consistent refinement, as mesh
    // smoothing would undo some of the requested coarsening/refinement.

    auto pack
    = [] (const typename DoFHandler<dim>::active_cell_iterator &cell) -> std::uint8_t
    {
      if (cell->refine_flag_set())
        return 1;
      if (cell->coarsen_flag_set())
        return 2;
      return 0;
    };
    auto unpack
    = [] (const typename DoFHandler<dim>::active_cell_iterator &cell, const std::uint8_t &flag) -> void
    {
      cell->clear_coarsen_flag();
      cell->clear_refine_flag();
      if (flag==1)
        cell->set_refine_flag();
      else if (flag==2)
        cell->set_coarsen_flag();
    };

    GridTools::exchange_cell_data_to_ghosts<std::uint8_t, DoFHandler<dim>>
    (dof_handler, pack, unpack);
  }



  template <int dim>
  void Simulator<dim>::maybe_refine_mesh (const double new_time_step,
                                          unsigned int &max_refinement_level)
  {
    /*
     * see if this is an additional refinement cycle. An additional refinement
     * cycle differs from a regular, because the maximal refinement level allowed
     * is increased by one from this time on.
     */
    if ((parameters.additional_refinement_times.size() > 0)
        &&
        (parameters.additional_refinement_times.front () < time+new_time_step))
      {
        // loop over as many times as this is necessary
        while ((parameters.additional_refinement_times.size() > 0)
               &&
               (parameters.additional_refinement_times.front () < time+new_time_step))
          {
            ++max_refinement_level;
            refine_mesh (max_refinement_level);

            parameters.additional_refinement_times
            .erase (parameters.additional_refinement_times.begin());
          }
      }
    // see if this is a time step where regular refinement is requested
    else if ((timestep_number > 0
              &&
              parameters.adaptive_refinement_interval > 0
              &&
              timestep_number % parameters.adaptive_refinement_interval == 0)
             ||
             (timestep_number == 0 && parameters.adaptive_refinement_interval == 1)
            )
      {
        refine_mesh (max_refinement_level);
      }
  }



  template <int dim>
  void Simulator<dim>::maybe_write_timing_output () const
  {
    bool write_timing_output = false;
    if (parameters.timing_output_frequency <= 1)
      write_timing_output = true;
    else if ((timestep_number > 0) &&
             (timestep_number % parameters.timing_output_frequency == 0))
      write_timing_output = true;

    // if requested output a summary of the current timing information
    if (write_timing_output)
      {
        computing_timer.print_summary ();
        pcout << "-- Total wallclock time elapsed including restarts: "
              << std::round(wall_timer.wall_time()+total_walltime_until_last_snapshot)
              << 's' << std::endl;
      }
  }



  template <int dim>
  bool Simulator<dim>::maybe_write_checkpoint (const time_t last_checkpoint_time,
                                               const bool force_writing_checkpoint)
  {
    // Do a checkpoint if this is the end of simulation,
    // and the termination criteria say to checkpoint at the end.
    bool write_checkpoint = force_writing_checkpoint;

    // If we base checkpoint frequency on timing, measure the time at process 0
    // This prevents race conditions where some processes will checkpoint and others won't
    if (!write_checkpoint && parameters.checkpoint_time_secs > 0)
      {
        int global_do_checkpoint = ((std::time(nullptr)-last_checkpoint_time) >=
                                    parameters.checkpoint_time_secs);
        const int ierr = MPI_Bcast(&global_do_checkpoint, 1, MPI_INT, 0, mpi_communicator);
        AssertThrowMPI(ierr);

        if (global_do_checkpoint == 1)
          write_checkpoint = true;
      }

    // If we base checkpoint frequency on steps, see if it's time for another checkpoint
    if (!write_checkpoint &&
        (parameters.checkpoint_time_secs == 0) &&
        (parameters.checkpoint_steps > 0) &&
        (timestep_number % parameters.checkpoint_steps == 0))
      write_checkpoint = true;

    // Do a checkpoint if indicated by checkpoint parameters
    if (write_checkpoint)
      {
        create_snapshot();
        // matrices will be regenerated after a resume, so do that here too
        // to be consistent. otherwise we would get different results
        // for a restarted computation than for one that ran straight
        // through
        rebuild_stokes_matrix =
          rebuild_stokes_preconditioner = true;
      }
    return write_checkpoint;
  }



  template <int dim>
  void Simulator<dim>::advance_time (const double step_size)
  {
    old_time_step = time_step;
    time_step = step_size;
    time += time_step;
    ++timestep_number;

    // prepare for the next time step by shifting solution vectors
    // by one time step. In timestep 0 (just increased in the
    // line above) initialize both old_solution
    // and old_old_solution with the currently computed solution.
    if (timestep_number == 1)
      {
        old_old_solution      = solution;
        old_solution          = solution;
      }
    else
      {
        old_old_solution      = old_solution;
        old_solution          = solution;
      }
  }



  template <int dim>
  std::pair<double,double>
  Simulator<dim>::
  get_extrapolated_advection_field_range (const AdvectionField &advection_field) const
  {
    const QIterated<dim> quadrature_formula (QTrapezoid<1>(),
                                             advection_field.polynomial_degree(introspection));

    const unsigned int n_q_points = quadrature_formula.size();

    const FEValuesExtractors::Scalar field = advection_field.scalar_extractor(introspection);

    FEValues<dim> fe_values (*mapping, finite_element, quadrature_formula,
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
           max_local_field = std::numeric_limits<double>::lowest();

    if (timestep_number > 1)
      {
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
        for (const auto &cell : dof_handler.active_cell_iterators())
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

    return std::make_pair(Utilities::MPI::min (min_local_field,
                                               mpi_communicator),
                          Utilities::MPI::max (max_local_field,
                                               mpi_communicator));
  }


  template <int dim>
  void Simulator<dim>::interpolate_onto_velocity_system(const TensorFunction<1,dim> &func,
                                                        LinearAlgebra::Vector &vec) const
  {
    Assert(introspection.block_indices.velocities == 0, ExcNotImplemented());

    const std::vector<Point<dim>> mesh_support_points = finite_element.base_element(introspection.base_elements.velocities).get_unit_support_points();
    const unsigned int n_velocity_dofs_per_cell = finite_element.base_element(introspection.base_elements.velocities).dofs_per_cell;

    FEValues<dim> mesh_points (*mapping, finite_element, mesh_support_points, update_quadrature_points);
    std::vector<types::global_dof_index> cell_dof_indices (finite_element.dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          mesh_points.reinit(cell);
          cell->get_dof_indices (cell_dof_indices);
          for (unsigned int j=0; j<n_velocity_dofs_per_cell; ++j)
            for (unsigned int dir=0; dir<dim; ++dir)
              {
                const unsigned int support_point_index
                  = finite_element.component_to_system_index(/*velocity component=*/ introspection.component_indices.velocities[dir],
                                                                                     /*dof index within component=*/ j);
                vec[cell_dof_indices[support_point_index]] = func.value(mesh_points.quadrature_point(j))[dir];
              }
        }

    vec.compress(VectorOperation::insert);

#if DEAL_II_VERSION_GTE(9,7,0)
    AffineConstraints<double> hanging_node_constraints(introspection.index_sets.system_relevant_set,
                                                       introspection.index_sets.system_relevant_set);
#else
    AffineConstraints<double> hanging_node_constraints(introspection.index_sets.system_relevant_set);
#endif

    DoFTools::make_hanging_node_constraints(dof_handler, hanging_node_constraints);
    hanging_node_constraints.close();

    // Create a view of all constraints that only pertains to the
    // Stokes subset of degrees of freedom. We can then use this later
    // to call constraints.distribute(), constraints.set_zero(), etc.,
    // on those block vectors that only have the Stokes components in
    // them.
    //
    // For the moment, assume that the Stokes degrees are first in the
    // overall vector, so that they form a contiguous range starting
    // at zero. The assertion checks this, but this could easily be
    // generalized if the Stokes block were not starting at zero.
#if DEAL_II_VERSION_GTE(9,6,0)
    Assert (introspection.block_indices.velocities == 0,
            ExcNotImplemented());
    if (parameters.use_direct_stokes_solver == false)
      Assert (introspection.block_indices.pressure == 1,
              ExcNotImplemented());

    IndexSet stokes_dofs (dof_handler.n_dofs());
    stokes_dofs.add_range (0, vec.size());
    const AffineConstraints<double> stokes_hanging_node_constraints
      = hanging_node_constraints.get_view (stokes_dofs);
#else
    const AffineConstraints<double> &stokes_hanging_node_constraints = hanging_node_constraints;
#endif

    stokes_hanging_node_constraints.distribute(vec);
  }



  template <int dim>
  double Simulator<dim>::normalize_pressure (LinearAlgebra::BlockVector &vector) const
  {
    if (parameters.pressure_normalization == "no")
      return 0;

    const FEValuesExtractors::Scalar &extractor_pressure =
      (parameters.include_melt_transport ?
       introspection.variable("fluid pressure").extractor_scalar()
       : introspection.extractors.pressure);

    double my_pressure = 0.0;
    double my_area = 0.0;
    if (parameters.pressure_normalization == "surface")
      {
        const types::boundary_id top_boundary_id = geometry_model->translate_symbolic_boundary_name_to_id("top");

        const Quadrature<dim-1> &quadrature = this->introspection.face_quadratures.pressure;
        const unsigned int n_q_points = quadrature.size();

        FEFaceValues<dim> fe_face_values (*mapping, finite_element,  quadrature,
                                          update_JxW_values | update_values);

        std::vector<double> pressure_values(n_q_points);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              for (const unsigned int face_no : cell->face_indices())
                {
                  const typename DoFHandler<dim>::face_iterator face = cell->face (face_no);
                  if (face->at_boundary() && face->boundary_id() == top_boundary_id)
                    {
                      fe_face_values.reinit (cell, face_no);
                      fe_face_values[extractor_pressure].get_function_values(vector,
                                                                             pressure_values);

                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          const double JxW = fe_face_values.JxW(q);
                          my_pressure += pressure_values[q]
                                         * JxW;
                          my_area += JxW;
                        }
                    }
                }
            }
      }
    else if (parameters.pressure_normalization == "volume")
      {
        const Quadrature<dim> &quadrature = this->introspection.quadratures.pressure;
        const unsigned int n_q_points = quadrature.size();

        FEValues<dim> fe_values (*mapping, finite_element,  quadrature,
                                 update_JxW_values | update_values);

        std::vector<double> pressure_values(n_q_points);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[extractor_pressure].get_function_values(vector,
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

    // sum up the integrals from each processor and compute the result we care about
    double pressure_adjustment = numbers::signaling_nan<double>();
    {
      const double my_temp[2] = {my_pressure, my_area};
      double temp[2];
      Utilities::MPI::sum (my_temp, mpi_communicator, temp);
      const double pressure = temp[0];
      const double area = temp[1];

      Assert (area > 0,
              ExcMessage("While computing the average pressure, the area/volume "
                         "to integrate over was found to be zero or negative. This "
                         "indicates that no appropriate surface faces were found, "
                         "which is typically the case if the geometry model is not "
                         "set up correctly."));

      if (parameters.pressure_normalization == "surface")
        pressure_adjustment = -pressure/area + parameters.surface_pressure;
      else if (parameters.pressure_normalization == "volume")
        pressure_adjustment = -pressure/area;
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
        if (introspection.block_indices.velocities != introspection.block_indices.pressure
            && !parameters.include_melt_transport)
          distributed_vector.block(introspection.block_indices.pressure).add(pressure_adjustment);
        else
          {
            // pressure is not in a separate block, so we have to modify the values manually
            const unsigned int pressure_component = (parameters.include_melt_transport ?
                                                     introspection.variable("fluid pressure").first_component_index
                                                     : introspection.component_indices.pressure);
            const unsigned int n_local_pressure_dofs = (parameters.include_melt_transport ?
                                                        finite_element.base_element(introspection.variable("fluid pressure").base_index).dofs_per_cell
                                                        : finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell);
            std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  cell->get_dof_indices (local_dof_indices);
                  for (unsigned int j=0; j<n_local_pressure_dofs; ++j)
                    {
                      unsigned int support_point_index
                        = finite_element.component_to_system_index(pressure_component,
                                                                   /*dof index within component=*/ j);

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
        Assert (dynamic_cast<const FE_DGP<dim>*>(&finite_element.base_element(introspection.base_elements.pressure)) != nullptr,
                ExcInternalError());
        const unsigned int pressure_component = (parameters.include_melt_transport ?
                                                 introspection.variable("fluid pressure").first_component_index
                                                 : introspection.component_indices.pressure);
        std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = finite_element.component_to_system_index (pressure_component, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (dof_handler.locally_owned_dofs().is_element(local_dof_indices[first_pressure_dof]),
                      ExcInternalError());

              // then adjust its value
              distributed_vector(local_dof_indices[first_pressure_dof]) = vector(local_dof_indices[first_pressure_dof])
                                                                          + pressure_adjustment;
            }
        distributed_vector.compress(VectorOperation::insert);
      }

    // now get back to the original vector and return the adjustment used
    // in the computations above
    vector = distributed_vector;

    return pressure_adjustment;
  }



  template <int dim>
  void
  Simulator<dim>::
  denormalize_pressure (const double                      pressure_adjustment,
                        LinearAlgebra::BlockVector       &vector) const
  {
    if (parameters.pressure_normalization == "no")
      return;

    if (parameters.use_locally_conservative_discretization == false)
      {
        if ((introspection.block_indices.velocities != introspection.block_indices.pressure)
            && !parameters.include_melt_transport)
          vector.block(introspection.block_indices.pressure).add(-1.0 * pressure_adjustment);
        else
          {
            // pressure is not in a separate block so we have to modify the values manually
            const unsigned int pressure_component = (parameters.include_melt_transport ?
                                                     introspection.variable("fluid pressure").first_component_index
                                                     : introspection.component_indices.pressure);
            const unsigned int n_local_pressure_dofs = (parameters.include_melt_transport ?
                                                        finite_element.base_element(introspection.variable("fluid pressure").base_index).dofs_per_cell
                                                        : finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell);

            // We may touch the same DoF multiple times, so we need to copy the
            // vector before modifying it to have access to the original value.
            LinearAlgebra::BlockVector vector_backup;
            vector_backup.reinit(vector, /* omit_zeroing_entries = */ true);
            const unsigned int pressure_block_index =
              parameters.include_melt_transport ?
              introspection.variable("fluid pressure").block_index
              :
              introspection.block_indices.pressure;
            vector_backup.block(pressure_block_index) = vector.block(pressure_block_index);

            std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  cell->get_dof_indices (local_dof_indices);
                  for (unsigned int j=0; j<n_local_pressure_dofs; ++j)
                    {
                      const unsigned int local_dof_index
                        = finite_element.component_to_system_index(pressure_component,
                                                                   /*dof index within component=*/ j);

                      // Then adjust its value. vector could be a vector with ghost elements
                      // or a fully distributed vector. In the latter case only access dofs that are
                      // locally owned. Note that because we end up touching
                      // entries more than once, we are not simply incrementing
                      // vector but copy from the vector_backup copy.
                      if (vector.has_ghost_elements() ||
                          dof_handler.locally_owned_dofs().is_element(local_dof_indices[local_dof_index]))
                        vector(local_dof_indices[local_dof_index])
                          = vector_backup(local_dof_indices[local_dof_index]) - pressure_adjustment;
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
        Assert (dynamic_cast<const FE_DGP<dim>*>(&finite_element.base_element(introspection.base_elements.pressure)) != nullptr,
                ExcInternalError());
        Assert(!parameters.include_melt_transport, ExcNotImplemented());
        const unsigned int pressure_component = introspection.component_indices.pressure;

        // We may touch the same DoF multiple times, so we need to copy the
        // vector before modifying it to have access to the original value.
        LinearAlgebra::BlockVector vector_backup;
        vector_backup.reinit(vector, /* omit_zeroing_entries = */ true);
        const unsigned int pressure_block_index =
          parameters.include_melt_transport ?
          introspection.variable("fluid pressure").block_index
          :
          introspection.block_indices.pressure;
        vector_backup.block(pressure_block_index) = vector.block(pressure_block_index);

        std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = finite_element.component_to_system_index (pressure_component, 0);

              // make sure that this DoF is really owned by the current processor
              Assert (dof_handler.locally_owned_dofs().is_element(local_dof_indices[first_pressure_dof]),
                      ExcInternalError());

              // then adjust its value
              vector (local_dof_indices[first_pressure_dof]) = vector_backup(local_dof_indices[first_pressure_dof])
                                                               - pressure_adjustment;
            }

        vector.compress(VectorOperation::insert);
      }
  }



  template <int dim>
  void
  Simulator<dim>::make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector)
  {
    // If the mass conservation is written as
    //   div u = f
    // make sure this is solvable by modifying f to ensure that
    // int_\Omega f = int_\Omega div u = 0
    //
    // We have to deal with several complications:
    // - we can have an FE_Q or an FE_DGP for the pressure
    // - we might use a direct solver, so pressure and velocity are in the same block
    // - we might have melt transport, where we need to operate only on p_f
    //
    // We ensure int_\Omega f = 0 by computing a correction factor
    //   c = \int f
    // and adjust pressure RHS to be
    //  fnew = f - c/|\Omega|
    // such that
    //   \int fnew = \int f - c/|\Omega| = -c + \int f = 0.
    //
    // We can compute
    //   c = \int f = (f, 1) = (f, \sum_i \phi_i) = \sum_i (f, \phi_i) = \sum_i F_i
    // which is just the sum over the RHS vector for FE_Q. For FE_DGP
    // we need to restrict to 0th shape functions on each cell,
    // because this is the shape function that is constant 1. (The
    // other shape functions have mean value zero.)
    //
    // To make the adjustment fnew = f - c/|\Omega|
    // note that
    // fnew_i = f_i - c/|\Omega| * (1, \phi_i)
    // and the same logic for FE_DGP applies


    if ((!parameters.use_locally_conservative_discretization)
        &&
        (!parameters.include_melt_transport)
        &&
        (introspection.block_indices.velocities != introspection.block_indices.pressure))
      {
        // Easy Case. We have an FE_Q in a separate block, so we can use
        // mean_value() and vector.block(p) += correction:
        const double mean = vector.block(introspection.block_indices.pressure).mean_value();
        Assert(std::isfinite(mean), ExcInternalError());
        const double int_rhs = mean * vector.block(introspection.block_indices.pressure).size();
        const double correction = -int_rhs / global_volume;
        Assert(global_volume > 0.0, ExcInternalError());

        vector.block(introspection.block_indices.pressure).add(correction, pressure_shape_function_integrals.block(introspection.block_indices.pressure));
      }
    else if (!parameters.use_locally_conservative_discretization)
      {
        // FE_Q but we can not access the pressure block separately (either
        // a direct solver or we have melt with p_f and p_c in the same block).
        // Luckily we don't need to go over DoFs on each cell, because we
        // have IndexSets to help us:

        // we need to operate only on p_f not on p_c
        const IndexSet &idxset = parameters.include_melt_transport ?
                                 introspection.index_sets.locally_owned_fluid_pressure_dofs
                                 :
                                 introspection.index_sets.locally_owned_pressure_dofs;
        double int_rhs = 0.0;

        for (unsigned int i=0; i < idxset.n_elements(); ++i)
          {
            types::global_dof_index idx = idxset.nth_index_in_set(i);
            int_rhs += vector(idx);
          }

        // We do not have to integrate over the normal velocity at the
        // boundaries with a prescribed velocity because the constraints
        // are already distributed to the right hand side in
        // current_constraints.distribute.
        const double global_int_rhs = Utilities::MPI::sum(int_rhs, mpi_communicator);
        const double correction = - global_int_rhs / global_volume;

        for (unsigned int i=0; i < idxset.n_elements(); ++i)
          {
            types::global_dof_index idx = idxset.nth_index_in_set(i);
            vector(idx) += correction * pressure_shape_function_integrals(idx);
          }

        vector.compress(VectorOperation::add);
      }
    else
      {
        // Locally conservative with or without direct solver and with or
        // without melt: grab a pickaxe and do everything by hand!
        AssertThrow(parameters.use_locally_conservative_discretization,
                    ExcInternalError());

        double int_rhs = 0.0;
        const unsigned int pressure_component = (parameters.include_melt_transport ?
                                                 introspection.variable("fluid pressure").first_component_index
                                                 : introspection.component_indices.pressure);
        std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = finite_element.component_to_system_index (pressure_component, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (dof_handler.locally_owned_dofs().is_element(local_dof_indices[first_pressure_dof]),
                      ExcInternalError());

              // compute integral:
              int_rhs += vector(local_dof_indices[first_pressure_dof]);
            }

        const double global_int_rhs = Utilities::MPI::sum(int_rhs, mpi_communicator);
        const double correction = - global_int_rhs / global_volume;

        // Now modify our RHS with the correction factor:
        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = finite_element.component_to_system_index (pressure_component, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (dof_handler.locally_owned_dofs().is_element(local_dof_indices[first_pressure_dof]),
                      ExcInternalError());

              // correct:
              types::global_dof_index idx = local_dof_indices[first_pressure_dof];
              vector(idx) += correction * pressure_shape_function_integrals(idx);
            }

        vector.compress(VectorOperation::add);
      }
  }


  template <int dim>
  double
  Simulator<dim>::compute_initial_stokes_residual()
  {
    LinearAlgebra::BlockVector linearized_stokes_variables (introspection.index_sets.stokes_partitioning, mpi_communicator);
    LinearAlgebra::BlockVector residual (introspection.index_sets.stokes_partitioning, mpi_communicator);
    const unsigned int pressure_block_index =
      parameters.include_melt_transport ?
      introspection.variable("fluid pressure").block_index
      :
      introspection.block_indices.pressure;

    // if velocity and pressure are in the same block, we have to copy the
    // pressure to the solution and RHS vector with a zero velocity
    if (pressure_block_index == introspection.block_indices.velocities)
      {
        const IndexSet &idxset = (parameters.include_melt_transport) ?
                                 introspection.index_sets.locally_owned_fluid_pressure_dofs
                                 :
                                 introspection.index_sets.locally_owned_pressure_dofs;

        for (unsigned int i=0; i < idxset.n_elements(); ++i)
          {
            types::global_dof_index idx = idxset.nth_index_in_set(i);
            linearized_stokes_variables(idx)        = current_linearization_point(idx);
          }
        linearized_stokes_variables.block(pressure_block_index).compress(VectorOperation::insert);
      }
    else
      linearized_stokes_variables.block (pressure_block_index) = current_linearization_point.block (pressure_block_index);

    denormalize_pressure (this->last_pressure_normalization_adjustment, linearized_stokes_variables);
    current_constraints.set_zero (linearized_stokes_variables);

    linearized_stokes_variables.block (pressure_block_index) /= pressure_scaling;

    // we calculate the velocity residual with a zero velocity,
    // computing only the part of the RHS not balanced by the static pressure
    if (pressure_block_index == introspection.block_indices.velocities)
      {
        // we can use the whole block here because we set the velocity to zero above
        return system_matrix.block(0,0).residual (residual.block(0),
                                                  linearized_stokes_variables.block(0),
                                                  system_rhs.block(0));
      }
    else
      {
        const double residual_u = system_matrix.block(0,1).residual (residual.block(0),
                                                                     linearized_stokes_variables.block(1),
                                                                     system_rhs.block(0));
        const double residual_p = system_rhs.block(pressure_block_index).l2_norm();
        return std::sqrt(residual_u*residual_u+residual_p*residual_p);
      }
  }



  template <int dim>
  bool
  Simulator<dim>::stokes_matrix_depends_on_solution() const
  {
    // Currently, the only coefficient that really appears on the
    // left hand side of the Stokes equation is the viscosity and possibly
    // the density in the case of the implicit reference density profile
    // approximation.
    // If melt transport is included in the simulation, we have an
    // additional equation with more coefficients on the left hand
    // side.

    return (material_model->get_model_dependence().viscosity != MaterialModel::NonlinearDependence::none)
           || (parameters.formulation_mass_conservation ==
               Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile)
           || parameters.include_melt_transport;
  }



  template <int dim>
  bool
  Simulator<dim>::stokes_A_block_is_symmetric() const
  {
    // The A block should be symmetric, unless there are free surface stabilization terms, or
    // the user has forced the use of a different solver.
    if (mesh_deformation && mesh_deformation->get_boundary_indicators_requiring_stabilization().empty() == false)
      return false;
    else if (parameters.force_nonsymmetric_A_block_solver)
      return false;
    else
      return true;
  }



  template <int dim>
  void Simulator<dim>::apply_limiter_to_dg_solutions (const AdvectionField &advection_field)
  {
    // TODO: Modify to more robust method
    // Skip if this composition field is being set from the volume_of_fluid handler
    if (!advection_field.is_temperature() &&
        parameters.volume_of_fluid_tracking_enabled)
      if (volume_of_fluid_handler->field_index_for_name(introspection.name_for_compositional_index(advection_field.compositional_variable))
          != volume_of_fluid_handler->get_n_fields())
        return;

    /*
     * First setup the quadrature points which are used to find the maximum and minimum solution values at those points.
     * A quadrature formula that combines all quadrature points constructed as all tensor products of
     * 1) one dimensional Gauss points; 2) one dimensional Gauss-Lobatto points.
     * We require that the Gauss-Lobatto points (2) appear in only one direction.
     * Therefore, possible combination
     * in 2d: the combinations are 21, 12
     * in 3d: the combinations are 211, 121, 112
     */
    const QGauss<1> quadrature_formula_1 (advection_field.polynomial_degree(introspection)+1);
    const QGaussLobatto<1> quadrature_formula_2 (advection_field.polynomial_degree(introspection)+1);

    const unsigned int n_q_points_1 = quadrature_formula_1.size();
    const unsigned int n_q_points_2 = quadrature_formula_2.size();
    const unsigned int n_q_points   = dim * n_q_points_2 * Utilities::fixed_power<dim-1>(n_q_points_1);

    std::vector<Point <dim>> quadrature_points;
    quadrature_points.reserve(n_q_points);

    switch (dim)
      {
        case 2:
        {
          // append quadrature points combination 12
          for (unsigned int i=0; i < n_q_points_1 ; ++i)
            {
              const double  x = quadrature_formula_1.point(i)(0);
              for (unsigned int j=0; j < n_q_points_2 ; ++j)
                {
                  const double y = quadrature_formula_2.point(j)(0);
                  quadrature_points.push_back(Point<dim>(x,y));
                }
            }
          // append quadrature points combination 21
          for (unsigned int i=0; i < n_q_points_2 ; ++i)
            {
              const double  x = quadrature_formula_2.point(i)(0);
              for (unsigned int j=0; j < n_q_points_1 ; ++j)
                {
                  const double y = quadrature_formula_1.point(j)(0);
                  quadrature_points.push_back(Point<dim>(x,y));
                }
            }
          break;
        }

        case 3:
        {
          // append quadrature points combination 121
          for ( unsigned int i=0; i < n_q_points_1 ; ++i)
            {
              const double x = quadrature_formula_1.point(i)(0);
              for ( unsigned int j=0; j < n_q_points_2 ; ++j)
                {
                  const double y = quadrature_formula_2.point(j)(0);
                  for ( unsigned int k=0; k < n_q_points_1 ; ++k)
                    {
                      const double z = quadrature_formula_1.point(k)(0);
                      quadrature_points.push_back(Point<dim>(x,y,z));
                    }
                }
            }
          // append quadrature points combination 112
          for (unsigned int i=0; i < n_q_points_1 ; ++i)
            {
              const double x = quadrature_formula_1.point(i)(0);
              for (unsigned int j=0; j < n_q_points_1 ; ++j)
                {
                  const double y = quadrature_formula_1.point(j)(0);
                  for (unsigned int k=0; k < n_q_points_2 ; ++k)
                    {
                      const double z = quadrature_formula_2.point(k)(0);
                      quadrature_points.push_back(Point<dim>(x,y,z));
                    }
                }
            }
          // append quadrature points combination 211
          for (unsigned int i=0; i < n_q_points_2 ; ++i)
            {
              const double x = quadrature_formula_2.point(i)(0);
              for ( unsigned int j=0; j < n_q_points_1 ; ++j)
                {
                  const double y = quadrature_formula_1.point(j)(0);
                  for ( unsigned int k=0; k < n_q_points_1 ; ++k)
                    {
                      const double z = quadrature_formula_1.point(k)(0);
                      quadrature_points.push_back(Point<dim>(x,y,z));
                    }
                }
            }
          break;
        }

        default:
          Assert (false, ExcNotImplemented());
      }

    Assert (quadrature_points.size() == n_q_points, ExcInternalError());
    const Quadrature<dim> quadrature_formula(quadrature_points);

    // Quadrature rules only used for the numerical integration for better accuracy
    const Quadrature<dim> &quadrature_formula_0
      = (advection_field.is_temperature() ?
         introspection.quadratures.temperature :
         introspection.quadratures.compositional_fields[advection_field.compositional_variable]);
    const unsigned int n_q_points_0 = quadrature_formula_0.size();

    // fe values for points evaluation
    FEValues<dim> fe_values (*mapping,
                             finite_element,
                             quadrature_formula,
                             update_values   |
                             update_quadrature_points);
    std::vector<double> values (n_q_points);
    // fe values for numerical integration, with a number of quadrature points
    // that is equal to 1/dim times the number of total points above
    FEValues<dim> fe_values_0 (*mapping,
                               finite_element,
                               quadrature_formula_0,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
    std::vector<double> values_0 (n_q_points_0);

    const FEValuesExtractors::Scalar field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    const double max_solution_exact_global = (advection_field.is_temperature()
                                              ?
                                              parameters.global_temperature_max_preset
                                              :
                                              parameters.global_composition_max_preset[advection_field.compositional_variable]
                                             );
    const double min_solution_exact_global = (advection_field.is_temperature()
                                              ?
                                              parameters.global_temperature_min_preset
                                              :
                                              parameters.global_composition_min_preset[advection_field.compositional_variable]
                                             );

    LinearAlgebra::BlockVector distributed_solution (introspection.index_sets.system_partitioning,
                                                     mpi_communicator);
    const unsigned int block_idx = advection_field.block_index(introspection);
    distributed_solution.block(block_idx) = solution.block(block_idx);

    std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (local_dof_indices);
          // used to find the maximum, minimum
          fe_values.reinit (cell);
          fe_values[field].get_function_values(solution, values);
          // used for the numerical integration
          fe_values_0.reinit (cell);
          fe_values_0[field].get_function_values(solution, values_0);

          // Find the local max and local min
          const double min_solution_local = *std::min_element (values.begin(), values.end());
          const double max_solution_local = *std::max_element (values.begin(), values.end());
          // Find the trouble cell
          if (min_solution_local < min_solution_exact_global
              || max_solution_local > max_solution_exact_global)
            {
              // Compute the cell area and cell solution average
              double local_area = 0.0;
              double local_solution_average = 0.0;
              for (unsigned int q = 0; q < n_q_points_0; ++q)
                {
                  local_area += fe_values_0.JxW(q);
                  local_solution_average += values_0[q]*fe_values_0.JxW(q);
                }
              local_solution_average /= local_area;

              /*
               * Define theta: a scaling constant used to correct the old solution by the formula
               *   new_value = theta * (old_value-old_solution_cell_average)+old_solution_cell_average
               * where theta \in [0,1] defined as below.
               * After the correction, the new solution does not exceed the user-given
               * exact global maximum/minimum values. Meanwhile, the new solution's cell average
               * equals to the old solution's cell average.
               */
              double theta = 1.0;
              if (std::abs(max_solution_local-local_solution_average) > std::numeric_limits<double>::min())
                {
                  theta = std::min(theta, std::abs((max_solution_exact_global-local_solution_average)
                                                   / (max_solution_local-local_solution_average)));
                }
              if (std::abs(min_solution_local-local_solution_average) > std::numeric_limits<double>::min())
                {
                  theta = std::min(theta, std::abs((min_solution_exact_global-local_solution_average)
                                                   / (min_solution_local-local_solution_average)));
                }

              /* Modify the advection degrees of freedom of the numerical solution.
               * Note that we are using DG elements, so every DoF on a locally owned cell is locally owned;
               * this means that we do not need to check whether the 'distributed_solution' vector actually
               * stores the element we read from/write to here.
               */
              for (unsigned int j = 0;
                   j < finite_element.base_element(advection_field.base_element(introspection)).dofs_per_cell;
                   ++j)
                {
                  const unsigned int support_point_index = finite_element.component_to_system_index(
                                                             (advection_field.is_temperature()
                                                              ?
                                                              introspection.component_indices.temperature
                                                              :
                                                              introspection.component_indices.compositional_fields[advection_field.compositional_variable]
                                                             ),
                                                             /*dof index within component=*/ j);
                  const double solution_value = solution(local_dof_indices[support_point_index]);
                  const double limited_solution_value = theta * (solution_value-local_solution_average) + local_solution_average;
                  distributed_solution(local_dof_indices[support_point_index]) = limited_solution_value;
                }
            }
        }

    distributed_solution.compress(VectorOperation::insert);
    // now get back to the original vector
    solution.block(block_idx) = distributed_solution.block(block_idx);
  }



  template <int dim>
  void Simulator<dim>::update_solution_vectors_with_reaction_results (const unsigned int block_index,
                                                                      const LinearAlgebra::BlockVector &distributed_vector,
                                                                      const LinearAlgebra::BlockVector &distributed_reaction_vector)
  {
    solution.block(block_index) = distributed_vector.block(block_index);

    // we have to update the old solution with our reaction update too
    // so that the advection scheme will have the correct time stepping in the next step
    LinearAlgebra::BlockVector tmp;
    tmp.reinit(distributed_vector, false);

    // What we really want to do is
    //     old_solution.block(block_index) += distributed_reaction_vector.block(block_index);
    // but because 'old_solution' is a ghosted vector, we can't write into it directly. Rather,
    // we have to go around with a completely distributed vector.
    tmp.block(block_index) = old_solution.block(block_index);
    tmp.block(block_index) +=  distributed_reaction_vector.block(block_index);
    old_solution.block(block_index) = tmp.block(block_index);

    // Same here with going through a distributed vector.
    tmp.block(block_index) = old_old_solution.block(block_index);
    tmp.block(block_index) +=  distributed_reaction_vector.block(block_index);
    old_old_solution.block(block_index) = tmp.block(block_index);

    operator_split_reaction_vector.block(block_index) = distributed_reaction_vector.block(block_index);
  }



  template <int dim>
  void Simulator<dim>::compute_unique_advection_support_points(const std::vector<AdvectionField> &advection_fields,
                                                               std::vector<Point<dim>> &unique_support_points,
                                                               std::vector<std::vector<unsigned int>> &support_point_index_by_field) const
  {
    const unsigned int n_fields = advection_fields.size();

    unique_support_points.clear();
    support_point_index_by_field.clear();
    support_point_index_by_field.resize(n_fields);

    // Reserve space to avoid most memory reallocations
    if (n_fields > 0)
      {
        const unsigned int likely_number_of_support_points = dof_handler.get_fe().base_element(advection_fields[0].base_element(introspection)).get_unit_support_points().size();
        unique_support_points.reserve(likely_number_of_support_points * n_fields);
        for (unsigned int i=0; i<n_fields; ++i)
          support_point_index_by_field[i].reserve(likely_number_of_support_points);
      }

    // Loop through all support points and fill the output data structures
    for (unsigned int i=0; i<advection_fields.size(); ++i)
      {
        const std::vector<Point<dim>> &support_points = dof_handler.get_fe().base_element(advection_fields[i].base_element(introspection)).get_unit_support_points();

        for (const auto &support_point: support_points)
          {
            // Naive n^2 algorithm to merge all points into the data structure, speed is likely irrelevant, because
            // n is small and this function is rarely called.
            const auto it = std::find(unique_support_points.begin(), unique_support_points.end(), support_point);
            if (it != unique_support_points.end())
              {
                // We already have this support point, record its number:
                support_point_index_by_field[i].push_back(std::distance(unique_support_points.begin(), it));
              }
            else
              {
                // This is a new point that needs to be added:
                unique_support_points.push_back(support_point);
                support_point_index_by_field[i].push_back(unique_support_points.size()-1);
              }
          }
      }
  }



  template <int dim>
  void Simulator<dim>::compute_reactions ()
  {
    // if the time step has a length of zero, there are no reactions
    if (time_step == 0)
      return;

    TimerOutput::Scope timer (computing_timer, "Solve composition reactions");

    // we need some temporary vectors to store our updates to composition and temperature in
    // while we do the time stepping, before we copy them over to the solution vector in the end
    LinearAlgebra::BlockVector distributed_vector (introspection.index_sets.system_partitioning,
                                                   mpi_communicator);

    LinearAlgebra::BlockVector distributed_reaction_vector (introspection.index_sets.system_partitioning,
                                                            mpi_communicator);

    // we use a different (potentially smaller) time step than in the advection scheme.
    // and for the fixed step scheme, we want all of our reaction time steps (within one advection step) to have the same size
    const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(time_step / parameters.reaction_time_step),
                                                           std::max(parameters.reaction_steps_per_advection_step,1U));

    const double reaction_time_step_size = time_step / static_cast<double>(number_of_reaction_steps);

    if (parameters.reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step)
      Assert (reaction_time_step_size > 0,
              ExcMessage("Reaction time step must be greater than 0."));

    pcout << "   Solving composition reactions... " << std::flush;


    // We want to compute reactions in each support point for all fields (compositional fields and temperature). The reaction
    // rate for an individual field depends on the values of all other fields, so we have to step them forward in time together.
    // The rates comes from the material and heating model, otherwise we have a simple ODE in each point on each cell to solve.
    //
    // So far so good. Except that fields can have different Finite Element discretizations (degree, continuous/discontinuous)
    // and will have different support points. We solve this by computing the union of all support points and evaluating all fields
    // in these points for every cell (if necessary by interpolation). Then we solve the ODEs on each cell together and
    // write back all values that correspond to support points of that particular field. Field values we computed that are
    // not actually support points, we just throw away.

    // First compute all unique support points (temperature and compositions):
    std::vector<Point<dim>> unique_support_points;
    std::vector<std::vector<unsigned int>> support_point_index_by_field;
    std::vector<AdvectionField> advection_fields;

    // First add the temperature field
    advection_fields.push_back(Simulator<dim>::AdvectionField::temperature());
    // Then add all compositional fields
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      advection_fields.push_back(Simulator<dim>::AdvectionField::composition(c));

    const unsigned int n_fields = advection_fields.size();
    compute_unique_advection_support_points(advection_fields, unique_support_points, support_point_index_by_field);

    const Quadrature<dim> combined_support_points(unique_support_points);
    FEValues<dim> fe_values (*mapping,
                             dof_handler.get_fe(),
                             combined_support_points,
                             update_quadrature_points | update_values | update_gradients);

    const unsigned int n_q_points = combined_support_points.size();
    std::vector<types::global_dof_index> local_dof_indices (dof_handler.get_fe().dofs_per_cell);
    MaterialModel::MaterialModelInputs<dim> in(n_q_points, introspection.n_compositional_fields);
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points, introspection.n_compositional_fields);
    HeatingModel::HeatingModelOutputs heating_model_outputs(n_q_points, introspection.n_compositional_fields);

    // add reaction rate outputs
    material_model->create_additional_named_outputs(out);

    MaterialModel::ReactionRateOutputs<dim> *reaction_rate_outputs
      = out.template get_additional_output<MaterialModel::ReactionRateOutputs<dim>>();

    AssertThrow(reaction_rate_outputs != nullptr,
                ExcMessage("You are trying to use the operator splitting solver scheme, "
                           "but the material model you use does not support operator splitting "
                           "(it does not create ReactionRateOutputs, which are required for this "
                           "solver scheme)."));

    // some heating models require the additional outputs
    heating_model_manager.create_additional_material_model_inputs_and_outputs(in, out);

    // We use SUNDIALs ARKode to compute the reactions. Set up the required parameters.
    // TODO: Should we change some of these based on the Reaction time step input parameter?
    using VectorType = Vector<double>;
    SUNDIALS::ARKode<VectorType>::AdditionalData data;
    data.initial_time = time;
    data.final_time = time + time_step;
    data.initial_step_size = 0.001 * time_step;
    data.output_period = time_step;
    data.minimum_step_size = 1.e-6 * time_step;

    // Both tolerances are added, but the composition might become 0.
    // We therefore set the absolute tolerance to a very small value.
    data.relative_tolerance = 1e-6;
    data.absolute_tolerance = 1e-10;

    SUNDIALS::ARKode<VectorType> ode(data);

    std::vector<std::vector<double>> initial_values_C (n_q_points, std::vector<double> (introspection.n_compositional_fields));
    std::vector<double> initial_values_T (n_q_points);

    // We have to store all values of the temperature and composition fields in one long vector.
    // Create the vector and functions for transferring values between this vector and the material model inputs object.
    VectorType fields (n_q_points * n_fields);

    auto copy_fields_into_one_vector = [n_q_points,n_fields](const MaterialModel::MaterialModelInputs<dim> &in,
                                                             VectorType &fields)
    {
      for (unsigned int j=0; j<n_q_points; ++j)
        for (unsigned int f=0; f<n_fields; ++f)
          if (f==0)
            fields[j*n_fields+f] = in.temperature[j];
          else
            fields[j*n_fields+f] = in.composition[j][f-1];
      return;
    };

    auto copy_fields_into_material_model_inputs = [n_q_points,n_fields](const VectorType &fields,
                                                                        MaterialModel::MaterialModelInputs<dim> &in)
    {
      for (unsigned int j=0; j<n_q_points; ++j)
        for (unsigned int f=0; f<n_fields; ++f)
          if (f==0)
            in.temperature[j]      = fields[j*n_fields+f];
          else
            in.composition[j][f-1] = fields[j*n_fields+f];
      return;
    };

    auto copy_rates_into_one_vector = [n_q_points,n_fields](const MaterialModel::ReactionRateOutputs<dim> *reaction_out,
                                                            const HeatingModel::HeatingModelOutputs &heating_out,
                                                            VectorType &rates)
    {
      for (unsigned int j=0; j<n_q_points; ++j)
        for (unsigned int f=0; f<n_fields; ++f)
          if (f==0)
            rates[j*n_fields+f] = heating_out.rates_of_temperature_change[j];
          else
            rates[j*n_fields+f] = reaction_out->reaction_rates[j][f-1];
      return;
    };

    unsigned int total_iteration_count = 0;
    unsigned int number_of_solves = 0;

    // Make a loop first over all cells, than over all reaction time steps, and then over
    // all degrees of freedom in each element to compute the reactions. This is possible
    // because the reactions only depend on the temperature and composition values at a given
    // degree of freedom (and are independent of the solution in other points).

    // Note that the values for some degrees of freedom are set more than once in the loop
    // below where we assign the new values to distributed_vector (if they are located on the
    // interface between cells), as we loop over all cells, and then over all degrees of freedom
    // on each cell. Although this means we do some additional work, the results are still
    // correct, as we never read from distributed_vector inside the loop over all cells.
    // We initialize the material model inputs object in using the solution vector
    // on every cell, compute the update, and then on every cell put the result into the
    // distributed_vector vector. Only after the loop over all cells do we copy distributed_vector
    // back onto the solution vector.
    // So even though we touch some DoF more than once, we always start from the same value, compute the
    // same value, and then overwrite the same value in distributed_vector.
    // TODO: make this even more efficient.
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          in.reinit(fe_values, cell, introspection, solution);

          std::vector<std::vector<double>> accumulated_reactions_C (n_q_points, std::vector<double> (introspection.n_compositional_fields));
          std::vector<double> accumulated_reactions_T (n_q_points);

          copy_fields_into_one_vector (in, fields);

          initial_values_C = in.composition;
          initial_values_T = in.temperature;

          if (parameters.reaction_solver_type == Parameters<dim>::ReactionSolverType::ARKode)
            {

              ode.explicit_function = [&] (const double /*time*/,
                                           const VectorType &y,
                                           VectorType &ydot)
              {
                copy_fields_into_material_model_inputs (y, in);
                material_model->fill_additional_material_model_inputs(in, solution, fe_values, introspection);
                material_model->evaluate(in, out);
                heating_model_manager.evaluate(in, out, heating_model_outputs);
                copy_rates_into_one_vector (reaction_rate_outputs, heating_model_outputs, ydot);
              };

              // Make the reaction time steps: We have to update the values of compositional fields and the temperature.
              // We can reuse the same material model inputs and outputs structure for each reaction time step.
              // We store the computed updates to temperature and composition in a separate (accumulated_reactions) vector,
              // so that we can later copy it over to the solution vector.
              const unsigned int iteration_count = ode.solve_ode(fields);

              total_iteration_count += iteration_count;
              number_of_solves += 1;

              for (unsigned int j=0; j<n_q_points; ++j)
                for (unsigned int f=0; f<n_fields; ++f)
                  {
                    if (f==0)
                      {
                        in.temperature[j]          = fields[j*n_fields+f];
                        accumulated_reactions_T[j] = in.temperature[j] - initial_values_T[j];
                      }
                    else
                      {
                        in.composition[j][f-1]          = fields[j*n_fields+f];
                        accumulated_reactions_C[j][f-1] = in.composition[j][f-1] - initial_values_C[j][f-1];
                      }
                  }
            }

          else if (parameters.reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step)
            {
              for (unsigned int i=0; i<number_of_reaction_steps; ++i)
                {
                  // Loop over composition element
                  material_model->fill_additional_material_model_inputs(in, solution, fe_values, introspection);

                  material_model->evaluate(in, out);
                  heating_model_manager.evaluate(in, out, heating_model_outputs);

                  for (unsigned int j=0; j<n_q_points; ++j)
                    {
                      for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
                        {
                          // simple forward euler
                          in.composition[j][c] = in.composition[j][c]
                                                 + reaction_time_step_size * reaction_rate_outputs->reaction_rates[j][c];
                          accumulated_reactions_C[j][c] += reaction_time_step_size * reaction_rate_outputs->reaction_rates[j][c];
                        }
                      in.temperature[j] = in.temperature[j]
                                          + reaction_time_step_size * heating_model_outputs.rates_of_temperature_change[j];
                      accumulated_reactions_T[j] += reaction_time_step_size * heating_model_outputs.rates_of_temperature_change[j];
                    }
                }
            }

          cell->get_dof_indices (local_dof_indices);
          const unsigned int component_idx_T = introspection.component_indices.temperature;

          for (unsigned int dof_idx = 0; dof_idx < local_dof_indices.size(); ++dof_idx)
            {
              const auto comp_pair = dof_handler.get_fe().system_to_component_index(dof_idx);
              const unsigned int component_idx = comp_pair.first;
              if (component_idx>=component_idx_T) // ignore velocity, pressure, etc.
                {
                  // We found a DoF that belongs to component component_idx, which is a temperature or compositional
                  // field. That means we want to find where this DoF in the computed reactions above to copy it
                  // back into the global solution vector.

                  // These two variables tell us the how-manyth shape function of which field (and therefore
                  // field) this DoF is:
                  const unsigned int index_within = comp_pair.second;
                  const unsigned int field_index = component_idx-component_idx_T;
                  // Now we can look up in the support_point_index_by_field data structure where this support
                  // point is in the list of unique_support_points (and in the Quadrature):
                  const unsigned int point_idx = support_point_index_by_field[field_index][index_within];

                  // The final step is grabbing the value from the reaction computation and write it into
                  // the global vector (if we own it, of course):
                  if (dof_handler.locally_owned_dofs().is_element(local_dof_indices[dof_idx]))
                    {
                      // temperatures and compositions are stored differently:
                      if (component_idx == component_idx_T)
                        {
                          distributed_vector(local_dof_indices[dof_idx]) = in.temperature[point_idx];
                          distributed_reaction_vector(local_dof_indices[dof_idx]) = accumulated_reactions_T[point_idx];
                        }
                      else
                        {
                          const unsigned int composition = field_index-1; // 0 is temperature...
                          distributed_vector(local_dof_indices[dof_idx]) = in.composition[point_idx][composition];
                          distributed_reaction_vector(local_dof_indices[dof_idx]) = accumulated_reactions_C[point_idx][composition];
                        }
                    }
                }
            }
        }

    distributed_vector.compress(VectorOperation::insert);
    distributed_reaction_vector.compress(VectorOperation::insert);

    // put the final values into the solution vector
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      update_solution_vectors_with_reaction_results(introspection.block_indices.compositional_fields[c],
                                                    distributed_vector,
                                                    distributed_reaction_vector);

    update_solution_vectors_with_reaction_results(introspection.block_indices.temperature,
                                                  distributed_vector,
                                                  distributed_reaction_vector);

    initialize_current_linearization_point();

    double average_iteration_count = number_of_reaction_steps;
    if (parameters.reaction_solver_type == Parameters<dim>::ReactionSolverType::ARKode)
      {
        if (number_of_solves > 0)
          average_iteration_count = total_iteration_count / number_of_solves;
        else
          average_iteration_count = total_iteration_count;
      }

    pcout << "in "
          << average_iteration_count
          << " substep(s)."
          << std::endl;
  }



  template <int dim>
  void Simulator<dim>::initialize_current_linearization_point ()
  {
    // Start with a simple copy of the last timestep
    current_linearization_point = old_solution;

    // If possible use an extrapolated solution from last and
    // previous to last timestep.
    if (timestep_number > 1)
      {
        // TODO: Trilinos sadd does not like ghost vectors even as input. Copy
        // into distributed vectors for now:
        LinearAlgebra::BlockVector distr_solution (system_rhs);
        distr_solution = old_solution;
        LinearAlgebra::BlockVector distr_old_solution (system_rhs);
        distr_old_solution = old_old_solution;
        distr_solution.sadd ((1 + time_step/old_time_step),
                             -time_step/old_time_step,
                             distr_old_solution);
        current_linearization_point = distr_solution;
      }
  }



  template <int dim>
  void Simulator<dim>::interpolate_material_output_into_advection_field (const std::vector<AdvectionField> &adv_fields)
  {
    // we need a temporary vector to store our updates to the advection fields in
    // before we copy them over to the solution vector in the end
    LinearAlgebra::BlockVector distributed_vector (introspection.index_sets.system_partitioning,
                                                   mpi_communicator);

    if (adv_fields.size() == 0)
      return;
    else if (adv_fields.size() == 1)
      {
        const AdvectionField &adv_field = adv_fields[0];
        if (adv_field.is_temperature())
          pcout << "   Copying properties into prescribed temperature field... "
                << std::flush;
        else
          {
            const std::string name_of_field = introspection.name_for_compositional_index(adv_field.compositional_variable);

            pcout << "   Copying properties into prescribed compositional field " + name_of_field + "... "
                  << std::flush;
          }
      }
    else
      {
        pcout << "   Copying properties into prescribed compositional fields... "
              << std::flush;
      }

    // Advection fields can have different Finite Element discretizations (degree, continuous/discontinuous)
    // and will have different support points. We solve this by computing the union of all support points and evaluating
    // the material output for all these support points for every cell.

    /// First compute all unique support points (temperature and compositions):
    const unsigned int n_fields = adv_fields.size();
    std::vector<Point<dim>> unique_support_points;
    std::vector<std::vector<unsigned int>> support_point_index_by_field;
    compute_unique_advection_support_points(adv_fields, unique_support_points, support_point_index_by_field);

    // Create an FEValues object that allows us to interpolate onto the solution
    // vector. To make this happen, we need to have a quadrature formula that
    // consists of the support points of all advection field finite elements
    const Quadrature<dim> quadrature(unique_support_points);

    FEValues<dim> fe_values (*mapping,
                             dof_handler.get_fe(),
                             quadrature,
                             update_quadrature_points | update_values | update_gradients);

    std::vector<types::global_dof_index> local_dof_indices (dof_handler.get_fe().dofs_per_cell);
    MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), introspection.n_compositional_fields);
    MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), introspection.n_compositional_fields);

    // add the prescribed field outputs that will be used for interpolating
    material_model->create_additional_named_outputs(out);

    MaterialModel::PrescribedFieldOutputs<dim> *prescribed_field_out
      = out.template get_additional_output<MaterialModel::PrescribedFieldOutputs<dim>>();
    MaterialModel::PrescribedTemperatureOutputs<dim> *prescribed_temperature_out
      = out.template get_additional_output<MaterialModel::PrescribedTemperatureOutputs<dim>>();

    // Make a loop first over all cells, and then over all degrees of freedom in each element
    // to interpolate material properties onto a solution vector.
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          cell->get_dof_indices (local_dof_indices);
          in.reinit(fe_values, cell, introspection, solution);
          material_model->evaluate(in, out);

          for (unsigned int i=0; i<n_fields; ++i)
            {
              const AdvectionField &adv_field = adv_fields[i];

              // Interpolate material properties onto the advection fields
              const unsigned int advection_dofs_per_cell =
                dof_handler.get_fe().base_element(adv_field.base_element(introspection)).dofs_per_cell;

              // Make sure data structures have the expected size
              Assert(advection_dofs_per_cell == support_point_index_by_field[i].size(), ExcInternalError());

              for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                {
                  const unsigned int dof_idx
                    = dof_handler.get_fe().component_to_system_index(adv_field.component_index(introspection),
                                                                     /*dof index within component=*/ j);

                  // Skip degrees of freedom that are not locally owned. These
                  // will eventually be handled by one of the other processors.
                  if (dof_handler.locally_owned_dofs().is_element(local_dof_indices[dof_idx]))
                    {
                      if (adv_field.is_temperature())
                        {
                          Assert(prescribed_temperature_out != nullptr,
                                 ExcMessage("You are trying to use a prescribed temperature field, "
                                            "but the material model you use does not support interpolating properties "
                                            "(it does not create PrescribedTemperatureOutputs, which is required for this "
                                            "temperature field type)."));
                          Assert(numbers::is_finite(prescribed_temperature_out->prescribed_temperature_outputs[j]),
                                 ExcMessage("You are trying to use a prescribed advection field, "
                                            "but the material model you use does not fill the PrescribedFieldOutputs "
                                            "for your prescribed field, which is required for this method."));

                          distributed_vector(local_dof_indices[dof_idx])
                            = prescribed_temperature_out->prescribed_temperature_outputs[support_point_index_by_field[i][j]];
                        }
                      else
                        {
                          Assert(prescribed_field_out != nullptr,
                                 ExcMessage("You are trying to use a prescribed advection field, "
                                            "but the material model you use does not support interpolating properties "
                                            "(it does not create PrescribedFieldOutputs, which is required for this "
                                            "advection field type)."));
                          Assert(numbers::is_finite(prescribed_field_out->prescribed_field_outputs[j][adv_field.compositional_variable]),
                                 ExcMessage("You are trying to use a prescribed advection field, "
                                            "but the material model you use does not fill the PrescribedFieldOutputs "
                                            "for your prescribed field, which is required for this method."));

                          distributed_vector(local_dof_indices[dof_idx])
                            = prescribed_field_out->prescribed_field_outputs[support_point_index_by_field[i][j]][adv_field.compositional_variable];
                        }
                    }
                }
            }
        }

    for (const auto &adv_field: adv_fields)
      {
        // Put the final values into the solution vector, also
        // updating the ghost elements of the solution vector.
        const unsigned int advection_block = adv_field.block_index(introspection);
        distributed_vector.block(advection_block).compress(VectorOperation::insert);

        // Apply boundary conditions and other constraints to all prescribed fields,
        // except for the density field (if it exists). See this PR for a justification:
        // https://github.com/geodynamics/aspect/pull/4450.
        if (adv_field.is_temperature() ||
            adv_field.compositional_variable != introspection.find_composition_type(CompositionalFieldDescription::density))
          current_constraints.distribute (distributed_vector);

        solution.block(advection_block) = distributed_vector.block(advection_block);
      }

    pcout << "done." << std::endl;
  }



  template <int dim>
  void
  Simulator<dim>::check_consistency_of_formulation()
  {
    // Replace Formulation::MassConservation::ask_material_model by the respective terms to avoid
    // complicated checks later on
    if (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::ask_material_model)
      {
        if (material_model->is_compressible() == true)
          parameters.formulation_mass_conservation = Parameters<dim>::Formulation::MassConservation::isentropic_compression;
        else
          parameters.formulation_mass_conservation = Parameters<dim>::Formulation::MassConservation::incompressible;
      }

    // Ensure the material model supports the selected formulation of the mass conservation equation
    if (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::incompressible)
      {
        AssertThrow(material_model->is_compressible() == false,
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The mass conservation equation was selected to be incompressible, "
                               "but the provided material model reports that it is compressible. "
                               "Please check the consistency of your material model and selected formulation."));
      }
    else if (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::isentropic_compression
             || parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::reference_density_profile
             || parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile
             || parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::projected_density_field)
      {
        AssertThrow(material_model->is_compressible() == true,
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The mass conservation equation was selected to be compressible, "
                               "but the provided material model reports that it is incompressible. "
                               "Please check the consistency of your material model and selected formulation."));
      }

    // Ensure that the correct heating terms have been selected for the chosen combined formulation
    // Note that if the combined formulation is 'custom' there is no check
    // (useful e.g. for smaller scale lithospheric models with shear heating but without adiabatic heating)
    if (parameters.formulation == Parameters<dim>::Formulation::isentropic_compression)
      {
        AssertThrow(heating_model_manager.adiabatic_heating_enabled(),
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The `isentropic compression' formulation expects adiabatic heating to be enabled, "
                               "but the `adiabatic heating' plugin has not been selected in the input file. "
                               "Please check the consistency of your input file."));

        AssertThrow(heating_model_manager.shear_heating_enabled(),
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The `isentropic compression' formulation expects shear heating to be enabled, "
                               "but the `shear heating' plugin has not been selected in the input file. "
                               "Please check the consistency of your input file."));
      }
    else if (parameters.formulation == Parameters<dim>::Formulation::boussinesq_approximation)
      {
        AssertThrow(!heating_model_manager.adiabatic_heating_enabled(),
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The 'Boussinesq approximation' formulation expects adiabatic heating to be disabled, "
                               "but the `adiabatic heating' plugin has been selected in the input file. "
                               "Please check the consistency of your input file."));

        AssertThrow(!heating_model_manager.shear_heating_enabled(),
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The 'Boussinesq approximation' formulation expects shear heating to be disabled, "
                               "but the `shear heating' plugin has been selected in the input file. "
                               "Please check the consistency of your input file."));
      }
    else if (parameters.formulation == Parameters<dim>::Formulation::anelastic_liquid_approximation)
      {
        AssertThrow(heating_model_manager.adiabatic_heating_enabled(),
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The `anelastic liquid approximation' formulation expects adiabatic heating to be enabled, "
                               "but the `adiabatic heating' plugin has not been selected in the input file. "
                               "Please check the consistency of your input file."));

        AssertThrow(heating_model_manager.shear_heating_enabled(),
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The `anelastic liquid approximation' formulation expects shear heating to be enabled, "
                               "but the `shear heating' plugin has not been selected in the input file. "
                               "Please check the consistency of your input file."));

        const bool use_simplified_adiabatic_heating =
          heating_model_manager.template get_matching_active_plugin<HeatingModel::AdiabaticHeating<dim>>()
          .use_simplified_adiabatic_heating();

        AssertThrow(use_simplified_adiabatic_heating == true,
                    ExcMessage("ASPECT detected an inconsistency in the provided input file. "
                               "The `anelastic liquid approximation' formulation expects adiabatic heating to use "
                               "a simplified heating term that neglects dynamic pressure influences, "
                               "but the adiabatic heating plugin does not report to simplify this term. "
                               "Please check the consistency of your input file."));
      }
  }



  template <int dim>
  void
  Simulator<dim>::replace_outflow_boundary_ids(const unsigned int offset)
  {
    const Quadrature<dim-1> &quadrature_formula = introspection.face_quadratures.temperature;

    FEFaceValues<dim> fe_face_values (*mapping,
                                      finite_element,
                                      quadrature_formula,
                                      update_values   | update_normal_vectors |
                                      update_quadrature_points | update_JxW_values);

    std::vector<Tensor<1,dim>> face_current_velocity_values (fe_face_values.n_quadrature_points);
    std::vector<Tensor<1,dim>> face_current_mesh_velocity_values (fe_face_values.n_quadrature_points);

    const auto &tangential_velocity_boundaries =
      boundary_velocity_manager.get_tangential_boundary_velocity_indicators();

    const auto &zero_velocity_boundaries =
      boundary_velocity_manager.get_zero_boundary_velocity_indicators();

    const auto &prescribed_velocity_boundaries =
      boundary_velocity_manager.get_prescribed_boundary_velocity_indicators();

    // Loop over all of the boundary faces, ...
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (!cell->is_artificial())
        for (const unsigned int face_number : cell->face_indices())
          {
            const typename DoFHandler<dim>::face_iterator face = cell->face(face_number);

            // If the face is at a boundary where we may want to replace its id
            if (face->at_boundary() &&
                tangential_velocity_boundaries.find(face->boundary_id()) == tangential_velocity_boundaries.end() &&
                zero_velocity_boundaries.find(face->boundary_id()) == zero_velocity_boundaries.end())
              {
                Assert(face->boundary_id() <= offset,
                       ExcMessage("If you do not 'Allow fixed temperature/composition on outflow boundaries', "
                                  "you are only allowed to use boundary ids between 0 and 128."));

                fe_face_values.reinit (cell, face_number);
                fe_face_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
                                                                                        face_current_velocity_values);
                // get the mesh velocity, as we need to subtract it off of the advection systems
                if (parameters.mesh_deformation_enabled)
                  fe_face_values[introspection.extractors.velocities].get_function_values(mesh_deformation->mesh_velocity,
                                                                                          face_current_mesh_velocity_values);

                // ... check if the face is an outflow boundary by integrating the normal velocities
                // (flux through the boundary) as: int u*n ds = Sum_q u(x_q)*n(x_q) JxW(x_q)...
                double integrated_flow = 0;

                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    Tensor<1,dim> boundary_velocity;
                    if (prescribed_velocity_boundaries.find(face->boundary_id()) == prescribed_velocity_boundaries.end())
                      boundary_velocity = face_current_velocity_values[q];
                    else
                      boundary_velocity = boundary_velocity_manager.boundary_velocity(face->boundary_id(),
                                                                                      fe_face_values.quadrature_point(q));

                    if (parameters.mesh_deformation_enabled)
                      boundary_velocity -= face_current_mesh_velocity_values[q];

                    integrated_flow += (boundary_velocity * fe_face_values.normal_vector(q)) *
                                       fe_face_values.JxW(q);
                  }

                // ... and change the boundary id of any outflow boundary faces.
                if (integrated_flow > 0)
                  face->set_boundary_id(face->boundary_id() + offset);
              }
          }
  }


  template <int dim>
  void
  Simulator<dim>::restore_outflow_boundary_ids(const unsigned int offset)
  {
    // Loop over all of the boundary faces...
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (!cell->is_artificial())
        for (const unsigned int face_number : cell->face_indices())
          {
            const typename DoFHandler<dim>::face_iterator face = cell->face(face_number);
            if (face->at_boundary())
              {
                // ... and reset all of the boundary ids we changed in replace_outflow_boundary_ids above.
                if (face->boundary_id() >= offset)
                  face->set_boundary_id(face->boundary_id() - offset);
              }
          }
  }


  namespace
  {
    /**
     * Return whether t is an element of the given container object.
     */
    template <typename Container>
    bool is_element (const typename Container::value_type &t,
                     const Container                      &container)
    {
      for (const auto &p : container)
        if (p == t)
          return true;

      return false;
    }
  }



  template <int dim>
  void
  Simulator<dim>::check_consistency_of_boundary_conditions() const
  {
    // a container for the indicators of all boundary conditions
    std::vector<std::set<types::boundary_id>> boundary_indicator_lists;
    boundary_indicator_lists.emplace_back(boundary_velocity_manager.get_zero_boundary_velocity_indicators());
    boundary_indicator_lists.emplace_back(boundary_velocity_manager.get_tangential_boundary_velocity_indicators());
    boundary_indicator_lists.emplace_back(boundary_velocity_manager.get_prescribed_boundary_velocity_indicators());
    boundary_indicator_lists.emplace_back(boundary_traction_manager.get_prescribed_boundary_traction_indicators());

    // Make sure that each combination of boundary velocity and boundary traction condition
    // either refers to different boundary indicators or to different components
    for (const auto velocity_boundary_id: boundary_indicator_lists[2])
      {
        bool found_compatible_duplicate_boundary_id = false;
        for (const auto traction_boundary_id: boundary_indicator_lists[3])
          {
            if (velocity_boundary_id == traction_boundary_id)
              {
                // if boundary ids are identical, make sure that the components are different
                AssertThrow((boundary_velocity_manager.get_component_mask(velocity_boundary_id) &
                             boundary_traction_manager.get_component_mask(traction_boundary_id)) ==
                            ComponentMask(introspection.n_components, false),
                            ExcMessage("Boundary indicator <"
                                       +
                                       Utilities::int_to_string(velocity_boundary_id)
                                       +
                                       "> with symbolic name <"
                                       +
                                       geometry_model->translate_id_to_symbol_name (velocity_boundary_id)
                                       +
                                       "> is listed as having both "
                                       "velocity and traction boundary conditions in the input file."));

                found_compatible_duplicate_boundary_id = true;
              }
          }
        // we have ensured the prescribed velocity and prescribed traction boundary conditions
        // for the current boundary id are compatible. In order to check them against the other
        // boundary conditions, we need to remove the boundary indicator from one of the lists
        // to make sure it only appears in one of them. We choose to remove the boundary
        // indicator from the traction list. We cannot do that in the loop above, because it
        // invalidates the range of the loop.
        if (found_compatible_duplicate_boundary_id)
          boundary_indicator_lists[3].erase(velocity_boundary_id);
      }

    // for each combination of velocity boundary indicator lists, make sure that the
    // intersection is empty
    for (unsigned int i=0; i<boundary_indicator_lists.size(); ++i)
      for (unsigned int j=i+1; j<boundary_indicator_lists.size(); ++j)
        {
          std::set<types::boundary_id> intersection;
          std::set_intersection (boundary_indicator_lists[i].begin(),
                                 boundary_indicator_lists[i].end(),
                                 boundary_indicator_lists[j].begin(),
                                 boundary_indicator_lists[j].end(),
                                 std::inserter(intersection, intersection.end()));

          // if the same indicators are specified for different boundary conditions, throw exception
          AssertThrow (intersection.empty(),
                       ExcMessage ("Boundary indicator <"
                                   +
                                   Utilities::int_to_string(*intersection.begin())
                                   +
                                   "> with symbolic name <"
                                   +
                                   geometry_model->translate_id_to_symbol_name (*intersection.begin())
                                   +
                                   "> is listed as having more "
                                   "than one type of velocity or traction boundary condition in the input file."));
        }

    // make sure temperature and heat flux boundary indicators don't appear in multiple lists
    // this is easier than for the velocity/traction, as there are no selectors
    boundary_indicator_lists.emplace_back(boundary_temperature_manager.get_fixed_temperature_boundary_indicators());
    boundary_indicator_lists.emplace_back(parameters.fixed_heat_flux_boundary_indicators);

    // are there any indicators that occur in both the prescribed temperature and heat flux list?
    std::set<types::boundary_id> T_intersection;
    std::set_intersection (boundary_temperature_manager.get_fixed_temperature_boundary_indicators().begin(),
                           boundary_temperature_manager.get_fixed_temperature_boundary_indicators().end(),
                           parameters.fixed_heat_flux_boundary_indicators.begin(),
                           parameters.fixed_heat_flux_boundary_indicators.end(),
                           std::inserter(T_intersection, T_intersection.end()));

    AssertThrow (T_intersection.empty(),
                 ExcMessage ("Boundary indicator <"
                             +
                             Utilities::int_to_string(*T_intersection.begin())
                             +
                             "> with symbolic name <"
                             +
                             geometry_model->translate_id_to_symbol_name (*T_intersection.begin())
                             +
                             "> is listed as having more "
                             "than one type of temperature or heat flux boundary condition in the input file."));

    boundary_indicator_lists.emplace_back(boundary_composition_manager.get_fixed_composition_boundary_indicators());

    // Check that the periodic boundaries do not have other boundary conditions set
    using periodic_boundary_set
      = std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>;

    periodic_boundary_set pbs = geometry_model->get_periodic_boundary_pairs();

    for (const auto &pb : pbs)
      for (const auto &boundary_indicators: boundary_indicator_lists)
        {
          AssertThrow(is_element(pb.first.first, boundary_indicators) == false,
                      ExcMessage ("Boundary indicator <"
                                  +
                                  Utilities::int_to_string(pb.first.first)
                                  +
                                  "> with symbolic name <"
                                  +
                                  geometry_model->translate_id_to_symbol_name (pb.first.first)
                                  +
                                  "> is listed as having a periodic boundary condition "
                                  "in the input file, but also has another type of boundary condition. "
                                  "Periodic boundaries cannot have other boundary conditions."));

          AssertThrow(is_element(pb.first.second, boundary_indicators) == false,
                      ExcMessage ("Boundary indicator <"
                                  +
                                  Utilities::int_to_string(pb.first.second)
                                  +
                                  "> with symbolic name <"
                                  +
                                  geometry_model->translate_id_to_symbol_name (pb.first.second)
                                  +
                                  "> is listed as having a periodic boundary condition "
                                  "in the input file, but also has another type of boundary condition. "
                                  "Periodic boundaries cannot have other boundary conditions."));
        }

    const std::set<types::boundary_id> all_boundary_indicators
      = geometry_model->get_used_boundary_indicators();

    // next make sure that all listed indicators are actually used by
    // this geometry
    for (const auto &list : boundary_indicator_lists)
      for (const auto &p : list)
        AssertThrow (all_boundary_indicators.find (p)
                     != all_boundary_indicators.end(),
                     ExcMessage ("Boundary indicator <"
                                 +
                                 Utilities::int_to_string(p)
                                 +
                                 "> is listed for a boundary condition, but is not used by the geometry model."));

    if (parameters.nonlinear_solver == NonlinearSolver::single_Advection_no_Stokes)
      {
        // make sure that there are no listed velocity boundary conditions
        for (unsigned int i=0; i<4; ++i)
          AssertThrow (boundary_indicator_lists[i].empty(),
                       ExcMessage ("With the solver scheme `single Advection, no Stokes', "
                                   "one cannot set boundary conditions for velocity or traction, "
                                   "but a boundary condition has been set."));
      }
  }



  template <int dim>
  double
  Simulator<dim>::compute_initial_newton_residual()
  {
    // Store the values of current_linearization_point to be able to restore it later.
    LinearAlgebra::BlockVector temp_linearization_point = current_linearization_point;

    // Set the velocity initial guess to zero.
    current_linearization_point.block(introspection.block_indices.velocities) = 0;

    // Rebuild the whole system to compute the rhs.
    assemble_newton_stokes_system = true;
    rebuild_stokes_preconditioner = false;

    // Technically we only need the rhs, but we have asserts in place that check if
    // the system is assembled correctly when boundary conditions are prescribed, so we assemble the whole system.
    // TODO: This is a waste of time in the first nonlinear iteration. Check if we can modify the asserts in the
    // assemble_stokes_system() function to only assemble the RHS.
    rebuild_stokes_matrix = boundary_velocity_manager.get_prescribed_boundary_velocity_indicators().size()!=0;
    assemble_newton_stokes_matrix = boundary_velocity_manager.get_prescribed_boundary_velocity_indicators().size()!=0;

    compute_current_constraints ();

    assemble_stokes_system();

    const double initial_newton_residual_vel = system_rhs.block(introspection.block_indices.velocities).l2_norm();
    const double initial_newton_residual_p = system_rhs.block(introspection.block_indices.pressure).l2_norm();
    const double initial_newton_residual = std::sqrt(initial_newton_residual_vel * initial_newton_residual_vel + initial_newton_residual_p * initial_newton_residual_p);

    current_linearization_point = temp_linearization_point;

    pcout << "   Initial Newton Stokes residual = " << initial_newton_residual << ", v = " << initial_newton_residual_vel << ", p = " << initial_newton_residual_p << std::endl << std::endl;
    return initial_newton_residual;
  }



  template <int dim>
  double
  Simulator<dim>::compute_Eisenstat_Walker_linear_tolerance(const bool EisenstatWalkerChoiceOne,
                                                            const double maximum_linear_stokes_solver_tolerance,
                                                            const double linear_stokes_solver_tolerance,
                                                            const double stokes_residual,
                                                            const double newton_residual,
                                                            const double newton_residual_old)
  {
    /**
       * The Eisenstat and Walker (1996) method is used for determining the linear tolerance of
       * the iteration after the first iteration. The paper gives two preferred choices of computing
       * this tolerance. Both choices are implemented here with the suggested parameter values and
       * safeguards.
     */
    double new_linear_stokes_solver_tolerance = linear_stokes_solver_tolerance;
    if (EisenstatWalkerChoiceOne)
      {
        // This is the preferred value for this parameter in the paper.
        // A value of 2 for the power-term might also work fine.
        const double powerterm = (1+std::sqrt(5))*0.5;
        if (std::pow(linear_stokes_solver_tolerance,powerterm) <= 0.1)
          {
            new_linear_stokes_solver_tolerance = std::min(maximum_linear_stokes_solver_tolerance,
                                                          std::fabs(newton_residual-stokes_residual)/(newton_residual_old));
          }
        else
          {
            new_linear_stokes_solver_tolerance = std::min(maximum_linear_stokes_solver_tolerance,
                                                          std::max(std::fabs(newton_residual-stokes_residual)/newton_residual_old,
                                                                   std::pow(linear_stokes_solver_tolerance,powerterm)));
          }
      }
    else
      {
        if (0.9*linear_stokes_solver_tolerance * linear_stokes_solver_tolerance <= 0.1)
          {
            new_linear_stokes_solver_tolerance =  std::min(maximum_linear_stokes_solver_tolerance,
                                                           0.9 * std::fabs(newton_residual * newton_residual) /
                                                           (newton_residual_old * newton_residual_old));
          }
        else
          {
            new_linear_stokes_solver_tolerance = std::min(newton_handler->parameters.maximum_linear_stokes_solver_tolerance,
                                                          std::max(0.9 * std::fabs(newton_residual*newton_residual)
                                                                   /
                                                                   (newton_residual_old*newton_residual_old),
                                                                   0.9*linear_stokes_solver_tolerance*linear_stokes_solver_tolerance));
          }
      }
    return new_linear_stokes_solver_tolerance;
  }



  template <int dim>
  void
  Simulator<dim>::select_default_solver_and_averaging()
  {
    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::default_solver)
      {
        // Catch all situations that are not supported by the GMG solver:
        //   - Melt transport
        //   - Ellipsoidal geometry
        //   - Locally conservative discretization
        //   - Implicit reference density profile
        //   - Periodic boundaries
        //   - Stokes velocity degree not 2 or 3
        //   - Material averaging explicitly disabled
        if (parameters.include_melt_transport == true ||
            dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*>(geometry_model.get()) != nullptr ||
            parameters.use_locally_conservative_discretization == true ||
            (material_model->is_compressible() == true && parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile) ||
            (geometry_model->get_periodic_boundary_pairs().size()) > 0 ||
            (parameters.stokes_velocity_degree < 2 || parameters.stokes_velocity_degree > 3) ||
            parameters.material_averaging == MaterialModel::MaterialAveraging::none)
          {
            // GMG is not supported (yet), by default fall back to AMG.
            parameters.stokes_solver_type = Parameters<dim>::StokesSolverType::block_amg;
          }
        else
          {
            // GMG is supported for all other cases
            parameters.stokes_solver_type = Parameters<dim>::StokesSolverType::block_gmg;
          }
      }

    // Now pick an appropriate material averaging for the chosen solver
    if (parameters.material_averaging == MaterialModel::MaterialAveraging::default_averaging)
      {
        if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
          {
            // project to Q1 is more accurate, but not supported if:
            //   - elasticity is enabled
            //   - the Newton solver is enabled
            if (parameters.enable_elasticity == true ||
                parameters.nonlinear_solver == NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
                parameters.nonlinear_solver == NonlinearSolver::single_Advection_iterated_Newton_Stokes)
              parameters.material_averaging = MaterialModel::MaterialAveraging::harmonic_average_only_viscosity;
            else
              parameters.material_averaging = MaterialModel::MaterialAveraging::project_to_Q1_only_viscosity;
          }
        else
          parameters.material_averaging = MaterialModel::MaterialAveraging::none;
      }
  }
}

// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Simulator<dim>::AdvectionField; \
  template double Simulator<dim>::normalize_pressure(LinearAlgebra::BlockVector &vector) const; \
  template void Simulator<dim>::denormalize_pressure(const double pressure_adjustment, \
                                                     LinearAlgebra::BlockVector &vector) const; \
  template double Simulator<dim>::compute_pressure_scaling_factor () const; \
  template double Simulator<dim>::get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const; \
  template std::pair<double,double> Simulator<dim>::get_extrapolated_advection_field_range (const AdvectionField &advection_field) const; \
  template void Simulator<dim>::maybe_write_timing_output () const; \
  template bool Simulator<dim>::maybe_write_checkpoint (const time_t, const bool); \
  template bool Simulator<dim>::maybe_do_initial_refinement (const unsigned int max_refinement_level); \
  template void Simulator<dim>::exchange_refinement_flags (); \
  template void Simulator<dim>::maybe_refine_mesh (const double new_time_step, unsigned int &max_refinement_level); \
  template void Simulator<dim>::advance_time (const double step_size); \
  template void Simulator<dim>::make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector); \
  template void Simulator<dim>::output_statistics(); \
  template void Simulator<dim>::write_plugin_graph(std::ostream &) const; \
  template double Simulator<dim>::compute_initial_stokes_residual(); \
  template bool Simulator<dim>::stokes_matrix_depends_on_solution() const; \
  template bool Simulator<dim>::stokes_A_block_is_symmetric() const; \
  template void Simulator<dim>::interpolate_onto_velocity_system(const TensorFunction<1,dim> &func, LinearAlgebra::Vector &vec) const;\
  template void Simulator<dim>::apply_limiter_to_dg_solutions(const AdvectionField &advection_field); \
  template void Simulator<dim>::compute_unique_advection_support_points(const std::vector<AdvectionField> &advection_fields, \
                                                                        std::vector<Point<dim>> &support_points, \
                                                                        std::vector<std::vector<unsigned int>> &support_point_index_by_field) const; \
  template void Simulator<dim>::compute_reactions(); \
  template void Simulator<dim>::initialize_current_linearization_point (); \
  template void Simulator<dim>::interpolate_material_output_into_advection_field(const std::vector<AdvectionField> &adv_field); \
  template void Simulator<dim>::check_consistency_of_formulation(); \
  template void Simulator<dim>::replace_outflow_boundary_ids(const unsigned int boundary_id_offset); \
  template void Simulator<dim>::restore_outflow_boundary_ids(const unsigned int boundary_id_offset); \
  template void Simulator<dim>::check_consistency_of_boundary_conditions() const; \
  template double Simulator<dim>::compute_initial_newton_residual(); \
  template double Simulator<dim>::compute_Eisenstat_Walker_linear_tolerance(const bool EisenstatWalkerChoiceOne, \
                                                                            const double maximum_linear_stokes_solver_tolerance, \
                                                                            const double linear_stokes_solver_tolerance, \
                                                                            const double stokes_residual, \
                                                                            const double newton_residual, \
                                                                            const double newton_residual_old); \
  template void Simulator<dim>::select_default_solver_and_averaging();


  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}

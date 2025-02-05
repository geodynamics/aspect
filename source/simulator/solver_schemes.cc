/*
  Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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
#include <aspect/global.h>
#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/newton.h>
#include <aspect/melt.h>

#include <deal.II/numerics/vector_tools.h>

#include <aspect/stokes_matrix_free.h>


namespace aspect
{

  namespace
  {
    /**
     * Converts a function with a certain number of components into a Function@<dim@>
     * with optionally having additional zero components.
     **/
    template <int dim>
    class VectorFunctionFromVectorFunctionObject : public Function<dim>
    {
      public:
        /**
         * Converts a function with @p n_object_components components into a Function@dim@
         * while optionally providing additional components that are set to zero.
         *
         * @param function_object The function that will form the components
         *     of the resulting Function object.
         * @param first_component The first component that should be
         *     filled.
         * @param n_object_components The number of components that should be
         *     filled from the first.
         * @param n_total_components The total number of vector components of the
         *     resulting Function object.
         **/

        VectorFunctionFromVectorFunctionObject (const std::function<void (const Point<dim> &,Vector<double> &)> &function_object,
                                                const unsigned int first_component,
                                                const unsigned int n_object_components,
                                                const unsigned int n_total_components)
          :
          Function<dim>(n_total_components),
          function_object (function_object),
          first_component (first_component),
          n_object_components (n_object_components)
        {
          Assert ((n_object_components > 0
                   &&
                   first_component+n_object_components <= n_total_components),
                  ExcMessage ("Number of objects components needs to be less than number of total components"));
        }



        double
        value (const Point<dim> &p,
               const unsigned int component) const override
        {
          Assert (component < this->n_components,
                  ExcIndexRange (component, 0, this->n_components));

          if (component < first_component)
            return 0;
          else if (component >= first_component + n_object_components)
            return 0;
          else
            {
              Vector<double> temp(n_object_components);
              function_object (p, temp);
              return temp(component - first_component);
            }
        }



        void
        vector_value (const Point<dim>   &p,
                      Vector<double>     &values) const override
        {
          AssertDimension(values.size(), this->n_components);

          // set everything to zero, and then the right components to their correct values
          values = 0;
          Vector<double> temp(n_object_components);
          function_object (p, temp);
          for (unsigned int i = 0; i < n_object_components; ++i)
            {
              values(first_component + i) = temp(i);
            }
        }



      private:
        /**
         * The function object which we call when this class's solution() function is called.
         **/
        const std::function<void (const Point<dim> &,Vector<double> &)> function_object;

        /**
         * The first vector component whose value is to be filled by the given
         * function.
         */
        const unsigned int first_component;
        /**
         * The number of vector components whose values are to be filled by the given
         * function.
         */
        const unsigned int n_object_components;

    };

  }



  template <int dim>
  double Simulator<dim>::assemble_and_solve_temperature (const double &initial_residual,
                                                         double *residual)
  {
    double current_residual = 0.0;

    switch (parameters.temperature_method)
      {
        case Parameters<dim>::AdvectionFieldMethod::fem_field:
        case Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion:
        {
          const AdvectionField adv_field (AdvectionField::temperature());

          // if this is a prescribed field with diffusion, we first have to copy the material model
          // outputs into the prescribed field before we assemble and solve the equation
          if (parameters.temperature_method == Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion)
            {
              TimerOutput::Scope timer (computing_timer, "Interpolate prescribed temperature");

              interpolate_material_output_into_advection_field(adv_field);

              // Also set the old_solution block to the prescribed field. The old
              // solution is the one that is used to assemble the diffusion system in
              // assemble_advection_system() for this solver scheme.
              old_solution.block(adv_field.block_index(introspection)) = solution.block(adv_field.block_index(introspection));
            }

          assemble_advection_system (adv_field);

          if (residual)
            *residual = system_rhs.block(introspection.block_indices.temperature).l2_norm();

          current_residual = solve_advection(adv_field);
          break;
        }

        case Parameters<dim>::AdvectionFieldMethod::prescribed_field:
        {
          const AdvectionField adv_field (AdvectionField::temperature());

          TimerOutput::Scope timer (computing_timer, "Interpolate prescribed temperature");

          interpolate_material_output_into_advection_field(adv_field);

          // Call the signal in case the user wants to do something with the variable:
          SolverControl dummy;
          signals.post_advection_solver(*this,
                                        adv_field.is_temperature(),
                                        adv_field.compositional_variable,
                                        dummy);
          break;
        }

        case Parameters<dim>::AdvectionFieldMethod::static_field:
        {
          const AdvectionField adv_field (AdvectionField::temperature());
          // Do nothing here, but at least call the signal in case the
          // user wants to do something with the variable:
          SolverControl dummy;
          signals.post_advection_solver(*this,
                                        adv_field.is_temperature(),
                                        adv_field.compositional_variable,
                                        dummy);
          break;
        }

        default:
          AssertThrow(false,ExcNotImplemented());
      }

    current_linearization_point.block(introspection.block_indices.temperature)
      = solution.block(introspection.block_indices.temperature);

    if (initial_residual > 0)
      return current_residual / initial_residual;

    return 0.0;
  }



  template <int dim>
  std::vector<double>
  Simulator<dim>::assemble_and_solve_composition (const std::vector<double> &initial_residual,
                                                  std::vector<double> *residual)
  {
    // Advect the particles before they are potentially used to
    // set up the compositional fields.
    for (auto &particle_manager : particle_managers)
      {
        // Do not advect the particles in the initial refinement stage
        const bool in_initial_refinement = (timestep_number == 0)
                                           && (pre_refinement_step < parameters.initial_adaptive_refinement);
        if (!in_initial_refinement)
          // Advance the particles in the manager to the current time
          particle_manager.advance_timestep();

        if (particle_manager.get_property_manager().need_update() == Particle::Property::update_output_step)
          particle_manager.update_particles();
      }

    std::vector<double> current_residual(introspection.n_compositional_fields,0.0);

    if (residual)
      Assert(residual->size() == introspection.n_compositional_fields, ExcInternalError());

    std::vector<AdvectionField> fields_advected_by_particles;

    for (unsigned int c=0; c < introspection.n_compositional_fields; ++c)
      {
        const AdvectionField adv_field (AdvectionField::composition(c));
        const typename Parameters<dim>::AdvectionFieldMethod::Kind method = adv_field.advection_method(introspection);
        switch (method)
          {
            case Parameters<dim>::AdvectionFieldMethod::fem_field:
            case Parameters<dim>::AdvectionFieldMethod::fem_melt_field:
            case Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion:
            case Parameters<dim>::AdvectionFieldMethod::fem_darcy_field:
            {
              // if this is a prescribed field with diffusion, we first have to copy the material model
              // outputs into the prescribed field before we assemble and solve the equation
              if (method == Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion)
                {
                  TimerOutput::Scope timer (computing_timer, "Interpolate prescribed composition");

                  interpolate_material_output_into_advection_field(adv_field);

                  // Also set the old_solution block to the prescribed field. The old
                  // solution is the one that is used to assemble the diffusion system in
                  // assemble_advection_system() for this solver scheme.
                  old_solution.block(adv_field.block_index(introspection)) = solution.block(adv_field.block_index(introspection));
                }

              assemble_advection_system (adv_field);

              if (residual)
                (*residual)[c] = system_rhs.block(introspection.block_indices.compositional_fields[c]).l2_norm();

              current_residual[c] = solve_advection(adv_field);

              // Release the contents of the matrix block we used again:
              const unsigned int block_idx = adv_field.block_index(introspection);
              if (adv_field.sparsity_pattern_block_index(introspection)!=block_idx)
                system_matrix.block(block_idx, block_idx).clear();

              // No need to call the post_advection_solver signal here: It is
              // automatically called from solve_advection() above.

              break;
            }

            case Parameters<dim>::AdvectionFieldMethod::particles:
            {
              // handle all particle fields together to increase efficiency
              fields_advected_by_particles.push_back(adv_field);
              break;
            }

            case Parameters<dim>::AdvectionFieldMethod::prescribed_field:
            {
              TimerOutput::Scope timer (computing_timer, "Interpolate prescribed composition");

              interpolate_material_output_into_advection_field(adv_field);

              // Call the signal in case the user wants to do something with the variable:
              SolverControl dummy;
              signals.post_advection_solver(*this,
                                            adv_field.is_temperature(),
                                            adv_field.compositional_variable,
                                            dummy);
              break;
            }

            case Parameters<dim>::AdvectionFieldMethod::volume_of_fluid:
            {
              volume_of_fluid_handler->do_volume_of_fluid_update(adv_field);
              break;
            }

            case Parameters<dim>::AdvectionFieldMethod::static_field:
            {
              // Do nothing here, but at least call the signal in case the
              // user wants to do something with the variable:
              SolverControl dummy;
              signals.post_advection_solver(*this,
                                            adv_field.is_temperature(),
                                            adv_field.compositional_variable,
                                            dummy);
              break;
            }

            default:
              AssertThrow(false,ExcNotImplemented());
          }
      }

    if (fields_advected_by_particles.size() > 0)
      interpolate_particle_properties(fields_advected_by_particles);


    // for consistency we update the current linearization point only after we have solved
    // all fields, so that we use the same point in time for every field when solving
    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      {
        current_linearization_point.block(introspection.block_indices.compositional_fields[c])
          = solution.block(introspection.block_indices.compositional_fields[c]);

        if ((initial_residual.size() > 0) && (initial_residual[c] > 0))
          current_residual[c] /= initial_residual[c];
        else
          current_residual[c] = 0.0;
      }

    return current_residual;
  }



  template <int dim>
  double
  Simulator<dim>::assemble_and_solve_stokes (const double &initial_nonlinear_residual,
                                             double *nonlinear_residual)
  {
    // If the Stokes matrix depends on the solution, or we have active
    // velocity boundary conditions, we need to re-assemble the system matrix
    // (and preconditioner) every time. If we have active boundary conditions,
    // they could a) depend on the solution, or b) be inhomogeneous. In both
    // cases, just assembling the RHS will be incorrect.  If no active
    // boundaries exist, we only have no-slip or free slip conditions, so we
    // don't need to force assembly of the matrix.
    if (stokes_matrix_depends_on_solution()
        ||
        (boundary_velocity_manager.get_prescribed_boundary_velocity_indicators().size() > 0)
        || parameters.mesh_deformation_enabled
       )
      rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

    // set constraints for p_c if porosity is below a threshold
    if (nonlinear_iteration == 0 && parameters.include_melt_transport)
      {
        compute_current_constraints();
        if (rebuild_sparsity_and_matrices)
          {
            setup_system_matrix (introspection.index_sets.system_partitioning);
            setup_system_preconditioner (introspection.index_sets.system_partitioning);

            rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;
          }
      }

    // Re-compute the pressure scaling factor for the Stokes assembly
    pressure_scaling = compute_pressure_scaling_factor();
    assemble_stokes_system ();

    // build the preconditioner
    if (stokes_matrix_free)
      stokes_matrix_free->build_preconditioner();
    else
      build_stokes_preconditioner();

    if (nonlinear_residual)
      *nonlinear_residual = compute_initial_stokes_residual();

    const double current_nonlinear_residual = solve_stokes(solution).first;

    current_linearization_point.block(introspection.block_indices.velocities)
      = solution.block(introspection.block_indices.velocities);

    if (introspection.block_indices.velocities != introspection.block_indices.pressure)
      current_linearization_point.block(introspection.block_indices.pressure)
        = solution.block(introspection.block_indices.pressure);

    if (parameters.include_melt_transport)
      {
        // Note that the compaction pressure is in the fluid pressure block
        // and will therefore be updated as well.
        const unsigned int fluid_velocity_block = introspection.variable("fluid velocity").block_index;
        const unsigned int fluid_pressure_block = introspection.variable("fluid pressure").block_index;
        current_linearization_point.block(fluid_velocity_block) = solution.block(fluid_velocity_block);
        current_linearization_point.block(fluid_pressure_block) = solution.block(fluid_pressure_block);
      }

    if (initial_nonlinear_residual > 0)
      return current_nonlinear_residual / initial_nonlinear_residual;
    else
      return 0.0;
  }


  template <int dim>
  void Simulator<dim>::do_one_defect_correction_Stokes_step(DefectCorrectionResiduals &dcr,
                                                            const bool use_picard)
  {
    // Many parts of the solver depend on the block layout (velocity = 0,
    // pressure = 1). For example the linearized_stokes_initial_guess vector or the StokesBlock matrix
    // wrapper. Let us make sure that this holds (and shorten their names):
    const unsigned int pressure_block_index = (parameters.include_melt_transport) ?
                                              introspection.variable("fluid pressure").block_index
                                              : introspection.block_indices.pressure;
    const unsigned int velocity_block_index = introspection.block_indices.velocities;
    Assert(velocity_block_index == 0, ExcNotImplemented());
    Assert(pressure_block_index == 1, ExcNotImplemented());
    (void) pressure_block_index;
    Assert(!parameters.include_melt_transport
           || introspection.variable("compaction pressure").block_index == 1, ExcNotImplemented());

    // Re-compute the pressure scaling factor for the Stokes assembly
    pressure_scaling = compute_pressure_scaling_factor();

    if (nonlinear_iteration == 0)
      {
        dcr.initial_residual = compute_initial_newton_residual();
        dcr.switch_initial_residual = dcr.initial_residual;
        dcr.residual_old = dcr.initial_residual;
        dcr.residual = dcr.initial_residual;
        assemble_newton_stokes_system = assemble_newton_stokes_matrix = false;
      }
    else
      assemble_newton_stokes_system = assemble_newton_stokes_matrix = true;

    // Make sure we use the correct assemblers for first and all other iterations
    // We solve the normal system in the first iteration, and the defect correction
    // system in all other iterations.
    if (nonlinear_iteration <= 1)
      {
        set_assemblers();
        compute_current_constraints ();
      }

    // If the Stokes matrix depends on the solution, or we have active
    // velocity boundary conditions, we need to re-assemble the system matrix
    // (and preconditioner) every time. If we have active boundary conditions,
    // they could a) depend on the solution, or b) be inhomogeneous. In both
    // cases, just assembling the RHS will be incorrect.  If no active
    // boundaries exist, we only have no-slip or free slip conditions, so we
    // don't need to force assembly of the matrix.
    if (stokes_matrix_depends_on_solution()
        ||
        nonlinear_iteration == 0)
      rebuild_stokes_matrix = rebuild_stokes_preconditioner = assemble_newton_stokes_matrix = true;
    else if (parameters.enable_prescribed_dilation)
      // The dilation requires the Stokes matrix (which is on the rhs
      // in the Newton solver) to be updated.
      rebuild_stokes_matrix = true;

    assemble_stokes_system();

    // recompute rhs
    dcr.velocity_residual = system_rhs.block(velocity_block_index).l2_norm();
    dcr.pressure_residual = system_rhs.block(pressure_block_index).l2_norm();
    dcr.residual = std::sqrt(dcr.velocity_residual * dcr.velocity_residual + dcr.pressure_residual * dcr.pressure_residual);

    // Test whether the rhs has dropped so much that we can assume that the iteration is done.
    if (dcr.residual < dcr.residual_old * 1e-8)
      {
        pcout << "   Nonlinear residual reduction has been very large ("
              << dcr.residual/dcr.residual_old << "); skipping Stokes solve"
              << std::endl;

        // Update old residual so that if the nonlinear solver is not done,
        // next time we will not skip the Stokes solve.
        dcr.residual_old = dcr.residual;
        return;
      }


    // Eisenstat Walker method for determining the linear solver tolerance
    if (nonlinear_iteration > 1)
      {
        if (!use_picard || newton_handler->parameters.use_Eisenstat_Walker_method_for_Picard_iterations)
          {
            const bool EisenstatWalkerChoiceOne = true;
            parameters.linear_stokes_solver_tolerance = compute_Eisenstat_Walker_linear_tolerance(EisenstatWalkerChoiceOne,
                                                        newton_handler->parameters.maximum_linear_stokes_solver_tolerance,
                                                        parameters.linear_stokes_solver_tolerance,
                                                        dcr.stokes_residuals.second,
                                                        dcr.residual,
                                                        dcr.residual_old);

            pcout << "   The linear solver tolerance is set to "
                  << parameters.linear_stokes_solver_tolerance
                  << ". ";
            if (!use_picard)
              {
                pcout << "Stabilization Preconditioner is "
                      << Newton::to_string(newton_handler->parameters.preconditioner_stabilization)
                      << " and A block is "
                      << Newton::to_string(newton_handler->parameters.velocity_block_stabilization)
                      << '.';
              }
            pcout << std::endl;
          }
      }

    if (stokes_matrix_free)
      stokes_matrix_free->build_preconditioner();
    else
      build_stokes_preconditioner();

    // The solution of the defect correction solver is an update direction
    LinearAlgebra::BlockVector search_direction;
    search_direction.reinit(solution, /* omit_zeroing_entries = */ true);

    try
      {
        dcr.stokes_residuals = solve_stokes(search_direction);
      }
    catch (const std::exception &exc)
      {
        // Test that we are trying to handle exceptions and that
        // the exception we got is one of the two documented by
        // throw_linear_solver_failure_exception(). If not, we have a genuine
        // problem here, and will need to get outta here right away:
        if (newton_handler->parameters.use_Newton_failsafe == false ||
            ((dynamic_cast<const ExcMessage *>(&exc)==nullptr) &&
             (dynamic_cast<const QuietException *>(&exc)==nullptr)))
          throw;

        // start the solve over again and try with a stabilized version
        pcout << "failed, trying again with stabilization" << std::endl;
        newton_handler->parameters.preconditioner_stabilization = Newton::Parameters::Stabilization::SPD;
        newton_handler->parameters.velocity_block_stabilization = Newton::Parameters::Stabilization::SPD;

        // If the Stokes matrix depends on the solution, or we have active
        // velocity boundary conditions, we need to re-assemble the system matrix
        // (and preconditioner) every time. If we have active boundary conditions,
        // they could a) depend on the solution, or b) be inhomogeneous. In both
        // cases, just assembling the RHS will be incorrect.  If no active
        // boundaries exist, we only have no-slip or free slip conditions, so we
        // don't need to force assembly of the matrix.
        if (stokes_matrix_depends_on_solution()
            ||
            (nonlinear_iteration == 0 && boundary_velocity_manager.get_prescribed_boundary_velocity_indicators().size() > 0))
          rebuild_stokes_matrix = rebuild_stokes_preconditioner = assemble_newton_stokes_matrix = true;
        else if (parameters.enable_prescribed_dilation)
          // The dilation requires the Stokes matrix (which is on the rhs
          // in the Newton solver) to be updated.
          rebuild_stokes_matrix = true;

        assemble_stokes_system();

        /**
         * Eisenstat Walker method for determining the tolerance
         */
        if (nonlinear_iteration > 1)
          {
            dcr.velocity_residual = system_rhs.block(velocity_block_index).l2_norm();
            dcr.pressure_residual = system_rhs.block(pressure_block_index).l2_norm();
            dcr.residual = std::sqrt(dcr.velocity_residual * dcr.velocity_residual + dcr.pressure_residual * dcr.pressure_residual);

            // Eisenstat Walker method for determining the linear solver tolerance
            if (!use_picard)
              {
                const bool EisenstatWalkerChoiceOne = true;
                parameters.linear_stokes_solver_tolerance = compute_Eisenstat_Walker_linear_tolerance(EisenstatWalkerChoiceOne,
                                                            newton_handler->parameters.maximum_linear_stokes_solver_tolerance,
                                                            parameters.linear_stokes_solver_tolerance,
                                                            dcr.stokes_residuals.second,
                                                            dcr.residual,
                                                            dcr.residual_old);

                pcout << "   The linear solver tolerance is set to " << parameters.linear_stokes_solver_tolerance << std::endl;
              }
          }

        if (stokes_matrix_free)
          stokes_matrix_free->build_preconditioner();
        else
          build_stokes_preconditioner();

        // Give this another try:
        dcr.stokes_residuals = solve_stokes(search_direction);
      }

    // Apply the solution or the update of the solution
    if (nonlinear_iteration == 0)
      {
        // In the first nonlinear iteration, the nonlinear residual and old
        // residual should be the same, both depend on the accuracy of
        // the initial guess.
        dcr.residual = dcr.stokes_residuals.first;
        dcr.residual_old = dcr.residual;

        // The first nonlinear iteration we are computing the whole system in a non-defect corrected Picard way,
        // to make sure that the boundary conditions are correct in combination with correct initial guesses.
        current_linearization_point.block(velocity_block_index) = search_direction.block(velocity_block_index);
        current_linearization_point.block(pressure_block_index) = search_direction.block(pressure_block_index);

        pcout << "      Newton system information: Norm of the rhs: " << dcr.stokes_residuals.first << std::endl;
      }
    else
      {
        // Backup the residual of the previous iteration, before we compute the new residual.
        dcr.residual_old = dcr.residual;

        // We have computed an update. Prepare for a line search.
        LinearAlgebra::BlockVector backup_linearization_point(introspection.index_sets.stokes_partitioning, mpi_communicator);

        double step_length_factor = 1.0;
        unsigned int line_search_iteration = 0;

        // Do the loop for the line search at least once with the full step length.
        // If line search is disabled we will exit the loop in the first iteration.
        do
          {
            if (line_search_iteration == 0)
              {
                // backup our starting point for the line search
                backup_linearization_point.block(pressure_block_index) = current_linearization_point.block(pressure_block_index);
                backup_linearization_point.block(velocity_block_index) = current_linearization_point.block(velocity_block_index);
              }
            else
              {
                // undo the last iteration and try again with smaller step length
                current_linearization_point.block(pressure_block_index) = backup_linearization_point.block(pressure_block_index);
                current_linearization_point.block(velocity_block_index) = backup_linearization_point.block(velocity_block_index);
                search_direction.block(pressure_block_index) *= step_length_factor;
                search_direction.block(velocity_block_index) *= step_length_factor;
              }

            // Update the current linearization point with the search direction
            current_linearization_point.block(pressure_block_index) += search_direction.block(pressure_block_index);
            current_linearization_point.block(velocity_block_index) += search_direction.block(velocity_block_index);

            // Rebuild the rhs to determine the new residual.
            assemble_newton_stokes_matrix = rebuild_stokes_preconditioner = false;
            rebuild_stokes_matrix = !boundary_velocity_manager.get_prescribed_boundary_velocity_indicators().empty();

            assemble_stokes_system();

            const double test_velocity_residual = system_rhs.block(velocity_block_index).l2_norm();
            const double test_pressure_residual = system_rhs.block(pressure_block_index).l2_norm();
            const double test_residual = std::sqrt(test_velocity_residual * test_velocity_residual
                                                   + test_pressure_residual * test_pressure_residual);

            // Determine if the residual has decreased sufficiently.
            const double alpha = 1e-4;
            if (test_residual < (1.0 - alpha * step_length_factor) * dcr.residual
                ||
                line_search_iteration >= newton_handler->parameters.max_newton_line_search_iterations
                ||
                use_picard)
              {
                pcout << "      Newton system information: Norm of the rhs: " << test_residual
                      << ", Derivative scaling factor: " << newton_handler->parameters.newton_derivative_scaling_factor << std::endl;
                dcr.residual = test_residual;
                break;
              }
            else
              {
                pcout << "   Line search iteration " << line_search_iteration << ", with norm of the rhs "
                      << test_residual << " and going to " << (1.0 - alpha * step_length_factor) * dcr.residual
                      << ", relative residual: " << test_residual/dcr.initial_residual << std::endl;

                // The current search direction has not decreased the residual
                // enough, so we take a smaller step and try again.
                // TODO: make a parameter out of this.
                step_length_factor = 2.0/3.0;
              }

            ++line_search_iteration;
            Assert(line_search_iteration <= newton_handler->parameters.max_newton_line_search_iterations,
                   ExcInternalError());
          }
        while (true);
      }


    if (use_picard == true)
      {
        // When we are using (defect corrected) Picard, keep the
        // newton_derivative_scaling_factor at zero. The newton_derivative_scaling_factor
        // is calculated above as std::max(0.0, (1.0-(
        // newton_residual_for_derivative_scaling_factor/switch_initial_residual)))
        dcr.switch_initial_residual = dcr.residual;
        dcr.newton_residual_for_derivative_scaling_factor = dcr.residual;
      }
    else
      {
        // This method allows to slowly introduce the derivatives based
        // on the improvement of the residual. This method was suggested
        // by Raid Hassani.
        if (newton_handler->parameters.use_newton_residual_scaling_method)
          dcr.newton_residual_for_derivative_scaling_factor = dcr.residual;
        else
          dcr.newton_residual_for_derivative_scaling_factor = 0;
      }

    if (nonlinear_iteration != 0)
      last_pressure_normalization_adjustment = normalize_pressure(current_linearization_point);

    // The rest of ASPECT assumes 'solution' contains the solution values, not the update
    // to the solution. Make sure the outside never sees the solution update we computed.
    solution.block(pressure_block_index) = current_linearization_point.block(pressure_block_index);
    solution.block(velocity_block_index) = current_linearization_point.block(velocity_block_index);
  }



  template <int dim>
  void Simulator<dim>::solve_single_advection_single_stokes ()
  {
    assemble_and_solve_temperature();
    assemble_and_solve_composition();
    assemble_and_solve_stokes();

    if (parameters.run_postprocessors_on_nonlinear_iterations)
      postprocess ();

    // Setup a nonlinear solver control that only allows a single iteration
    SolverControl nonlinear_solver_control(1,1.0);
    // Announce that we did a single iteration, and assume we have converged
    nonlinear_solver_control.check(1,0.0);
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_no_advection_iterated_stokes ()
  {
    double initial_stokes_residual = 0.0;

    const unsigned int max_nonlinear_iterations =
      (pre_refinement_step < parameters.initial_adaptive_refinement)
      ?
      std::min(parameters.max_nonlinear_iterations,
               parameters.max_nonlinear_iterations_in_prerefinement)
      :
      parameters.max_nonlinear_iterations;

    SolverControl nonlinear_solver_control(max_nonlinear_iterations,
                                           parameters.nonlinear_tolerance);

    double relative_residual = std::numeric_limits<double>::max();
    nonlinear_iteration = 0;
    do
      {
        relative_residual =
          assemble_and_solve_stokes(initial_stokes_residual,
                                    nonlinear_iteration == 0 ? &initial_stokes_residual : nullptr);

        pcout << "      Relative nonlinear residual (Stokes system) after nonlinear iteration " << nonlinear_iteration+1
              << ": " << relative_residual
              << std::endl
              << std::endl;

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        ++nonlinear_iteration;
      }
    while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);

    AssertThrow(nonlinear_solver_control.last_check() != SolverControl::failure, ExcNonlinearSolverNoConvergence());
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_no_advection_iterated_defect_correction_stokes ()
  {
    // Now store the linear_tolerance we started out with, because we might change
    // it within this timestep.
    const double begin_linear_tolerance = parameters.linear_stokes_solver_tolerance;

    DefectCorrectionResiduals dcr;
    dcr.initial_residual = 1;

    dcr.velocity_residual = 0;
    dcr.pressure_residual = 0;
    dcr.residual = 1;
    dcr.residual_old = 1;

    dcr.switch_initial_residual = 1;
    dcr.newton_residual_for_derivative_scaling_factor = 1;

    const unsigned int max_nonlinear_iterations =
      (pre_refinement_step < parameters.initial_adaptive_refinement)
      ?
      std::min(parameters.max_nonlinear_iterations,
               parameters.max_nonlinear_iterations_in_prerefinement)
      :
      parameters.max_nonlinear_iterations;

    SolverControl nonlinear_solver_control(max_nonlinear_iterations,
                                           parameters.nonlinear_tolerance);

    // Now iterate out the nonlinearities.
    dcr.stokes_residuals = std::pair<double,double>  (numbers::signaling_nan<double>(),
                                                      numbers::signaling_nan<double>());

    double relative_residual = std::numeric_limits<double>::max();
    nonlinear_iteration = 0;
    do
      {
        do_one_defect_correction_Stokes_step(dcr, /* use_picard= */true);

        relative_residual = dcr.residual/dcr.initial_residual;

        pcout << "      Relative nonlinear residual (Stokes system) after nonlinear iteration " << nonlinear_iteration+1
              << ": " << relative_residual
              << std::endl
              << std::endl;

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        ++nonlinear_iteration;
      }
    while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);

    AssertThrow(nonlinear_solver_control.last_check() != SolverControl::failure, ExcNonlinearSolverNoConvergence());

    // Reset the linear tolerance to what it was at the beginning of the time step.
    parameters.linear_stokes_solver_tolerance = begin_linear_tolerance;

    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_single_advection_iterated_defect_correction_stokes ()
  {
    // The defect correction solver is just the Newton solver without derivatives.
    solve_single_advection_and_iterated_newton_stokes(/*use_newton_iterations = */ false);
  }

  template <int dim>
  void Simulator<dim>::solve_iterated_advection_and_defect_correction_stokes ()
  {
    // The defect correction solver is just the Newton solver without derivatives.
    solve_iterated_advection_and_newton_stokes(/*use_newton_iterations = */ false);
  }



  template <int dim>
  void Simulator<dim>::solve_no_advection_single_stokes ()
  {
    assemble_and_solve_stokes();

    if (parameters.run_postprocessors_on_nonlinear_iterations)
      postprocess ();

    // Setup a nonlinear solver control that only allows a single iteration
    SolverControl nonlinear_solver_control(1,1.0);
    // Announce that we did a single iteration, and assume we have converged
    nonlinear_solver_control.check(1,0.0);
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }

  template <int dim>
  void Simulator<dim>::solve_first_timestep_only_single_stokes ()
  {
    if (timestep_number == 0)
      assemble_and_solve_stokes();

    // Setup a nonlinear solver control that only allows a single iteration
    SolverControl nonlinear_solver_control(1,1.0);
    // Announce that we did a single iteration, and assume we have converged
    nonlinear_solver_control.check(1,0.0);
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_iterated_advection_and_stokes ()
  {
    double initial_temperature_residual = 0;
    double initial_stokes_residual      = 0;
    std::vector<double> initial_composition_residual (introspection.n_compositional_fields,0);

    const unsigned int max_nonlinear_iterations =
      (pre_refinement_step < parameters.initial_adaptive_refinement)
      ?
      std::min(parameters.max_nonlinear_iterations,
               parameters.max_nonlinear_iterations_in_prerefinement)
      :
      parameters.max_nonlinear_iterations;

    SolverControl nonlinear_solver_control(max_nonlinear_iterations,
                                           parameters.nonlinear_tolerance);

    double relative_residual = std::numeric_limits<double>::max();
    nonlinear_iteration = 0;

    do
      {
        // Restore particles through stored copy of particle handler,
        // but only if they have already been displaced in a nonlinear
        // iteration (in the assemble_and_solve_composition call).

        if (nonlinear_iteration > 0)
          for (auto &particle_manager : particle_managers)
            particle_manager.restore_particles();

        const double relative_temperature_residual =
          assemble_and_solve_temperature(initial_temperature_residual,
                                         nonlinear_iteration == 0 ? &initial_temperature_residual : nullptr);

        const std::vector<double>  relative_composition_residual =
          assemble_and_solve_composition(initial_composition_residual,
                                         nonlinear_iteration == 0 ? &initial_composition_residual : nullptr);

        const double relative_nonlinear_stokes_residual =
          assemble_and_solve_stokes(initial_stokes_residual,
                                    nonlinear_iteration == 0 ? &initial_stokes_residual : nullptr);

        // write the residual output in the same order as the solutions
        pcout << "      Relative nonlinear residuals (temperature"
              << (introspection.n_compositional_fields > 0 ? ", compositional fields" : "")
              << ", Stokes system): " << relative_temperature_residual;
        for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
          pcout << ", " << relative_composition_residual[c];
        pcout << ", " << relative_nonlinear_stokes_residual;
        pcout << std::endl;

        // Find the maximum residual of the individual equations
        relative_residual = relative_temperature_residual;
        for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
          {
            // in models with melt migration the melt advection equation includes the divergence of the velocity
            // and can not be expected to converge to a smaller value than the residual of the Stokes equation.
            // thus, we set a threshold for the initial composition residual.
            // this only plays a role if the right-hand side of the advection equation is very small.
            const double threshold = (parameters.include_melt_transport && c == introspection.compositional_index_for_name("porosity")
                                      ?
                                      parameters.linear_stokes_solver_tolerance * time_step
                                      :
                                      0.0);
            if (initial_composition_residual[c]>threshold)
              relative_residual = std::max(relative_composition_residual[c],relative_residual);
          }
        relative_residual = std::max(relative_nonlinear_stokes_residual, relative_residual);

        pcout << "      Relative nonlinear residual (total system) after nonlinear iteration " << nonlinear_iteration+1
              << ": " << relative_residual
              << std::endl
              << std::endl;

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        ++nonlinear_iteration;
      }
    while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);

    AssertThrow(nonlinear_solver_control.last_check() != SolverControl::failure, ExcNonlinearSolverNoConvergence());
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_single_advection_iterated_stokes ()
  {
    // solve the temperature and composition systems once...
    assemble_and_solve_temperature();

    assemble_and_solve_composition();

    // ...and then iterate the solution of the Stokes system
    double initial_stokes_residual = 0;

    const unsigned int max_nonlinear_iterations =
      (pre_refinement_step < parameters.initial_adaptive_refinement)
      ?
      std::min(parameters.max_nonlinear_iterations,
               parameters.max_nonlinear_iterations_in_prerefinement)
      :
      parameters.max_nonlinear_iterations;

    SolverControl nonlinear_solver_control(max_nonlinear_iterations,
                                           parameters.nonlinear_tolerance);

    double relative_residual = std::numeric_limits<double>::max();
    nonlinear_iteration = 0;
    do
      {
        relative_residual =
          assemble_and_solve_stokes(initial_stokes_residual,
                                    nonlinear_iteration == 0 ? &initial_stokes_residual : nullptr);

        pcout << "      Relative nonlinear residual (Stokes system) after nonlinear iteration " << nonlinear_iteration+1
              << ": " << relative_residual
              << std::endl
              << std::endl;

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        ++nonlinear_iteration;
      }
    while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);

    AssertThrow(nonlinear_solver_control.last_check() != SolverControl::failure, ExcNonlinearSolverNoConvergence());
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_iterated_advection_and_newton_stokes (bool use_newton_iterations)
  {
    // Now store the linear_tolerance we started out with, because we might change
    // it within this timestep.
    const double begin_linear_tolerance = parameters.linear_stokes_solver_tolerance;
    double initial_temperature_residual = 0;
    std::vector<double> initial_composition_residual (parameters.n_compositional_fields,0);

    DefectCorrectionResiduals dcr;
    dcr.initial_residual = 1;

    dcr.velocity_residual = 0;
    dcr.pressure_residual = 0;
    dcr.residual = 1;
    dcr.residual_old = 1;

    dcr.switch_initial_residual = 1;
    dcr.newton_residual_for_derivative_scaling_factor = 1;

    bool use_picard = true;

    const Newton::Parameters::Stabilization starting_preconditioner_stabilization = newton_handler->parameters.preconditioner_stabilization;
    const Newton::Parameters::Stabilization starting_velocity_block_stabilization = newton_handler->parameters.velocity_block_stabilization;

    const unsigned int max_nonlinear_iterations =
      (pre_refinement_step < parameters.initial_adaptive_refinement)
      ?
      std::min(parameters.max_nonlinear_iterations,
               parameters.max_nonlinear_iterations_in_prerefinement)
      :
      parameters.max_nonlinear_iterations;

    // Now iterate out the nonlinearities.
    dcr.stokes_residuals = std::pair<double,double>  (numbers::signaling_nan<double>(),
                                                      numbers::signaling_nan<double>());

    SolverControl nonlinear_solver_control(max_nonlinear_iterations,
                                           parameters.nonlinear_tolerance);

    SolverControl nonlinear_solver_control_picard(newton_handler->parameters.max_pre_newton_nonlinear_iterations,
                                                  newton_handler->parameters.nonlinear_switch_tolerance);

    double relative_residual = std::numeric_limits<double>::max();
    nonlinear_iteration = 0;

    do
      {
        // Restore particles through stored copy of particle handler,
        // but only if they have already been displaced in a nonlinear
        // iteration (in the assemble_and_solve_composition call).

        if (nonlinear_iteration > 0)
          for (auto &particle_manager : particle_managers)
            particle_manager.restore_particles();

        const double relative_temperature_residual =
          assemble_and_solve_temperature(initial_temperature_residual,
                                         nonlinear_iteration == 0 ? &initial_temperature_residual : nullptr);

        const std::vector<double>  relative_composition_residual =
          assemble_and_solve_composition(initial_composition_residual,
                                         nonlinear_iteration == 0 ? &initial_composition_residual : nullptr);

        if (use_newton_iterations && use_picard &&
            nonlinear_solver_control_picard.check(nonlinear_iteration, relative_residual) != SolverControl::iterate)
          {
            use_picard = false;
            pcout << "   Switching from defect correction form of Picard to the Newton solver scheme." << std::endl;

            /**
             * This method allows to slowly introduce the derivatives based
             * on the improvement of the residual. If we do not use it, we
             * just set it so the newton_derivative_scaling_factor goes from
             * zero to one when switching on the Newton solver.
             */
            if (!newton_handler->parameters.use_newton_residual_scaling_method)
              dcr.newton_residual_for_derivative_scaling_factor = 0;
          }

        newton_handler->parameters.newton_derivative_scaling_factor
          = (std::max(0.0,
                      (1.0-(dcr.newton_residual_for_derivative_scaling_factor/dcr.switch_initial_residual))));

        do_one_defect_correction_Stokes_step(dcr, use_picard);

        double max_residual = relative_temperature_residual;
        for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
          {
            // in models with melt migration the melt advection equation includes the divergence of the velocity
            // and can not be expected to converge to a smaller value than the residual of the Stokes equation.
            // thus, we set a threshold for the initial composition residual.
            // this only plays a role if the right-hand side of the advection equation is very small.
            const double threshold = (parameters.include_melt_transport && c == introspection.compositional_index_for_name("porosity")
                                      ?
                                      parameters.linear_stokes_solver_tolerance * time_step
                                      :
                                      0.0);
            if (initial_composition_residual[c]>threshold)
              max_residual = std::max(relative_composition_residual[c],max_residual);
          }

        // write the residual output in the same order as the solutions
        pcout << "      Relative nonlinear residuals (temperature"
              << (introspection.n_compositional_fields > 0 ? ", compositional fields" : "")
              << ", Stokes system): " << relative_temperature_residual;
        for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
          pcout << ", " << relative_composition_residual[c];
        pcout << ", " << dcr.residual/dcr.initial_residual << std::endl;

        relative_residual = std::max(dcr.residual/dcr.initial_residual, max_residual);
        pcout << "      Relative nonlinear residual (total system) after nonlinear iteration " << nonlinear_iteration+1
              << ": " << relative_residual
              << std::endl
              << std::endl;

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        ++nonlinear_iteration;
      }
    while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);

    // Reset the Newton stabilization at the end of the timestep.
    newton_handler->parameters.preconditioner_stabilization = starting_preconditioner_stabilization;
    newton_handler->parameters.velocity_block_stabilization = starting_velocity_block_stabilization;

    // Reset the linear tolerance to what it was at the beginning of the time step.
    parameters.linear_stokes_solver_tolerance = begin_linear_tolerance;

    AssertThrow(nonlinear_solver_control.last_check() != SolverControl::failure, ExcNonlinearSolverNoConvergence());

    signals.post_nonlinear_solver(nonlinear_solver_control);
  }


  template <int dim>
  void Simulator<dim>::solve_single_advection_and_iterated_newton_stokes (bool use_newton_iterations)
  {
    // First assemble and solve the temperature and compositional fields
    assemble_and_solve_temperature();
    assemble_and_solve_composition();

    // Now store the linear_tolerance we started out with, because we might change
    // it within this timestep.
    double begin_linear_tolerance = parameters.linear_stokes_solver_tolerance;

    std::vector<double> initial_composition_residual (parameters.n_compositional_fields,0);

    DefectCorrectionResiduals dcr;
    dcr.initial_residual = 1;

    dcr.velocity_residual = 0;
    dcr.pressure_residual = 0;
    dcr.residual = 1;
    dcr.residual_old = 1;

    dcr.switch_initial_residual = 1;
    dcr.newton_residual_for_derivative_scaling_factor = 1;

    bool use_picard = true;

    const Newton::Parameters::Stabilization starting_preconditioner_stabilization = newton_handler->parameters.preconditioner_stabilization;
    const Newton::Parameters::Stabilization starting_velocity_block_stabilization = newton_handler->parameters.velocity_block_stabilization;

    const unsigned int max_nonlinear_iterations =
      (pre_refinement_step < parameters.initial_adaptive_refinement)
      ?
      std::min(parameters.max_nonlinear_iterations,
               parameters.max_nonlinear_iterations_in_prerefinement)
      :
      parameters.max_nonlinear_iterations;

    SolverControl nonlinear_solver_control(max_nonlinear_iterations,
                                           parameters.nonlinear_tolerance);

    SolverControl nonlinear_solver_control_picard(newton_handler->parameters.max_pre_newton_nonlinear_iterations,
                                                  newton_handler->parameters.nonlinear_switch_tolerance);

    // Now iterate out the nonlinearities.
    dcr.stokes_residuals = {numbers::signaling_nan<double>(), numbers::signaling_nan<double>()};

    double relative_residual = std::numeric_limits<double>::max();
    nonlinear_iteration = 0;
    do
      {
        // If we are in the Picard phase, check if we can switch to Newton
        if (use_newton_iterations && use_picard &&
            nonlinear_solver_control_picard.check(nonlinear_iteration, relative_residual) != SolverControl::iterate)
          {
            use_picard = false;
            pcout << "   Switching from defect correction form of Picard to the Newton solver scheme." << std::endl;

            /**
             * This method allows to slowly introduce the derivatives based
             * on the improvement of the residual. If we do not use it, we
             * just set it so the newton_derivative_scaling_factor goes from
             * zero to one when switching on the Newton solver.
             */
            if (!newton_handler->parameters.use_newton_residual_scaling_method)
              dcr.newton_residual_for_derivative_scaling_factor = 0;
          }

        newton_handler->parameters.newton_derivative_scaling_factor
          = (std::max(0.0,
                      (1.0-(dcr.newton_residual_for_derivative_scaling_factor/dcr.switch_initial_residual))));

        do_one_defect_correction_Stokes_step(dcr, use_picard);

        relative_residual = dcr.residual/dcr.initial_residual;

        pcout << "      Relative nonlinear residual (Stokes system) after nonlinear iteration " << nonlinear_iteration+1
              << ": " << relative_residual
              << std::endl
              << std::endl;

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess ();

        ++nonlinear_iteration;
      }
    while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);

    // Reset the Newton stabilization at the end of the timestep.
    newton_handler->parameters.preconditioner_stabilization = starting_preconditioner_stabilization;
    newton_handler->parameters.velocity_block_stabilization = starting_velocity_block_stabilization;

    // Reset the linear tolerance to what it was at the beginning of the time step.
    parameters.linear_stokes_solver_tolerance = begin_linear_tolerance;

    AssertThrow(nonlinear_solver_control.last_check() != SolverControl::failure, ExcNonlinearSolverNoConvergence());

    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_single_advection_no_stokes ()
  {
    assemble_and_solve_temperature();
    assemble_and_solve_composition();

    {
      TimerOutput::Scope timer (computing_timer, "Interpolate Stokes solution");

      // Assign Stokes solution
      LinearAlgebra::BlockVector distributed_stokes_solution (introspection.index_sets.system_partitioning, mpi_communicator);

      VectorFunctionFromVectorFunctionObject<dim> func(
        [&](const Point<dim> &p, Vector<double> &result)
      {
        prescribed_stokes_solution->stokes_solution(p, result);
      },
      0,
      parameters.include_melt_transport ? 2*dim+3 : dim+1, // velocity and pressure
      introspection.n_components);

      VectorTools::interpolate (*mapping, dof_handler, func, distributed_stokes_solution);

      // distribute hanging node and other constraints
      current_constraints.distribute (distributed_stokes_solution);

      solution.block(introspection.block_indices.velocities) =
        distributed_stokes_solution.block(introspection.block_indices.velocities);
      solution.block(introspection.block_indices.pressure) =
        distributed_stokes_solution.block(introspection.block_indices.pressure);

      if (parameters.include_melt_transport)
        {
          const unsigned int block_u_f = introspection.variable("fluid velocity").block_index;
          const unsigned int block_p_f = introspection.variable("fluid pressure").block_index;
          solution.block(block_u_f) = distributed_stokes_solution.block(block_u_f);
          solution.block(block_p_f) = distributed_stokes_solution.block(block_p_f);
        }

    }

    if (parameters.run_postprocessors_on_nonlinear_iterations)
      postprocess ();

    // Setup a nonlinear solver control that only allows a single iteration
    SolverControl nonlinear_solver_control(1,1.0);
    // Announce that we did a single iteration, and assume we have converged
    nonlinear_solver_control.check(1,0.0);
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }



  template <int dim>
  void Simulator<dim>::solve_no_advection_no_stokes ()
  {
    if (parameters.run_postprocessors_on_nonlinear_iterations)
      postprocess ();

    // Setup a nonlinear solver control that only allows a single iteration
    SolverControl nonlinear_solver_control(1,1.0);
    // Announce that we did a single iteration, and assume we have converged
    nonlinear_solver_control.check(1,0.0);
    signals.post_nonlinear_solver(nonlinear_solver_control);
  }
}

// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template double Simulator<dim>::assemble_and_solve_temperature(const double &, double*); \
  template std::vector<double> Simulator<dim>::assemble_and_solve_composition(const std::vector<double> &, std::vector<double> *); \
  template double Simulator<dim>::assemble_and_solve_stokes(const double &, double*); \
  template void Simulator<dim>::solve_single_advection_single_stokes(); \
  template void Simulator<dim>::solve_no_advection_iterated_stokes(); \
  template void Simulator<dim>::solve_no_advection_single_stokes(); \
  template void Simulator<dim>::solve_iterated_advection_and_stokes(); \
  template void Simulator<dim>::solve_single_advection_iterated_stokes(); \
  template void Simulator<dim>::solve_no_advection_iterated_defect_correction_stokes(); \
  template void Simulator<dim>::solve_single_advection_iterated_defect_correction_stokes(); \
  template void Simulator<dim>::solve_iterated_advection_and_defect_correction_stokes(); \
  template void Simulator<dim>::solve_iterated_advection_and_newton_stokes(bool); \
  template void Simulator<dim>::solve_single_advection_and_iterated_newton_stokes(bool); \
  template void Simulator<dim>::solve_single_advection_no_stokes(); \
  template void Simulator<dim>::solve_first_timestep_only_single_stokes(); \
  template void Simulator<dim>::solve_no_advection_no_stokes();

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}

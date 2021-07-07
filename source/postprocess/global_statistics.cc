/*
  Copyright (C) 2017 - 2020 by the authors of the ASPECT code.

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


#include <aspect/postprocess/global_statistics.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    GlobalStatistics<dim>::initialize()
    {
      this->get_signals().post_stokes_solver.connect(
        [&](const SimulatorAccess<dim> &/*simulator_access*/,
            const unsigned int number_S_iterations,
            const unsigned int number_A_iterations,
            const SolverControl &solver_control_cheap,
            const SolverControl &solver_control_expensive)
      {
        this->store_stokes_solver_history(number_S_iterations,
                                          number_A_iterations,
                                          solver_control_cheap,
                                          solver_control_expensive);
      });

      this->get_signals().post_advection_solver.connect(
        [&](const SimulatorAccess<dim> &/*simulator_access*/,
            const bool solved_temperature_field,
            const unsigned int compositional_index,
            const SolverControl &solver_control)
      {
        this->store_advection_solver_history(solved_temperature_field,
                                             compositional_index,
                                             solver_control);
      });

      // delete the data after the initial refinement steps, to not mix it up
      // with the first time step
      if (!this->get_parameters().run_postprocessors_on_initial_refinement)
        this->get_signals().post_set_initial_state.connect(
          [&] (const SimulatorAccess<dim> &)
        {
          this->clear_data();
        });
    }



    template <int dim>
    void
    GlobalStatistics<dim>::clear_data()
    {
      list_of_S_iterations.clear();
      list_of_A_iterations.clear();
      stokes_iterations_cheap.clear();
      stokes_iterations_expensive.clear();
      advection_iterations.clear();
    }



    template <int dim>
    void
    GlobalStatistics<dim>::store_stokes_solver_history(const unsigned int number_S_iterations,
                                                       const unsigned int number_A_iterations,
                                                       const SolverControl &solver_control_cheap,
                                                       const SolverControl &solver_control_expensive)
    {
      list_of_S_iterations.push_back(number_S_iterations);
      list_of_A_iterations.push_back(number_A_iterations);
      stokes_iterations_cheap.push_back(solver_control_cheap.last_step());
      stokes_iterations_expensive.push_back(solver_control_expensive.last_step());
    }



    template <int dim>
    void
    GlobalStatistics<dim>::store_advection_solver_history(const bool solved_temperature_field,
                                                          const unsigned int compositional_index,
                                                          const SolverControl &solver_control)
    {
      const std::string column_name = (solved_temperature_field) ?
                                      "Iterations for temperature solver"
                                      :
                                      "Iterations for composition solver " + Utilities::int_to_string(compositional_index+1);

      unsigned int column_position = numbers::invalid_unsigned_int;

      for (unsigned int i=0; i<advection_iterations.size(); ++i)
        if (column_name == advection_iterations[i].first)
          column_position = i;

      if (column_position == numbers::invalid_unsigned_int)
        advection_iterations.emplace_back(column_name,std::vector<unsigned int>(1,solver_control.last_step()));
      else
        advection_iterations[column_position].second.push_back(solver_control.last_step());
    }



    template <int dim>
    std::pair<std::string,std::string>
    GlobalStatistics<dim>::execute (TableHandler &statistics)
    {
      // We would like to know how many nonlinear iterations were performed since
      // the last call to execute(). Unfortunately some solver schemes iterate only the
      // Stokes equation, some Stokes and advection, and some do not even call the
      // Stokes or advection solver. Therefore, the best estimate for the number
      // of nonlinear iterations since the last call is the maximum size of the
      // vectors in which we store our data.
      unsigned int nonlinear_iterations = stokes_iterations_cheap.size();

      for (unsigned int column=0; column<advection_iterations.size(); ++column)
        nonlinear_iterations = std::max(nonlinear_iterations, static_cast<unsigned int> (advection_iterations[column].second.size()));

      if (one_line_per_iteration)
        for (unsigned int iteration = 0; iteration < nonlinear_iterations; ++iteration)
          {
            generate_global_statistics(statistics);

            statistics.add_value("Nonlinear iteration number",
                                 iteration);

            for (unsigned int column=0; column<advection_iterations.size(); ++column)
              if (iteration < advection_iterations[column].second.size())
                statistics.add_value(advection_iterations[column].first,
                                     advection_iterations[column].second[iteration]);

            if (iteration < stokes_iterations_cheap.size())
              {
                statistics.add_value("Iterations for Stokes solver",
                                     stokes_iterations_cheap[iteration] + stokes_iterations_expensive[iteration]);
                statistics.add_value("Velocity iterations in Stokes preconditioner",
                                     list_of_A_iterations[iteration]);
                statistics.add_value("Schur complement iterations in Stokes preconditioner",
                                     list_of_S_iterations[iteration]);
              }

          }
      else
        {
          generate_global_statistics(statistics);

          unsigned int A_iterations = 0;
          unsigned int S_iterations = 0;
          unsigned int Stokes_outer_iterations = 0;
          std::vector<unsigned int> advection_outer_iterations(advection_iterations.size(),0);

          for (unsigned int iteration = 0; iteration < nonlinear_iterations; ++iteration)
            {
              for (unsigned int column=0; column<advection_iterations.size(); ++column)
                if (iteration < advection_iterations[column].second.size())
                  advection_outer_iterations[column] += advection_iterations[column].second[iteration];

              if (iteration < stokes_iterations_cheap.size())
                {
                  A_iterations += list_of_A_iterations[iteration];
                  S_iterations += list_of_S_iterations[iteration];
                  Stokes_outer_iterations += stokes_iterations_cheap[iteration] +
                                             stokes_iterations_expensive[iteration];
                }
            }

          // only output the number of nonlinear iterations if we actually
          // use a nonlinear solver scheme
          if (!(this->get_parameters().nonlinear_solver == Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                || this->get_parameters().nonlinear_solver == Parameters<dim>::NonlinearSolver::single_Advection_no_Stokes))
            statistics.add_value("Number of nonlinear iterations",
                                 nonlinear_iterations);

          // Only output statistics columns if the solver actually signaled at least one
          // successful solve. Some solver schemes might need no advection or Stokes solver
          for (unsigned int column=0; column<advection_iterations.size(); ++column)
            statistics.add_value(advection_iterations[column].first,
                                 advection_outer_iterations[column]);

          // Note that even if no cheap solver iterations were done, the solver control
          // object is still stored, so this line works in that case as well.
          if (stokes_iterations_cheap.size() > 0)
            {
              statistics.add_value("Iterations for Stokes solver",
                                   Stokes_outer_iterations);
              statistics.add_value("Velocity iterations in Stokes preconditioner",
                                   A_iterations);
              statistics.add_value("Schur complement iterations in Stokes preconditioner",
                                   S_iterations);
            }
        }

      clear_data();

      return std::make_pair (std::string(),std::string());
    }



    template <int dim>
    void
    GlobalStatistics<dim>::generate_global_statistics(TableHandler &statistics)
    {
      // set global statistics about this time step
      statistics.add_value("Time step number", this->get_timestep_number());

      if (this->convert_output_to_years() == true)
        {
          statistics.add_value("Time (years)", this->get_time() / year_in_seconds);
          statistics.set_precision("Time (years)", 12);
          statistics.set_scientific("Time (years)", true);

          statistics.add_value("Time step size (years)", this->get_timestep() / year_in_seconds);
          statistics.set_precision("Time step size (years)", 12);
          statistics.set_scientific("Time step size (years)", true);
        }
      else
        {
          statistics.add_value("Time (seconds)", this->get_time());
          statistics.set_precision("Time (seconds)", 12);
          statistics.set_scientific("Time (seconds)", true);

          statistics.add_value("Time step size (seconds)", this->get_timestep());
          statistics.set_precision("Time step size (seconds)", 12);
          statistics.set_scientific("Time step size (seconds)", true);
        }

      // set global statistics about the mesh and problem size
      statistics.add_value("Number of mesh cells",
                           this->get_triangulation().n_global_active_cells());

      types::global_dof_index n_stokes_dofs = this->introspection().system_dofs_per_block[0];
      if (this->introspection().block_indices.velocities != this->introspection().block_indices.pressure)
        n_stokes_dofs += this->introspection().system_dofs_per_block[this->introspection().block_indices.pressure];

      statistics.add_value("Number of Stokes degrees of freedom", n_stokes_dofs);
      statistics.add_value("Number of temperature degrees of freedom",
                           this->introspection().system_dofs_per_block[this->introspection().block_indices.temperature]);
      if (this->n_compositional_fields() > 0)
        statistics.add_value("Number of degrees of freedom for all compositions",
                             static_cast<types::global_dof_index>(this->n_compositional_fields())
                             * this->introspection().system_dofs_per_block[this->introspection().block_indices.compositional_fields[0]]);
    }



    template <int dim>
    void
    GlobalStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Global statistics");
        {
          prm.declare_entry ("Write statistics for each nonlinear iteration", "false",
                             Patterns::Bool (),
                             "Whether to put every nonlinear iteration into a separate "
                             "line in the statistics file (if true), or to output only "
                             "one line per time step that contains the total number of "
                             "iterations of the Stokes and advection linear system solver.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    GlobalStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Global statistics");
        {
          one_line_per_iteration = prm.get_bool("Write statistics for each nonlinear iteration");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(GlobalStatistics,
                                  "global statistics",
                                  "A postprocessor that outputs all the global statistics "
                                  "information, e.g. the time of the simulation, the timestep "
                                  "number, number of degrees of freedom and solver iterations "
                                  "for each timestep. The postprocessor can output different "
                                  "formats, the first printing one line in the statistics file "
                                  "per nonlinear solver iteration (if a nonlinear solver scheme "
                                  "is selected). The second prints one line per timestep, "
                                  "summing the information about all nonlinear iterations in "
                                  "this line. Note that this postprocessor is always active "
                                  "independent on whether or not it is selected in the "
                                  "parameter file.")
  }
}

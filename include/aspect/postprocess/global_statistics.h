/*
  Copyright (C) 2017 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_global_statistics_h
#define _aspect_postprocess_global_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/lac/solver_control.h>

namespace aspect
{

  namespace Postprocess
  {
    /**
     * A postprocessor that outputs all the global statistics
     * information, e.g. the time of the simulation, the timestep
     * number, number of degrees of freedom and solver iterations
     * for each timestep. The postprocessor can output different
     * formats, the first printing one line in the statistics file
     * per nonlinear solver iteration (if a nonlinear solver scheme
     * is selected). The second prints one line per timestep,
     * summing the information about all nonlinear iterations in
     * this line. Note that this postprocessor is always active
     * independent on whether or not it is selected in the
     * parameter file.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class GlobalStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Connect the callback functions to the respective signals.
         */
        void initialize() override;

        /**
         * Write all global statistics columns into the statistics object.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;


      private:
        /**
         * Write global statistics such as the time step number and the
         * number of degrees of freedom into the statistics object.
         */
        void
        generate_global_statistics(TableHandler &statistics);

        /**
         * This function clears all the saved solver information that is
         * stored in this class. It is executed every time the output is
         * written into the statistics object, but also after each initial
         * refinement step, in case it is not written (to avoid summing
         * information across different refinement steps).
         */
        void
        clear_data();

        /**
         * Callback function that is connected to the post_stokes_solver
         * signal to store the solver history.
         */
        void
        store_stokes_solver_history(const unsigned int number_S_iterations,
                                    const unsigned int number_A_iterations,
                                    const SolverControl &solver_control_cheap,
                                    const SolverControl &solver_control_expensive);

        /**
         * Callback function that is connected to the post_advection_solver
         * signal to store the solver history.
         */
        void
        store_advection_solver_history(const bool solved_temperature_field,
                                       const unsigned int compositional_index,
                                       const SolverControl &solver_control);

        /**
         * Variables that store the Stokes solver history of the current
         * timestep, until they are written into the statistics object
         * upon the call to execute(). They are cleared after writing
         * the content.
         */
        std::vector<unsigned int> list_of_S_iterations;
        std::vector<unsigned int> list_of_A_iterations;
        std::vector<unsigned int> stokes_iterations_cheap;
        std::vector<unsigned int> stokes_iterations_expensive;

        /**
         * A container that stores the advection solver history of the current
         * timestep, until it is written into the statistics object
         * upon the call to execute(). It is cleared after writing
         * the content. The vector contains pairs, which consist
         * of a column name (for the temperature or one of the
         * compositional fields), and a vector of SolverControl
         * objects (one per nonlinear iteration for this particular
         * field). This layout allows storing varying numbers of
         * nonlinear iterations for temperature and compositional fields
         * (if any nonlinear solver scheme would implement that at
         * some point).
         */
        std::vector<std::pair<std::string, std::vector<unsigned int>>> advection_iterations;

        /**
         * Whether to put every nonlinear iteration into a separate
         * line in the statistics file or to only output one line
         * per time step.
         */
        bool one_line_per_iteration;
    };
  }
}


#endif

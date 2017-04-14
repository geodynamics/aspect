/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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
     *
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class GlobalStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void initialize();

        /**
         * Evaluate the solution for some general statistics.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


      private:
        /**
         * Write global statistics such as the time step number and the
         * number of degrees of freedom into the statistics object.
         */
        void
        generate_global_statistics(TableHandler &statistics);

        void
        store_stokes_solver_history(const unsigned int number_S_iterations,
                                    const unsigned int number_A_iterations,
                                    const SolverControl &solver_control_cheap,
                                    const SolverControl &solver_control_expensive);

        void
        store_advection_solver_history(const bool solved_temperature_field,
                                       const unsigned int compositional_index,
                                       const SolverControl &solver_control);

        std::vector<unsigned int> list_of_S_iterations;
        std::vector<unsigned int> list_of_A_iterations;
        std::vector<SolverControl> solver_controls_cheap;
        std::vector<SolverControl> solver_controls_expensive;

        std::vector<std::pair<std::string, std::vector<SolverControl> > > advection_solver_controls;

        /**
         * Whether to put every nonlinear iteration into a separate
         * line in the statistics file or to only output one line
         * per time step.
         */
        double one_line_per_iteration;
    };
  }
}


#endif

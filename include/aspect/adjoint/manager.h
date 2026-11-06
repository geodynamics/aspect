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

#ifndef _aspect_adjoint_manager_h
#define _aspect_adjoint_manager_h

#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/adjoint/objective_functional.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parameter_handler.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace aspect
{
  template <int dim> class Simulator;

  namespace Adjoint
  {
    /**
     * A class that manages instantaneous adjoint workflows. It collects the
     * active objective functionals, assembles and solves the corresponding
     * adjoint systems, and calculates physical-property kernels.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Manager();

        /**
         * Declare the parameters used by the adjoint manager and its
         * subobjects.
         *
         * @param prm The parameter handler in which parameters are declared.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Initialize this manager and all objective plugins with access to the
         * simulator object. Unlike most SimulatorAccess users, the adjoint
         * manager needs mutable access to the simulator: the instantaneous
         * adjoint solve temporarily replaces the Stokes right hand side and
         * solution vectors, calls simulator solve routines, and then restores
         * the forward state.
         *
         * This function is therefore the only valid initialization path for
         * this class. It stores the mutable simulator pointer used by the
         * adjoint workflow, initializes the SimulatorAccess base class for
         * standard read-only simulator queries, and forwards simulator access
         * to the contained objective manager.
         *
         * @param simulator The simulator object that owns the forward model
         * and adjoint workflow.
         */
        void
        initialize_simulator(Simulator<dim> &simulator);

        /**
         * Override the generic SimulatorAccess initialization interface, but
         * reject it deliberately. The base class interface only provides a
         * const simulator reference, which would leave this manager without
         * the mutable simulator access required by solve_instantaneous_stokes().
         *
         * Defining this override prevents accidental initialization through a
         * SimulatorAccess pointer from silently setting only the base-class
         * simulator pointer. Call initialize_simulator(Simulator<dim> &) on
         * the concrete adjoint manager instead.
         */
        void
        initialize_simulator(const Simulator<dim> &simulator) override;

        /**
         * Read adjoint parameters from the parameter file and create the
         * selected objective plugins.
         *
         * @param prm The parameter handler from which adjoint parameters are
         * read.
         */
        void
        parse_parameters(ParameterHandler &prm);

        /**
         * Initialize the adjoint manager after parameters have been parsed.
         */
        void
        initialize();

        /**
         * Run the instantaneous Stokes adjoint workflow for the current model.
         */
        void
        solve_instantaneous_stokes();

        /**
         * Return the most recently computed physical-property kernels.
         *
         * @return The kernel repository.
         */
        const KernelRepository<dim> &
        get_kernels() const;

        /**
         * Return the most recently evaluated objective values.
         *
         * @return A map from objective functional names to objective values.
         */
        std::map<std::string, double>
        get_objective_values() const;

      private:
        void
        validate_instantaneous_stokes_setup() const;

        void
        prepare_required_forward_postprocessors();

        void
        solve_adjoint_states();

        void
        calculate_kernels();

        /**
         * Mutable simulator pointer used only by the adjoint workflow for
         * operations that cannot be expressed through SimulatorAccess, such as
         * temporarily replacing the system right hand side or current solution
         * before an adjoint Stokes solve. The SimulatorAccess base class is
         * still initialized separately and remains the preferred route for
         * ordinary read-only simulator queries.
         */
        Simulator<dim> *simulator;
        ObjectiveManager<dim> objective_manager;
        std::vector<std::unique_ptr<ObjectiveResult<dim>>> objective_results;
        std::vector<std::unique_ptr<AdjointState<dim>>> adjoint_states;
        KernelRepository<dim> kernels;
    };
  }
}

#endif

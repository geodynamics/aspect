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

#ifndef _aspect_adjoint_kernel_calculator_h
#define _aspect_adjoint_kernel_calculator_h

#include <aspect/adjoint/state.h>
#include <aspect/adjoint/types.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/lac/vector.h>

#include <map>
#include <memory>
#include <vector>

namespace aspect
{
  template <int dim> class Simulator;

  namespace Adjoint
  {
    /**
     * A data structure that stores cellwise adjoint kernel contributions.
     * Contributions are kept separated by objective, physical property, and
     * physics term so that output plugins can inspect or combine them
     * explicitly.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    class KernelRepository
    {
      public:
        /**
         * Add one kernel contribution to the repository.
         *
         * @param key The key that identifies objective, physics term, and
         * physical property of the contribution.
         * @param values The cellwise values of the contribution.
         */
        void
        add_contribution(const KernelContributionKey &key,
                         const Vector<double> &values);

        /**
         * Store the cell volumes associated with the cellwise kernel values.
         *
         * @param values The volume of every active cell, indexed by active
         * cell index.
         */
        void
        set_cell_volumes(const Vector<double> &values);

        /**
         * Return the cell volumes associated with the stored kernel values.
         *
         * @return The cell volume vector.
         */
        const Vector<double> &
        cell_volumes() const;

        /**
         * Return whether the repository contains no kernel contributions.
         *
         * @return true if no kernel contributions are stored.
         */
        bool
        empty() const;

        /**
         * Return the number of kernel contributions stored in the repository.
         *
         * @return The number of kernel contributions.
         */
        unsigned int
        n_contributions() const;

        /**
         * Return all stored kernel contributions.
         *
         * @return A map from contribution keys to cellwise contribution
         * values.
         */
        const std::map<KernelContributionKey, Vector<double>> &
        contributions() const;

      private:
        std::map<KernelContributionKey, Vector<double>> contribution_map;
        Vector<double> cell_volume_values;
    };



    /**
     * A class that computes physical-property adjoint kernels from the
     * current forward solution and one adjoint state per objective.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    class KernelCalculator : public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         *
         * @param simulator The simulator object that provides access to the
         * forward model state and finite element data.
         */
        explicit KernelCalculator(const Simulator<dim> &simulator);

        /**
         * Calculate physical-property kernels from the current forward solution and the
         * adjoint states associated with the objective results.
         *
         * @param objective_results The objective values and adjoint right hand
         * sides associated with each objective.
         * @param adjoint_states The adjoint solutions associated with each
         * objective in @p objective_results.
         *
         * @return A repository containing the computed physical-property
         * kernels.
         */
        KernelRepository<dim>
        calculate(const std::vector<std::unique_ptr<ObjectiveResult<dim>>> &objective_results,
                  const std::vector<std::unique_ptr<AdjointState<dim>>> &adjoint_states) const;
    };
  }
}

#endif

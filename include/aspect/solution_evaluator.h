/*
 Copyright (C) 2024 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_solution_evaluator_h
#define _aspect_solution_evaluator_h

#include <aspect/global.h>


#include <aspect/simulator_access.h>

#include <deal.II/base/array_view.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/non_matching/mapping_info.h>
#include <deal.II/dofs/dof_handler.h>

namespace aspect
{
  using namespace dealii;

  /**
   * This class evaluates the solution vector at arbitrary positions inside a cell.
   * This base class only provides the interface for SolutionEvaluatorImplementation.
   * See there for more details.
   */
  template <int dim>
  class SolutionEvaluator
  {
    public:
      /**
       * virtual Destructor.
       */
      virtual ~SolutionEvaluator() = default;

      /**
       * Reinitialize all variables to evaluate the given solution for the given cell
       * and the given positions. The update flags control if only the solution or
       * also the gradients should be evaluated.
       * If other flags are set an assertion is triggered.
       */
      virtual void
      reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const ArrayView<Point<dim>> &positions,
             const ArrayView<double> &solution_values,
             const UpdateFlags update_flags) = 0;

      /**
       * Fill @p solution with all solution components at the given @p evaluation_point. Note
       * that this function only works after a successful call to reinit(),
       * because this function only returns the results of the computation that
       * happened in reinit().
       */
      virtual void get_solution(const unsigned int evaluation_point,
                                const ArrayView<double> &solution) = 0;

      /**
       * Fill @p gradients with all solution gradients at the given @p evaluation_point. Note
       * that this function only works after a successful call to reinit(),
       * because this function only returns the results of the computation that
       * happened in reinit().
       */
      virtual void get_gradients(const unsigned int evaluation_point,
                                 const ArrayView<Tensor<1, dim>> &gradients) = 0;

      /**
       * Return the evaluator for velocity or fluid velocity. This is the only
       * information necessary for advecting particles.
       */
      virtual FEPointEvaluation<dim, dim> &
      get_velocity_or_fluid_velocity_evaluator(const bool use_fluid_velocity) = 0;

      /**
       * Return the cached mapping information.
       */
      virtual NonMatching::MappingInfo<dim> &
      get_mapping_info() = 0;
  };

  /**
   * A function to create a pointer to a SolutionEvaluator object.
   */
  template <int dim>
  std::unique_ptr<SolutionEvaluator<dim>>
  construct_solution_evaluator(const SimulatorAccess<dim> &simulator_access,
                               const UpdateFlags update_flags);
}

#endif

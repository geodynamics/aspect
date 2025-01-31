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
  namespace internal
  {
    /**
     * Wrapper around dealii::FEPointEvaluation to choose number of components dynamically.
     * This class is only an abstract base class that can be used to implement derived classes.
     */
    template <int dim>
    class DynamicFEPointEvaluation
    {
      public:
        /**
         * Constructor which allows to select at run time the @p first_component and the @p n_components.
         */
        DynamicFEPointEvaluation(const unsigned int first_component, const unsigned int n_components)
          : first_component (first_component),
            n_components (n_components)
        {}

        /**
         * Destructor.
         */
        virtual ~DynamicFEPointEvaluation() = default;

        /**
         * Evaluate the solution at the given positions.
         *
         * @p solution_values contains the values of the degrees of freedom.
         * @p flags controls which values should be computed.
         */
        virtual void evaluate(const ArrayView<double> &solution_values,
                              const EvaluationFlags::EvaluationFlags flags) = 0;

        /**
         * Return the gradient of the solution at the given @p evaluation_point.
         */
        virtual
        small_vector<Tensor<1,dim>>
        get_gradient(const unsigned int evaluation_point) const = 0;

        /**
         * Copy the gradient of the solution at the given @p evaluation_point
         * into the array @p gradients.
         */
        virtual
        void
        get_gradient(const unsigned int evaluation_point,
                     const ArrayView<Tensor<1,dim>> &gradients) const = 0;

        /**
         * Return the value of the solution at the given @p evaluation_point.
         */
        virtual
        small_vector<double>
        get_value(const unsigned int evaluation_point) const = 0;

        /**
         * Copy the value of the solution at the given @p evaluation_point.
         * into the array @p solution.
         */
        virtual
        void
        get_value(const unsigned int evaluation_point,
                  const ArrayView<double> &solution) const = 0;

        /**
         * Return the first component of the solution vector.
         */
        virtual
        unsigned int
        get_first_component() const final
        {
          return first_component;
        }

        /**
         * Return the number of components of the solution vector.
         */
        virtual
        unsigned int
        get_n_components() const final
        {
          return n_components;
        }

      private:
        /**
         * The first component of the solution vector.
         */
        unsigned int first_component;

        /**
         * The number of components of the solution vector.
         */
        unsigned int n_components;
    };
  }

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
       * Constructor. Create the member variables given a simulator access @p simulator and a set of
       * update flags @p update_flags. The @p update_flags control if only the solution or also the gradients
       * should be evaluated.
       */
      SolutionEvaluator(const SimulatorAccess<dim> &simulator,
                        const UpdateFlags update_flags);

      /**
       * Reinitialize all variables to prepare for evaluation for the given @p cell
       * and at the given @p positions in the reference coordinate system of that cell.
       */
      void
      reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const ArrayView<Point<dim>> &positions);

      /**
       * Evaluate all variables in the cell and at the positions controlled by a previous
       * call to reinit().
       *
       * @p solution_values contains the values of the degrees of freedom.
       * @p evaluation_flags controls if nothing, the solution values, and/or the gradients should be computed.
       * The size of @p solution_flags has to be equal to the number of components as returned by
       * n_components().
       */
      void
      evaluate(const ArrayView<double> &solution_values,
               const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags);

      /**
       * Fill @p solution with all solution components at the given @p evaluation_point. Note
       * that this function only works after a successful call to reinit() and evaluate()
       * because this function only returns the results of the computation that
       * happened in those functions.
       *
       * @param evaluation_point The index of the evaluation point in the positions array.
       * @param solution The array to fill with the solution values. This array has to be
       *                 of size n_components().
       * @param evaluation_flags The flags that indicate which values should be copied into the
       *                        solution array. This vector has to be of size n_components().
       */
      void get_solution(const unsigned int evaluation_point,
                        const ArrayView<double> &solution,
                        const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags) const;

      /**
       * Fill @p gradients with all solution gradients at the given @p evaluation_point. Note
       * that this function only works after a successful call to reinit() and evaluate()
       * because this function only returns the results of the computation that
       * happened in those functions.
       *
       * @param evaluation_point The index of the evaluation point in the positions array.
       * @param gradients The array to fill with the solution gradients. This array has to be
       *                  of size n_components().
       * @param evaluation_flags The flags that indicate which gradients should be copied into the
       *                        solution array. This vector has to be of size n_components().
       */
      void get_gradients(const unsigned int evaluation_point,
                         const ArrayView<Tensor<1, dim>> &gradients,
                         const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags) const;

      /**
       * Return the evaluator for velocity or fluid velocity. This is the only
       * information necessary for advecting particles.
       */
      FEPointEvaluation<dim, dim> &
      get_velocity_or_fluid_velocity_evaluator(const bool use_fluid_velocity);

      /**
       * Return the cached mapping information.
       */
      NonMatching::MappingInfo<dim> &
      get_mapping_info();

      /**
       * Return the number of components in the solution vector.
       */
      unsigned int
      n_components() const;

    private:
      /**
       * MappingInfo object for the FEPointEvaluation objects
       */
      NonMatching::MappingInfo<dim> mapping_info;

      /**
       * FEPointEvaluation objects for all common
       * components of ASPECT's finite element solution.
       * These objects are used inside of the member functions of this class.
       */
      FEPointEvaluation<dim, dim> velocity;
      std::unique_ptr<FEPointEvaluation<1, dim>> pressure;
      FEPointEvaluation<1, dim> temperature;

      /**
       * We group compositions by type (FiniteElement) and evaluate
       * them together.
       */
      std::vector<std::unique_ptr<internal::DynamicFEPointEvaluation<dim>>> compositions;

      /**
       * Pointers to FEPointEvaluation objects for all melt
       * components of ASPECT's finite element solution, which only
       * point to valid objects in case we use melt transport. Other
       * documentation like for the objects directly above.
       */
      std::unique_ptr<FEPointEvaluation<dim, dim>> fluid_velocity;
      std::unique_ptr<FEPointEvaluation<1, dim>> compaction_pressure;
      std::unique_ptr<FEPointEvaluation<1, dim>> fluid_pressure;

      /**
       * The component indices for the three melt formulation
       * variables fluid velocity, compaction pressure, and
       * fluid pressure (in this order). They are cached
       * to avoid repeated expensive lookups.
       */
      std::array<unsigned int, 3> melt_component_indices;

      /**
       * Reference to the active simulator access object. Provides
       * access to the general simulation variables.
       */
      const SimulatorAccess<dim> &simulator_access;
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

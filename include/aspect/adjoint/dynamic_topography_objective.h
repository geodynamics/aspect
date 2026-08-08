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

#ifndef _aspect_adjoint_dynamic_topography_objective_h
#define _aspect_adjoint_dynamic_topography_objective_h

#include <aspect/adjoint/objective_functional.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * A class that implements the instantaneous dynamic topography objective
     * functional. It evaluates a surface dynamic topography misfit, assembles
     * the corresponding top-boundary adjoint right hand side, and contributes
     * direct density and viscosity kernel terms associated with the objective.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    class DynamicTopographyObjective : public ObjectiveFunctional<dim>
    {
      public:
        /**
         * Constructor.
         */
        DynamicTopographyObjective ();

        /**
         * Declare the parameters used by this objective functional.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters used by this objective functional.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Evaluate one half of the squared surface dynamic topography misfit.
         *
         * @return The scalar value of the objective functional.
         */
        double
        evaluate() const override;

        /**
         * Assemble the top-boundary adjoint right hand side associated with
         * the dynamic topography objective.
         *
         * @param[in,out] rhs The adjoint right hand side vector to which this
         * objective adds its contribution.
         */
        void
        assemble_adjoint_rhs(LinearAlgebra::BlockVector &rhs) const override;

        /**
         * Add direct density and viscosity kernel terms associated with the
         * dynamic topography objective.
         *
         * @param objective_name The objective functional name used to identify
         * the contributions in @p kernels.
         * @param[in,out] kernels The repository to which this objective adds
         * direct kernel contributions.
         */
        void
        add_direct_kernel_contributions(const std::string &objective_name,
                                        KernelRepository<dim> &kernels) const override;

        /**
         * Return the forward postprocessors required by this objective.
         *
         * @return A set containing the dynamic topography postprocessor name.
         */
        std::set<std::string>
        required_forward_postprocessors() const override;

      private:
        /**
         * Return the prescribed target dynamic topography at @p position.
         */
        double
        target_topography (const Point<dim> &position) const;

        /**
         * The currently selected target topography model.
         */
        std::string target_topography_model;

        /**
         * A function object representing the target topography.
         */
        Functions::ParsedFunction<dim> target_topography_function;

        /**
         * The coordinate representation in which the target function is
         * evaluated.
         */
        Utilities::Coordinates::CoordinateSystem target_coordinate_system;
    };
  }
}

#endif

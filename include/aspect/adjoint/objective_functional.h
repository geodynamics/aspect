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

#ifndef _aspect_adjoint_objective_functional_h
#define _aspect_adjoint_objective_functional_h

#include <aspect/adjoint/state.h>
#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <memory>
#include <set>
#include <string>

namespace aspect
{
  /**
   * A namespace for everything to do with adjoint methods.
   *
   * @ingroup Adjoint
   */
  namespace Adjoint
  {
    template <int dim> class KernelRepository;

    /**
     * A base class for adjoint objective functionals. Implementations
     * evaluate an objective, assemble its contribution to the adjoint right
     * hand side, and may add direct physical-property kernel contributions.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    class ObjectiveFunctional : public Plugins::InterfaceBase,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor.
         */
        virtual
        ~ObjectiveFunctional() override = default;

        /**
         * Evaluate the objective functional for the current forward solution.
         *
         * @return The scalar value of the objective functional.
         */
        virtual
        double
        evaluate() const = 0;

        /**
         * Assemble the derivative of the objective functional with respect to
         * the Stokes solution into the adjoint right hand side for the current
         * forward solution.
         *
         * @param[in,out] rhs The adjoint right hand side vector. Implementations
         * add their contribution to this vector. Implementations may leave this
         * object unchanged if the objective has no dependence on the
         * instantaneous Stokes state.
         */
        virtual
        void
        assemble_adjoint_rhs(LinearAlgebra::BlockVector &rhs) const = 0;

        /**
         * Add direct physical-property kernel contributions that do not pass
         * through the adjoint solution.
         *
         * @param objective_name The objective functional name used to identify
         * the contributions in @p kernels.
         * @param[in,out] kernels The repository to which implementations add
         * direct kernel contributions. Implementations may leave this object
         * unchanged if the objective has no direct physical-property
         * dependence.
         */
        virtual
        void
        add_direct_kernel_contributions(const std::string &objective_name,
                                        KernelRepository<dim> &kernels) const = 0;

        /**
         * Return the names of forward postprocessors that must be executed
         * before this objective functional is evaluated or assembled. The
         * names are the same strings used in parameter files. The default
         * implementation returns an empty set.
         *
         * @return A set of required forward postprocessor names.
         */
        virtual
        std::set<std::string>
        required_forward_postprocessors() const;
    };



    /**
     * A class that manages all active adjoint objective functionals. It reads
     * the list of selected objectives from the parameter file and creates the
     * corresponding objective plugins. It also provides functions that call
     * all active objective plugins to collect required postprocessors,
     * assemble objective right hand sides, and add direct kernel
     * contributions.
     *
     * @ingroup Adjoint
     */
    template <int dim>
    class ObjectiveManager : public Plugins::ManagerBase<ObjectiveFunctional<dim>>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Register an adjoint objective functional so that it can be selected
         * from the parameter file. This allows objective implementations to
         * register themselves without modifying the manager whenever a new
         * objective functional is added.
         *
         * @param name A string that identifies the objective functional.
         * @param description A text description of what this objective
         * functional does and that will be listed in the documentation of the
         * parameter file.
         * @param declare_parameters_function A pointer to a function that can
         * be used to declare the parameters that this objective functional
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this objective functional.
         */
        static
        void
        register_objective_functional(const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      std::unique_ptr<ObjectiveFunctional<dim>> (*factory_function) ());

        /**
         * Declare the parameters used by the objective manager and all
         * registered objective functional plugins.
         *
         * @param prm The parameter handler in which parameters are declared.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the list of objective functionals from the parameter file and
         * create the corresponding objective plugin objects.
         *
         * @param prm The parameter handler from which the objective list is
         * read.
         */
        void
        parse_parameters(ParameterHandler &prm) override;

        /**
         * Return the names of all forward postprocessors required by active
         * objective functionals.
         *
         * @return A set of required forward postprocessor names.
         */
        std::set<std::string>
        required_forward_postprocessors() const;

        /**
         * Evaluate all active objective functionals for the current forward
         * solution and assemble their adjoint right hand side contributions.
         *
         * @return Objective values and right hand sides for all active
         * objective functionals.
         */
        std::vector<std::unique_ptr<ObjectiveResult<dim>>>
        evaluate_and_assemble_right_hand_sides() const;

        /**
         * Add direct physical-property kernel contributions from all active
         * objective functionals for the current forward solution.
         *
         * @param[in,out] kernels The repository to which active objectives add
         * direct kernel contributions.
         */
        void
        add_direct_kernel_contributions(KernelRepository<dim> &kernels) const;
    };



#define ASPECT_REGISTER_ADJOINT_OBJECTIVE(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_ADJOINT_OBJECTIVE_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Adjoint::ObjectiveFunctional<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::Adjoint::ObjectiveManager<2>::register_objective_functional, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Adjoint::ObjectiveFunctional<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Adjoint::ObjectiveManager<3>::register_objective_functional, \
                                name, description); \
  }
  }
}

#endif

// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2013-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_composition_initial_composition_h
#define _aspect_boundary_composition_initial_composition_h

#include <aspect/initial_composition/interface.h>
#include <aspect/boundary_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a
     * spherical shell geometry in which the composition at the inner and
     * outer surfaces (i.e. at the core-mantle and the
     * mantle-lithosphere/atmosphere boundaries) are constant.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class InitialComposition : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run.
         *
         * This specific function makes sure that the objects that describe
         * initial conditions remain available throughout the run of the
         * program.
         */
        void
        initialize () override;

        /**
         * This function returns the boundary compositions that are defined
         * by the initial conditions.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;

        /**
         * Return the minimal composition on that part of the boundary on
         * which Dirichlet conditions are posed.
         */
        virtual
        double minimal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Return the maximal composition on that part of the boundary on
         * which Dirichlet conditions are posed.
         */
        virtual
        double maximal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary compositions.
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
         * Compositions at the inner and outer boundaries.
         *
         * This variable is read from the parameter file through a parameter called 'Minimal composition'.
         */
        double min_composition;

        /**
         * Compositions at the inner and outer boundaries.
         *
         * This variable is read from the parameter file through a parameter called 'Maximal composition'.
         */
        double max_composition;

        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition;
    };
  }
}


#endif

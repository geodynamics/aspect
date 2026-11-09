// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_initial_composition_world_builder_h
#define _aspect_initial_composition_world_builder_h

#include <aspect/global.h>

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>

namespace WorldBuilder
{
  class World;
}


namespace aspect
{

#ifdef ASPECT_WITH_WORLD_BUILDER

  namespace InitialComposition
  {
    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file through the
     * World builder.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class WorldBuilder : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

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
         * A vector that specifies for each compositional field index if the world builder
         * should be evaluated for this compositional field.
         */
        std::vector<bool> relevant_compositions;

        /**
         * A pointer to the WorldBuilder object. Keeping this pointer ensures
         * that the object doesn't go away while we still need it.
         */
        std::shared_ptr<const ::WorldBuilder::World> world_builder;
    };
  }
#endif

}

#endif

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

#ifndef _aspect_initial_temperature_world_builder_h
#define _aspect_initial_temperature_world_builder_h

#include <aspect/global.h>

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace WorldBuilder
{
  class World;
}



namespace aspect
{

#ifdef ASPECT_WITH_WORLD_BUILDER

  namespace InitialTemperature
  {
    /**
     * A class that implements temperature initial conditions based on a
     * functional description provided in the input file through the
     * World builder.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class WorldBuilder : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        WorldBuilder ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

      private:
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

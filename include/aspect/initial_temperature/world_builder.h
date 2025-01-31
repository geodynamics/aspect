/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_initial_temperature_world_builder_h
#define _aspect_initial_temperature_world_builder_h

#include <aspect/global.h>

#ifdef ASPECT_WITH_WORLD_BUILDER

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace WorldBuilder
{
  class World;
}



namespace aspect
{
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
}

#endif
#endif

// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_geometry_model_initial_topography_model_world_builder_h
#define _aspect_geometry_model_initial_topography_model_world_builder_h

#include <aspect/global.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#  include <world_builder/config.h>
#endif

namespace WorldBuilder
{
  class World;
}


namespace aspect
{
#ifdef ASPECT_WITH_WORLD_BUILDER
#  if WORLD_BUILDER_VERSION_GTE(1,1,1)
  namespace InitialTopographyModel
  {
    /**
     * An initial topography model that queries the Geodynamic World Builder.
     *
     * World Builder property 6 represents topography. The World Builder file
     * is shared with the other World Builder initial-condition plugins and is
     * selected through the top-level "World builder file" parameter.
     */
    template <int dim>
    class WorldBuilder : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the shared World Builder object.
         */
        void
        initialize () override;

        /**
         * Return the World Builder topography at the given surface point.
         */
        double
        value (const Point<dim-1> &surface_point) const override;

        /**
         * Return World Builder's upper bound for the topography.
         */
        double
        max_topography () const override;

      private:
        std::shared_ptr<const ::WorldBuilder::World> world_builder;
        double maximum_topography;
    };
  }
#  endif
#endif
}

#endif

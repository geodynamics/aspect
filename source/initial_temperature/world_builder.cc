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

#include <aspect/global.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#  include <world_builder/config.h>
#  include <world_builder/world.h>

#  include <aspect/initial_temperature/world_builder.h>
#  include <aspect/geometry_model/interface.h>
#  include <aspect/gravity_model/interface.h>
#endif


namespace aspect
{
#ifdef ASPECT_WITH_WORLD_BUILDER
  namespace InitialTemperature
  {
    template <int dim>
    WorldBuilder<dim>::WorldBuilder ()
      = default;

    template <int dim>
    void
    WorldBuilder<dim>::
    initialize()
    {
      CitationInfo::add("GWB");
      world_builder = this->get_world_builder_pointer();
    }


    template <int dim>
    double
    WorldBuilder<dim>::
    initial_temperature (const Point<dim> &position) const
    {
#if WORLD_BUILDER_VERSION_MAJOR > 0 || WORLD_BUILDER_VERSION_MINOR >= 5
      return world_builder->temperature(Utilities::convert_point_to_array(position),
                                        -this->get_geometry_model().height_above_reference_surface(position));
#else

      return world_builder->temperature(Utilities::convert_point_to_array(position),
                                        -this->get_geometry_model().height_above_reference_surface(position),
                                        this->get_gravity_model().gravity_vector(position).norm());
#endif
    }

  }
#endif
}

// explicit instantiations
namespace aspect
{
#ifdef ASPECT_WITH_WORLD_BUILDER
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(WorldBuilder,
                                              "world builder",
                                              "Specify the initial temperature through the World Builder. "
                                              "More information on the World Builder can be found at "
                                              "\\url{https://geodynamicworldbuilder.github.io}. "
                                              "Make sure to specify the location of the World Builder file "
                                              "in the parameter 'World builder file'.")
  }
#endif
}

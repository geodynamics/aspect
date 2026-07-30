/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#include <aspect/global.h>

#ifdef ASPECT_WITH_WORLD_BUILDER
#  include <world_builder/config.h>
#  if WORLD_BUILDER_VERSION_GTE(1,1,1)
#    include <aspect/geometry_model/initial_topography_model/world_builder.h>
#    include <aspect/geometry_model/interface.h>
#    include <aspect/geometry_model/chunk.h>
#    include <aspect/geometry_model/spherical_shell.h>
#    include <aspect/geometry_model/two_merged_chunks.h>
#    include <aspect/utilities.h>

#    include <world_builder/world.h>
#  endif
#endif


namespace aspect
{
#ifdef ASPECT_WITH_WORLD_BUILDER
#  if WORLD_BUILDER_VERSION_GTE(1,1,1)
  namespace InitialTopographyModel
  {
    namespace
    {
      /**
       * Return the radius of the undeformed reference surface for supported
       * spherical geometry models.
       */
      template <int dim>
      double
      reference_surface_radius (const GeometryModel::Interface<dim> &geometry_model)
      {
        if (const auto *geometry =
              dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&geometry_model))
          return geometry->outer_radius();

        if (const auto *geometry =
              dynamic_cast<const GeometryModel::Chunk<dim> *>(&geometry_model))
          return geometry->outer_radius();

        if (const auto *geometry =
              dynamic_cast<const GeometryModel::TwoMergedChunks<dim> *>(&geometry_model))
          return geometry->outer_radius();

        AssertThrow(false,
                    ExcMessage("The World Builder initial topography model does "
                               "not know how to determine the reference-surface "
                               "radius for the selected spherical geometry model."));
        return 0.0;
      }

    }



    template <int dim>
    void
    WorldBuilder<dim>::initialize ()
    {
      CitationInfo::add("GWB");
      world_builder = this->get_world_builder_pointer();

      maximum_topography = world_builder->maximum_topography();
    }



    template <int dim>
    double
    WorldBuilder<dim>::value (const Point<dim-1> &surface_point) const
    {
      Assert(world_builder != nullptr, ExcInternalError());

      const GeometryModel::Interface<dim> &geometry_model =
        this->get_geometry_model();

      Point<dim> position;
      switch (geometry_model.natural_coordinate_system())
        {
          case Utilities::Coordinates::CoordinateSystem::cartesian:
          {
            // The last coordinate is vertical and does not affect a World
            // Builder topography query. The first dim-1 coordinates locate
            // the point along the reference surface.
            for (unsigned int d = 0; d < dim-1; ++d)
              position[d] = surface_point[d];
            position[dim-1] = 0.0;
            break;
          }

          case Utilities::Coordinates::CoordinateSystem::spherical:
          {
            // ASPECT surface points contain longitude (and colatitude in
            // 3d). Reconstruct a point on the undeformed outer sphere and
            // let World Builder convert the Cartesian point into its own
            // natural coordinates.
            std::array<double,dim> spherical_position;
            spherical_position[0] =
              reference_surface_radius<dim>(geometry_model);
            for (unsigned int d = 0; d < dim-1; ++d)
              spherical_position[d+1] = surface_point[d];

            position =
              geometry_model.natural_to_cartesian_coordinates(spherical_position);
            break;
          }

          default:
            AssertThrow(false,
                        ExcMessage("The World Builder initial topography model "
                                   "currently supports Cartesian and spherical "
                                   "geometry models."));
        }

      return world_builder->properties(Utilities::convert_point_to_array(position),
                                       0.0,
      {{{6,0,0}}})[0];
    }



    template <int dim>
    double
    WorldBuilder<dim>::max_topography () const
    {
      return maximum_topography;
    }



  }
#  endif
#endif
}


namespace aspect
{
#ifdef ASPECT_WITH_WORLD_BUILDER
#  if WORLD_BUILDER_VERSION_GTE(1,1,1)
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(
      WorldBuilder,
      "world builder",
      "Specify the initial topography through the Geodynamic World Builder. "
      "World Builder topography models are queried from the file selected by "
      "the top-level parameter `World builder file'. ASPECT obtains the "
      "maximum expected topography directly from World Builder.")
  }
#  endif
#endif
}

/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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



#ifndef _aspect_world_feature_features_continental_plate_h
#define _aspect_world_feature_features_continental_plate_h

#include <aspect/world_builder/coordinate_system/interface.h>
#include <aspect/world_builder/utilities.h>
#include <aspect/world_builder/coordinate_systems.h>


namespace aspect
{
  namespace WorldBuilder
  {
    namespace CoordinateSystem
    {

      class Cartesian : public Interface
      {
        public:
          /**
           * constructor
           */
          Cartesian();

          /**
           * Destructor
           */
          ~Cartesian();

          /**
           * Read in the world builder file
           */
          virtual
          void read(ptree &property_tree);

          /**
           * Returns what the natural coordinate system for this geometry model is.
           */
          virtual
          aspect::WorldBuilder::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const;

          /**
           * Takes the Cartesian points (x,z or x,y,z) and returns standardized
           * coordinates which are most 'natural' to the geometry model. For a box
           * this will  be (x,z) in 2d or (x,y,z) in 3d, and for a spheroid geometry
           * model it  will be (radius, longitude) in 2d and (radius, longitude,
           * latitude) in 3d.
           */
          virtual
          std::array<double,3> cartesian_to_natural_coordinates(const std::array<double,3> &position) const;

          /**
           * Undoes the action of cartesian_to_natural_coordinates, and turns the
           * coordinate system which is most 'natural' to the geometry model into
           * Cartesian coordinates.
           */
          virtual
          std::array<double,3> natural_to_cartesian_coordinates(const std::array<double,3> &position) const;


        private:

      };
    }
  }
}

#endif

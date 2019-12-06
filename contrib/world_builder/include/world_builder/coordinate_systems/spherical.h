/*
  Copyright (C) 2018 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef _world_builder_coordinate_systems_spherical_h
#define _world_builder_coordinate_systems_spherical_h

#include <world_builder/utilities.h>
#include <world_builder/coordinate_systems/interface.h>


namespace WorldBuilder
{
  namespace CoordinateSystems
  {
    /**
     * Register header file
     */
    //WB_REGISTER_COORDINATE_SYSTEM_HEADER(Spherical)

    /**
     * This implements a Spherical geometry model.The Cartesian geometry model
     * doesn't do anything with the coordinates, but is needed to have a common
     * interface for all the geometry models.
     */
    class Spherical : public Interface
    {
      public:
        /**
         * constructor
         */
        Spherical(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~Spherical();

        /**
         * declare and read in the world builder file into the parameters class
         */
        static
        void declare_entries(Parameters &prm, const std::string &parent_name = "");

        /**
         * declare and read in the world builder file into the parameters class
         */
        virtual
        void parse_entries(Parameters &prm);

        /**
         * Returns what the natural coordinate system for this Coordinate System is.
         */
        virtual
        CoordinateSystem natural_coordinate_system() const;

        /**
         * Returns what method should be used to go down with an angle into
         * the domain.
         * \sa DepthMethod
         */
        virtual
        DepthMethod depth_method() const;

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


        /**
         * Computes the distance between two points which are on the same depth.
         * The input is two 2d points at that depth. The distance is returned
         * in meters. To compute this distance from spherical coordinates the
         * radius is multiplied with the central angle. The central angle is
         * compute in the most accurate way given on wikipedia
         * (https://en.wikipedia.org/wiki/Great-circle_distance), through a
         * special case of the Vincenty formula for an ellipsoid with equal
         * major and minor axes (https://doi.org/10.1179/sre.1975.23.176.88).
         */
        virtual
        double distance_between_points_at_same_depth(const Point<3> &point_1, const Point<3> &point_3) const;

        /**
         * What depth method the spherical coordinates use.
         */
        DepthMethod used_depth_method;

      private:

    };
  }
}

#endif

/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#ifndef WORLD_BUILDER_COORDINATE_SYSTEMS_CARTESIAN_H
#define WORLD_BUILDER_COORDINATE_SYSTEMS_CARTESIAN_H

#include "world_builder/coordinate_systems/interface.h"

namespace WorldBuilder
{

  namespace CoordinateSystems
  {
    /**
     * Register header file
     */
    //WB_REGISTER_COORDINATE_SYSTEM_HEADER(Cartesian)


    /**
     * This implements a Cartesian geometry model.The Cartesian geometry model
     * doesn't do anything with the coordinates, but is needed to have a common
     * interface for all the geometry models.
     */
    class Cartesian final : public Interface
    {
      public:
        /**
         * constructor
         */
        Cartesian(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~Cartesian() override final;

        /**
         * declare and read in the world builder file into the parameters class
         */
        static
        void declare_entries(Parameters &prm, const std::string &parent_name = "");

        /**
         * declare and read in the world builder file into the parameters class
         */
        void parse_entries(Parameters &prm) override final;


        /**
         * Returns what the natural coordinate system for this Coordinate System is.
         */
        CoordinateSystem natural_coordinate_system() const override final;

        /**
         * Returns what method should be used to go down with an angle into
         * the domain.
         * \sa DepthMethod
         */
        DepthMethod depth_method() const override final;

        /**
         * Takes the Cartesian points (x,z or x,y,z) and returns standardized
         * coordinates which are most 'natural' to the geometry model. For a box
         * this will  be (x,z) in 2d or (x,y,z) in 3d, and for a spheroid geometry
         * model it  will be (radius, longitude) in 2d and (radius, longitude,
         * latitude) in 3d.
         */
        std::array<double,3> cartesian_to_natural_coordinates(const std::array<double,3> &position) const override final;

        /**
         * Undoes the action of cartesian_to_natural_coordinates, and turns the
         * coordinate system which is most 'natural' to the geometry model into
         * Cartesian coordinates.
         */
        std::array<double,3> natural_to_cartesian_coordinates(const std::array<double,3> &position) const override final;


        /**
         * Computes the distance between two points which are on the same depth.
         * The input is two 3d points at that depth. It is implemented as the
         */
        double distance_between_points_at_same_depth(const Point<3> &point_1, const Point<3> &point_2) const override final;

        /**
         * Returns the max model depth. This always returns infinity for Cartesian
         * models.
         */
        virtual
        double max_model_depth() const override final;

      private:

    };
  } // namespace CoordinateSystems
} // namespace WorldBuilder

#endif

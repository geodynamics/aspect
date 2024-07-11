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

#ifndef WORLD_BUILDER_COORDINATE_SYSTEMS_UNIFORM_H
#define WORLD_BUILDER_COORDINATE_SYSTEMS_UNIFORM_H

#include "world_builder/gravity_model/interface.h"

#include "world_builder/utilities.h"


namespace WorldBuilder
{

  namespace GravityModel
  {
    /**
     * Register header file
     */
    //WB_REGISTER_COORDINATE_SYSTEM_HEADER(Uniform)


    /**
     * This implements a Uniform geometry model.The Uniform geometry model
     * doesn't do anything with the coordinates, but is needed to have a common
     * interface for all the geometry models.
     */
    class Uniform final : public Interface
    {
      public:
        /**
         * constructor
         */
        Uniform(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~Uniform() override final;

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
         * Returns the gravity vector in a Cartesian coordinate system at a given position,
         * which has a constant magitude for the whole domain. The vector points down in
         * cartesian coordinates and to the center of the sphere in spherical coordinates.
         */
        virtual
        Point<3> gravity_vector(Point<3> point) const override final;

        /**
         * Returns the norm of the gravity at a given position, which is a constant
         * number for the whole domain.
         */
        virtual
        double gravity_norm(Point<3> point) const override final;

      private:

        /**
         * The uniform gravity.
         */
        double gravity_magnitude;

    };
  } // namespace CoordinateSystems
} // namespace WorldBuilder

#endif

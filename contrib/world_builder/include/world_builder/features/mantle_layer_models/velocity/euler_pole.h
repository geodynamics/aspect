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

#ifndef WORLD_BUILDER_FEATURES_MANTLE_LAYER_MODELS_VELOCITY_EULER_POLE_H
#define WORLD_BUILDER_FEATURES_MANTLE_LAYER_MODELS_VELOCITY_EULER_POLE_H


#include "world_builder/features/mantle_layer_models/velocity/interface.h"
#include "world_builder/features/feature_utilities.h"
#include "world_builder/objects/surface.h"


namespace WorldBuilder
{
  namespace Features
  {
    using namespace FeatureUtilities;
    namespace MantleLayerModels
    {
      namespace Velocity
      {
        /**
         * This class represents a mantle layer and can implement submodules
         * for velocity and composition. These submodules determine what
         * the returned velocity or composition of the velocity and composition
         * functions of this class will be.
         */
        class EulerPole final: public Interface
        {
          public:
            /**
             * constructor
             */
            EulerPole(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~EulerPole() override final;

            /**
             * declare and read in the world builder file into the parameters class
             */
            static
            void declare_entries(Parameters &prm, const std::string &parent_name = "");

            /**
             * declare and read in the world builder file into the parameters class
             */
            void parse_entries(Parameters &prm, const std::vector<Point<2>> &coordinates) override final;


            /**
             * Returns a velocity based on the given position, depth in the model,
             * gravity and current velocity.
             */
            std::array<double,3> get_velocity(const Point<3> &position,
                                              const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                              const double depth,
                                              const double gravity,
                                              std::array<double,3> velocity,
                                              const double feature_min_depth,
                                              const double feature_max_depth) const override final;


          private:
            // euler pole velocity submodule parameters
            double min_depth;
            Objects::Surface min_depth_surface;
            double max_depth;
            Objects::Surface max_depth_surface;
            Point<3> euler_pole;
            Operations operation;

        };
      } // namespace Velocity
    } // namespace MantleLayerModels
  } // namespace Features
} // namespace WorldBuilder

#endif

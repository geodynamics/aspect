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

#ifndef WORLD_BUILDER_FEATURES_OCEANIC_PLATE_MODELS_TEMPERATURE_PLATE_MODEL_CONSTANT_AGE_H
#define WORLD_BUILDER_FEATURES_OCEANIC_PLATE_MODELS_TEMPERATURE_PLATE_MODEL_CONSTANT_AGE_H


#include "world_builder/features/oceanic_plate_models/temperature/interface.h"
#include "world_builder/features/feature_utilities.h"
#include "world_builder/objects/surface.h"


namespace WorldBuilder
{
  class Parameters;
  class World;

  namespace Features
  {
    using namespace FeatureUtilities;
    namespace OceanicPlateModels
    {
      namespace Temperature
      {
        /**
         * This class represents a oceanic plate of constant age and can implement submodules
         * for temperature and composition. These submodules determine what
         * the returned temperature or composition of the temperature and composition
         * functions of this class will be.
         * In this plugin, the temperature of the plate is derived from the plate model,
         * featuring cooling of a plate from a ridge by vertical heat conduction,
         * while the effects of horizontal heat conduction are ignored.
         */
        class PlateModelConstantAge final: public Interface
        {
          public:
            /**
             * constructor
             */
            PlateModelConstantAge(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~PlateModelConstantAge() override final;

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
             * Returns a temperature based on the given position, depth in the model,
             * gravity and current temperature.
             */
            double get_temperature(const Point<3> &position,
                                   const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                   const double depth,
                                   const double gravity,
                                   double temperature,
                                   const double feature_min_depth,
                                   const double feature_max_depth) const override final;


          private:
            // plate model constant age temperature submodule parameters
            double min_depth;
            Objects::Surface min_depth_surface;
            double max_depth;
            Objects::Surface max_depth_surface;
            double top_temperature;
            double bottom_temperature;
            double plate_age;
            Operations operation;

        };
      } // namespace Temperature
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder

#endif

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

#ifndef WORLD_BUILDER_FEATURES_OCEANIC_PLATE_MODELS_TEMPERATURE_PLATE_MODEL_H
#define WORLD_BUILDER_FEATURES_OCEANIC_PLATE_MODELS_TEMPERATURE_PLATE_MODEL_H


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
         * This class represents a oceanic plate and can implement submodules
         * for temperature and composition. These submodules determine what
         * the returned temperature or composition of the temperature and composition
         * functions of this class will be.
         * In this plugin, the temperature of the plate is derived from the plate model,
         * featuring continuous cooling of a plate from a spreading ridge by both vertical
         * and horizontal heat conduction.
         */
        class PlateModel final: public Interface
        {
          public:
            /**
             * constructor
             */
            PlateModel(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~PlateModel() override final;

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
            // plate model temperature submodule parameters
            double min_depth;
            Objects::Surface min_depth_surface;
            double max_depth;
            Objects::Surface max_depth_surface;
            double top_temperature;
            double bottom_temperature;
            std::pair<std::vector<double>,std::vector<double>> spreading_velocities;
            std::vector<std::vector<Point<2> > > mid_oceanic_ridges;
            std::vector<std::vector<double>> spreading_velocities_at_each_ridge_point;
            Operations operation;

        };
      } // namespace Temperature
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder

#endif

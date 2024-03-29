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

#ifndef WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_MODELS_TEMPERATURE_ADIABATIC_H
#define WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_MODELS_TEMPERATURE_ADIABATIC_H


#include "world_builder/features/subducting_plate_models/temperature/interface.h"
#include "world_builder/features/feature_utilities.h"


namespace WorldBuilder
{

  namespace Features
  {
    namespace SubductingPlateModels
    {
      namespace Temperature
      {
        /**
         * This class represents a subducting plate and can implement submodules
         * for temperature and composition. These submodules determine what
         * the returned temperature or composition of the temperature and composition
         * functions of this class will be.
         */
        class Adiabatic final: public Interface
        {
          public:
            /**
             * constructor
             */
            Adiabatic(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~Adiabatic() override final;

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
             * Returns a temperature based on the given position, depth in the model,
             * gravity and current temperature.
             */
            double get_temperature(const Point<3> &position,
                                   const double depth,
                                   const double gravity,
                                   double temperature,
                                   const double feature_min_depth,
                                   const double feature_max_depth,
                                   const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_planes,
                                   const AdditionalParameters &additional_parameters) const override final;


          private:
            // adiabatic temperature submodule parameters
            double min_depth;
            double max_depth;
            /**
             * Todo
             */
            double potential_mantle_temperature;


            /**
             * Todo
             */
            double thermal_expansion_coefficient;

            /**
             * Todo
             */
            double specific_heat;

            Operations operation;

        };
      } // namespace Temperature
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder

#endif

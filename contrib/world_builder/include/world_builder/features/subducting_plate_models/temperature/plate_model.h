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

#ifndef _world_builder_features_subducting_plate_temperature_plate_model_h
#define _world_builder_features_subducting_plate_temperature_plate_model_h

#include <world_builder/features/subducting_plate_models/temperature/interface.h>
#include <world_builder/world.h>


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
        class PlateModel : public Interface
        {
          public:
            /**
             * constructor
             */
            PlateModel(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~PlateModel();

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
             * Returns a temperature based on the given position, depth in the model,
             * gravity and current temperature.
             */
            virtual
            double get_temperature(const Point<3> &position,
                                   const double depth,
                                   const double gravity,
                                   double temperature,
                                   const double feature_min_depth,
                                   const double feature_max_depth,
                                   const std::map<std::string,double> &distance_from_planes) const;


          private:
            // plate model temperature submodule parameters
            double min_depth;
            double max_depth;
            double density;
            double plate_velocity;
            double thermal_conductivity;
            double thermal_expansion_coefficient;
            double specific_heat;
            double potential_mantle_temperature;
            double surface_temperature;
            bool adiabatic_heating;
            std::string operation;

        };
      }
    }
  }
}

#endif

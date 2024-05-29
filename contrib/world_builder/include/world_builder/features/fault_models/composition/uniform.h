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

#ifndef WORLD_BUILDER_FEATURES_FAULT_MODELS_COMPOSITION_UNIFORM_H
#define WORLD_BUILDER_FEATURES_FAULT_MODELS_COMPOSITION_UNIFORM_H


#include "world_builder/features/fault_models/composition/interface.h"


namespace WorldBuilder
{

  namespace Features
  {
    namespace FaultModels
    {
      namespace Composition
      {
        /**
         * This class represents a subducting plate and can implement submodules
         * for temperature and composition. These submodules determine what
         * the returned temperature or composition of the temperature and composition
         * functions of this class will be.
         */
        class Uniform final: public Interface
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
             * Returns a composition based on the given position, depth in the model,
             * gravity and current composition.
             */
            double get_composition(const Point<3> &position,
                                   const double depth,
                                   const unsigned int composition_number,
                                   double composition,
                                   const double feature_min_depth,
                                   const double feature_max_depth,
                                   const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_planes,
                                   const AdditionalParameters &additional_parameters) const override final;


          private:
            // uniform composition submodule parameters
            double min_depth;
            double max_depth;
            std::vector<unsigned int> compositions;
            std::vector<double> fractions;
            Operations operation;

        };
      } // namespace Composition
    } // namespace FaultModels
  } // namespace Features
} // namespace WorldBuilder

#endif

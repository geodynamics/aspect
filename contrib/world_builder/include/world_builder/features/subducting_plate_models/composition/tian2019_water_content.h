/*
  Copyright (C) 2018 - 2024 by the authors of the World Builder code.

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

#ifndef WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_MODELS_COMPOSITION_WATER_CONTENT_H
#define WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_MODELS_COMPOSITION_WATER_CONTENT_H


#include "world_builder/features/subducting_plate_models/composition/interface.h"
#include <mutex>


namespace WorldBuilder
{

  namespace Features
  {
    namespace SubductingPlateModels
    {
      namespace Composition
      {
        /**
         * This class represents the bound water content in a subducting plate and can implement
         * submodules for computing this bound water content. These submodules determine what
         * the returned water content will be based on the the temperature and pressure at a point
         * within the world. Currently, the bound water content can be determined for four different
         * lithologies: sediments, mid-ocean ridge basalts, gabbros, and peridotites, using
         * parameterized phase diagrams from Tian et al., 2019 (https://doi.org/10.1029/2019GC008488).
         */
        class TianWaterContent final: public Interface
        {
          public:
            /**
             * constructor
             */
            TianWaterContent(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~TianWaterContent() override final;

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
             *  Calculates what the maximum bound water content is at a point given a pressure and
             * temperature. This can be done for 4 different lithologies, "sediment", "gabbro",
             * "MORB", and "peridotite"
             */
            double calculate_water_content(double pressure,
                                           double temperature) const;
            /**
             * Returns a value for the bound water contend based on the given position, depth in the model,
             * gravity, and the temperature at that point.
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
            // TianWaterContent composition submodule parameters
            double min_depth;
            double max_depth;
            double density; // Density used to calculate the lithostatic pressure
            std::vector<unsigned int> compositions;
            Operations operation;
            std::string lithology_str;
            double max_water_content;

            enum LithologyName
            {
              peridotite,
              gabbro,
              MORB,
              sediment
            };
            LithologyName lithology_type;

            // Define the coefficients for the polynomials for 3 quantities: LR which represents the
            // enthalpy change of the dehydration reactions, c_sat which represents the volatile saturation
            // content, and Td which represents the onset temperature of dehydration. The first row is for
            // 'peridotite', the second row is for 'gabbro', the third row is for 'MORB', and the fourth row
            // is for 'sediment'.
            std::vector<std::vector<double>> LR_poly =
            {
              {-19.0609, 168.983, -630.032, 1281.84, -1543.14, 1111.88, -459.142, 95.4143, 1.97246},
              {-1.81745, 7.67198, -10.8507, 5.09329, 8.14519},
              {-1.78177, 7.50871, -10.4840, 5.19725, 7.96365},
              {-2.03283, 10.8186, -21.2119, 18.3351, -6.48711, 8.32459}
            };

            std::vector<std::vector<double>> c_sat_poly =
            {
              {0.00115628, 2.42179},
              {-0.0176673, 0.0893044, 1.52732},
              {0.0102725, -0.115390, 0.324452, 1.41588},
              {-0.150662, 0.301807, 1.01867}
            };

            std::vector<std::vector<double>> Td_poly =
            {
              {-15.4627, 94.9716, 636.603},
              {-1.72277, 20.5898, 637.517},
              {-3.81280, 22.7809, 638.049},
              {2.83277, -24.7593, 85.9090, 524.898}
            };

            // Maximum pressure for the lithologies (Peridotite, Gabbro, MORB, Sediment). These are required because
            // Above these pressures, the parameterized phase diagrams break down and the solubility goes to infinity.
            double cutoff_pressure;
        };
      } // namespace Composition
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder

#endif

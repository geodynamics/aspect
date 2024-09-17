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

#ifndef WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_MODELS_TEMPERATURE_MASS_CONSERVING_H
#define WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_MODELS_TEMPERATURE_MASS_CONSERVING_H


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
         * This class represents a subducting plate temperature model. The temperature
         * model uses the heat content (proportional to to thermal mass anomaly) to
         * define a smooth temperature profile that conserves mass along the slab length.
         * An empirical (linear) model is used to define how the minimum temperature
         * increases with depth and how the location of the minimum temperature shifts
         * into the slab interior. The slab is divided in to top and bottom parts,
         * which meet at the location where the minimum temperature occurs in the slab.
         * For the bottom slab the temperature is defined by a half-space cooling model.
         * For the top of the slab the temperature is defined by one side of a 1D infinite
         * space cooling model. The age of the overriding plate is used so the slab temperature
         * at shallow depth smoothly transitions to the temperature of the overriding plate:
         * this is not perfect, and is affected by the value of "top truncation" parameter
         * subducting plate. Also note that the parameter "thickness" for the subducting plate
         * segments needs to be defined but is not used.
         * Note that the empirical model used to define how Tmin increases with depth
         * and how the position of Tmin shift with depth is expected to change somewhat
         * after better calibrating with further tests.
         */
        class MassConserving final: public Interface
        {
          public:
            /**
             * constructor
             */
            MassConserving(WorldBuilder::World *world);

            /**
             * Destructor
             */
            ~MassConserving() override final;

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

            /**
             * Returns a temperature based on the given heat content, temperatures, effective plate age,
             * and adjusted distance to the coldest point in the slab. The temperature is formulated by
             * the analytical solutions.
             */
            double get_temperature_analytic(const double top_heat_content,
                                            const double min_temperature,
                                            const double background_temperature,
                                            const double temperature_,
                                            const double plate_velocity,
                                            const double effective_plate_age,
                                            const double adjusted_distance) const;


          private:
            //  temperature submodule parameters
            double min_depth;
            double max_depth;
            double density;
            std::vector<std::vector<double>> subducting_velocities;
            std::pair<std::vector<double>,std::vector<double>> ridge_spreading_velocities;
            std::vector<std::vector<double>> ridge_spreading_velocities_at_each_ridge_point;
            double mantle_coupling_depth;
            double forearc_cooling_factor;
            double thermal_conductivity;
            double thermal_expansion_coefficient;
            double specific_heat;
            double thermal_diffusivity;
            double potential_mantle_temperature;
            double surface_temperature;
            double taper_distance;
            bool adiabatic_heating;
            std::vector<std::vector<Point<2>>> mid_oceanic_ridges;
            Operations operation;
            enum ReferenceModelName
            {
              half_space_model,
              plate_model
            };
            ReferenceModelName reference_model_name;
            const int plate_model_summation_number = 100; // for the plate model
            bool apply_spline;
            unsigned int spline_n_points;
        };
      } // namespace Temperature
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder

#endif

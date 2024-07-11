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

#include "world_builder/features/plume_models/temperature/gaussian.h"


#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/array.h"
#include "world_builder/world.h"

#include <algorithm>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace PlumeModels
    {
      namespace Temperature
      {
        Gaussian::Gaussian(WorldBuilder::World *world_)
          :
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "gaussian";
        }

        Gaussian::~Gaussian()
          = default;

        void
        Gaussian::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `max distance fault` center to the required parameters.
          prm.declare_entry("", Types::Object({"centerline temperatures"}),
                            "Gaussian temperature model. The temperature is interpolated between the plume center "
                            "and margin (as defined by the plume feature) using a Gaussian function: "
                            "T(r) = T_center(z) exp(-r^2/(2 sigma^2). "
                            "The temperature at the plume centerline T_center can be changed with depth by defining "
                            "an array of depths and centerline temperatures, and temperature is interpolated linearly "
                            "with depth. Similarly, the sigma of the Gaussian function (relative to the width of "
                            "the plume as given by the plume feature) can be changed with depth. "
                            "Temperature is always interpolated in a horizonzal/radial plane, except for the plume "
                            "head: If the first depth of the plume centerline and the minimum depth of the plume "
                            "feature are different, an ellipsoidal plume head is created in this depth range. "
                            "Within this plume head, temperature is interpolated radially, i.e., depending on the "
                            "distance from the center of the ellipsoid.");

          prm.declare_entry("depths", Types::Array(Types::Double(0)),
                            "The list of depths where both the temperature in the center of the plume "
                            "and the width of the temperature anomaly in terms of the sigma of a Gaussian "
                            "function can be provided. Temperature is interpolated linearly in vertical "
                            "direction between these depths. Units: m.");
          prm.declare_entry("centerline temperatures", Types::Array(Types::Double(0)),
                            "The temperature at the center of this feature in degree Kelvin."
                            "If the value is below zero, then an adiabatic temperature is used.");
          prm.declare_entry("gaussian sigmas", Types::Array(Types::Double(0.3)),
                            "The sigma (standard deviation) of the Gaussian function used to compute the "
                            "temperature distribution within the plume. This sigma is non-dimensional, i.e. "
                            "it is defined relative to the distance between the plume center and margin as "
                            "defined by the plume feature. Choosing a sigma of 1 therefore means that the "
                            "temperature at the plume margin is set to a fraction of 1/sqrt(e) (approx. 0.61) "
                            "of the centerline temperature. To achieve a smoother transition between the "
                            "plume temperature and the outside temperature a smaller values has to be chosen "
                            "for the gaussian sigmas.");
        }

        void
        Gaussian::parse_entries(Parameters &prm)
        {
          operation = string_operations_to_enum(prm.get<std::string>("operation"));

          depths = prm.get_vector<double>("depths");
          center_temperatures = prm.get_vector<double>("centerline temperatures");
          gaussian_sigmas = prm.get_vector<double>("gaussian sigmas");

          WBAssert(center_temperatures.size() == depths.size() && gaussian_sigmas.size() == depths.size(),
                   "The depths, center_temperatures and gaussian_sigmas arrays need to have the same number of entries. "
                   "At the moment there are: "
                   << depths.size() << " depth entries, "
                   << center_temperatures.size() << " centerline temperature entries, and "
                   << gaussian_sigmas.size() << " gaussian sigma entries!");
        }


        double
        Gaussian::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                  const Objects::NaturalCoordinate & /*position_in_natural_coordinates*/,
                                  const double depth,
                                  const double gravity_norm,
                                  double temperature_,
                                  const double feature_min_depth,
                                  const double feature_max_depth,
                                  const double relative_distance_from_center) const
        {
          if (depth <= feature_max_depth && depth >= feature_min_depth && relative_distance_from_center <= 1.)
            {
              // Figure out if the point is within the plume
              auto upper = std::upper_bound(depths.begin(), depths.end(), depth);

              double center_temperature_local;
              double gaussian_sigma;

              if (upper - depths.begin() == 0)
                {
                  center_temperature_local = center_temperatures.front();
                  gaussian_sigma = gaussian_sigmas.front();
                }
              else if (upper - depths.end() == 0)
                {
                  center_temperature_local = center_temperatures.back();
                  gaussian_sigma = gaussian_sigmas.back();
                }
              else
                {
                  const unsigned int index = static_cast<unsigned int>(std::distance(depths.begin(), upper));
                  const double fraction = (depth - depths[index-1]) / (depths[index] - depths[index-1]);

                  center_temperature_local = (1-fraction) * center_temperatures[index-1] + fraction * center_temperatures[index];
                  gaussian_sigma = (1-fraction) * gaussian_sigmas[index-1] + fraction * gaussian_sigmas[index];
                }

              if (center_temperature_local < 0)
                {
                  center_temperature_local =  this->world->potential_mantle_temperature *
                                              std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                        this->world->specific_heat) * depth);
                }

              const double new_temperature =   center_temperature_local * std::exp(-relative_distance_from_center/(2.*std::pow(gaussian_sigma, 2)));

              return apply_operation(operation,temperature_,new_temperature);

            }

          return temperature_;
        }

        WB_REGISTER_FEATURE_PLUME_TEMPERATURE_MODEL(Gaussian, gaussian)
      } // namespace Temperature
    } // namespace PlumeModels
  } // namespace Features
} // namespace WorldBuilder


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

#include "world_builder/features/plume.h"


#include "world_builder/features/plume_models/composition/interface.h"
#include "world_builder/features/plume_models/grains/interface.h"
#include "world_builder/features/plume_models/temperature/interface.h"
#include "world_builder/features/feature_utilities.h"
#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/plugin_system.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"

#include <iostream>
#include <algorithm>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    Plume::Plume(WorldBuilder::World *world_)
      :
      min_depth(NaN::DSNAN),
      max_depth(NaN::DSNAN)
    {
      this->world = world_;
      this->name = "plume";
    }

    Plume::~Plume()
      = default;



    void Plume::make_snippet(Parameters &prm)
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;

      const std::string path = prm.get_full_json_path();

      Pointer((path + "/body").c_str()).Set(declarations,"object");
      Pointer((path + "/body/model").c_str()).Set(declarations,"plume");
      Pointer((path + "/body/name").c_str()).Set(declarations,"${1:My Plume}");
      Pointer((path + "/body/coordinates").c_str()).Create(declarations).SetArray();
      Pointer((path + "/body/temperature models").c_str()).Create(declarations).SetArray();
      Pointer((path + "/body/composition models").c_str()).Create(declarations).SetArray();
    }



    void
    Plume::declare_entries(Parameters &prm,
                           const std::string & /*unused*/,
                           const std::vector<std::string> &required_entries)
    {
      prm.declare_entry("", Types::Object(required_entries), "Plume object. Requires properties `model` and `coordinates`.");

      prm.declare_entry("min depth", Types::Double(0),
                        "The depth from which this feature is present, in other words, the "
                        "depth of the tip of the plume. If the first entry in the cross "
                        "section depths has a greater depth, an ellipsoidal plume head will "
                        "be added in between. Units: m.");
      prm.declare_entry("max depth", Types::Double(std::numeric_limits<double>::max()),
                        "The depth to which this feature is present. Units: m.");
      prm.declare_entry("cross section depths", Types::Array(Types::Double(0)),
                        "The depths of the elliptic cross section of the plume. Units: m.");
      prm.declare_entry("semi-major axis", Types::Array(Types::Double(100.e3)),
                        "The lengths of the semi-major axes of the elliptic cross sections of the plume. "
                        "In spherical coordinates, this is in degrees, otherwise in meters.");
      prm.declare_entry("eccentricity", Types::Array(Types::Double(0)),
                        "The eccentricities of the cross sections.");
      prm.declare_entry("rotation angles", Types::Array(Types::Double(0)),
                        "The directions that the semi-major axis of the elliptic cross-sections "
                        "are pointing to, in degrees. This direction is expressed as the angle from "
                        "geographic North in spherical coordinates, or as the angle from the Y axis "
                        "(clockwise) in Cartesian coordinates. "
                        "The angle should be between 0 and 360 degrees.");

      prm.declare_entry("temperature models",
                        Types::PluginSystem("", Features::PlumeModels::Temperature::Interface::declare_entries, {"model"}),
                        "A list of temperature models.");
      prm.declare_entry("composition models",
                        Types::PluginSystem("", Features::PlumeModels::Composition::Interface::declare_entries, {"model"}),
                        "A list of composition models.");
      prm.declare_entry("grains models",
                        Types::PluginSystem("", Features::PlumeModels::Grains::Interface::declare_entries, {"model"}),
                        "A list of grains models.");
    }

    void
    Plume::parse_entries(Parameters &prm)
    {
      const CoordinateSystem coordinate_system = prm.coordinate_system->natural_coordinate_system();

      this->name = prm.get<std::string>("name");

      std::string tag = prm.get<std::string>("tag");
      if (tag == "")
        {
          tag = "plume";
        }
      this->tag_index = FeatureUtilities::add_vector_unique(this->world->feature_tags,tag);

      this->get_coordinates("coordinates", prm, coordinate_system);

      min_depth = prm.get<double>("min depth");
      max_depth = prm.get<double>("max depth");

      depths = prm.get_vector<double>("cross section depths");
      semi_major_axis_lengths = prm.get_vector<double>("semi-major axis");
      eccentricities = prm.get_vector<double>("eccentricity");
      rotation_angles = prm.get_vector<double>("rotation angles");

      for (unsigned int i = 0; i < depths.size()-1; ++i)
        WBAssert(depths[i] < depths[i+1],
                 "The depths of the elliptic cross sections of the plume need to be listed in ascending order.");

      WBAssert(depths.size() == coordinates.size(),
               "The cross section depths array needs to have the same number of entries as there are coordinates. At the moment there are: "
               << depths.size()
               << " depth entries but "
               << coordinates.size()
               << " coordinates!");

      WBAssert(semi_major_axis_lengths.size() == coordinates.size(),
               "The semi-major axis array needs to have the same number of entries as there are coordinates. At the moment there are: "
               << semi_major_axis_lengths.size()
               << " semi-major axis entries but "
               << coordinates.size()
               << " coordinates!");

      WBAssert(eccentricities.size() == coordinates.size(),
               "The eccentricity array needs to have the same number of entries as there are coordinates. At the moment there are: "
               << eccentricities.size()
               << " eccentricity entries but "
               << coordinates.size()
               << " coordinates!");

      WBAssert(rotation_angles.size() == coordinates.size(),
               "The rotation angles array needs to have the same number of entries as there are coordinates. At the moment there are: "
               << rotation_angles.size()
               << " rotation angle entries but "
               << coordinates.size()
               << " coordinates!");


      // Convert degrees to radians, convert from geographical to mathematical
      for (double &rotation_angle : rotation_angles)
        rotation_angle = Consts::PI/2. - rotation_angle * Consts::PI/180.;

      // convert semi_major_axis_lengths to radians if we are in spherical coordinates
      if (world->parameters.coordinate_system->natural_coordinate_system() == CoordinateSystem::spherical)
        for (double &semi_major_axis_length : semi_major_axis_lengths)
          {
            semi_major_axis_length *= Consts::PI/180.;
          }

      prm.get_unique_pointers<Features::PlumeModels::Temperature::Interface>("temperature models", temperature_models);

      prm.enter_subsection("temperature models");
      {
        for (unsigned int i = 0; i < temperature_models.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              temperature_models[i]->parse_entries(prm);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();


      prm.get_unique_pointers<Features::PlumeModels::Composition::Interface>("composition models", composition_models);

      prm.enter_subsection("composition models");
      {
        for (unsigned int i = 0; i < composition_models.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              composition_models[i]->parse_entries(prm);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();


      prm.get_unique_pointers<Features::PlumeModels::Grains::Interface>("grains models", grains_models);

      prm.enter_subsection("grains models");
      {
        for (unsigned int i = 0; i < grains_models.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              grains_models[i]->parse_entries(prm);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();

    }



    void
    Plume::properties(const Point<3> &position_in_cartesian_coordinates,
                      const Objects::NaturalCoordinate &position_in_natural_coordinates,
                      const double depth,
                      const std::vector<std::array<unsigned int,3>> &properties,
                      const double gravity_norm,
                      const std::vector<size_t> &entry_in_output,
                      std::vector<double> &output) const
    {
      // Figure out if the point is within the plume
      auto upper = std::upper_bound(depths.begin(), depths.end(), depth);

      Point<2> plume_center(coordinates[0]);
      double semi_major_axis_length;
      double eccentricity;
      double rotation_angle;

      if (depth < min_depth)
        return;
      else if (upper - depths.begin() == 0)
        {
          // interpolate to make the top of the plume spherical
          plume_center = coordinates.front();
          eccentricity = eccentricities.front();
          rotation_angle = rotation_angles.front();

          const double fraction = (depth - min_depth) / (depths.front() - min_depth);

          const double a = depths.front() - min_depth;
          const double b = semi_major_axis_lengths.front();
          const double y = (1. - fraction) * a;

          // use ellipse equation:
          semi_major_axis_length = std::sqrt((1 - std::pow(y/a,2)) * b*b);
        }
      else if (upper - depths.end() == 0)
        {
          plume_center = coordinates.back();
          semi_major_axis_length = semi_major_axis_lengths.back();
          eccentricity = eccentricities.back();
          rotation_angle = rotation_angles.back();
        }
      else
        {
          const unsigned int index = static_cast<unsigned int>(std::distance(depths.begin(), upper));
          const double fraction = (depth - depths[index-1]) / (depths[index] - depths[index-1]);

          plume_center[0] = (1-fraction) * coordinates[index-1][0] + fraction * (coordinates[index][0]);
          plume_center[1] = (1-fraction) * coordinates[index-1][1] + fraction * (coordinates[index][1]);

          semi_major_axis_length = (1-fraction) * semi_major_axis_lengths[index-1] + fraction * semi_major_axis_lengths[index];
          eccentricity = (1-fraction) * eccentricities[index-1] + fraction * eccentricities[index];

          // For the angles, we only want to go between zero and pi, and we have to make sure we
          // interpolate the values close to zero/pi correctly:
          rotation_angle = interpolate_angle_across_zero(rotation_angles[index-1], rotation_angles[index], fraction);
        }

      const Point<2> surface_point(position_in_natural_coordinates.get_surface_coordinates(),
                                   world->parameters.coordinate_system->natural_coordinate_system());

      double relative_distance_from_center =
        WorldBuilder::Utilities::fraction_from_ellipse_center(plume_center,
                                                              semi_major_axis_length,
                                                              eccentricity,
                                                              rotation_angle,
                                                              surface_point);

      // If we are in the tip, we have to compute the difference diffently:
      if (depth >= min_depth && depth < depths.front())
        {
          const double a = semi_major_axis_lengths.front();
          const double b = a * std::sqrt(1 - std::pow(eccentricity, 2));
          const double c = depths.front() - min_depth;

          const double x = (surface_point[0] - plume_center[0]) * std::cos(rotation_angle) + (surface_point[1] - plume_center[1])* std::sin(rotation_angle);
          const double y = -(surface_point[0] - plume_center[0]) * std::sin(rotation_angle) + (surface_point[1] - plume_center[1])* std::cos(rotation_angle);
          const double z = depths.front() - depth;

          // use ellipsoid equation:
          relative_distance_from_center = (x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c);
        }

      if (depth <= max_depth && depth >= min_depth && relative_distance_from_center <= 1.)
        {

          for (unsigned int i_property = 0; i_property < properties.size(); ++i_property)
            {
              switch (properties[i_property][0])
                {
                  case 1:  // temperature
                  {
                    for (const auto &temperature_model: temperature_models)
                      {
                        output[entry_in_output[i_property]] = temperature_model->get_temperature(position_in_cartesian_coordinates,
                                                                                                 position_in_natural_coordinates,
                                                                                                 depth,
                                                                                                 gravity_norm,
                                                                                                 output[entry_in_output[i_property]],
                                                                                                 min_depth,
                                                                                                 max_depth,
                                                                                                 relative_distance_from_center);

                        WBAssert(!std::isnan(output[entry_in_output[i_property]]), "Temperature is not a number: " << output[entry_in_output[i_property]]
                                 << ", based on a temperature model with the name " << temperature_model->get_name() << ", in feature " << this->name);
                        WBAssert(std::isfinite(output[entry_in_output[i_property]]), "Temperature is not finite: " << output[entry_in_output[i_property]]
                                 << ", based on a temperature model with the name " << temperature_model->get_name() << ", in feature " << this->name);

                      }
                    break;
                    case 2: // composition

                      for (const auto &composition_model: composition_models)
                        {
                          output[entry_in_output[i_property]] = composition_model->get_composition(position_in_cartesian_coordinates,
                                                                                                   position_in_natural_coordinates,
                                                                                                   depth,
                                                                                                   properties[i_property][1],
                                                                                                   output[entry_in_output[i_property]],
                                                                                                   min_depth,
                                                                                                   max_depth);

                          WBAssert(!std::isnan(output[entry_in_output[i_property]]), "Composition is not a number: " << output[entry_in_output[i_property]]
                                   << ", based on a composition model with the name " << composition_model->get_name() << ", in feature " << this->name);
                          WBAssert(std::isfinite(output[entry_in_output[i_property]]), "Composition is not finite: " << output[entry_in_output[i_property]]
                                   << ", based on a composition model with the name " << composition_model->get_name() << ", in feature " << this->name);

                        }

                      break;
                    }
                  case 3: // grains
                  {
                    WorldBuilder::grains  grains(output,properties[i_property][2],entry_in_output[i_property]);
                    for (const auto &grains_model: grains_models)
                      {
                        grains = grains_model->get_grains(position_in_cartesian_coordinates,
                                                          position_in_natural_coordinates,
                                                          depth,
                                                          properties[i_property][1],
                                                          grains,
                                                          min_depth,
                                                          max_depth);

                      }
                    grains.unroll_into(output,entry_in_output[i_property]);
                    break;
                  }
                  case 4:
                  {
                    output[entry_in_output[i_property]] = static_cast<double>(tag_index);
                    break;
                  }
                  default:
                  {
                    WBAssertThrow(false,
                                  "Internal error: Unimplemented property provided. " <<
                                  "Only temperature (1), composition (2), grains (3) or tag (4) are allowed. "
                                  "Provided property number was: " << properties[i_property][0]);
                  }
                }
            }
        }
    }

    WB_REGISTER_FEATURE(Plume, plume)

  } // namespace Features
} // namespace WorldBuilder
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

#include "world_builder/features/oceanic_plate.h"

#include "world_builder/features/oceanic_plate_models/composition/interface.h"
#include "world_builder/features/oceanic_plate_models/grains/interface.h"
#include "world_builder/features/oceanic_plate_models/temperature/interface.h"
#include "world_builder/features/feature_utilities.h"
#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/plugin_system.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    OceanicPlate::OceanicPlate(WorldBuilder::World *world_)
      :
      min_depth(NaN::DSNAN),
      max_depth(NaN::DSNAN)
    {
      this->world = world_;
      this->name = "oceanic plate";
    }

    OceanicPlate::~OceanicPlate()
      = default;



    void OceanicPlate::make_snippet(Parameters &prm)
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;

      const std::string path = prm.get_full_json_path();

      /*
      ideally:
       {
        "model": "fault",
        "name": "${1:default name}",
        "dip point":[0.0,0.0],
        "coordinates": [[0.0,0.0]],
        "segments": [],
        "sections": [],
        "temperature models":[{"model":"uniform", "temperature":600.0}],
        "composition models":[{"model":"uniform", "compositions": [0], "fractions":[1.0]}]
       }
       */

      Pointer((path + "/body").c_str()).Set(declarations,"object");
      Pointer((path + "/body/model").c_str()).Set(declarations,"oceanic plate");
      Pointer((path + "/body/name").c_str()).Set(declarations,"${1:My Oceanic Plate}");
      Pointer((path + "/body/coordinates").c_str()).Create(declarations).SetArray();
      Pointer((path + "/body/temperature models").c_str()).Create(declarations).SetArray();
      Pointer((path + "/body/composition models").c_str()).Create(declarations).SetArray();
    }



    void
    OceanicPlate::declare_entries(Parameters &prm,
                                  const std::string & /*unused*/,
                                  const std::vector<std::string> &required_entries)
    {
      prm.declare_entry("", Types::Object(required_entries), "Oceanic plate object. Requires properties `model` and `coordinates`.");

      prm.declare_entry("min depth", Types::OneOf(Types::Double(0),Types::Array(Types::ValueAtPoints(0., 2.))),
                        "The depth from which this feature is present");
      prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(), 2.))),
                        "The depth to which this feature is present");
      prm.declare_entry("temperature models",
                        Types::PluginSystem("", Features::OceanicPlateModels::Temperature::Interface::declare_entries, {"model"}),
                        "A list of temperature models.");
      prm.declare_entry("composition models",
                        Types::PluginSystem("", Features::OceanicPlateModels::Composition::Interface::declare_entries, {"model"}),
                        "A list of composition models.");
      prm.declare_entry("grains models",
                        Types::PluginSystem("", Features::OceanicPlateModels::Grains::Interface::declare_entries, {"model"}),
                        "A list of grains models.");
    }

    void
    OceanicPlate::parse_entries(Parameters &prm)
    {

      const CoordinateSystem coordinate_system = prm.coordinate_system->natural_coordinate_system();

      this->name = prm.get<std::string>("name");

      std::string tag = prm.get<std::string>("tag");
      if (tag == "")
        {
          tag = "oceanic plate";
        }
      this->tag_index = FeatureUtilities::add_vector_unique(this->world->feature_tags,tag);

      this->get_coordinates("coordinates", prm, coordinate_system);

      min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
      min_depth = min_depth_surface.minimum;

      max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
      max_depth = max_depth_surface.maximum;


      prm.get_unique_pointers<Features::OceanicPlateModels::Temperature::Interface>("temperature models", temperature_models);

      prm.enter_subsection("temperature models");
      {
        for (unsigned int i = 0; i < temperature_models.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              temperature_models[i]->parse_entries(prm,coordinates);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();


      prm.get_unique_pointers<Features::OceanicPlateModels::Composition::Interface>("composition models", composition_models);

      prm.enter_subsection("composition models");
      {
        for (unsigned int i = 0; i < composition_models.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              composition_models[i]->parse_entries(prm,coordinates);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();


      prm.get_unique_pointers<Features::OceanicPlateModels::Grains::Interface>("grains models", grains_models);

      prm.enter_subsection("grains models");
      {
        for (unsigned int i = 0; i < grains_models.size(); ++i)
          {
            prm.enter_subsection(std::to_string(i));
            {
              grains_models[i]->parse_entries(prm,coordinates);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();
    }


    void
    OceanicPlate::properties(const Point<3> &position_in_cartesian_coordinates,
                             const Objects::NaturalCoordinate &position_in_natural_coordinates,
                             const double depth,
                             const std::vector<std::array<unsigned int,3>> &properties,
                             const double gravity_norm,
                             const std::vector<size_t> &entry_in_output,
                             std::vector<double> &output) const
    {
      if (depth <= max_depth && depth >= min_depth &&
          WorldBuilder::Utilities::polygon_contains_point(coordinates, Point<2>(position_in_natural_coordinates.get_surface_coordinates(),
                                                                                world->parameters.coordinate_system->natural_coordinate_system())))
        {
          const double min_depth_local = min_depth_surface.constant_value ? min_depth : min_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
          const double max_depth_local = max_depth_surface.constant_value ? max_depth : max_depth_surface.local_value(position_in_natural_coordinates.get_surface_point()).interpolated_value;
          if (depth <= max_depth_local &&  depth >= min_depth_local)
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
                                                                                                     min_depth_local,
                                                                                                     max_depth_local);

                            WBAssert(!std::isnan(output[entry_in_output[i_property]]), "Temperature is not a number: " << output[entry_in_output[i_property]]
                                     << ", based on a temperature model with the name " << temperature_model->get_name() << ", in feature " << this->name);
                            WBAssert(std::isfinite(output[entry_in_output[i_property]]), "Temperature is not a finite: " << output[entry_in_output[i_property]]
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
                                                                                                       min_depth_local,
                                                                                                       max_depth_local);

                              WBAssert(!std::isnan(output[entry_in_output[i_property]]), "Composition is not a number: " << output[entry_in_output[i_property]]
                                       << ", based on a composition model with the name " << composition_model->get_name() << ", in feature " << this->name);
                              WBAssert(std::isfinite(output[entry_in_output[i_property]]), "Composition is not a finite: " << output[entry_in_output[i_property]]
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
                                                              min_depth_local,
                                                              max_depth_local);

                          }
                        grains.unroll_into(output,entry_in_output[i_property]);
                        break;
                      }
                      case 4:
                      {
                        output[entry_in_output[i_property]] = static_cast<double>(tag_index);
                        break;
                      }
                      break;
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
    }

    /**
     * Register plugin
     */
    WB_REGISTER_FEATURE(OceanicPlate, oceanic plate)
  } // namespace Features
} // namespace WorldBuilder


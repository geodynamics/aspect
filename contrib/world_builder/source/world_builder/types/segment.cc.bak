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

#include "world_builder/types/segment.h"


namespace WorldBuilder
{
  namespace Features
  {
    namespace FaultModels
    {
      namespace Composition
      {
        class Interface;
      }  // namespace Composition
      namespace Grains
      {
        class Interface;
      }  // namespace Grains
      namespace Temperature
      {
        class Interface;
      }  // namespace Temperature
    }  // namespace FaultModels
    namespace SubductingPlateModels
    {
      namespace Composition
      {
        class Interface;
      }  // namespace Composition
      namespace Grains
      {
        class Interface;
      }  // namespace Grains
      namespace Temperature
      {
        class Interface;
      }  // namespace Temperature
    }  // namespace SubductingPlateModels
  }  // namespace Features

  namespace Types
  {
    Segment::Segment(const double default_length_,
                     const WorldBuilder::Point<2> &default_thickness_,
                     const WorldBuilder::Point<2> &default_top_truncation_,
                     const WorldBuilder::Point<2> &default_angle_,
                     const Types::Interface &temperature_plugin_system_,
                     const Types::Interface &composition_plugin_system_,
                     const Types::Interface &grains_plugin_system_)
      :
      value_length(default_length_),
      default_length(default_length_),
      value_thickness(default_thickness_),
      default_thickness(default_thickness_),
      default_top_truncation(default_top_truncation_),
      value_angle(default_angle_),
      default_angle(default_angle_),
      temperature_plugin_system(temperature_plugin_system_.clone()),
      composition_plugin_system(composition_plugin_system_.clone()),
      grains_plugin_system(grains_plugin_system_.clone())
    {
      this->type_name = Types::type::Segment;
    }


    Segment::Segment(Segment const &other)
      :
      value_length(other.default_length),
      default_length(other.default_length),
      value_thickness(other.default_thickness),
      default_thickness(other.default_thickness),
      default_top_truncation(other.default_top_truncation),
      value_angle(other.default_angle),
      default_angle(other.default_angle),
      temperature_plugin_system(other.temperature_plugin_system->clone()),
      composition_plugin_system(other.composition_plugin_system->clone()),
      grains_plugin_system(other.grains_plugin_system->clone())
    {
      this->type_name = Types::type::Segment;
    }


    Segment::~Segment ()
      = default;



    void
    Segment::write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const
    {
      using namespace rapidjson;
      prm.enter_subsection(name);
      {
        Document &declarations = prm.declarations;
        std::string base = prm.get_full_json_path();

        Pointer((base + "/type").c_str()).Set(declarations,"object");
        Pointer((base + "/additionalProperties").c_str()).Set(declarations,false);
        Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());
        std::vector<std::string> restricted_values = {"length", "thickness", "angle"};
        for (unsigned int i = 0; i < restricted_values.size(); ++i)
          {
            if (!restricted_values[i].empty())
              {
                if (i == 0 && Pointer((base + "/required").c_str()).Get(declarations) == nullptr)
                  {
                    // The enum array doesn't exist yet, so we create it and fill it.
                    Pointer((base + "/required/0").c_str()).Create(declarations);
                    Pointer((base + "/required/0").c_str()).Set(declarations, restricted_values[i].c_str());
                  }
                else
                  {
                    // The enum array already exist yet, so we add an element to the end.
                    Pointer((base + "/required/-").c_str()).Set(declarations, restricted_values[i].c_str());
                  }
              }
          }

        prm.enter_subsection("properties");
        {
          base = prm.get_full_json_path();
          Pointer((base + "/length/type").c_str()).Set(declarations,"number");

          Pointer((base + "/thickness/type").c_str()).Set(declarations,"array");
          Pointer((base + "/thickness/minItems").c_str()).Set(declarations,1);
          Pointer((base + "/thickness/maxItems").c_str()).Set(declarations,2);
          Pointer((base + "/thickness/items/type").c_str()).Set(declarations,"number");

          Pointer((base + "/top truncation/type").c_str()).Set(declarations,"array");
          Pointer((base + "/top truncation/minItems").c_str()).Set(declarations,1);
          Pointer((base + "/top truncation/maxItems").c_str()).Set(declarations,2);
          Pointer((base + "/top truncation/items/type").c_str()).Set(declarations,"number");

          Pointer((base + "/angle/type").c_str()).Set(declarations,"array");
          Pointer((base + "/angle/minItems").c_str()).Set(declarations,1);
          Pointer((base + "/angle/maxItems").c_str()).Set(declarations,2);
          Pointer((base + "/angle/items/type").c_str()).Set(declarations,"number");

          temperature_plugin_system->write_schema(prm, "temperature models", "");
          composition_plugin_system->write_schema(prm, "composition models", "");
          grains_plugin_system->write_schema(prm, "grains models", "");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }
  } // namespace Types
} // namespace WorldBuilder


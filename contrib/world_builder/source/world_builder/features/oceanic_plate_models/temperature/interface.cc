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

#include "world_builder/features/oceanic_plate_models/temperature/interface.h"

#include <algorithm>

#include "world_builder/types/string.h"

namespace WorldBuilder
{
  namespace Features
  {
    namespace OceanicPlateModels
    {
      namespace Temperature
      {
        Interface::Interface()
          = default;

        Interface::~Interface ()
          = default;

        void
        Interface::declare_entries(Parameters &prm,
                                   const std::string &parent_name,
                                   const std::vector<std::string> &required_entries)
        {
          // The type needs to be stored in a separate value, otherwise there are memory issues
          const Types::String type = Types::String("replace", std::vector<std::string> {"replace", "add", "subtract"});

          prm.declare_model_entries("temperature",parent_name, get_declare_map(),required_entries,
          {
            std::tuple<std::string,const WorldBuilder::Types::Interface &, std::string>
            {
              "operation", type,
              "Whether the value should replace any value previously defined at this location (replace), "
              "add the value to the previously define value (add) or subtract the value to the previously "
              "define value (subtract)."
            }
          });
        }


        void
        Interface::registerType(const std::string &name,
                                void ( *declare_entries)(Parameters &, const std::string &),
                                ObjectFactory *factory)
        {
          get_factory_map()[name] = factory;
          get_declare_map()[name] = declare_entries;
        }

        std::unique_ptr<Interface>
        Interface::create(const std::string &name, WorldBuilder::World *world)
        {
          std::string lower_case_name;
          std::transform(name.begin(),
                         name.end(),
                         std::back_inserter(lower_case_name),
                         ::tolower);;

          // Have a nice assert message to check whether a plugin exists in the case
          // of a debug compilation.
          WBAssertThrow(get_factory_map().find(lower_case_name) != get_factory_map().end(),
                        "Internal error: Plugin with name '" << lower_case_name << "' is not found. "
                        "The size of factories is " << get_factory_map().size() << '.');

          // Using at() because the [] will just insert values
          // which is undesirable in this case. An exception is
          // thrown when the name is not present.
          return get_factory_map().at(lower_case_name)->create(world);
        }
      } // namespace Temperature
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder


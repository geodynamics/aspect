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

#include "world_builder/features/oceanic_plate_models/grains/interface.h"

#include <algorithm>

#include "world_builder/types/object.h"
#include "world_builder/types/string.h"

namespace WorldBuilder
{
  namespace Features
  {
    namespace OceanicPlateModels
    {
      namespace Grains
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
          unsigned int counter = 0;
          for (auto &it : get_declare_map())
            {
              // prevent infinite recursion
              if (it.first != parent_name)
                {
                  prm.enter_subsection("oneOf");
                  {
                    prm.enter_subsection(std::to_string(counter));
                    {
                      prm.enter_subsection("properties");
                      {
                        prm.declare_entry("", Types::Object(required_entries), "grains object");

                        prm.declare_entry("model", Types::String("",it.first),
                                          "The name of the grains model.");

                        it.second(prm, parent_name);
                      }
                      prm.leave_subsection();
                    }
                    prm.leave_subsection();
                  }
                  prm.leave_subsection();

                  counter++;

                }
            }
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
      } // namespace Grains
    } // namespace OceanicPlateModels
  } // namespace Features
} // namespace WorldBuilder


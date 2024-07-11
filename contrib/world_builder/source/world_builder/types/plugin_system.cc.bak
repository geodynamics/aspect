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
#include "world_builder/types/plugin_system.h"



namespace WorldBuilder
{
  namespace Types
  {
    PluginSystem::PluginSystem(std::string default_value_,
                               void ( *declare_entries_)(Parameters &, const std::string &, const std::vector<std::string> &),
                               std::vector<std::string> required_entries_,
                               const bool allow_multiple_)
      :
      default_value(std::move(default_value_)),
      declare_entries(declare_entries_),
      required_entries(std::move(required_entries_)),
      allow_multiple(allow_multiple_)
    {
      this->type_name = Types::type::PluginSystem;

      WBAssert(declare_entries_ != nullptr, "declare entries may not be a null pointer.");
    }


    PluginSystem::PluginSystem(PluginSystem const &plugin_system)
      :
      default_value(plugin_system.default_value),
      declare_entries(plugin_system.declare_entries),
      required_entries(plugin_system.required_entries),
      allow_multiple(plugin_system.allow_multiple)
    {
      this->type_name = Types::type::PluginSystem;
    }

    PluginSystem::~PluginSystem ()
      = default;

    void
    PluginSystem::write_schema(Parameters &prm,
                               const std::string &name,
                               const std::string &documentation) const
    {
      using namespace rapidjson;

      prm.enter_subsection(name);
      {
        const std::string path = prm.get_full_json_path();
        Pointer((path + "/description").c_str()).Set(prm.declarations,documentation.c_str());
        Pointer((path + "/default value").c_str()).Set(prm.declarations,default_value.c_str());

        if (allow_multiple)
          {
            Pointer((path + "/type").c_str()).Set(prm.declarations,"array");

            prm.enter_subsection("items");
            {
              WBAssert(this->declare_entries != nullptr, "No declare entries given.");

              this->declare_entries(prm, name, required_entries);
            }
            prm.leave_subsection();
          }
        else
          {
            Pointer((path + "/type").c_str()).Set(prm.declarations,"object");

            WBAssert(this->declare_entries != nullptr, "No declare entries given.");
            this->declare_entries(prm, name, required_entries);


          }
      }
      prm.leave_subsection();
    }
  } // namespace Types
} // namespace WorldBuilder


/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

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

#ifndef _world_feature_types_plugin_system_h
#define _world_feature_types_plugin_system_h

#include <world_builder/types/interface.h>
#include <world_builder/features/interface.h>


namespace WorldBuilder
{
  namespace Types
  {

    /**
     * This class represents a plate tectonic feature class, such as the
     * continental plate class, oceanic plate class and subduction zone class.
     */
    class PluginSystem : public Interface
    {
      public:
        /**
         * constructor
         */
        PluginSystem(std::string default_value_,
                     void ( *declare_entries)(Parameters &, const std::string &, const std::vector<std::string> &),
                     std::vector<std::string> required_entries,
                     const bool allow_multiple = true);


        /**
         * Copy constructor
         */
        PluginSystem(PluginSystem const &plugin_system);

        /**
         * Destructor
         */
        ~PluginSystem();

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;


        std::string default_value;
        void( *declare_entries)(Parameters &, const std::string &, const std::vector<std::string> &);
        std::vector<std::string> required_entries;
        bool allow_multiple;

      protected:
        PluginSystem *clone_impl() const override final
        {
          return new PluginSystem(*this);
        };

      private:

    };
  }
}

#endif

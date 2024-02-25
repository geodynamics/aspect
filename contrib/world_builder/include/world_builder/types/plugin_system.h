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

#ifndef WORLD_BUILDER_TYPES_PLUGIN_SYSTEM_H
#define WORLD_BUILDER_TYPES_PLUGIN_SYSTEM_H


#include "world_builder/features/interface.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {

    /**
     * This class represents a plate tectonic feature class, such as the
     * continental plate class, oceanic plate class and subduction zone class.
     */
    class PluginSystem final: public Interface
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
        ~PluginSystem() override final;

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
  } // namespace Types
} // namespace WorldBuilder

#endif

/*
  Copyright (C) 2018 by the authors of the World Builder code.

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

#ifndef _world_feature_types_string_h
#define _world_feature_types_string_h

#include <world_builder/types/interface.h>


namespace WorldBuilder
{
  namespace Types
  {

    /**
     * This class represents a continental plate and can implement submodules
     * for temperature and composition. These submodules determine what
     * the returned temperature or composition of the temperature and composition
     * functions of this class will be.
     */
    class String : public Interface
    {
      public:
        /**
         * constructor
         */
        String(const std::string default_value);

        /**
         * constructor
         */
        String(const std::string default_value, const std::string restricted_value);

        /**
         * constructor
         */
        String(const std::string default_value, const std::vector<std::string> &restricted_values);

        /**
         * constructor
         */
        //String(std::string default_value, std::string description);


        /**
         * A constructor for the clone and set_entry function
         */
        String(std::string value, std::string default_value, std::string description);

        /**
         * Destructor
         */
        ~String();

        /**
         * Todo
         */
        virtual
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const;

        /**
         * clone
         */
        virtual
        std::unique_ptr<Interface> clone() const;


        std::string value;
        std::string default_value;
        std::string description;
        std::vector<std::string> restricted_values;
    };
  }
}

#endif

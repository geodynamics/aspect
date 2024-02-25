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

#ifndef WORLD_BUILDER_TYPES_STRING_H
#define WORLD_BUILDER_TYPES_STRING_H

#include <vector>

#include "world_builder/types/interface.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {

    /**
     * This class represents a continental plate and can implement submodules
     * for temperature and composition. These submodules determine what
     * the returned temperature or composition of the temperature and composition
     * functions of this class will be.
     */
    class String final: public Interface
    {
      public:
        /**
         * constructor
         */
        String(std::string default_value);

        /**
         * constructor
         */
        String(std::string default_value, const std::string &restricted_value);

        /**
         * constructor
         */
        String(std::string default_value, std::vector<std::string> restricted_values);

        /**
         * constructor
         */
        //String(std::string default_value, std::string description);


        /**
         * A constructor for the clone and set_entry function
         */
        String(std::string value, std::string default_value, std::string description);

        /**
         * Copy constructor
         */
        String(String const &other);

        /**
         * Destructor
         */
        ~String() final;

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;


        std::string value;
        std::string default_value;
        std::string description;
        std::vector<std::string> restricted_values;


      protected:
        String *clone_impl() const override final
        {
          return new String(*this);
        };
    };
  } // namespace Types
} // namespace WorldBuilder

#endif

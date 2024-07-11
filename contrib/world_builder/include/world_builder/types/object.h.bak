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

#ifndef WORLD_BUILDER_TYPES_OBJECT_H
#define WORLD_BUILDER_TYPES_OBJECT_H

#include <vector>

#include "world_builder/types/interface.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {
    /**
     * This class represents an object of values which can be of any Types. It
     * stores the type of values it is holding and an vector of indices where
     * the values are actually stored in the parameters vector which hold all
     * the values of that type. This type can also hold a unique pointer to
     * the type is should hold. This is used for declaring types.The difference
    * between an object and a list is that the object holds the values retrievable
    * by index, and a list holds the values strictly ordered, only accessible
    * an iterator. An other difference it that lists have a name.
     */
    class Object final: public Interface
    {
      public:
        /**
         * Constructor for the declaration
         */
        Object(std::vector<std::string> required = std::vector<std::string>(),
               const bool additional_properties = false);

        /**
         * Copy constructor
         */
        Object(Object const &other);

        /**
         * Destructor
         */
        ~Object() final;

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;


        /**
         * Todo
         */
        std::vector<std::string> required;
        bool additional_properties;

      protected:
        Object *clone_impl() const override final
        {
          return new Object(*this);
        };

    };
  } // namespace Types
} // namespace WorldBuilder

#endif

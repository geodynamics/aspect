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

#ifndef WORLD_BUILDER_TYPES_ARRAY_H
#define WORLD_BUILDER_TYPES_ARRAY_H


#include "world_builder/types/interface.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {
    /**
     * This class represents an array of values which can be of any Types. It
     * stores the type of values it is holding and an vector of indices where
     * the values are actually stored in the parameters vector which hold all
     * the values of that type. This type can also hold a unique pointer to
     * the type is should hold. This is used for declaring types.The difference
    * between an array and a list is that the array holds the values retrievable
    * by index, and a list holds the values strictly ordered, only accessible
    * an iterator. An other difference it that lists have a name.
     */
    class Array final: public Interface
    {
      public:
        /**
         * Constructor for the declaration
         */
        Array(const Interface &type,
              const unsigned int min_items = 0,
              const unsigned int max_items = std::numeric_limits<unsigned int>::max(),
              const bool unique_items = false);

        /**
         * Constructor for cloning an array.
         */
        Array(Array const &other);


        /**
         * Destructor
         */
        ~Array() override final;

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;

        /**
         * An enum of the type which this class points to
         * @see Types::type
         */
        Types::type inner_type;

        /**
         * This class is sometimes responsible for the object it points to, but
         * sometimes it is not responsible for the object is points to.
         * When it is responsible the unique_inner_type points to it and the
         * inner_type should have size zero. When it is not responsible,
         * unique_inner_type should point to the nullptr and inner_type should
         * have a size larger then zero.
         * @see inner_type_index
         */
        std::unique_ptr<Interface> inner_type_ptr;


        /**
         * Todo
         */
        bool required;

        /**
         * Todo
         */
        unsigned int min_items;

        /**
         * Todo
         */
        unsigned int max_items;

        /**
         * Todo
         */
        bool unique_items;

      protected:
        Array *clone_impl() const override final
        {
          return new Array(*this);
        };
    };
  } // namespace Types
} // namespace WorldBuilder

#endif

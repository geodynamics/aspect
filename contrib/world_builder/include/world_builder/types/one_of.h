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

#ifndef WORLD_BUILDER_TYPES_ONE_OF_H
#define WORLD_BUILDER_TYPES_ONE_OF_H


#include "world_builder/types/interface.h"
#include <vector>


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {
    /**
     * This class represents an a single choice between options (one of).
     */
    class OneOf final: public Interface
    {
      public:
        /**
         * Constructor for the declaration
         */
        OneOf(const Interface &type_1,
              const Interface &type_2);

        /**
         * Constructor for cloning an array.
         */
        OneOf(OneOf const &other);


        /**
         * Destructor
         */
        ~OneOf() override final;

        /**
         * Write schema
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
        std::vector<std::unique_ptr<Interface>> inner_types_ptr;

      protected:
        OneOf *clone_impl() const override final
        {
          return new OneOf(*this);
        };
    };
  } // namespace Types
} // namespace WorldBuilder

#endif

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

#ifndef _world_feature_types_bool_h
#define _world_feature_types_bool_h

#include <world_builder/types/interface.h>


namespace WorldBuilder
{
  namespace Types
  {

    /**
     * This class represents a bool value with documentation
     */
    class Bool : public Interface
    {
      public:
        /**
         * A constructor for the load_entry function
         */
        Bool(const bool default_value);

        /**
         * Copy constructor
         */
        Bool(Bool const &other);

        /**
         * Destructor
         */
        ~Bool();


        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;

        bool default_value;

      protected:
        Bool *clone_impl() const override final
        {
          return new Bool(*this);
        };

      private:

    };
  }
}

#endif

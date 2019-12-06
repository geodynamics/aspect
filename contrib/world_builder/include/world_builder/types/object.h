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

#ifndef _world_feature_types_object_h
#define _world_feature_types_object_h

#include <world_builder/types/interface.h>


namespace WorldBuilder
{
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
    class Object : public Interface
    {
      public:
        /**
         * Constructor for the declaration
         */
        Object(const std::vector<std::string> required = std::vector<std::string>(),
               const bool additional_properties = false);

        /**
         * Destructor
         */
        ~Object();


        /**
         * Clone. The caller will be responsible for the liftime of this
         * object, return a unique pointer. This clone can only be used
         * when inner_type.size() == 0.
         */
        virtual
        std::unique_ptr<Interface>   clone() const;

        /**
         * Todo
         */
        virtual
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const;


        /**
         * Todo
         */
        std::vector<std::string> required;
        bool additional_properties;


    };
  }
}

#endif

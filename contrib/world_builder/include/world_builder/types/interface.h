/*
  Copyright (C) 2018 - 2020 by the authors of the World Builder code.

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

#ifndef _world_builder_types_interface_h
#define _world_builder_types_interface_h

#include <string>
#include <vector>
#include <memory>
#include <map>

#include <rapidjson/document.h>
#include <rapidjson/pointer.h>

namespace WorldBuilder
{
  class Parameters;
  /**
   * This class is an interface for the specific plate tectonic feature classes,
   * such as continental plate, oceanic plate and subduction zone.
   */
  namespace Types
  {

    enum class type
    {
      None,Bool,String,Double,Int,UnsignedInt,Array,Object,List,Point2D,Point3D,CoordinateSystem,PluginSystem,Segment,ConstantLayer
    };

    class Interface
    {
      public:
        /**
         * Constructor
         */
        Interface();

        /**
         * Destructor
         */
        virtual
        ~Interface();

        /**
         * Todo
         */
        virtual
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const = 0;

        /**
         * clone
         */
        virtual
        std::unique_ptr<Interface> clone() const
        {
          return std::unique_ptr<Interface>(clone_impl());
        }

        /**
         * read in the world builder file
         */

        type get_type() const;


      protected:
        type type_name;


        virtual
        Interface *clone_impl() const = 0;
    };
  }
}

#endif

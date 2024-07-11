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

#ifndef WORLD_BUILDER_TYPES_POINT_H
#define WORLD_BUILDER_TYPES_POINT_H


#include "world_builder/point.h"
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
    template <unsigned int dim>
    class Point final: public Interface
    {
      public:

        /**
         * A constructor used for the load_entry function
         */
        Point();

        /**
         * A constructor used for the load_entry function
         */
        Point(const WorldBuilder::Point<dim> &default_value, std::string description);

        /**
         * A constructor used for cloning and the set_entry function
         */
        Point(const WorldBuilder::Point<dim> &value, const WorldBuilder::Point<dim> &default_value, std::string description);

        /**
         * Copy constructor
         */
        Point(Point const &other);

        /**
         * Destructor
         */
        ~Point() override final;

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;

        /**
         * dot product
         */
        double operator*(const Point<dim> &point) const;


        /**
         * Multiply the vector with a scalar
         */
        WorldBuilder::Point<dim> operator*(const double scalar) const;

        /**
         * add two points
         */
        WorldBuilder::Point<dim> operator+(const Point<dim> &point) const;


        /**
         * Subtract two points
         */
        WorldBuilder::Point<dim> operator-(const Point<dim> &point) const;

        /**
         * access index (const)
         */
        const double &operator[](const unsigned int index) const;


        /**
         * access index
         */
        double &operator[](const unsigned int index);



        WorldBuilder::Point<dim> value;
        WorldBuilder::Point<dim> default_value;
        std::string description;

      protected:
        Point *clone_impl() const override final
        {
          return new Point(*this);
        };

      private:

    };

    template<unsigned int dim>
    WorldBuilder::Point<dim> operator*(const double scalar, const Point<dim> &point);
  } // namespace Types
} // namespace WorldBuilder

#endif

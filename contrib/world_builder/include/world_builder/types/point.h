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

#ifndef _world_feature_types_point_h
#define _world_feature_types_point_h

#include <world_builder/types/interface.h>
#include <world_builder/point.h>
#include <world_builder/assert.h>
#include <world_builder/utilities.h>
#include <world_builder/parameters.h>

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
    template <int dim>
    class Point : public Interface
    {
      public:

        /**
         * A constructor used for the load_entry function
         */
        Point();

        /**
         * A constructor used for the load_entry function
         */
        Point(const WorldBuilder::Point<dim> &default_value, const std::string &description);

        /**
         * A constructor used for cloning and the set_entry function
         */
        Point(const WorldBuilder::Point<dim> &value, const WorldBuilder::Point<dim> &default_value, const std::string &description);

        /**
         * Copy constructor
         */
        Point(Point const &other);

        /**
         * Destructor
         */
        ~Point();

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
         * Substract two points
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

    template <int dim>
    inline
    Point<dim>::Point()
      :
      value(WorldBuilder::Point<dim>(CoordinateSystem::cartesian)),
      default_value(WorldBuilder::Point<dim>(CoordinateSystem::cartesian)),
      description("")
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <int dim>
    inline
    Point<dim>::Point(const WorldBuilder::Point<dim> &default_value_,
                      const std::string &description_)
      :
      value(default_value_),
      default_value(default_value_),
      description(description_)
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <int dim>
    inline
    Point<dim>::Point(const WorldBuilder::Point<dim> &value_,
                      const WorldBuilder::Point<dim> &default_value_,
                      const std::string &description_)
      :
      value(value_),
      default_value(default_value_),
      description(description_)
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <int dim>
    inline
    Point<dim>::Point(Point const &other)
      :
      value(other.value),
      default_value(other.default_value),
      description(other.description)
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <int dim>
    inline
    Point<dim>::~Point ()
    {}

    template<int dim>
    inline void
    Point<dim>::write_schema(Parameters &prm,
                             const std::string &name,
                             const std::string &documentation) const
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;
      const std::string base = prm.get_full_json_path() + "/" + name;

      Pointer((base + "/type").c_str()).Set(declarations,"array");
      Pointer((base + "/minItems").c_str()).Set(declarations,dim);
      Pointer((base + "/maxItems").c_str()).Set(declarations,dim);
      Pointer((base + "/documentation").c_str()).Set(declarations,documentation.c_str());
      Pointer((base + "/items/type").c_str()).Set(declarations,"number");
      // todo: default value
    }


    template<int dim>
    inline double 
    Point<dim>::operator*(const Point<dim> &point_) const
    {
      const std::array<double,dim> array = point_.value.get_array();
      double dot_product = 0;
      for (unsigned int i = 0; i < dim; ++i)
        dot_product += value[i] * array[i];
      return dot_product;
    }


    template<int dim>
    inline WorldBuilder::Point<dim>
    Point<dim>::operator*(const double scalar) const
    {
      // initialize the array to zero.
      std::array<double,dim> array = WorldBuilder::Point<dim>(value.get_coordinate_system()).get_array();
      for (unsigned int i = 0; i < dim; ++i)
        array[i] += value[i] * scalar;
      return WorldBuilder::Point<dim>(array,value.get_coordinate_system());
    }

    template<int dim>
    inline WorldBuilder::Point<dim>
    Point<dim>::operator+(const Point<dim> &point_) const
    {
      WorldBuilder::Point<dim> point_tmp(value);
      point_tmp += point_.value;

      return point_tmp;
    }

    template<int dim>
    inline WorldBuilder::Point<dim>
    Point<dim>::operator-(const Point<dim> &point_) const
    {
      WorldBuilder::Point<dim> point_tmp(value);
      point_tmp -= point_.value;

      return point_tmp;
    }

    /**
     * access index
     */
    template<int dim>
    inline const double &
    Point<dim>::operator[](const unsigned int index) const
    {
      return value[index];
    }


    /**
     * access index
     */
    template<int dim>
    inline double &
    Point<dim>::operator[](const unsigned int index)
    {
      return value[index];
    }

    /**
     * Multiplies a Types::Point<dim> with a scalr and returns a
     * WorldBuilder::Point<dim>.
     */
    template<int dim>
    inline WorldBuilder::Point<dim>
    operator*(const double scalar, const Point<dim> &point)
    {
      return point.value*scalar;
    }
  }
}

#endif

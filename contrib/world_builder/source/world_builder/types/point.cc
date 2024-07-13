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
#include "world_builder/types/point.h"


#include "world_builder/parameters.h"

namespace WorldBuilder
{
  namespace Types
  {

    template <unsigned int dim>
    Point<dim>::Point()
      :
      value(WorldBuilder::Point<dim>(CoordinateSystem::cartesian)),
      default_value(WorldBuilder::Point<dim>(CoordinateSystem::cartesian))

    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <unsigned int dim>
    Point<dim>::Point(const WorldBuilder::Point<dim> &default_value_,
                      std::string description_)
      :
      value(default_value_),
      default_value(default_value_),
      description(std::move(description_))
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <unsigned int dim>
    Point<dim>::Point(const WorldBuilder::Point<dim> &value_,
                      const WorldBuilder::Point<dim> &default_value_,
                      std::string description_)
      :
      value(value_),
      default_value(default_value_),
      description(std::move(description_))
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <unsigned int dim>
    Point<dim>::Point(Point const &other)
      :
      value(other.value),
      default_value(other.default_value),
      description(other.description)
    {
      this->type_name = dim == 2 ? Types::type::Point2D : Types::type::Point3D;
    }

    template <unsigned int dim>
    Point<dim>::~Point ()
      = default;

    template<unsigned int dim>
    void
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
      Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());
      Pointer((base + "/items/type").c_str()).Set(declarations,"number");
      // todo: default value
    }


    template<unsigned int dim>
    double Point<dim>::operator*(const Point<dim> &point_) const
    {
      const std::array<double,dim> array = point_.value.get_array();
      double dot_product = 0;
      for (unsigned int i = 0; i < dim; ++i)
        dot_product += value[i] * array[i];
      return dot_product;
    }


    template<unsigned int dim>
    WorldBuilder::Point<dim>
    Point<dim>::operator*(const double scalar) const
    {
      // initialize the array to zero.
      std::array<double,dim> array = WorldBuilder::Point<dim>(value.get_coordinate_system()).get_array();
      for (unsigned int i = 0; i < dim; ++i)
        array[i] += value[i] * scalar;
      return WorldBuilder::Point<dim>(array,value.get_coordinate_system());
    }

    template<unsigned int dim>
    WorldBuilder::Point<dim>
    Point<dim>::operator+(const Point<dim> &point_) const
    {
      WorldBuilder::Point<dim> point_tmp(value);
      point_tmp += point_.value;

      return point_tmp;
    }

    template<unsigned int dim>
    WorldBuilder::Point<dim>
    Point<dim>::operator-(const Point<dim> &point_) const
    {
      WorldBuilder::Point<dim> point_tmp(value);
      point_tmp -= point_.value;

      return point_tmp;
    }

    /**
     * access index
     */
    template<unsigned int dim>
    double &
    Point<dim>::operator[](const unsigned int index)
    {
      return value[index];
    }

    /**
     * access index
     */
    template<unsigned int dim>
    const double &
    Point<dim>::operator[](const unsigned int index) const
    {
      return value[index];
    }

    /**
     * Multiplies a Types::Point<dim> with a scalr and returns a
     * WorldBuilder::Point<dim>.
     */
    template<unsigned int dim>
    WorldBuilder::Point<dim>
    operator*(const double scalar, const Point<dim> &point)
    {
      return point.value*scalar;
    }

    template class Point<2>;
    template class Point<3>;

    /**
     * Multiplies a Types::Point<2> with a scalr and returns a
     * WorldBuilder::Point<2>.
     */
    template WorldBuilder::Point<2> operator*(const double scalar, const Point<2> &point);


    /**
     * Multiplies a Types::Point<3> with a scalr and returns a
     * WorldBuilder::Point<3>.
     */
    template WorldBuilder::Point<3> operator*(const double scalar, const Point<3> &point);
  } // namespace Types
} // namespace WorldBuilder


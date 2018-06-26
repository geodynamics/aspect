/*
  Copyright (C) 2018 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <limits>
#include <iostream>

#include <world_builder/point.h>
#include <world_builder/assert.h>

namespace WorldBuilder
{
template<>
Point<3>::Point()
:
point({0,0,0})
{}

template<>
Point<2>::Point()
:
point({0,0}),
coordinate_system(CoordinateSystem::cartesian)
{}

template<int dim>
Point<dim>::Point(const std::array<double,dim>& array, CoordinateSystem coordinate_system_)
:
point(array),
coordinate_system(coordinate_system_)
{}

template<int dim>
Point<dim>::Point(const Point<dim>& point_, CoordinateSystem coordinate_system_)
:
point(point_.get_array()),
coordinate_system(coordinate_system_)
{}


template<>
Point<2>::Point(const double x, const double y, CoordinateSystem coordinate_system_)
:
point({x,y}),
coordinate_system(coordinate_system_)
{}

template<>
Point<3>::Point(const double /*x*/, const double /*y*/, CoordinateSystem coordinate_system_)
:
point({std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}),
coordinate_system(coordinate_system_)
{
	AssertThrow(false,"Can't use the 2d constructor in 3d.");
}


template<>
Point<2>::Point(const double /*x*/, const double /*y*/, const double /*z*/, CoordinateSystem coordinate_system_)
:
point({std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN()}),
coordinate_system(coordinate_system_)
{
	AssertThrow(false,"Can't use the 3d constructor in 2d.");
}


template<>
Point<3>::Point(const double x, const double y, const double z, CoordinateSystem coordinate_system_)
:
point({x,y,z}),
coordinate_system(coordinate_system_)
{}


template<int dim>
Point<dim>::~Point()
{}

template<int dim>
Point<dim> Point<dim>::operator=(const Point<dim>& point_)
{
	return Point<dim>(point_);
}

template<int dim>
double Point<dim>::operator*(const Point<dim>& point_) const
{
	const std::array<double,dim> array = point_.get_array();
	double dot_product = 0;
	for(unsigned int i = 0; i < dim; ++i)
		dot_product += point[i] * array[i];
	return dot_product;
}


template<int dim>
Point<dim> Point<dim>::operator*(const double scalar) const
{
	// initialize the array to zero.
	std::array<double,dim> array = Point<dim>().get_array();
	for(unsigned int i = 0; i < dim; ++i)
		array[i] += point[i] * scalar;
	return Point<dim>(array);
}

template<int dim>
Point<dim> Point<dim>::operator+(const Point<dim>& point_) const
{
	std::array<double,dim> array = point_.get_array();
	for(unsigned int i = 0; i < dim; ++i)
		array[i] += point[i];
	return Point<dim>(array);
}

template<int dim>
Point<dim> Point<dim>::operator-(const Point<dim>& point_) const
{
	std::array<double,dim> array = point_.get_array();
	for(unsigned int i = 0; i < dim; ++i)
		array[i] -= point[i];
	return Point<dim>(array);
}


template<int dim>
Point<dim>&
Point<dim>::operator*=(const double scalar)
{
	for(unsigned int i = 0; i < dim; ++i)
		point[i] *= scalar;
	return *this;
}

template<int dim>
Point<dim>&
Point<dim>::operator+=(const Point<dim>& point_)
{
	for(unsigned int i = 0; i < dim; ++i)
		point[i] += point_[i];
	return *this;
}


template<int dim>
Point<dim>&
Point<dim>::operator-=(const Point<dim>& point_)
{
	for(unsigned int i = 0; i < dim; ++i)
		point[i] -= point_[i];
	return *this;
}


/**
 * access index
 */
template<int dim>
const double&
Point<dim>::operator[](const unsigned int index) const
{
	return point[index];
}


/**
 * access index
 */
template<int dim>
double&
Point<dim>::operator[](const unsigned int index)
{
	return point[index];
}

/**
 * access index
 */
template<int dim>
const double&
Point<dim>::operator()(const unsigned int index) const
{
	return point[index];
}

/**
 * access index
 */
template<int dim>
double&
Point<dim>::operator()(const unsigned int index)
{
	return point[index];
}

template<int dim>
std::array<double,dim>
Point<dim>::get_array() const
{
	return point;
}


template<int dim>
CoordinateSystem
Point<dim>::get_coordinate_system() const
{
	return coordinate_system;
}


template<int dim>
double
Point<dim>::norm() const
{
	return std::sqrt(this->norm_square());
}


template<>
double
Point<2>::norm_square() const
{
	return point[0] * point[0] + point[1] * point[1];
}

template<>
double
Point<3>::norm_square() const
{
	return point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
}


template<int dim>
Point<dim>
operator*(const double scalar, const Point<dim>& point)
{
	return point*scalar;
}


template class Point<2>;
template class Point<3>;
template Point<2> operator*(const double scalar, const Point<2>& point);
template Point<3> operator*(const double scalar, const Point<3>& point);
}

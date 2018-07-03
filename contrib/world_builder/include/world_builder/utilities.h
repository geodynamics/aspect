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

#ifndef _world_builder_utilities_h
#define _world_builder_utilities_h

#include <vector>

#include <boost/property_tree/json_parser.hpp>

#include <world_builder/point.h>
#include <world_builder/coordinate_system.h>
#include <world_builder/coordinate_systems/interface.h>

using boost::property_tree::ptree;

namespace WorldBuilder
{

  namespace CoordinateSystems
  {
    class Interface;
  }
  namespace Utilities
  {

    /**
     * Given a 2d point and a list of points which form a polygon, computes if the point
     * falls within the polygon.
     */
    bool
    polygon_contains_point(const std::vector<Point<2> > &point_list,
                           const Point<2> &point);

    /**
     * Given a 2d point and a list of points which form a polygon, compute the smallest
     * distance of the point to the polygon. The sign is negative for points outside of
     * the polygon and positive for points inside the polygon.
     */
    double
    signed_distance_to_polygon(const std::vector<Point<2> > &point_list_,
                               const Point<2> &point_);


    /*
    * A class that represents a point in a chosen coordinate system.
    */
    class NaturalCoordinate
    {
      public:
        /**
         * Constructor based on providing the geometry model as a pointer
         */
        NaturalCoordinate(const std::array<double,3> &position,
                          const ::WorldBuilder::CoordinateSystems::Interface &coordinate_system);

        /**
         * Constructor based on providing the geometry model as a pointer
         */
        NaturalCoordinate(const Point<3> &position,
                          const ::WorldBuilder::CoordinateSystems::Interface &coordinate_system);

        /**
         * Returns the coordinates in the given coordinate system, which may
         * not be Cartesian.
         */
        const std::array<double,3> &get_coordinates();

        /**
         * The coordinate that represents the 'surface' directions in the
         * chosen coordinate system.
         */
        const std::array<double,2> get_surface_coordinates() const;

        /**
         * The coordinate that represents the 'depth' direction in the chosen
         * coordinate system.
         */
        double get_depth_coordinate() const;

      private:
        /**
         * An enum which stores the the coordinate system of this natural
         * point
         */
        CoordinateSystem coordinate_system;

        /**
         * An array which stores the coordinates in the coordinates system
         */
        std::array<double,3> coordinates;
    };

    /**
     * Returns spherical coordinates of a Cartesian point. The returned array
     * is filled with radius, phi and theta (polar angle). If the dimension is
     * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
     *
     */
    std::array<double,3>
    cartesian_to_spherical_coordinates(const Point<3> &position);

    /**
     * Return the Cartesian point of a spherical position defined by radius,
     * phi and theta (polar angle). If the dimension is set to 2 theta is
     * omitted.
     */
    Point<3>
    spherical_to_cartesian_coordinates(const std::array<double,3> &scoord);

    /**
     * Returns ellipsoidal coordinates of a Cartesian point. The returned array
     * is filled with phi, theta and radius.
     *
     */
    std::array<double,3>
    cartesian_to_ellipsoidal_coordinates(const Point<3> &position,
                                         const double semi_major_axis_a,
                                         const double eccentricity);

    /**
     * Return the Cartesian point of a ellipsoidal position defined by phi,
     * phi and radius.
     */
    Point<3>
    ellipsoidal_to_cartesian_coordinates(const std::array<double,3> &phi_theta_d,
                                         const double semi_major_axis_a,
                                         const double eccentricity);

    /**
     * A function that takes a string representation of the name of a
     * coordinate system (as represented by the CoordinateSystem enum)
     * and returns the corresponding value.
     */
    CoordinateSystem
    string_to_coordinate_system (const std::string &);


    /**
     * Convert point to array
     */
    template<int dim>
    const std::array<double,dim> convert_point_to_array(Point<dim> &point);

    /**
     * Converts a string to a double
     */
    double
    string_to_double(const std::string &string);

    /**
     * Converts a string to a int
     */
    double
    string_to_int(const std::string &string);


    /**
     * Converts a string to a unsigned int
     */
    double
    string_to_unsigned_int(const std::string &string);

    /**
     * Returns a value from the property tree and asserts with
     * the path and value in the error message when the value
     * was not present.
     */
    std::string
    get_from_ptree(const ptree &tree,
                   const std::string &path,
                   const std::string &key,
                   const std::string &path_separator = ".");
  }
}


#endif

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

#ifndef _world_builder_utilities_h
#define _world_builder_utilities_h

#include <vector>
#include <map>

#include <world_builder/point.h>
#include <world_builder/coordinate_system.h>
#include <world_builder/coordinate_systems/interface.h>


namespace WorldBuilder
{

  namespace CoordinateSystems
  {
    class Interface;
  }
  namespace Utilities
  {

    // define pi. Directly defining it seems to be
    // the safest option.
    constexpr double const_pi = 3.141592653589793238462643383279502884;

    /**
     * Given a 2d point and a list of points which form a polygon, computes if
     * the point falls within the polygon. For spherical coordinates it will
     * return true if the point or the point where the longitude is shifted
     * by 2 * PI is inside on of the polygons. It calls
     * polygon_contains_point_implementation to do the real work.
     */
    bool
    polygon_contains_point(const std::vector<Point<2> > &point_list,
                           const Point<2> &point);

    /**
     * Given a 2d point and a list of points which form a polygon, computes if the point
     * falls within the polygon.
     */
    bool
    polygon_contains_point_implementation(const std::vector<Point<2> > &point_list,
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

        /**
         * get the coordinate system type of this coordinate.
         */
        CoordinateSystem get_coordinate_system() const;

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
    const std::array<double,dim> convert_point_to_array(const Point<dim> &point);

    /**
     * Converts a string to a double
     */
    double
    string_to_double(const std::string &string);

    /**
     * Converts a string to a int
     */
    int
    string_to_int(const std::string &string);


    /**
     * Converts a string to a unsigned int
     */
    unsigned int
    string_to_unsigned_int(const std::string &string);


    /**
     * Cross product between two 3d points.
     */
    Point<3> cross_product(const Point<3> &a, const Point<3> &b);

    /**
     * Computes the distance of a point to a curved plane.
     * TODO: add more info on how this works/is implemented.
     * \param point This is the cartesian point of which we want to know the
     * distance to the curved planes
     * \param reference_point This is a 2d point in natural coordinates at the
     * surface which the curved planes dip towards. Natural coordinates are in
     * cartesian (x,y,z) in meters and in spherical radius in meters and longitude
     * and latitude in radians.
     * \param point_list This is a vector of 2d Points in natural coordinates at the
     * surface which define the line along the surface at which the curved planes
     * start. Natural coordinates are in cartesian (x,y,z) in meters and in spherical
     * radius in meters and longitude and latitude in radians.
     * \param plane_segment_lengths This is a vector of vectors of doubles. It contains
     * the length of every segment at point in the point_list (in the same order as
     * the point_list.
     * \param plane_segment_angles This is a vector of vectors of 2d points. It contains
     * the begin and end angle of every segment at point in the point_list (in the same
     * order as the point_list.
     * \param start_depth This value contains the depth at which the plane starts. This
     * means that the start_depth effectively becomes the surface for this slab.
     * \param coordinate_system This is a reference to the coordinate system of the
     * World Builder. This is used to convert cartesian to natural coordinates and back.
     * \param only_positive This value deterines whether only the the part below the
     * plane should count as distance or both sides of the plane. It is called only_positive
     * because the area below the plane, the distance is positve, and above the plane the
     * distance is negative.
     * \param global_x_list This is a list of one dimensional coorindates, with zero or the
     * amount of coordinates entries, used for interpolation. An empty list is interpretated
     * as a list filled with {0,1,2,...,number of coordinates}. Filling this list with other
     * values changes the returned section fraction. It allows for, for example, adding
     * extra coordinates automatically, and still reference the user provided coordinates by
     * the original number. Note that no whole numbers may be skiped. So for a list of 4 points,
     * {0,0.5,1,2} is allowed, but {0,2,3,4} is not.
     */
    std::map<std::string,double> distance_point_from_curved_planes(const Point<3> &point,
                                                                   const Point<2> &reference_point,
                                                                   const std::vector<Point<2> > &point_list,
                                                                   const std::vector<std::vector<double> > &plane_segment_lengths,
                                                                   const std::vector<std::vector<Point<2> > > &plane_segment_angles,
                                                                   const double start_depth,
                                                                   const std::unique_ptr<CoordinateSystems::Interface> &coordinate_system,
                                                                   const bool only_positive,
                                                                   std::vector<double> global_x_list = {});

    /**
     * Class for linear and monotone spline interpolation
     */
    class interpolation
    {
      public:
        /**
         * Initialize the spline.
         *
         * @param x X coordinates of interpolation points.
         * @param y Values in the interpolation points.
         * @param monotone_spline Whether to construct a monotone cubic spline or just do linear interpolation.
         */
        void set_points(const std::vector<double> &x,
                        const std::vector<double> &y,
                        const bool monotone_spline = false);
        /**
         * Evaluate at point @p x.
         */
        double operator() (double x) const;

      private:
        /**
         * x coordinates of points
         */
        std::vector<double> m_x;

        /**
         * interpolation parameters
         * \[
         * f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
         * \]
         */
        std::vector<double> m_a, m_b, m_c, m_y;
    };


    /**
     * Ensure angle is between 0 and 360 degrees
     */
    double wrap_angle(const double angle);


    /**
     * Transorm a rotation matrix into euler angles
     */
    std::array<double,3>
    euler_angles_from_rotation_matrix(const std::array<std::array<double,3>,3> &rotation_matrix);

    /**
     * Transform euler angles into a rotation matrix
     */
    std::array<std::array<double,3>,3>
    euler_angles_to_rotation_matrix(double phi1, double theta, double phi2);
  }
}


#endif

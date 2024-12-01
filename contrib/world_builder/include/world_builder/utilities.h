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

#ifndef WORLD_BUILDER_UTILITIES_H
#define WORLD_BUILDER_UTILITIES_H


#include "world_builder/nan.h"
#include "world_builder/coordinate_systems/interface.h"
#include "world_builder/objects/natural_coordinate.h"
#include "world_builder/objects/bezier_curve.h"
#include <iostream>


namespace WorldBuilder
{

  namespace CoordinateSystems
  {
    class Interface;
  } // namespace CoordinateSystems
  namespace Utilities
  {

    /**
     * provide a short way to test if two doubles are equal.
     * Based on https://stackoverflow.com/a/4010279.
     * Removed a==b test since it triggers warnings. If used in
     * performance critical parts where this could matter, a fast
     * version could be added.
     */
    inline bool approx(double a, double b, double error_factor=1e4)
    {
      return std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*
             error_factor;
    }

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
     * Given a 2d point, a semi-major axis, and an eccentricity, computes where
     * the point falls within the ellipse. If the fraction is larger than 1, the
     * point is outside the ellipse.
     */
    double
    fraction_from_ellipse_center (const Point<2> &ellipse_center,
                                  const double semi_major_axis,
                                  const double eccentricity,
                                  const double rotation_angle,
                                  const Point<2> &point);


    /**
     * Given a 2d point and a list of points which form a polygon, compute the smallest
     * distance of the point to the polygon. The sign is negative for points outside of
     * the polygon and positive for points inside the polygon.
     */
    double
    signed_distance_to_polygon(const std::vector<Point<2> > &point_list_,
                               const Point<2> &point_);


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
    string_to_coordinate_system (const std::string & /*coordinate_system*/);


    /**
     * Convert point to array
     */
    template<unsigned int dim>
    std::array<double,dim> convert_point_to_array(const Point<dim> &point);

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
     * Enum class for interolation type
     */
    enum class InterpolationType
    {
      None,
      Linear,
      MonotoneSpline,
      ContinuousMonotoneSpline,
      Invalid,
    };

    /**
     * Class for linear and monotone spline interpolation
     */
    class interpolation
    {
      public:
        /**
         * Initialize the spline. This function assumes that all y points are spaced 1 in the x direction.
         *
         * @param y Values in the interpolation points.
         */
        void set_points(const std::vector<double> &y);


        /**
         * Evaluate at point @p x.
         */
        inline
        double operator() (const double x) const
        {
          if (x >= 0. && x <= static_cast<double>(mx_size_min))
            {
              const size_t idx = static_cast<size_t>(x);
              const double h = x-static_cast<double>(idx);
              return ((m[idx][0]*h + m[idx][1])*h + m[idx][2])*h + m[idx][3];
            }
          const size_t idx = std::min(static_cast<size_t>(std::max( static_cast<int>(x), static_cast<int>(0))),mx_size_min);
          const double h = x-static_cast<double>(idx);
          return (m[idx][1]*h + m[idx][2])*h + m[idx][3];
        }


        inline
        double operator() (const double x, const size_t idx, const double h) const
        {
          return (x >= 0. && x <= static_cast<double>(mx_size_min))
                 ?
                 ((m[idx][0]*h + m[idx][1])*h + m[idx][2])*h + m[idx][3]
                 :
                 (m[idx][1]*h + m[idx][2])*h + m[idx][3];
        }


        /**
         * Evaluate at point @p x. assumes x is between 0 and mx_size_min.
         * assume size_t idx = (size_t)x and h = x-idx.
         */
        inline
        double value_inside (const size_t idx, const double h) const
        {
          WBAssert(idx <= mx_size_min, "Internal error: using value_inside outside the range of 0 to " << mx_size_min << ", but value was outside of this range: " << idx << ".");
          WBAssert(h >= 0 && h <= 1., "Internal error: using value_inside outside the range of 0 to " << mx_size_min << ", but value was outside of this range: " << h << ".");
          return ((m[idx][0]*h + m[idx][1])*h + m[idx][2])*h + m[idx][3];
        }


        /**
         * Evaluate at point @p x. assumes x is between 0 and mx_size_min.
         * assume size_t idx = (size_t)x and h = x-idx.
         */
        inline
        double value_outside (const size_t idx, const double h) const
        {
          WBAssert(idx <= mx_size_min, "Internal error: using value_inside outside the range of 0 to " << mx_size_min << ", but value was outside of this range: " << idx << ".");
          WBAssert(!(static_cast<double>(idx) + h >= 0 && static_cast<double>(idx) + h <= 1.), "Internal error: using value_inside outside the range of 0 to " << mx_size_min << ", but value was outside of this range: " << static_cast<double>(idx) + h << " (h=" << h << ", idx = " << idx << ").");
          return (m[idx][1]*h + m[idx][2])*h + m[idx][3];
        }


        /**
         * number of x coordinates of points
         */
        size_t mx_size_min;

        /**
         * interpolation parameters
         * \[
         * f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
         * \]
         */
        std::vector<std::array<double,4>> m; //m_a, m_b, m_c, m_y;

      private:
    };

    /**
     * A struct that is used to hold the return values of the function
     * distance_point_from_curved_planes(). See there for a documentation
     * of the meaning of the member variables. The variables are describing
     * a where a point is with respect to a plane/surface. The surface is
     * meshed by a grid. The axis parallel to the surface are formed by
     * sections, and the axis perpendicuar to the surface are formed segments.
     * Both sections and elements represent a whole cell in the grid, which is
     * an integer. This structure also provides the fraction in each direction
     * the closest point on the plane is along these two axes (sections and
     * segments). These variables are called fractions.
     *
     * This structure furthermore provides the distance the provided point is
     * from the closest point on the surface and the average angle. The variable
     * local_thickness will not be automatically filled by the distance_point_from_curved_planes
     * function.
     */
    struct PointDistanceFromCurvedPlanes
    {
      /**
       * Constructor
       */
      PointDistanceFromCurvedPlanes(CoordinateSystem coordinate_system)
        :
        distance_from_plane(NaN::DSNAN),
        distance_along_plane(NaN::DSNAN),
        fraction_of_section(NaN::DSNAN),
        fraction_of_segment(NaN::DSNAN),
        section(NaN::ISNAN),
        segment(NaN::ISNAN),
        average_angle(NaN::DSNAN),
        depth_reference_surface(NaN::DSNAN),
        closest_trench_point(Point<3>(coordinate_system))
      {}

      /**
       * The shortest distance between point and plane.
       */
      double distance_from_plane;

      /**
       * The distance between the start of the first segment (usually at
       * the surface) to the provided point following the (curved) plane.
       */
      double distance_along_plane;

      /**
       * The fraction of the section that lies before the projected point
       * on the plane (when looking from the start point of the section).
       */
      double fraction_of_section;

      /**
       * The fraction of the segment that lies before the projected point
       * on the plane (when looking from the start point of the segment).
       */
      double fraction_of_segment;

      /**
       * The number of the section that is closest to the point
       */
      size_t section;

      /**
       * The number of the segment that is closest to the point.
       */
      size_t segment;

      /**
       * The average dip angle of the plane at the location where the
       * point is projected onto the plane.
       */
      double average_angle;

      /**
       * The depth of the closest point on reference surface.
       */
      double depth_reference_surface;

      /**
       * The closest point on the trench line in cartesian coordinates.
       */
      Point<3> closest_trench_point;
    };

    /**
     * Computes the distance of a point to a curved plane.
     * TODO: add more info on how this works/is implemented.
     * \param check_point This is the cartesian point of which we want to know the
     * distance to the curved planes
     * \param check_point_natural the check_point in the natural coordinates of the
     * current coordinate system.
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
     * \param start_radius This value contains the radius or height from bottom of the box
     * at which the plane starts. This means that the start_radius effectively becomes the
     * surface for this slab.
     * \param coordinate_system This is a reference to the coordinate system of the
     * World Builder. This is used to convert cartesian to natural coordinates and back.
     * \param only_positive This value determines whether only the part below the
     * plane should count as distance or both sides of the plane. It is called only_positive
     * because the area below the plane, the distance is positive, and above the plane the
     * distance is negative.
     * \param interpolation_type This value determines what interpolation type should be used
     * when determining the location with respect to the curved plane.
     * \param spline_x the spline representing the x coordinate.
     * \param spline_y the spline representing the y coordinate.
     * \param global_x_list This is a list of one dimensional coorindates, with zero or the
     * amount of coordinates entries, used for interpolation. An empty list is interpreted
     * as a list filled with {0,1,2,...,number of coordinates}. Filling this list with other
     * values changes the returned section fraction. It allows for, for example, adding
     * extra coordinates automatically, and still reference the user provided coordinates by
     * the original number. Note that no whole numbers may be skipped. So for a list of 4 points,
     * {0,0.5,1,2} is allowed, but {0,2,3,4} is not.
     *
     * The function returns a struct that contains which segment and section of the curved
     * planes the point is closest to, what fraction of those segment and section lies before
     * the point (looking from the start of segment/section), the distance
     * of the point from the plane and the distance of the point along the plane,
     * and the average angle of the closest segment/section.
     */
    PointDistanceFromCurvedPlanes distance_point_from_curved_planes(const Point<3> &check_point,
                                                                    const Objects::NaturalCoordinate &check_point_natural,
                                                                    const Point<2> &reference_point,
                                                                    const std::vector<Point<2> > &point_list,
                                                                    const std::vector<std::vector<double> > &plane_segment_lengths,
                                                                    const std::vector<std::vector<Point<2> > > &plane_segment_angles,
                                                                    const double start_radius,
                                                                    const std::unique_ptr<CoordinateSystems::Interface> &coordinate_system,
                                                                    const bool only_positive,
                                                                    const Objects::BezierCurve &bezier_curve);



    /**
     * Ensure angle is between 0 and 360 degrees
     */
    double wrap_angle(const double angle);

    /**
     * Interpolate between two angles (angle1 and angle2),
     * with fraction defining the weighting between the two,
     * taking into account we might cross over from 360 to 0 degrees.
     */
    double interpolate_angle_across_zero(const double angle_1,
                                         const double angle_2,
                                         const double fraction);

    /**
     * Transform a rotation matrix into euler angles
     */
    std::array<double,3>
    euler_angles_from_rotation_matrix(const std::array<std::array<double,3>,3> &rotation_matrix);

    /**
     * Transform euler angles into a rotation matrix
     */
    std::array<std::array<double,3>,3>
    euler_angles_to_rotation_matrix(double phi1, double theta, double phi2);

    /**
     * Read a file and distribute the content over all MPI processes.
     * If WB_WITH_MPI is not defined, this function will just read the file.
     *
     * @param filename The name of the file to read.
     * @return The content of the file.
    */
    std::string
    read_and_distribute_file_content(const std::string &filename);

    /**
     * Calculate the distance of a point from a mid oceanic ridge, and also calculate
     * the spreading velocity of the ridge at this point.
     * TODO: make the spreading velocity spatially/temporally variable
     *
     * @param mid_oceanic_ridges The coordinates of the mid oceanic ridges
     * @param mid_oceanic_spreading_velocities The spreading rate of the mid oceanic ridges at each ridge coordinate
     * @param coordinate_system The coordinate system
     * @param position_in_natural_coordinates_at_min_depth the current position in natural_coordinates
     * @param subducting_plate_velocities the subducting plate velocities, currently this is only an optional parameter
     * that can be used in the mass conserving subducting plate temperature model. This parameter allows the user to
     * track the effect of ridge migration on the slab thermal structure
     * @param ridge_migration_times the times that the corresponding section of the ridge has been moving, in years. This
     * is used in combination with subducting_plate_velocities, and mid_oceanic_spreading_velocities to compute the distance
     * that the spreading center has migrated. This vector is obtained from the input parameter "spreading velocity" in the
     * mass conserving model when "spreading velocity" has the form: [ [t1,[[v11, v12, ...]], [t2,[[v21, v22, ...]], ... ].
     * where tn is the time that ridge section n has been moving.
     * @return The content of the file.
    */
    std::vector<double>
    calculate_ridge_distance_and_spreading(std::vector<std::vector<Point<2>>> mid_oceanic_ridges,
                                           std::vector<std::vector<double>> mid_oceanic_spreading_velocities,
                                           const std::unique_ptr<WorldBuilder::CoordinateSystems::Interface> &coordinate_system,
                                           const Objects::NaturalCoordinate &position_in_natural_coordinates_at_min_depth,
                                           const std::vector<std::vector<double>> &subducting_plate_velocities,
                                           const std::vector<double> &ridge_migration_times);

    // todo_effective
    /**
     * Calculate the effective plate ages of a point on the slab surface, and also calculates
     * the effective trench ages at the start of subduction.
     * @param ridge_parameters The distance and spreading velocity relative to a mid ocean ridge
     * @param distance_along_plane The distance along the slab surface plane
     * @return The effective plate age and the trench age
    */
    std::vector<double>
    calculate_effective_trench_and_plate_ages(std::vector<double> ridge_parameters, double distance_along_plane);

    /*
     * Returns the result of the multiplication of two 3*3 matrix,
     * used in applying the random uniform distribution rotation matrix
     * to a given orientation (rotation matrix)
     */
    std::array<std::array<double,3>,3>
    multiply_3x3_matrices(const std::array<std::array<double,3>,3> mat1, std::array<std::array<double,3>,3> const mat2);

  } // namespace Utilities
} // namespace WorldBuilder


#endif

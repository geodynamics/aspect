

#include <aspect/world_builder/utilities.h>
#include <aspect/utilities.h>
#include <aspect/world_builder/coordinate_systems.h>
#include <aspect/world_builder/coordinate_system/interface.h>
#include <deal.II/base/exceptions.h>

namespace aspect
{
  namespace WorldBuilder
  {
    namespace Utilities
    {
      bool
      polygon_contains_point(const std::vector<std::array<double,2> > &point_list_,
                             const std::array<double,2> &point_)
      {
    	  std::vector<Point<2> > point_list(point_list_.size());
    	  for(unsigned int i = 0; i < point_list_.size(); ++i)
    		  point_list[i] = aspect::Utilities::convert_array_to_point<2>(point_list_[i]);

    	  const Point<2> point = aspect::Utilities::convert_array_to_point<2>(point_);

        /**
         * This code has been based on http://geomalgorithms.com/a03-_inclusion.html,
         * and therefore requires the following copyright notice:
         *
         * Copyright 2000 softSurfer, 2012 Dan Sunday
         * This code may be freely used and modified for any purpose
         * providing that this copyright notice is included with it.
         * SoftSurfer makes no warranty for this code, and cannot be held
         * liable for any real or imagined damage resulting from its use.
         * Users of this code must verify correctness for their application.
         *
         * The main functional difference between the original code and this
         * code is that all the boundaries are considered to be inside the
         * polygon. One should of course realize that with floating point
         * arithmetic no guarantees can be made for the borders, but for
         * exact arithmetic this algorithm would work (also see polygon
         * in point test).
         */
        int pointNo = point_list.size();
        int    wn = 0;    // the  winding number counter
        int   j=pointNo-1;

        // loop through all edges of the polygon
        for (int i=0; i<pointNo; i++)
          {
            // edge from V[i] to  V[i+1]
            if (point_list[j][1] <= point[1])
              {
                // start y <= P.y
                if (point_list[i][1] >= point[1])      // an upward crossing
                  {
                    const double is_left = (point_list[i][0] - point_list[j][0]) * (point[1] - point_list[j][1])
                                           - (point[0] -  point_list[j][0]) * (point_list[i][1] - point_list[j][1]);

                    if ( is_left > 0 && point_list[i][1] > point[1])
                      {
                        // P left of  edge
                        ++wn;            // have  a valid up intersect
                      }
                    else if ( is_left == 0)
                      {
                        // The point is exactly on the infinite line.
                        // determine if it is on the segment
                        const double dot_product = (point - point_list[j])*(point_list[i] - point_list[j]);

                        if (dot_product >= 0)
                          {
                            const double squaredlength = (point_list[i] - point_list[j]).norm_square();

                            if (dot_product <= squaredlength)
                              {
                                return true;
                              }
                          }
                      }
                  }
              }
            else
              {
                // start y > P.y (no test needed)
                if (point_list[i][1]  <= point[1])     // a downward crossing
                  {
                    const double is_left = (point_list[i][0] - point_list[j][0]) * (point[1] - point_list[j][1])
                                           - (point[0] -  point_list[j][0]) * (point_list[i][1] - point_list[j][1]);

                    if ( is_left < 0)
                      {
                        // P right of  edge
                        --wn;            // have  a valid down intersect
                      }
                    else if ( is_left == 0)
                      {
                        // This code is to make sure that the boundaries are included in the polygon.
                        // The point is exactly on the infinite line.
                        // determine if it is on the segment
                        const double dot_product = (point - point_list[j])*(point_list[i] - point_list[j]);

                        if (dot_product >= 0)
                          {
                            const double squaredlength = (point_list[i] - point_list[j]).norm_square();

                            if (dot_product <= squaredlength)
                              {
                                return true;
                              }
                          }
                      }
                  }
              }
            j=i;
          }

        return (wn != 0);
      }

      double
      signed_distance_to_polygon(const std::vector<std::array<double,2> > &point_list_,
                                 const std::array<double,2> &point_)
      {
        // If the point lies outside polygon, we give it a negative sign,
        // inside a positive sign.
        const double sign = polygon_contains_point(point_list_, point_) ? 1.0 : -1.0;

  	  std::vector<Point<2> > point_list(point_list_.size());
  	  for(unsigned int i = 0; i < point_list_.size(); ++i)
  		  point_list[i] = aspect::Utilities::convert_array_to_point<2>(point_list_[i]);

  	  const Point<2> point = aspect::Utilities::convert_array_to_point<2>(point_);

        /**
         * This code is based on http://geomalgorithms.com/a02-_lines.html#Distance-to-Infinite-Line,
         * and therefore requires the following copyright notice:
         *
         * Copyright 2000 softSurfer, 2012 Dan Sunday
         * This code may be freely used and modified for any purpose
         * providing that this copyright notice is included with it.
         * SoftSurfer makes no warranty for this code, and cannot be held
         * liable for any real or imagined damage resulting from its use.
         * Users of this code must verify correctness for their application.
         *
         */

        const unsigned int n_poly_points = point_list.size();
        AssertThrow(n_poly_points >= 3, ExcMessage("Not enough polygon points were specified."));

        // Initialize a vector of distances for each point of the polygon with a very large distance
        std::vector<double> distances(n_poly_points, 1e23);

        // Create another polygon but with all points shifted 1 position to the right
        std::vector<Point<2> > shifted_point_list(n_poly_points);
        shifted_point_list[0] = point_list[n_poly_points-1];

        for (unsigned int i = 0; i < n_poly_points-1; ++i)
          shifted_point_list[i+1] = point_list[i];

        for (unsigned int i = 0; i < n_poly_points; ++i)
          {
            // Create vector along the polygon line segment
            Tensor<1,2> vector_segment = shifted_point_list[i] - point_list[i];
            // Create vector from point to the second segment point
            Tensor<1,2> vector_point_segment = point - point_list[i];

            // Compute dot products to get angles
            const double c1 = vector_point_segment * vector_segment;
            const double c2 = vector_segment * vector_segment;

            // point lies closer to not-shifted polygon point, but perpendicular base line lies outside segment
            if (c1 <= 0.0)
              distances[i] = (Tensor<1,2> (point_list[i] - point)).norm();
            // point lies closer to shifted polygon point, but perpendicular base line lies outside segment
            else if (c2 <= c1)
              distances[i] = (Tensor<1,2> (shifted_point_list[i] - point)).norm();
            // perpendicular base line lies on segment
            else
              {
                const Point<2> point_on_segment = point_list[i] + (c1/c2) * vector_segment;
                distances[i] = (Tensor<1,2> (point - point_on_segment)).norm();
              }
          }

        // Return the minimum of the distances of the point to all polygon segments
        return *std::min_element(distances.begin(),distances.end()) * sign;
      }

      //typedef class NaturalCoordinate;
      //typedef class Coordinates;

      NaturalCoordinate::NaturalCoordinate(const std::array<double,3> &position,
                                           const CoordinateSystem::Interface &coordinate_system_)
      {
        coordinate_system = coordinate_system_.natural_coordinate_system();
        coordinates = coordinate_system_.cartesian_to_natural_coordinates(position);
      }

      std::array<double,3> &NaturalCoordinate::get_coordinates()
      {
        return coordinates;
      }


      /*std::array<double,1> NaturalCoordinate::get_surface_coordinates() const
      {
        std::array<double,1> coordinate;

        switch (coordinate_system)
          {
            case Coordinates::CoordinateSystem::cartesian:
              coordinate[0] = coordinates[0];
              break;

            case Coordinates::CoordinateSystem::spherical:
              coordinate[0] = coordinates[1];
              break;

            case Coordinates::CoordinateSystem::ellipsoidal:
              coordinate[0] = coordinates[1];
              break;

            default:
              coordinate[0] = 0;
              Assert (false, ExcNotImplemented());
              break;
          }

        return coordinate;
      }*/


      std::array<double,2> NaturalCoordinate::get_surface_coordinates() const
      {
        std::array<double,2> coordinate;

        switch (coordinate_system)
          {
            case Coordinates::CoordinateSystem::cartesian:
              coordinate[0] = coordinates[0];
              coordinate[1] = coordinates[1];
              break;

            case Coordinates::CoordinateSystem::spherical:
              coordinate[0] = coordinates[1];
              coordinate[1] = coordinates[2];
              break;

            case Coordinates::CoordinateSystem::ellipsoidal:
              coordinate[0] = coordinates[1];
              coordinate[1] = coordinates[2];
              break;

            default:
              Assert (false, ExcNotImplemented());
          }

        return coordinate;
      }


      double NaturalCoordinate::get_depth_coordinate() const
      {
        switch (coordinate_system)
          {
            case Coordinates::CoordinateSystem::cartesian:
              return coordinates[2];

            case Coordinates::CoordinateSystem::spherical:
              return coordinates[0];

            case Coordinates::CoordinateSystem::ellipsoidal:
              return coordinates[0];

            default:
              Assert (false, ExcNotImplemented());
          }

        return 0;
      }


      std::array<double,3>
      cartesian_to_spherical_coordinates(const Point<3> &position)
      {
        std::array<double,3> scoord;

        scoord[0] = position.norm(); // R
        scoord[1] = std::atan2(position(1),position(0)); // Phi
        if (scoord[1] < 0.0)
          scoord[1] += 2.0*numbers::PI; // correct phi to [0,2*pi]

        if (scoord[0] > std::numeric_limits<double>::min())
          scoord[2] = std::acos(position(2)/scoord[0]);
        else
          scoord[2] = 0.0;

        return scoord;
      }

      Point<3>
      spherical_to_cartesian_coordinates(const std::array<double,3> &scoord)
      {
        Point<3> ccoord;

        ccoord[0] = scoord[0] * std::sin(scoord[2]) * std::cos(scoord[1]); // X
        ccoord[1] = scoord[0] * std::sin(scoord[2]) * std::sin(scoord[1]); // Y
        ccoord[2] = scoord[0] * std::cos(scoord[2]); // Z


        return ccoord;
      }

      std::array<double,3>
      cartesian_to_ellipsoidal_coordinates(const Point<3> &x,
                                           const double semi_major_axis_a,
                                           const double eccentricity)
      {
        const double R    = semi_major_axis_a;
        const double b      = std::sqrt(R * R * (1 - eccentricity * eccentricity));
        const double ep     = std::sqrt((R * R - b * b) / (b * b));
        const double p      = std::sqrt(x(0) * x(0) + x(1) * x(1));
        const double th     = std::atan2(R * x(2), b * p);
        const double phi    = std::atan2(x(1), x(0));
        const double theta  = std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th),3),
                                         (p - (eccentricity * eccentricity * R  * std::pow(std::cos(th),3))));
        const double R_bar = R / (std::sqrt(1 - eccentricity * eccentricity * std::sin(theta) * std::sin(theta)));
        const double R_plus_d = p / std::cos(theta);

        std::array<double,3> phi_theta_d;
        phi_theta_d[0] = phi;

        phi_theta_d[1] = theta;
        phi_theta_d[2] = R_plus_d - R_bar;
        return phi_theta_d;
      }

      Point<3>
      ellipsoidal_to_cartesian_coordinates(const std::array<double,3> &phi_theta_d,
                                           const double semi_major_axis_a,
                                           const double eccentricity)
      {
        const double phi   = phi_theta_d[0];
        const double theta = phi_theta_d[1];
        const double d     = phi_theta_d[2];

        const double R_bar = semi_major_axis_a / std::sqrt(1 - (eccentricity * eccentricity *
                                                                std::sin(theta) * std::sin(theta)));

        return Point<3> ((R_bar + d) * std::cos(phi) * std::cos(theta),
                         (R_bar + d) * std::sin(phi) * std::cos(theta),
                         ((1 - eccentricity * eccentricity) * R_bar + d) * std::sin(theta));

      }



      Coordinates::CoordinateSystem
      string_to_coordinate_system(const std::string &coordinate_system)
      {
        if (coordinate_system == "cartesian")
          return Coordinates::CoordinateSystem::cartesian;
        else if (coordinate_system == "spherical")
          return Coordinates::CoordinateSystem::spherical;
        else if (coordinate_system == "ellipsoidal")
          return Coordinates::CoordinateSystem::ellipsoidal;
        else
          AssertThrow(false, ExcNotImplemented());

        return Coordinates::invalid;
      }
    }
  }
}



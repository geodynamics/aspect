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
#include "world_builder/assert.h"
#include "world_builder/nan.h"
#include "world_builder/objects/bezier_curve.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace WorldBuilder;

namespace WorldBuilder
{
  namespace Objects
  {
    BezierCurve::BezierCurve(const std::vector<Point<2> > &p, const std::vector<double> &angle_constrains_input)
    {
      points = p;
      const size_t n_points = p.size();
      control_points.resize(n_points-1, {{p[0],p[0]}});
      lengths.resize(n_points-1,NaN::DSNAN);
      angles.resize(n_points,NaN::DSNAN);
      std::vector<double> angle_constrains = angle_constrains_input;
      angle_constrains.resize(n_points,NaN::DQNAN);

      // if no angle is provided, compute the angle as the average angle between the previous and next point.
      // The first angle points at the second point and the last angle points at the second to last point.
      // The check points are set at a distance of 1/10th the line length from the point in the direction of the angle.
      if (std::isnan(angle_constrains[0]))
        {
          Point<2> P1P2 = points[1]-points[0];
          angles[0] = atan2(P1P2[1],P1P2[0]);
        }
      else
        {
          angles[0] = angle_constrains[0];
        }

      for (size_t p_i = 1; p_i < n_points-1; ++p_i)
        {
          // first determine the angle
          if (std::isnan(angle_constrains[p_i]))
            {
              // get the average angle
              const Point<2> P1P2 = points[p_i-1]-points[p_i];
              const Point<2> P3P2 = points[p_i+1]-points[p_i];

              const double angle_p1p2 = atan2(P1P2[1],P1P2[0]);
              const double angle_p3p1 = atan2(P3P2[1],P3P2[0]);
              const double average_angle = (angle_p1p2 + angle_p3p1)*0.5;
              angles[p_i] = average_angle;
              angles[p_i] -= Consts::PI*0.5;
            }
          else
            {
              angles[p_i] = angle_constrains[p_i];
            }

        }

      if (std::isnan(angle_constrains[n_points-1]))
        {
          Point<2> P1P2 = points[n_points-2]-points[n_points-1];
          angles[n_points-1] =  atan2(P1P2[1],P1P2[0]);
        }
      else
        {
          angles[n_points-1] = angle_constrains[n_points-1];
        }

      if (points.size() > 2)
        {
          // next determine the location of the control points
          // the location of the control point is 1/10th p1p2 distance in the direction of the angle.
          // make sure the angle is pointing away from the next point, e.g.
          // the check point is on the other side of the of the line p1p2 compared to p3.
          const double fraction_of_length = 0.2;
          {
            const Point<2> &p1 = points[0];
            const Point<2> &p2 = points[1];
            const Point<2> &p3 = points[2];
            const double length = (points[0]-points[1]).norm(); // can be squared
            control_points[0][0][0] = cos(angles[0])*length*fraction_of_length+p1[0];
            control_points[0][0][1] = sin(angles[0])*length*fraction_of_length+p1[1];
            control_points[0][1][0] = cos(angles[1])*length*fraction_of_length+p2[0];
            control_points[0][1][1] = sin(angles[1])*length*fraction_of_length+p2[1];
            {
              const bool side_of_line_1 =  (p1[0] - p2[0]) * (control_points[0][1][1] - p1[1])
                                           - (p1[1] - p2[1]) * (control_points[0][1][0] - p1[0])
                                           < 0;
              const bool side_of_line_2 =  (p1[0] - p2[0]) * (p3[1] - p1[1])
                                           - (p1[1] - p2[1]) * (p3[0] - p1[0])
                                           < 0;
              if (side_of_line_1 == side_of_line_2)
                {
                  // use a 180 degree rotated angle to create this control_point
                  control_points[0][1][0] = cos(angles[1]+Consts::PI)*length*fraction_of_length+p2[0];
                  control_points[0][1][1] = sin(angles[1]+Consts::PI)*length*fraction_of_length+p2[1];
                }
            }
          }

          for (size_t p_i = 1; p_i < n_points-1; ++p_i)
            {
              const Point<2> &p1 = points[p_i];
              const Point<2> &p2 = points[p_i+1];
              const double length = (points[p_i]-points[p_i+1]).norm(); // can be squared
              control_points[p_i][0][0] = cos(angles[p_i])*length*fraction_of_length+p1[0];
              control_points[p_i][0][1] = sin(angles[p_i])*length*fraction_of_length+p1[1];

              {
                const bool side_of_line_1 =  (p1[0] - p2[0]) * (control_points[p_i-1][1][1] - p1[1])
                                             - (p1[1] - p2[1]) * (control_points[p_i-1][1][0] - p1[0])
                                             < 0;
                const bool side_of_line_2 =  (p1[0] - p2[0]) * (control_points[p_i][0][1] - p1[1])
                                             - (p1[1] - p2[1]) * (control_points[p_i][0][0] - p1[0])
                                             < 0;
                if (side_of_line_1 == side_of_line_2)
                  {
                    // use a 180 degree rotated angle to create this control_point
                    control_points[p_i][0][0] = cos(angles[p_i]+Consts::PI)*length*fraction_of_length+p1[0];
                    control_points[p_i][0][1] = sin(angles[p_i]+Consts::PI)*length*fraction_of_length+p1[1];
                  }
              }

              control_points[p_i][1][0] = cos(angles[p_i+1])*length*fraction_of_length+points[p_i+1][0];
              control_points[p_i][1][1] = sin(angles[p_i+1])*length*fraction_of_length+points[p_i+1][1];

              if (p_i+1 < n_points-1)
                {
                  const Point<2> &p3 = points[p_i+2];
                  const bool side_of_line_1 =  (p1[0] - p2[0]) * (control_points[p_i][1][1] - p1[1])
                                               - (p1[1] - p2[1]) * (control_points[p_i][1][0] - p1[0])
                                               < 0;
                  const bool side_of_line_2 =  (p1[0] - p2[0]) * (p3[1] - p1[1])
                                               - (p1[1] - p2[1]) * (p3[0] - p1[0])
                                               < 0;
                  if (side_of_line_1 == side_of_line_2)
                    {
                      // use a 180 degree rotated angle to create this control_point
                      control_points[p_i][1][0] = cos(angles[p_i+1]+Consts::PI)*length*fraction_of_length+p2[0];
                      control_points[p_i][1][1] = sin(angles[p_i+1]+Consts::PI)*length*fraction_of_length+p2[1];
                    }
                }
            }
        }
    }


    Point<2>
    BezierCurve::operator()(const size_t i, const double t) const
    {
      WBAssert(i < points.size()-1 && i < control_points.size() ,
               "Trying to access index " << i << ", but points.size() = " << points.size() << ", and control_points = " << control_points.size() << ".");
      return (1-t)*(1-t)*(1-t)*points[i] + 3*(1-t)*(1-t)*t*control_points[i][0] + 3.*(1-t)*t*t*control_points[i][1]+t*t*t*points[i+1];
    }


    ClosestPointOnCurve
    BezierCurve::closest_point_on_curve_segment(const Point<2> &check_point,
                                                const bool verbose) const
    {
      ClosestPointOnCurve closest_point_on_curve;
      const Point<2> &cp = check_point;
      double min_squared_distance = std::numeric_limits<double>::infinity();
      if (check_point.get_coordinate_system() == CoordinateSystem::cartesian)
        {
          for ( size_t cp_i = 0; cp_i < control_points.size(); ++cp_i)
            {
#ifndef NDEBUG
              std::stringstream output;
#endif
              const Point<2> &p1 = points[cp_i];
              const Point<2> &p2 = points[cp_i+1];
              //min_squared_distance = std::min(std::min(min_squared_distance,(check_point-p1).norm_square()),(check_point-p1).norm_square());

              // Getting an estimate for where the closest point is with a linear approximation
              const Point<2> P1P2 = p2-p1;
              const Point<2> P1Pc = check_point-p1;

              const double P2P2_dot = P1P2*P1P2;

              double est =  P2P2_dot > 0.0 ? std::min(1.,std::max(0.,(P1Pc*P1P2) / P2P2_dot)) : 1.0; // est=estimate of solution
              bool found = false;

              // based on https://stackoverflow.com/questions/2742610/closest-point-on-a-cubic-bezier-curve
              const double a_0 = 3.*control_points[cp_i][0][0]-3.*control_points[cp_i][1][0]+points[cp_i+1][0]-points[cp_i][0];
              const double a_1 = 3.*control_points[cp_i][0][1]-3.*control_points[cp_i][1][1]+points[cp_i+1][1]-points[cp_i][1];
              const double b_0 = 3.*points[cp_i][0] - 6.*control_points[cp_i][0][0]+3.*control_points[cp_i][1][0];
              const double b_1 = 3.*points[cp_i][1] - 6.*control_points[cp_i][0][1]+3.*control_points[cp_i][1][1];
              const double c_0 = -3.*points[cp_i][0] + 3.*control_points[cp_i][0][0];
              const double c_1 = -3.*points[cp_i][1] + 3.*control_points[cp_i][0][1];
              const double d_0 = points[cp_i][0];
              const double d_1 = points[cp_i][1];

              const double d_min_cp_0 = d_0-cp[0];
              const double d_min_cp_1 = d_1-cp[1];

#ifndef NDEBUG
              double estimate_point_min_cp_0_dg = a_0*est*est*est+b_0*est*est+c_0*est+d_min_cp_0;
              double estimate_point_min_cp_1_dg = a_1*est*est*est+b_1*est*est+c_1*est+d_min_cp_1;
              double min_squared_distance_cartesian_temp_dg = (estimate_point_min_cp_0_dg*estimate_point_min_cp_0_dg)+(estimate_point_min_cp_1_dg*estimate_point_min_cp_1_dg);
#endif

              for (size_t newton_i = 0; newton_i < 150; newton_i++)
                {
#ifndef NDEBUG
                  output << "  wolfram alpha: (" << a_0 << "*x^3+" << b_0 << "*x^2+"<< c_0 << "*x+" << d_0 << "-" << cp[0] << ")^2+(" << a_1 << "*x^3+" << b_1 << "*x^2+"<< c_1 << "*x+" << d_1 << "-" << cp[1] << ")^2 with x=" << est << std::endl;
#endif
                  const double est_sq = est*est;
                  const double estimate_point_min_cp_0 = a_0*est_sq*est+b_0*est_sq+c_0*est+d_min_cp_0;
                  const double estimate_point_min_cp_1 = a_1*est_sq*est+b_1*est_sq+c_1*est+d_min_cp_1;

                  const double deriv_0 = 3.0*a_0*est_sq+2.0*b_0*est+c_0;
                  const double deriv_1 = 3.0*a_1*est_sq+2.0*b_1*est+c_1;
                  const double squared_distance_cartesian = (estimate_point_min_cp_0*estimate_point_min_cp_0)+(estimate_point_min_cp_1*estimate_point_min_cp_1);

                  const double squared_distance_cartesian_derivative = 2.0*(deriv_0*estimate_point_min_cp_0 + deriv_1*estimate_point_min_cp_1);
                  const double squared_distance_cartesian_second_derivative_abs  = std::fabs(2.0*((6.0*a_0*est+2.0*b_0)*estimate_point_min_cp_0 + deriv_0*deriv_0
                                                                                                  + (6.0*a_1*est+2.0*b_1)*estimate_point_min_cp_1 + deriv_1*deriv_1));

                  if (squared_distance_cartesian_second_derivative_abs <= 0.0)
                    {
                      found = true;
                      break;
                    }

                  // the local minimum is where  squared_distance_cartesian_derivative=0 and squared_distance_cartesian_derivative>=0
                  const double update = std::min(0.5,std::max(-0.5,squared_distance_cartesian_derivative/squared_distance_cartesian_second_derivative_abs));
                  double line_search = 1.;

                  if (std::fabs(update) > 1e-1)
                    {
                      double est_test = est-update*line_search;
                      double squared_distance_cartesian_test = squared_distance_cartesian;
                      double squared_distance_cartesian_test_previous = squared_distance_cartesian;

                      for (unsigned int i = 0; i < 10; i++)
                        {
                          est_test = est-update*line_search;
                          const double est_test_sq = est_test*est_test;
                          const double estimate_point_min_cp_test_0 = a_0*est_test_sq*est_test+b_0*est_test_sq+c_0*est_test+d_min_cp_0;
                          const double estimate_point_min_cp_test_1 = a_1*est_test_sq*est_test+b_1*est_test_sq+c_1*est_test+d_min_cp_1;

                          squared_distance_cartesian_test = (estimate_point_min_cp_test_0*estimate_point_min_cp_test_0)+(estimate_point_min_cp_test_1*estimate_point_min_cp_test_1);

#ifndef NDEBUG
                          const Point<2> a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
                          const Point<2> b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
                          const Point<2> c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
                          const Point<2> d = points[cp_i];
                          const double squared_distance_cartesian_derivative_test = 2.0*(3.0*a_0*est_test*est_test+2.0*b_0*est_test+c_0)*(a_0*est_test*est_test*est_test+b_0*est_test*est_test+c_0*est_test+d_0-cp[0])
                                                                                    + 2.0*(3.0*a_1*est_test*est_test+2.0*b_1*est_test+c_1)*(a_1*est_test*est_test*est_test+b_1*est_test*est_test+c_1*est_test+d_1-cp[1]);
                          const double squared_distance_cartesian_second_derivative_test = 2.0*(6.0*a_0*est_test + 2.0*b_0)*(a_0*est_test*est_test*est_test+b_0*est_test*est_test+c_0*est_test+d_0-cp[0])
                                                                                           + 2.0*(3.0*a_0*est_test*est_test + 2.0*b_0*est_test + c_0)*(3.0*a_0*est_test*est_test + 2.0*b_0*est_test + c_0)
                                                                                           + 2.0*(6.0*a_1*est_test + 2.0*b_1)*(a_1*est_test*est_test*est_test+b_1*est_test*est_test+c_1*est_test+d_1-cp[1])
                                                                                           + 2.0*(3.0*a_1*est_test*est_test + 2.0*b_1*est_test + c_1)*(3.0*a_1*est_test*est_test + 2.0*b_1*est_test + c_1) ;
                          output << "    i: " << cp_i << ", ni: " << newton_i<< ", lsi: " << i << ", line_search_step=" << 2./3. << ": squared_distance_cartesian_test = " << squared_distance_cartesian_test << ", diff= " << squared_distance_cartesian_test-squared_distance_cartesian
                                 << ", tests: " << (squared_distance_cartesian_test_previous < squared_distance_cartesian ? "true" : "false") << ":" << (squared_distance_cartesian_test > squared_distance_cartesian_test_previous ? "true" : "false") << ", est_test=" << est_test
                                 << ", update=" << update << ", ls=" << line_search << ", up*ls=" << update *line_search << ", test deriv =" << squared_distance_cartesian_derivative_test  << ", test update=" << squared_distance_cartesian_derivative_test/fabs(squared_distance_cartesian_second_derivative_test)
                                 << ", p1=" << p1 << ", p2= " << p2 << ", poc= " << a *est_test *est_test *est_test + b *est_test *est_test+c *est_test+d << ", cp= " <<  check_point << ", ds:" << ((a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d)-check_point).norm_square() << ":" << min_squared_distance_cartesian_temp_dg
                                 << ", diff = " << squared_distance_cartesian_test-min_squared_distance_cartesian_temp_dg << std::endl;
#endif
                          if (i > 0 && (squared_distance_cartesian_test > squared_distance_cartesian_test_previous))
                            {
                              if (squared_distance_cartesian_test_previous-squared_distance_cartesian < 0)
                                {
                                  line_search *= 3./2.;
                                  break;
                                }
                            }
                          squared_distance_cartesian_test_previous = squared_distance_cartesian_test;

                          line_search *=  2./3.;
                        }
                    }

                  est -= update*line_search;

                  if (std::fabs(update) < 1e-4 || est < -0.1 || est > 1.1)
                    {
                      found = true;
                      break;
                    }
                }
#ifndef NDEBUG
              WBAssertThrow(found, "Could not find a good solution. " << output.str());
#else
              WBAssertThrow(found, "Could not find a good solution. Enable debug mode for more info.");
#endif

              const double est_min_cp_end_0 = a_0*est*est*est+b_0*est*est+c_0*est+d_min_cp_0;
              const double est_min_cp_end_1 = a_1*est*est*est+b_1*est*est+c_1*est+d_min_cp_1;
              const double min_squared_distance_temp = (est_min_cp_end_0*est_min_cp_end_0)+(est_min_cp_end_1*est_min_cp_end_1);
              if (min_squared_distance_temp < min_squared_distance)
                {
                  if (est >= -1e-8 && static_cast<double>(cp_i)+est > 0 && est-1. <= 1e-8 && est-1. < static_cast<double>(cp_i))
                    {
                      min_squared_distance = min_squared_distance_temp;
                      const Point<2> point_on_curve = Point<2>(a_0*est*est*est+b_0*est*est+c_0*est+d_0,a_1*est*est*est+b_1*est*est+c_1*est+d_1,cp.get_coordinate_system());
                      WBAssert(!std::isnan(point_on_curve[0]) && !std::isnan(point_on_curve[1]), "Point on curve has NAN entries: " << point_on_curve);
                      // the sign is rotating the derivative by 90 degrees.
                      // When moving in the direction of increasing t, left is positive and right is negative.
                      // https://www.wolframalpha.com/input?i=d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                      const Point<2> derivative_point = points[cp_i]*((6.-3.*est)*est-3.) + control_points[cp_i][0]*(est*(9*est-12)+3)
                                                        + control_points[cp_i][1]*(6.-9.*est)*est + points[cp_i+1]*3.*est*est;

                      Point<2> tangent_point = derivative_point - point_on_curve;
                      // if angle between check point and tangent point is larger than 90 degrees, return a negative distance
                      const double dot_product = (tangent_point*(check_point-point_on_curve));
                      const double sign = dot_product < 0. ? -1. : 1.;
                      tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                      closest_point_on_curve.distance = sign*std::sqrt(min_squared_distance);
                      closest_point_on_curve.parametric_fraction = est;
                      closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                      closest_point_on_curve.index = cp_i;
                      closest_point_on_curve.point = point_on_curve;
                      WBAssert(!std::isnan(point_on_curve[0]) && !std::isnan(point_on_curve[1]), "Point on curve has NAN entries: " << point_on_curve);
                      Point<2> normal = point_on_curve;
                      {
                        Point<2> derivative = Point<2>(a_0*est*est+b_0*est+c_0,a_1*est*est+b_1*est+c_1,cp.get_coordinate_system());
                        normal=derivative;
                        const double normal_size = derivative.norm();
                        if (normal_size > 0.)
                          {
                            normal[0] = derivative[1]/normal_size;
                            normal[1] = -derivative[0]/normal_size;
                          }
                      }
                      closest_point_on_curve.normal = normal;
                    }
                }
            }
        }
      else
        {
          const double cos_cp_lat = cos(cp[1]);
          for ( size_t cp_i = 0; cp_i < control_points.size(); ++cp_i)
            {
              const Point<2> &p1 = points[cp_i];
              const Point<2> &p2 = points[cp_i+1];
              // Getting an estimate for where the closest point is with a linear approximation
              const Point<2> P1P2 = p2-p1;
              const Point<2> P1Pc = check_point-p1;

              const double P2P2_dot = P1P2*P1P2;

              double est =  P2P2_dot > 0.0 ? std::min(1.,std::max(0.,(P1Pc*P1P2) / P2P2_dot)) : 1.0; // est=estimate of solution
              bool found = false;

              // only used if verbose is true
              std::stringstream output;

              Point<2> a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
              Point<2> b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
              Point<2> c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
              const Point<2> d = points[cp_i];

              Point<2> estimate_point = a*est*est*est+b*est*est+c*est+d;

              double cos_lat_dg = NaN::DSNAN;
              double sin_d_long_h_dg = NaN::DSNAN;
              double sin_d_lat_h_dg = NaN::DSNAN;
              double min_squared_distance_cartesian_temp_dg = NaN::DSNAN;
              if (verbose == true)
                {
                  cos_lat_dg = cos(estimate_point[1]);
                  sin_d_long_h_dg = sin((estimate_point[0]-cp[0])*0.5);
                  sin_d_lat_h_dg = sin((estimate_point[1]-cp[1])*0.5);
                  min_squared_distance_cartesian_temp_dg = sin_d_lat_h_dg*sin_d_lat_h_dg+sin_d_long_h_dg*sin_d_long_h_dg*cos_cp_lat*cos_lat_dg;
                  output << "cp_i=" << cp_i << ", init est = " << est << ", min_squared_distance = " << min_squared_distance << ", min_squared_distance_cartesian_temp_dg: " << min_squared_distance_cartesian_temp_dg << ", p1: " << p1 << ", p2: " << p2 << std::endl;
                  output  << std::setprecision(6) << "  wolfram: sin((" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ")*.5)^2+sin((" << a[0] << "*x^3+" << b[0] << "*x^2+"<< c[0] << "*x+" << d[0] << "-" << cp[0] << ")*.5)^2*cos(" << cp[1] << ")*cos(" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ") with x=" << est << std::endl;
                  output  << std::setprecision(10) << "  python: y=np.sin((" << a[1] << "*x**3+" << b[1] << "*x**2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << ")*.5)**2+np.sin((" << a[0] << "*x**3+" << b[0] << "*x**2+"<< c[0] << "*x+" << d[0] << "-" << cp[0] << ")*.5)**2*np.cos(" << cp[1] << ")*np.cos(" << a[1] << "*x**3+" << b[1] << "*x**2+"<< c[1] << "*x+" << d[1] << "-" << cp[1] << "); x=" << est << std::endl;
                }
              for (size_t newton_i = 0; newton_i < 150; newton_i++)
                {
                  // based on https://stackoverflow.com/questions/2742610/closest-point-on-a-cubic-bezier-curve
                  estimate_point = a*est*est*est+b*est*est+c*est+d;

                  double sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                  double sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
                  const double cos_d_lat = cos(estimate_point[1]-cp[1]);
                  const double squared_distance_cartesian = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos_d_lat;

                  double sin_dlat = sin(estimate_point[1]-cp[1]);
                  const double cos_dlong_h = cos(0.5*(estimate_point[0]-cp[0]));
                  double cos_dlat_h = cos(0.5*(estimate_point[1]-cp[1]));
                  double deriv_long = (3.0*a[0]*est*est+2.0*b[0]*est+c[0]);
                  double deriv_lat = (3.0*a[1]*est*est+2.0*b[1]*est+c[1]);

                  const double squared_distance_cartesian_derivative = cos_cp_lat*(-deriv_lat)*sin_d_long_h*sin_d_long_h*sin_dlat+cos_cp_lat*deriv_long*sin_d_long_h*cos_dlong_h*cos_d_lat+deriv_lat*sin_d_lat_h*cos_dlat_h;
                  double update = NaN::DSNAN;
                  if (std::fabs(squared_distance_cartesian_derivative) > 1e-15)
                    {
                      const double squared_distance_cartesian_second_derivative = cos_cp_lat*cos_d_lat*(-0.5*deriv_long*deriv_long*sin_d_long_h*sin_d_long_h+0.5*deriv_long*deriv_long*cos_dlong_h*cos_dlong_h+(6.0*a[0]*est+2.0*b[0])*sin_d_long_h*cos_dlong_h)+cos_cp_lat*sin_d_long_h*sin_d_long_h*(deriv_lat*deriv_lat*(-cos_d_lat)-(6.0*a[1]*est+2.0*b[1])*sin_dlat)-2.0*cos_cp_lat*deriv_long*deriv_lat*sin_d_long_h*cos_dlong_h*sin_dlat-0.5*deriv_lat*deriv_lat*sin_d_lat_h*sin_d_lat_h+0.5*deriv_lat*deriv_lat*cos_dlat_h*cos_dlat_h+(6.0*a[1]*est+2.0*b[1])*sin_d_lat_h*cos_dlat_h;

                      if (verbose == true)
                        {
                          const double squared_distance_cartesian_full = sin((a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*0.5)*sin((a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*0.5)+sin((a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])*0.5)*sin((a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])*0.5)*cos(cp[1])*cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]);
                          const double squared_distance_cartesian_derivative_full = cos(cp[1])*(-(3.0*a[1]*est*est+2.0*b[1]*est+c[1]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])+cos(cp[1])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])+(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]));
                          const double squared_distance_cartesian_second_derivative_full = cos(cp[1])*cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])*(-0.5*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))+0.5*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))+(6.0*a[0]*est+2.0*b[0])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0])))+cos(cp[1])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*((3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(-cos(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))-(6.0*a[1]*est+2.0*b[1])*sin(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))-2.0*cos(cp[1])*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*cos(0.5*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-cp[0]))*sin(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1])-0.5*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))+0.5*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))+(6.0*a[1]*est+2.0*b[1])*sin(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]))*cos(0.5*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-cp[1]));
                          output <<"sqd = " << squared_distance_cartesian <<":" << squared_distance_cartesian_full << ", diff=" << squared_distance_cartesian-squared_distance_cartesian_full << ", sqdd: " << squared_distance_cartesian_derivative <<":" << squared_distance_cartesian_derivative_full << ", diff="<< squared_distance_cartesian_derivative-squared_distance_cartesian_derivative_full << ", sqdd: " << squared_distance_cartesian_second_derivative << ":" << squared_distance_cartesian_second_derivative_full << ", diff= " << squared_distance_cartesian_second_derivative-squared_distance_cartesian_second_derivative_full << ", est: " << est << std::endl;
                        }
                      // the local minimum is where  squared_distance_cartesian_derivative=0 and squared_distance_cartesian_derivative>=0
                      update = std::min(0.5,std::max(-0.5,squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative)));
                      double line_search = 1.;
                      double est_test = est-update*line_search;
                      double squared_distance_cartesian_test = squared_distance_cartesian;
                      double squared_distance_cartesian_test_previous = squared_distance_cartesian;
                      double line_search_step = 2./3.;

                      for (unsigned int i = 0; i < 10; i++)
                        {
                          est_test = est-update*line_search;
                          estimate_point = a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d;

                          sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                          sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
                          squared_distance_cartesian_test = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos(estimate_point[1]-cp[1]);

                          if (verbose == true)
                            {
                              sin_dlat = sin(estimate_point[1]-cp[1]);
                              deriv_long = (3.0*a[0]*est_test*est_test+2.0*b[0]*est_test+c[0]);
                              deriv_lat = (3.0*a[1]*est_test*est_test+2.0*b[1]*est_test+c[1]);
                              const double squared_distance_cartesian_derivative_test = cos_cp_lat*(-deriv_lat)*sin_d_long_h*sin_d_long_h*sin_dlat+cos_cp_lat*deriv_long*sin_d_long_h*cos_dlong_h*cos_d_lat+deriv_lat*sin_d_lat_h*cos_dlat_h;
                              const double squared_distance_cartesian_second_derivative_test = cos_cp_lat*cos_d_lat*(-0.5*deriv_long*deriv_long*sin_d_long_h*sin_d_long_h+0.5*deriv_long*deriv_long*cos_dlong_h*cos_dlong_h+(6.0*a[0]*est_test+2.0*b[0])*sin_d_long_h*cos_dlong_h)+cos_cp_lat*sin_d_long_h*sin_d_long_h*(deriv_lat*deriv_lat*(-cos_d_lat)-(6.0*a[1]*est_test+2.0*b[1])*sin_dlat)-2.0*cos_cp_lat*deriv_long*deriv_lat*sin_d_long_h*cos_dlong_h*sin_dlat-0.5*deriv_lat*deriv_lat*sin_d_lat_h*sin_d_lat_h+0.5*deriv_lat*deriv_lat*cos_dlat_h*cos_dlat_h+(6.0*a[1]*est_test+2.0*b[1])*sin_d_lat_h*cos_dlat_h;
                              output << "    i: " << cp_i << ", ni: " << newton_i<< ", lsi: " << i << ", line_search_step=" << line_search_step << ": squared_distance_cartesian_test = " << squared_distance_cartesian_test << ", diff= " << squared_distance_cartesian_test-squared_distance_cartesian << ", tests: " << (squared_distance_cartesian_test_previous < squared_distance_cartesian ? "true" : "false") << ":" << (squared_distance_cartesian_test > squared_distance_cartesian_test_previous ? "true" : "false") << ", est_test=" << est_test << ", update=" << update << ", ls=" << line_search << ", up*ls=" << update *line_search << ", test deriv =" << squared_distance_cartesian_derivative_test  << ", test update=" << squared_distance_cartesian_derivative_test/fabs(squared_distance_cartesian_second_derivative_test) << ", p1=" << p1 << ", p2= " << p2 << ", poc= " << a *est_test *est_test *est_test+b *est_test *est_test+c *est_test+d << ", cp= " <<  check_point << ", ds:" << ((a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d)-check_point).norm_square() << ":" << min_squared_distance_cartesian_temp_dg << ", diff = " << squared_distance_cartesian_test-min_squared_distance_cartesian_temp_dg<< std::endl;
                            }
                          if (i > 0 && (squared_distance_cartesian_test > squared_distance_cartesian_test_previous))
                            {
                              if (squared_distance_cartesian_test_previous-squared_distance_cartesian < 0)
                                {
                                  line_search *= 1/line_search_step;
                                  break;
                                }
                              if (i> 1)
                                {
                                  line_search *= (1/line_search_step)*(1/line_search_step);
                                  est_test = est-update*line_search;
                                  estimate_point = a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d;

                                  sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
                                  sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);
                                  squared_distance_cartesian_test_previous = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos(estimate_point[1]-cp[1]);
                                  line_search_step = std::min(line_search_step*(11./10.),0.95);
                                  continue;
                                }
                            }
                          squared_distance_cartesian_test_previous = squared_distance_cartesian_test;

                          line_search *= line_search_step;
                        }
                      if (verbose == true)
                        {
                          output << "    i: " << cp_i << ", ni: " << newton_i<< ", est= " << est-update *line_search << ", ls:" << line_search << ": squared_distance_cartesian_test = " << squared_distance_cartesian_test << ", diff= " << squared_distance_cartesian_test-squared_distance_cartesian << std::endl;
                        }
                      est -= update*line_search;
                    }

                  if (std::fabs(squared_distance_cartesian_derivative) <= 1e-15 || std::fabs(update) < 1e-4 || est < -0.1 || est > 1.1)
                    {
                      found = true;
                      break;
                    }
                }

              if (verbose == true && found == false)
                {
                  // report the error and print the output
                  WBAssertThrow(found, "Could not find a good solution. " << output.str());
                }
              else if (verbose == false && found == false)
                {
                  // redo the iteration with verbose=true to be able to report the error
                  return closest_point_on_curve_segment(check_point, true);
                }

              estimate_point = a*est*est*est+b*est*est+c*est+d;

              const double sin_d_long_h = sin((estimate_point[0]-cp[0])*0.5);
              const double sin_d_lat_h = sin((estimate_point[1]-cp[1])*0.5);

              const double min_squared_distance_cartesian_temp = sin_d_lat_h*sin_d_lat_h+sin_d_long_h*sin_d_long_h*cos_cp_lat*cos(estimate_point[1]-cp[1]);

              if (min_squared_distance_cartesian_temp < min_squared_distance)
                {
                  if (est >= -1e-8 && static_cast<double>(cp_i)+est > 0 && est-1. <= 1e-8 && est-1. < static_cast<double>(cp_i))
                    {
                      min_squared_distance = min_squared_distance_cartesian_temp;
                      const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;

                      // the sign is rotating the derivative by 90 degrees.
                      // When moving in the direction of increasing t, left is positive and right is negative.
                      // https://www.wolframalpha.com/input?i=d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                      const Point<2> derivative_point = points[cp_i]*((6.-3.*est)*est-3.) + control_points[cp_i][0]*(est*(9*est-12)+3)
                                                        + control_points[cp_i][1]*(6.-9.*est)*est + points[cp_i+1]*3.*est*est;

                      Point<2> tangent_point = derivative_point - point_on_curve;
                      // if angle between check point and tangent point is larger than 90 degrees, return a negative distance
                      const double dot_product = (tangent_point*(check_point-point_on_curve));
                      const double sign = dot_product < 0. ? -1. : 1.;
                      tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                      closest_point_on_curve.distance = sign*std::sqrt(min_squared_distance);
                      closest_point_on_curve.parametric_fraction = est;
                      closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                      closest_point_on_curve.index = cp_i;
                      closest_point_on_curve.point = point_on_curve;
                      Point<2> normal = point_on_curve;
                      {
                        Point<2> derivative = a*est*est+b*est+c;
                        normal=derivative;
                        const double normal_size = derivative.norm();
                        if (normal_size > 0.)
                          {
                            normal[0] = derivative[1]/normal_size;
                            normal[1] = -derivative[0]/normal_size;
                          }
                      }
                      closest_point_on_curve.normal = normal;
                    }
                }
            }

        }
      return closest_point_on_curve;
    }
  } // namespace Objects
} // namespace WorldBuilder

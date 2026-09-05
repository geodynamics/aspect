/*
  Copyright (C) 2017 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include "lithosphere_rift_IC.h"
#include "lithosphere_rift_IT.h"
#include <aspect/geometry_model/box.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    LithosphereRift<dim>::
    initialize ()
    {
      // Check that the corresponding initial temperature model is used.
      const std::vector<std::string> active_initial_temperature_models = this->get_initial_temperature_manager().get_active_initial_temperature_names();
      AssertThrow(find(active_initial_temperature_models.begin(),active_initial_temperature_models.end(), "lithosphere with rift") != active_initial_temperature_models.end(),
                  ExcMessage("The lithosphere with rift initial composition plugin requires the lithosphere with rift initial temperature plugin."));
    }

    template <int dim>
    double
    LithosphereRift<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      // Retrieve the indices of the fields that represent the lithospheric layers.
      // We assume a 3-layer system with an upper crust, lower crust and lithospheric mantle.
      const unsigned int id_upper = this->introspection().compositional_index_for_name("upper");
      const unsigned int id_lower = this->introspection().compositional_index_for_name("lower");
      const unsigned int id_mantle_L = this->introspection().compositional_index_for_name("mantle_L");

      // Determine coordinate system
      const bool cartesian_geometry = Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model());

      // Get the surface coordinates of the point under consideration
      const Point<dim-1> surface_point = surface_position(position, cartesian_geometry);

      // Compute the local thickness of the upper crust, lower crust and mantle part of the lithosphere
      // (in this exact order) based on the distance from the rift axis and the polygons.
      const std::vector<double> local_thicknesses = compute_local_thicknesses(surface_point);

      // Get depth with respect to the surface.
      const double depth = this->get_geometry_model().depth(position);

      // Check which layer the current point lies in and return a value of 1 if the field corresponds to the layer.
      if (depth <= local_thicknesses[0] && compositional_index == id_upper)
        return 1.;
      else if (depth > local_thicknesses[0] && depth <= local_thicknesses[0] + local_thicknesses[1]
               && compositional_index == id_lower)
        return 1.;
      else if (depth > local_thicknesses[0] + local_thicknesses[1] && depth <= local_thicknesses[0] + local_thicknesses[1] + local_thicknesses[2]
               && compositional_index == id_mantle_L)
        return 1.;
      else
        return 0.;
    }

    template <int dim>
    double
    LithosphereRift<dim>::
    distance_to_rift (const Point<dim-1> &surface_position) const
    {
      // Initiate distance with large value
      double distance_to_rift_axis = 1e23;
      double temp_distance = 0;

      // Loop over all line segments
      for (unsigned int i_segments = 0; i_segments < rift_point_list.size(); ++i_segments)
        {
          // The Utilities function only works in 3d, so compute distance here in 2d.
          temp_distance = (dim == 2) ? std::abs(surface_position[0]-rift_point_list[i_segments][0][0])
                          : std::abs(Utilities::distance_to_line(rift_point_list[i_segments], Point<2>(surface_position[0],surface_position[dim-2])));

          // Get the minimum distance
          distance_to_rift_axis = std::min(distance_to_rift_axis, temp_distance);
        }

      return distance_to_rift_axis;
    }

    template <int dim>
    Point<dim-1>
    LithosphereRift<dim>::
    surface_position (const Point<dim> &position,
                      const bool cartesian_geometry) const
    {
      Point<dim-1> surface_point;
      if (cartesian_geometry)
        {
          for (unsigned int d=0; d<dim-1; ++d)
            surface_point[d]=position[d];
        }
      // chunk or spherical geometries
      else
        {
          // Get spherical coordinates in radius [m], lon [rad], colat [rad] format
          const std::array<double,dim> spherical_point = Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
          // return lon [degrees], lat [degrees]
          for (unsigned int d=0; d<dim-1; ++d)
            surface_point[d] = spherical_point[d+1]*180./numbers::PI;
          if (dim == 3)
            surface_point[dim-2] = 90. - surface_point[dim-2];
        }

      return surface_point;
    }

    template <int dim>
    std::pair<double,unsigned int>
    LithosphereRift<dim>::
    distance_to_polygon (const Point<dim-1> &surface_position) const
    {
      // Inside the polygon is positive, outside negative. We assume the different polygons do not overlap.
      // The Utilities function only works in 3d, so compute distance here in 2d.
      double max_distance = -1e24;
      unsigned int max_distance_polygon = 0;
      for (unsigned int n = 0; n<polygon_point_list.size(); ++n)
        {
          double temp_distance = 0;
          if (dim == 2)
            {
              double sign = -1.;
              if (surface_position[0]>polygon_point_list[n][0][0] && surface_position[0]<polygon_point_list[n][1][0])
                sign = 1.;
              temp_distance = sign * std::min(std::abs(polygon_point_list[n][1][0] - surface_position[0]), std::abs(surface_position[0] - polygon_point_list[n][0][0]));
            }
          else
            {
              temp_distance = Utilities::signed_distance_to_polygon<dim>(polygon_point_list[n], Point<2>(surface_position[0],surface_position[dim-2]));
            }

          if (temp_distance > max_distance)
            {
              max_distance = temp_distance;
              max_distance_polygon = n;
            }
        }
      return std::pair<double, unsigned int> (max_distance, max_distance_polygon);
    }

    template <int dim>
    std::vector<double>
    LithosphereRift<dim>::
    compute_local_thicknesses(const Point<dim - 1> &surface_point) const
    {
      // Get the distance to the rift segments along a path parallel to the surface
      const double distance_to_rift_axis = distance_to_rift(surface_point);

      // Get the signed distance to potential polygons of different lithospheric thicknesses
      const std::pair<double, unsigned int> distance_to_L_polygon = distance_to_polygon(surface_point);

      // Compute the local thickness of the upper crust, lower crust and mantle part of the lithosphere
      // (in this exact order) based on the distance from the rift axis and the polygons.
      // The transition from reference and/or rift perturbation lithosphere to polygon lithosphere
      // is smoothed by a hyperbolic tangent.
      const double polygon_contribution = (0.5 + 0.5 * std::tanh(distance_to_L_polygon.first / sigma_polygon));
      const double rift_contribution = (0.5 - 0.5 * std::tanh(distance_to_L_polygon.first / sigma_polygon));

      // The rift perturbation follows a Gaussian contribution.
      std::vector<double> local_thicknesses(3);
      for (unsigned int i = 0; i < 3; ++i)
        local_thicknesses[i] = (1.0 - A_rift[i] * std::exp((-std::pow(distance_to_rift_axis, 2) / (2.0 * std::pow(sigma_rift, 2)))))
                               * reference_thicknesses[i] * rift_contribution + polygon_contribution * polygon_thicknesses[distance_to_L_polygon.second][i];

      return local_thicknesses;
    }

    template <int dim>
    void
    LithosphereRift<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          prm.declare_entry("Standard deviation of Gaussian rift geometry", "20000",
                            Patterns::Double(0),
                            "The standard deviation of the Gaussian distribution of the thinning/thickening "
                            "of the lithosphere thicknesses. This parameter is taken to be the same "
                            "for all rift segments. Units: \\si{\\meter} or degrees.");
          prm.declare_entry("Amplitude of Gaussian rift geometry", "0.2",
                            Patterns::List(Patterns::Double(-1, 1)),
                            "The amplitude of the Gaussian distribution of the thinning/thickening of the. "
                            "lithosphere thicknesses. This parameter is taken to be the same for all rift segments, "
                            "but can be set to vary per lithosphere layer. "
                            "Units: none.");
          prm.declare_entry("Half width of polygon smoothing", "20000",
                            Patterns::Double(0),
                            "The half width of the hyperbolic tangent smoothing used to transition to the "
                            "lithospheric thicknesses of the polygon. This parameter is taken to be the same for all polygons. "
                            "Units: \\si{\\meter} or degrees.");
          prm.declare_entry("Layer thicknesses", "30000.",
                            Patterns::List(Patterns::Double(0)),
                            "List of reference lithospheric layer thicknesses, i.e., the thicknesses of "
                            "the upper crust, lower crust and lithospheric mantle layers. "
                            "Units: \\si{\\meter}");
          prm.declare_entry("Rift axis line segments",
                            "",
                            Patterns::Anything(),
                            "The line segments that represent the rift axis separated by a semi-colon. "
                            "Each segment is made up of two points that represent horizontal coordinates (x,y) or (lon,lat). "
                            "The exact format for the point list describing the segments is "
                            "\"x1,y1>x2,y2;x2,y2>x3,y3;x4,y4>x5,y5\". Note that the segments can be connected "
                            "or isolated. The units of the coordinates are "
                            "dependent on the geometry model. In the box model they are in meters, in the "
                            "chunks they are in degrees. Units: \\si{\\meter} or degrees.");
          prm.declare_entry("Lithospheric polygons",
                            "",
                            Patterns::List(Patterns::Anything()),
                            "The points making up polygons that represent an area of different lithospheric thickness. "
                            "The polygons are separated by semicolons. Each polygon is a list of "
                            "points that represent horizontal coordinates (x,y) or (lon,lat). "
                            "The exact format for the point list describing a polygon is "
                            "\"x1,y1>x2,y2>x3,y3>x4,y4>x5,y5\". Note that the polygon is assumed to be closed. "
                            "The units of the coordinates are dependent on the geometry model. "
                            "In the box model they are in meters, in the chunks they are in degrees. Units: \\si{\\meter} or degrees.");
          prm.declare_entry("Lithospheric polygon layer thicknesses", "30000.",
                            Patterns::List(Patterns::List(Patterns::Double(0), 0, 3, ","), 0, 10, ";"),
                            "List of thicknesses of the lithospheric layers for each polygon."
                            "For each polygon, a total of 3 thicknesses should be given (upper crust, lower crust, mantle lithosphere)."
                            "If only one value is given, then all layers are assigned the same value. Units: \\si{\\meter}");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    LithosphereRift<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Check that the required compositional fields exist.
      AssertThrow(this->introspection().compositional_name_exists("upper"),
                  ExcMessage("We need a compositional field called 'upper' representing the upper crust."));
      AssertThrow(this->introspection().compositional_name_exists("lower"),
                  ExcMessage("We need a compositional field called 'lower' representing the lower crust."));
      AssertThrow(this->introspection().compositional_name_exists("mantle_L"),
                  ExcMessage("We need a compositional field called 'mantle_L' representing the lithospheric part of the mantle."));

      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          sigma_rift             = prm.get_double ("Standard deviation of Gaussian rift geometry");
          A_rift                 = Utilities::possibly_extend_from_1_to_N(Utilities::string_to_double(Utilities::split_string_list(prm.get("Amplitude of Gaussian rift geometry"))),
                                                                          3,
                                                                          "Amplitude of Gaussian rift geometry");
          sigma_polygon          = prm.get_double ("Half width of polygon smoothing");
          reference_thicknesses  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Layer thicknesses"))),
                                                                           3,
                                                                           "Layer thicknesses");
          // Read in the string of segments
          const std::string temp_all_segments = prm.get("Rift axis line segments");
          // Split the string into segment strings
          const std::vector<std::string> temp_segments = Utilities::split_string_list(temp_all_segments,';');
          const unsigned int n_temp_segments = temp_segments.size();
          rift_point_list.resize(n_temp_segments);
          // Loop over the segments to extract the points
          for (unsigned int i_segment = 0; i_segment < n_temp_segments; ++i_segment)
            {
              // In 3d a line segment consists of 2 points,
              // in 2d of only 1 (the horizontal ridge axis runs orthogonal to x and y into the screen)
              const std::vector<std::string> temp_segment = Utilities::split_string_list(temp_segments[i_segment],'>');

              if (dim == 3)
                {
                  AssertThrow(temp_segment.size() == 2,ExcMessage ("The given segment '" + temp_segment[i_segment] + "' is not correct. "
                                                                   "It should only contain 2 parts: "
                                                                   "the two points of the segment, separated by a '>'."));
                }
              else
                {
                  AssertThrow(temp_segment.size() == 1,ExcMessage ("The given segment '" + temp_segment[i_segment] + "' is not correct. "
                                                                   "It should only contain 1 part: "
                                                                   "the point representing the rift axis."));
                }

              // Loop over the dim-1 points of each segment
              for (unsigned int i_points = 0; i_points < dim-1; ++i_points)
                {
                  const std::vector<double> temp_point = Utilities::string_to_double(Utilities::split_string_list(temp_segment[i_points],','));

                  if (dim == 3)
                    {
                      AssertThrow(temp_point.size() == 2,ExcMessage ("The given coordinates of segment '" + temp_segment[i_points] + "' are not correct. "
                                                                     "It should only contain 2 parts: "
                                                                     "the two coordinates of the segment end point, separated by a ','."));
                    }
                  else
                    {
                      AssertThrow(temp_point.size() == 1,ExcMessage ("The given coordinates of segment '" + temp_segment[i_points] + "' are not correct. "
                                                                     "It should only contain 1 part: "
                                                                     "the one coordinate of the segment end point."));
                    }

                  // Add the point to the list of points for this segment. Two points are set
                  // even in 2d, although in 2d they are the same.
                  rift_point_list[i_segment][i_points][0] = temp_point[0];
                  rift_point_list[i_segment][i_points][1] = temp_point[dim-2];
                }
            }

          // Split the polygon string into the separate polygons
          const std::vector<std::string> temp_polygons = Utilities::split_string_list(prm.get("Lithospheric polygons"),';');
          const std::vector<std::string> temp_thicknesses = Utilities::split_string_list(prm.get("Lithospheric polygon layer thicknesses"),';');
          const unsigned int n_polygons = temp_polygons.size();
          AssertThrow(temp_thicknesses.size() == n_polygons,
                      ExcMessage("The number of polygons (" + Utilities::int_to_string(n_polygons) +
                                 ") does not correspond to the number of polygons for which a thickness is prescribed (" +
                                 Utilities::int_to_string(temp_thicknesses.size()) + ")."));

          polygon_point_list.resize(n_polygons);
          polygon_thicknesses.resize(n_polygons);

          for (unsigned int i_polygons = 0; i_polygons < n_polygons; ++i_polygons)
            {
              polygon_thicknesses[i_polygons] = Utilities::string_to_double(Utilities::split_string_list(temp_thicknesses[i_polygons],','));
              AssertThrow(polygon_thicknesses[i_polygons].size() == 3,
                          ExcMessage ("The number of layer thicknesses should be equal to 3 for polygon: " + Utilities::int_to_string(i_polygons) +
                                      " but it is " + Utilities::int_to_string(polygon_thicknesses[i_polygons].size())));

              // Split the polygon string into point strings
              const std::vector<std::string> temp_points = Utilities::split_string_list(temp_polygons[i_polygons],'>');
              const unsigned int n_temp_points = temp_points.size();
              if (dim == 3)
                {
                  AssertThrow(n_temp_points>=3, ExcMessage ("The number of polygon points should be equal to or larger than 3 in 3d."));
                }
              else
                {
                  AssertThrow(n_temp_points==2, ExcMessage ("The number of polygon points should be equal to 2 in 2d."));
                }

              polygon_point_list[i_polygons].resize(n_temp_points);

              // Loop over the points of the polygon. Each point coordinate should consist of 2 (in 3d) or 1 (in 2d) values.
              for (unsigned int i_points = 0; i_points < n_temp_points; ++i_points)
                {
                  const std::vector<double> temp_point = Utilities::string_to_double(Utilities::split_string_list(temp_points[i_points],','));
                  AssertThrow(temp_point.size() == dim-1,ExcMessage ("The given coordinates of point '" + temp_points[i_points] + "' are not correct. "
                                                                     "It should only contain 1 (2d) or 2 (in 3d) parts: "
                                                                     "the longitude/x (and latitude/y in 3d) coordinate (separated by a ',')."));

                  // Add the point to the list of points for this segment
                  polygon_point_list[i_polygons][i_points][0] = temp_point[0];
                  polygon_point_list[i_polygons][i_points][1] = temp_point[dim-2];
                }

              // For simplicity later on, we want the polygon coordinates of each point in 2d to be increasing.
              if  (dim == 2)
                AssertThrow(polygon_point_list[i_polygons][0][0] < polygon_point_list[i_polygons][1][0], ExcMessage("The order of the x coordinates of the 2 points "
                            "of each 2d polygon should be ascending. "));

            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(LithosphereRift,
                                              "lithosphere with rift",
                                              "A class that implements initial conditions for a continental rift "
                                              "by computing the variable thickness of the upper crust, lower crust and mantle "
                                              "lithosphere layers. Thicknesses can change towards the axis of the "
                                              "rift according to a Gaussian distribution whose amplitude and standard "
                                              "deviation can be specified from the input file. "
                                              "Additional variations in layer thicknesses can be introduced through the "
                                              "addition of polygons. Note that the different polygons "
                                              "are assumed not to overlap. The user can decide how overlapping rifts "
                                              "and polygons are treated. ")
  }
}

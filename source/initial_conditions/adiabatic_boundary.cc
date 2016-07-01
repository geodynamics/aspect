/*
   Copyright (C) 2016 by the authors of the ASPECT code.

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
   along with ASPECT; see the file doc/COPYING.  If not see
   <http://www.gnu.org/licenses/>.
*/


#include <aspect/initial_conditions/adiabatic_boundary.h>
#include <fstream>
#include <iostream>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    template <int dim>
    std::pair<double, double>
    AdiabaticBoundary<dim>::lat_long_from_xyz_WGS84(const Point<3> &pos) const
    {
      /* Define WGS84 ellipsoid constants. */
      const double radius = 6378137;
      const double ellipticity = 8.1819190842622e-2;

      const double b = std::sqrt(radius * radius
                                 * (1 - ellipticity * ellipticity));
      const double ep = std::sqrt((radius * radius - b * b) / (b * b));
      const double p = std::sqrt(pos(0) * pos(0) + pos(1) * pos(1));
      const double th = std::atan2(radius * pos(2), b * p);
      const double lon = std::atan2(pos(1), pos(0));
      const double lat = std::atan2((pos(2) + ep * ep * b * std::sin(th)
                                     * std::sin(th) * std::sin(th)),
                                    (p - (ellipticity
                                          * ellipticity
                                          * radius
                                          * (std::cos(th) * std::cos(th)
                                             * std::cos(th)))));

      /* Convert to degrees. */
      const double lon_degrees = lon * (180. / numbers::PI);
      const double lat_degrees = lat * (180. / numbers::PI);

      /* Set all longitudes between [0,360]. */
      if (lon_degrees < 0.)
        return std::make_pair(lat_degrees, lon_degrees + 360.);
      else if (lon_degrees > 360.)
        return std::make_pair(lat_degrees, lon_degrees - 360.);
      else
        return std::make_pair(lat_degrees, lon_degrees);

    }

    template <>
    double
    AdiabaticBoundary<2>::get_isotherm_depth(const double,
                                             const double) const
    {
      AssertThrow (false, ExcMessage ("The 'adiabatic boundary' initial temperature plugin is only implemented for 3d cases."));
      return 0;
    }

    template <>
    double
    AdiabaticBoundary<3>::get_isotherm_depth(const double latitude,
                                             const double longitude) const
    {
      /**
       * Loop over the entire array and see if we find a point
       * that's within delta of what we're looking for. If this
       * point is not within the delta then test if it is within twice
       * as much as delta.
       */
      for (unsigned int i = 0; i <= depths_iso.size();)
        {
          if (std::fabs(latitude - latitudes_iso[i]) <= delta && std::fabs(longitude - longitudes_iso[i]) <= delta)
            {
              return depths_iso[i];
            }
          else
            i++;
        }
      Assert (false, ExcInternalError());
      return 0;
    }

    template <int dim>
    double
    AdiabaticBoundary<dim>::radius_WGS84(const double theta) const
    {
      const double eccentricity    = 8.1819190842622e-2;
      const double semi_major_axis = 6378137.0;
      return semi_major_axis/std::sqrt(1- eccentricity * eccentricity
                                       * std::sin(numbers::PI*theta/180)*std::sin(numbers::PI*theta/180));
    }

    template <>
    double
    AdiabaticBoundary<2>::initial_temperature (const Point<2> &) const
    {
      AssertThrow (false, ExcMessage ("The 'adiabatic boundary' initial temperature plugin is only implemented for 3d cases."));
      return 0;
    }

    template <int dim>
    double
    AdiabaticBoundary<dim>::initial_temperature (const Point<dim> &position) const
    {
      const std::pair<double, double> lat_long = lat_long_from_xyz_WGS84(position);
      const double depth                       = radius_WGS84(lat_long.first) - position.norm();
      const double isotherm_depth              = get_isotherm_depth(lat_long.first, lat_long.second);
      if (depth > isotherm_depth)
        return isotherm_temperature + (depth - isotherm_depth) * temperature_gradient;
      else
        return surface_temperature + (depth/isotherm_depth)*(isotherm_temperature - surface_temperature);
    }


    template <int dim>
    void
    AdiabaticBoundary<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Adiabatic boundary");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-conditions/adiabatic-boundary/",
                             Patterns::DirectoryName (),
                             "The path to the isotherm depth data file");
          prm.declare_entry ("Isotherm depth filename",
                             "adiabatic_boundary.txt",
                             Patterns::FileName (),
                             "File from which the isotherm depth data is read.");
          prm.declare_entry ("Isotherm temperature", "1673.15",
                             Patterns::Double (0),
                             "The value of the isothermal boundary temperature. Units: Kelvin.");
          prm.declare_entry ("Surface temperature", "273.15",
                             Patterns::Double (0),
                             "The value of the suface temperature. Units: Kelvin.");
          prm.declare_entry ("Adiabatic temperature gradient", "0.0005",
                             Patterns::Double (0),
                             "The value of the adiabatic temperature gradient. Units: $K m^{-1}$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    AdiabaticBoundary<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Adiabatic boundary");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));

          isotherm_file_name = prm.get("Isotherm depth filename");
          isotherm_temperature = prm.get_double("Isotherm temperature");
          surface_temperature  = prm.get_double("Surface temperature");
          temperature_gradient = prm.get_double("Adiabatic temperature gradient");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      AssertThrow ((dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*>
                    (&this->get_geometry_model()) != 0),
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is an ellipsoidal chunk."));

      const std::string filename = data_directory+isotherm_file_name;

      /**
       * Read data from disk and distribute among processes
       */
      std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

      /**
       * Reading data lines.
       */
      double latitude_iso, longitude_iso, depth_iso;
      while (in >> latitude_iso >> longitude_iso >> depth_iso)
        {
          latitudes_iso.push_back(latitude_iso);
          longitudes_iso.push_back(longitude_iso);
          depths_iso.push_back(depth_iso*1000.);
        }

      /**
       * Find first 2 numbers that are different to use in
       * calculating half the difference between each position as delta.
       */
      if (std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9)
        {
          /**
           * Calculate delta as half the latitude distance.
          */
          delta = std::fabs((0.5)*(latitudes_iso[0] - latitudes_iso[1]));
        }
      else
        {
          /**
           * Calculate delta as half the longitude distance.
          */
          delta = std::fabs((0.5)*(longitudes_iso[0] - longitudes_iso[1]));
        }

    }

  }
}


namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(AdiabaticBoundary,
                                       "adiabatic boundary",
                                       "An initial temperature condition that allows for discretizing "
                                       "the model domain into two layers separated by a user-defined "
                                       "isothermal boundary using a table look-up approach. The user includes an "
                                       "input ascii data file that is formatted as 3 columns of 'latitude', "
                                       "'longitude', and 'depth', where 'depth' is in kilometers and "
                                       "represents the depth of an initial temperature of 1673.15 K (by default). "
                                       "The temperature is defined from the surface (273.15 K) to the isotherm "
                                       "as a linear gradient. Below the isotherm the temperature increases "
                                       "approximately adiabatically (0.0005 K per meter). This initial temperature condition "
                                       "is designed specifically for the ellipsoidal chunk geometry model.")
  }
}

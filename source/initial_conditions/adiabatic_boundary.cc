/*
   Copyright (C) 2015 by the authors of the ASPECT code.
 
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


#include <aspect/initial_conditions/interface.h>
#include <aspect/initial_conditions/adiabatic_boundary.h>
#include <fstream>
#include <iostream>

/**
 * A class that computes an initial temperature field based
 * on a user-defined adiabatic boundary for the 3D ellipsoid
 * chunk geometry model. It discretizes the model domain into
 * two regions separated by an isotherm boundary below which
 * the temperature increases adiabatically. Above the
 * user-defined isotherm boundary the temperature linearly
 * increases from a surface temperature (273.15 K or 0 degree C)
 * to the isotherm (1673.15 K or 1600 degree C). The user
 * defines the location of the isotherm boundary with an ascii
 * data file with the format defined in the ASPECT manual. Note
 * that the latitudinal and longitudinal bounds of the ascii input
 * data file need to be at least 1 degree wider than the bounds
 * you use to define the region of your ellipsoid chunk geometry.
 *
 * This plugin is developed by Tahiry Rajaonarison, D. Sarah Stamps,
 * and Wolfgang Bangerth.
 */


namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * Function that returns latitude and longitude from ECEF Cartesian
     * coordinates that account for ellipsoidal shape of the Earth
     * with WGS84 parameters.
     */

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
    	    const double lon_degrees = lon * (180 / numbers::PI);
    	    const double lat_degrees = lat * (180 / numbers::PI);

    	    /* Set all longitudes between [0,360]. */
    	    if (lon_degrees < 0)
    	      return std::make_pair(lat_degrees, lon_degrees + 360);
    	    else if (lon_degrees > 360)
    	      return std::make_pair(lat_degrees, lon_degrees - 360);
    	    else
    	      return std::make_pair(lat_degrees, lon_degrees);

    }

    template <>
    double
    AdiabaticBoundary<2>::get_isotherm_depth(const double,
                                              const double) const
    {
    	    AssertThrow(false, ExcInternalError());
            return 0;
    } 

    template <int dim>
    double
    AdiabaticBoundary<dim>::get_isotherm_depth(const double latitude,
                                                  const double longitude) const
    {
          /**
           * Loop over the entire array and see if we find a point
           * that's within delta of what we're looking for. If this
           * point is within the delta then test if it is within twice
           * as much as delta.
           */
        for (unsigned int i = 0; i <= depths_iso.size();)
	    if (std::fabs(latitude - latitudes_iso[i]) <= delta && std::fabs(longitude - longitudes_iso[i]) < delta)
            {
                  return -depths_iso[i]*1000;
            }
        else if (std::fabs(latitude - latitudes_iso[i]) <= delta*2 && std::fabs(longitude - longitudes_iso[i]) < delta*2)         
            {
                 return -depths_iso[i]*1000;
            }
        else 
                 i++;
        Assert(false, ExcInternalError());
        return 0;
     }
     
    /**
     * Calculate distance from the center of the Earth to the surface of WGS84.
     */

    template <int dim>
    double
    AdiabaticBoundary<dim>::radius_WGS84(const double &theta) const
    {
          const double eccentricity    = 8.1819190842622e-2;
          const double semi_major_axis = 6378137.0;
          return semi_major_axis/std::sqrt(1- eccentricity * eccentricity
        		                          * std::sin(numbers::PI*theta/180)*std::sin(numbers::PI*theta/180));
    }

    template <int dim>
    double
    AdiabaticBoundary<dim>::initial_temperature (const Point<dim> &) const
    {
          Assert (false, ExcNotImplemented());
          return 0;
    }

    /**
     * Get the depth of the adiabatic boundary isotherm for the current lat/long.
     * If above the isotherm use a linear interpolation to the surface.
     * If below the isotherm use an increase of .5 degrees per kilometer.
     */

    template <>
    double
    AdiabaticBoundary<3>::initial_temperature (const Point<3> &position) const
    {
          const std::pair<double, double> lat_long = lat_long_from_xyz_WGS84(position);
          const double depth                       = position.norm() - radius_WGS84(lat_long.first);
          const double isotherm_depth              = get_isotherm_depth(lat_long.first, lat_long.second);
          if (depth < isotherm_depth)
          return isotherm_temperature - (depth - isotherm_depth) * 0.0005;
          else
          return surface_temperature + ( depth/isotherm_depth)*(isotherm_temperature - surface_temperature);
     } 


    template <int dim>
    void
    AdiabaticBoundary<dim>::declare_parameters(ParameterHandler &prm)
    {
         prm.enter_subsection("Initial conditions");
         {
           prm.enter_subsection("Adiabatic boundary");
           {
             prm.declare_entry ("Isotherm depth filename",
                                "adiabatic_boundary.txt",
                                Patterns::FileName(),
                                "File from which the base of the lithosphere depth data is read. "
                                "Note that full path (/$DIRECTORY_PATH/filename) of the location "
                                "of the file must be provided. ");
             prm.declare_entry ("Isotherm temperature", "1673.15",
            		            Patterns::Double (0),
                                "The value of the isothermal boundary temperature. Units: Kelvin. ");
             prm.declare_entry ("Surface temperature", "273.15",
                         		Patterns::Double (0),
                                "The value of the suface temperature. Units: Kelvin. ");

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
            isotherm_file_name = prm.get("Isotherm depth filename");
            isotherm_temperature = prm.get_double("Isotherm temperature");
            surface_temperature  = prm.get_double("Surface temperature");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        std::ifstream input1(isotherm_file_name.c_str());
        AssertThrow (input1.is_open(),
                     ExcMessage (std::string("Can't read from file <") + isotherm_file_name + ">"));


        /**
         * Reading the data file.
         */
        while (true)
             {
               double latitude_iso, longitude_iso, depth_iso;
               input1 >> latitude_iso >>longitude_iso >> depth_iso;
               if (input1.eof())
               break;
               latitudes_iso.push_back(latitude_iso);
               longitudes_iso.push_back(longitude_iso);
               depths_iso.push_back(depth_iso);
             }
        /**
         * Find first 2 numbers that are different to use in
         * calculating half the difference between each position as delta.
         */

        if (std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9)
             {
               /**
                * Calculate delta as half the distance between points.
		        */
               delta = std::fabs((0.5)*(latitudes_iso[0] - latitudes_iso[1]));
             }
        else
             {
               /**
                * Calculate delta as half the distance between points.
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
                                       "isotherm boundary using a table look-up approach. The user includes an "
                                       "input ascii data file that is formatted as 3 columns of 'latitude', "
									   "'longitude', and 'depth', where 'depth' is in kilometers and "
									   "represents that depth of an initial temperature of 1673.15 K (by default). "
									   "The temperature is defined from the surface (273.15 K) to the isotherm "
									   "as a linear gradient. Below the isotherm the temperature increases " 
							           "adiabatically (0.0005 K per meter). This initial temperature condition "
									   "is designed specifically for the ellipsoidal chunk geometry model.")



  }
}

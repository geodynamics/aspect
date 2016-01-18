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

/** This initial condition is designed for the 3D ellipsoid chunk geometry
 * model. It discretizes the model domain into two regions separated by an
 * isotherm below with the temperature increases adiabatically. The user
 * defines the location of the thermal isotherm with a data file with the
 * format defined in the ASPECT manual. Note that the latitudinal and
 * longitudinal bounds of the ascii input data file needs to be at least 1
 * degree wider than the bounds you use to define the ellipsoid chunk geometry
 * model.
 * This plugin is developed by Tahiry Rajaonarison, D. Sarah Stamps, and Wolfgang Bangerth.”
 */


namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * Function that return latitude and longitude from cartesian
     * coordinates that account for ellipsoidal shape of the Earth
     */

    template <int dim>
    std::pair<double, double>
    AdiabaticBoundary<dim>::lat_long_from_xyz_wgs84(const Point<3> &pos) const
    {
    	    /* WGS84 ellipsoid constants */
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

    	    /* convert to degrees */
    	    const double lon_degrees = lon * (180 / numbers::PI);
    	    const double lat_degrees = lat * (180 / numbers::PI);

    	    /* set all longitudes between [0,360] */
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
    {       abort();
            return 0;
    } 

    template <int dim>
    double
    AdiabaticBoundary<dim>::get_isotherm_depth(const double latitude,
                                                  const double longitude) const
    {
          /* Loop over the entire array and see if we find a point
           * that's within delta of what we're looking for. The data
           * is arranged in a way that keeps the latitude constant
           * while running over all longitudes. When we're finished
           * with that change the latitude by one step. So if the
           * latitude is wrong, we can definitely skip ahead a whole
           * set of longitudes. The number of values to skip is calculated.
           */

        for (unsigned int i = 0; i <= latitudes_iso.size();)
	    if (std::fabs(latitude - latitudes_iso[i]) <= delta)
              {
                if (std::fabs(longitude - longitudes_iso[i]) <= delta)
                  return -depths_iso[i]*1000;
                else
                  ++i;
              }
         else
              i += number_coords_depth;

          Assert(false, ExcInternalError());
          return 0;
     }
     
    /**
     * Calculate distance from the center of the Earth to the surface of WGS84.
     */

    template <int dim>
    double
    AdiabaticBoundary<dim>::radius_wgs84(const double &theta) const
    {
          const double eccentricity    = 8.1819190842622e-2;
          const double semi_major_axis = 6378137.0;
          return semi_major_axis/std::sqrt(1- eccentricity * eccentricity * std::sin(numbers::PI*theta/180)*std::sin(numbers::PI*theta/180));
    }

    template <int dim>
    double
    AdiabaticBoundary<dim>::initial_temperature (const Point<dim> &) const
    {
          Assert (false, ExcNotImplemented());
          return 0;
    }

    /** Get the depth of the adiabatic boundary isotherm for the current lat/long.
     * If above the isotherm use a linear interpolation to the surface.
     * If below the isotherm use an increase of .5 degrees per kilometer.
     */

    template <>
    double
    AdiabaticBoundary<3>::initial_temperature (const Point<3> &position) const
    {
          const std::pair<double, double> lat_long = lat_long_from_xyz_wgs84(position);
          const double depth                       = position.norm() - radius_wgs84(lat_long.first);
          const double isotherm_temp               = 1673.15;
          const double surface_temp                = 273.15;
          const double isotherm_depth              = get_isotherm_depth(lat_long.first, lat_long.second);
          if (depth < isotherm_depth)
          return isotherm_temp - (depth - isotherm_depth) * 0.0005;
          else
          return surface_temp + ( depth/isotherm_depth)*(isotherm_temp - surface_temp);
     } 


    template <int dim>
    void
    AdiabaticBoundary<dim>::declare_parameters(ParameterHandler &prm)
    {
         prm.enter_subsection("Initial conditions");
         {
           prm.enter_subsection("Adiabatic boundary isotherm");
           {
             prm.declare_entry("Adiabatic Boundary Filename",
                               "adiabatic_boundary.txt",
                               Patterns::FileName(),
                               "Surface coordinates and depths to the 1673.15 K isotherm. Units: degrees and kilometers.");
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
          prm.enter_subsection("Adiabatic boundary isotherm");
          {
            adiabatic_boundary_file_name = prm.get("Adiabatic Boundary Filename");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        std::ifstream input1(adiabatic_boundary_file_name.c_str());
        AssertThrow (input1.is_open(),
                     ExcMessage (std::string("Can't read from file <") + adiabatic_boundary_file_name + ">"));


        /** Start loop with int count to calculate delta below for adiabatic boundary
         */
	int count = 0;
        while (true)
          { 
            double latitude_iso, longitude_iso, depth_iso;
            input1 >> latitude_iso >>longitude_iso >> depth_iso;
            if (input1.eof())
            break;
            latitudes_iso.push_back(latitude_iso);
            longitudes_iso.push_back(longitude_iso);
            depths_iso.push_back(depth_iso);

            /** Find first 2 numbers that are different to use in
             * calculating half the difference between each position as delta.
             */
            if (count < 2 )
              {
                count ++;
              }
            if (count == 2)
              {
                if ((std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9) && (std::fabs(longitudes_iso[0] - longitudes_iso[1]) > 1e-9))
                  {
                    /** Stop program if file formatted incorrectly.
		             */
                    std::cout << ""<< std::endl;
                    throw std::ios_base::failure("Adiabatic boundary file not formatted correctly. " + adiabatic_boundary_file_name + "Make sure you have lat, lon, value with lat. or lon. varying.");
                  }
                if ((std::fabs(latitudes_iso[0] - latitudes_iso[1]) < 1e-9) && (std::fabs(longitudes_iso[0] - longitudes_iso[1]) < 1e-9))
                  {
                    /** Stop program if file formatted incorrectly.
		             */
                    std::cout << ""<< std::endl;
                    throw std::ios_base::failure("Adiabatic boundary file not formatted correctly. " + adiabatic_boundary_file_name + "Make sure you have lat, lon, value with lat. or lon. varying.");
                  }

                if (std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9)
                  {
                    /** Calculate delta as half the distance between points.
		             */
                    delta = std::fabs((0.5)*(latitudes_iso[0] - latitudes_iso[1]));
                    /** If flag is 0 then longitudes grouped and we calculate delta from latitudes
		             */
                    data_flag = 0;
                  }
                else
                  {
                    /** Calculate delta as half the distance between points.
		             */
                    delta = std::fabs((0.5)*(longitudes_iso[0] - longitudes_iso[1]));
                    /** If flag is 1 then latitudes are grouped and we calculate delta from longitudes
                     */
		    data_flag = 1;
                  }
                std::cout << ""<< std::endl;
                std::cout<<"Adiabatic boundary data interaval  delta = "<< delta << std::endl;
                std::cout << ""<< std::endl;
                std::cout<<"Resolution of input data in meters is approximately = "<< delta*111*2 << std::endl;
                std::cout << ""<< std::endl;
                count++;
              }
          }

    /** Calculate the number of unique longitudes or latitudes from the adiabatic boundary file.
     */
	double c,d;
        int count3;
        count3 = 0;
        c = 0;
        d = 0;

        if ( data_flag == 1 )
          {
            c = latitudes_iso[0];
            d = latitudes_iso[1];
            count3 = 2;

            while (c-d < 1e-9)
              {
                c = d;
                d = latitudes_iso[count3];
                count3++;
              }
          }

        if ( data_flag == 0 )
          {
            c = longitudes_iso[0];
            d = longitudes_iso[1];
            count3 = 2;

            while (c-d < 1e-9)
              {
                c = d;
                d = longitudes_iso[count3];
                count3++;
              }
           }

        std::cout << ""<< std::endl;
        std::cout<<"number of unique latitudes or longitudes in Adiabatic boundary file= "<< count3 - 1 << std::endl;
        number_coords_depth = count3-1;
        std::cout << ""<< std::endl;
      }

  }
}


namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(AdiabaticBoundary,
                                       "adiabatic boundary",
                                       "In subsection Adiabatic boundary isotherm we define "
                                       "the extent of the conductive heat "
                                       "equation is assumed to be 1400 C (1673.15 K) "
                                       "as previously used by Bird et al., 2008 "
                                       "and (Stamps et al., 2015). This assumption "
                                       "is consistent with the Schubert 2001 "
                                       "definition of the mechanical lithosphere.")
  }
}

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

    /** Here is class that compute initial temperature field based on a 
     *  user-defined adiabatic boundary. Below the adiabatic boundary the
     *  temperature increases adiabatically while above the adiabatic 
     *  boundary the temperature linearly increases from a surface temperature
     *  (273.15 K or 0 degree C) to an isotherm (1673.15 K or 1600 degree C)
     *  that is the adiabatic boundary.
     */
     template <int dim>
     class AdiabaticBoundary : public Interface<dim>
     {
       public:
    	  /** 
	   * Return the initial temperature as a function of position
    	   */
    	  virtual
          double initial_temperature (const Point<dim> &position) const;

    	  /**
    	   * declare the parameters that this class needs
    	   */
          static
          void declare_parameters (ParameterHandler &prm);
          /**
           * Read the parameters above from the parameter file
           */
          virtual
	  void parse_parameters (ParameterHandler &prm);

       private:
          std::vector<double>  latitudes_iso;
          std::vector<double>  longitudes_iso;
          std::vector<double>  depths_iso;
	  std::string adiabatic_boundary_file_name;
           std::string line;
          double delta;
          int data_flag;
          int number_coords_depth;

          /**
           * A function that reads adiabatic boundary depth form ascii data file and return the value of the depth for each position
           */
          double
          get_isotherm_depth (const double latitude,
                                    const double longitude) const;

          /**
           * Return latitude and longitude from cartesian x,y and z (wgs84)
           */
          std::pair<double, double>
          lat_long_from_xyz_wgs84(const Point<3> &pos) const;

         /**
	  * Return distance from center of the WGS84 to a point on the surface
	  */

	  double
          radius_wgs84(const double &theta) const;
     };
  }
}


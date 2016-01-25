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
#include <aspect/simulator.h>


namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

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
     template <int dim>
     class AdiabaticBoundary : public Interface<dim>
     {
       public:
    	  /** 
	       * Return the initial temperature as a function of position.
    	   */
    	  virtual
          double initial_temperature (const Point<dim> &position) const;

    	  /**
    	   * Declare the parameters that this class needs.
    	   */
          static
          void declare_parameters (ParameterHandler &prm);
          /**
           * Read the parameters above from the parameter file.
           */
          virtual
	      void parse_parameters (ParameterHandler &prm);

       private:
          std::vector<double>  latitudes_iso;
          std::vector<double>  longitudes_iso;
          std::vector<double>  depths_iso;
	      std::string isotherm_file_name;
          double isotherm_temperature;
          double surface_temperature;
          double delta;

          /**
           * A function that reads the depth of the user-defined
           * adiabatic boundary from an ascii data file and returns
           * the value of the depth for each position.
           */
          double
          get_isotherm_depth (const double latitude,
                              const double longitude) const;

          /**
           * Return latitude and longitude from x,y and z in the WGS84
           * Earth-Center Earth-Fixed (ECEF) reference frame.
           */
          std::pair<double, double>
          lat_long_from_xyz_WGS84(const Point<3> &pos) const;

          /**
	       * Return distance from the Earth's center to a point on the
	       * surface in the WGS84 ECEF reference frane.
	       */
	      double
          radius_WGS84(const double &theta) const;
     };
  }
}


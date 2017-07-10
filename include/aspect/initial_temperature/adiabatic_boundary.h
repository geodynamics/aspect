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
   along with ASPECT; see the file LICENSE.  If not see
   <http://www.gnu.org/licenses/>.
 */


#include <aspect/initial_temperature/interface.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
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
     * This plugin is developed by Tahiry Rajaonarison and D. Sarah Stamps.
     */
    template <int dim>
    class AdiabaticBoundary : public Interface<dim>, public SimulatorAccess<dim>
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
        std::string data_directory;
        double isotherm_temperature;
        double surface_temperature;
        double temperature_gradient;
        double delta;

        /**
         * A function that returns the isotherm depth for a given position.
         */
        double
        get_isotherm_depth (const double latitude,
                            const double longitude) const;
    };
  }
}


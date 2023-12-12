/*
  Copyright (C) 2014 - 2016 by the authors of the ASPECT code.

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



#ifndef _aspect_initial_temperature_lithosphere_rift_h
#define _aspect_initial_temperature_lithosphere_rift_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{


  namespace InitialTemperature
  {

    /**
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class LithosphereRift : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        LithosphereRift ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Return the initial temperature as a function of depth and
         * the local layer thicknesses.
         */
        virtual
        double temperature (const double depth,
                            const std::vector<double> thicknesses) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Surface temperature
         */
        double T0;

        /**
         * LAB isotherm temperature
         */
        double LAB_isotherm;

        /**
         * Wheter or not to take the polygon thicknesses as dominant, or to smooth them
         * gradually into rift areas.
         */
        bool blend_rift_and_polygon;

        /**
         * Vector for field heat production rates.
         */
        std::vector<double> heat_productivities;

        /**
         * Vector for thermal conductivities.
         */
        std::vector<double> conductivities;

        /**
         * Vector for field densities.
         */
        std::vector<double> densities;

        /**
         * Vector for the reference field thicknesses.
         */
        std::vector<double> thicknesses;

        /**
         * Vector for the reference field thicknesses away from the rift.
         */
        std::vector<std::vector<double> > polygon_thicknesses;

        /**
         * The standard deviation of the Gaussian amplitude of the lithospheric thicknesses.
         */
        double sigma_rift;
        double sigma_polygon;

        /**
         * The maximum amplitude of the Gaussian distribution of the lithospheric thicknesses
         * with distance from the rift axis. It should have values between -1 and 1, where positive
         * numbers represent a reduction in thickness and negative numbers an increase.
         */
        std::vector<double> A;

        /**
         * Whether or not to use a compensation depth for the temperature
         * up to which the LAB isotherm is prescribed. If yes, temperatures
         * in the sub-lithospheric mantle are set to the LAB isotherm up
         * to the depth held by the parameter compensation_depth.
         */
        bool use_compensation_depth;
        double compensation_depth;
    };
  }
}


#endif

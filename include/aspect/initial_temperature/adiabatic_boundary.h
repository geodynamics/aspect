/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_initial_temperature_adiabatic_boundary_h
#define _aspect_initial_temperature_adiabatic_boundary_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that computes an initial temperature field based
     * on a user-defined adiabatic boundary. It discretizes the model domain into
     * two regions separated by an isotherm boundary below which
     * the temperature increases adiabatically. Above the
     * user-defined isotherm boundary the temperature linearly
     * increases from a surface temperature (273.15 K or 0 degree C)
     * to the isotherm (1673.15 K or 1400 degree C). The user
     * defines the location of the isotherm boundary with an ascii
     * data file with the format defined in the ASPECT manual.
     *
     * This plugin is developed by Tahiry Rajaonarison, Emmanuel Njinju, and D. Sarah Stamps.
     */
    template <int dim>
    class AdiabaticBoundary : public Interface<dim>, public Utilities::AsciiDataBoundary<dim>
    {
      public:

        /**
         * Constructor.
         */
        AdiabaticBoundary ();

        void initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

        /**
         * Declare the parameters that this class needs.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters above from the parameter file.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:

        double isotherm_temperature;
        double surface_temperature;
        double temperature_gradient;
        types::boundary_id surface_boundary_id;

    };
  }
}
#endif

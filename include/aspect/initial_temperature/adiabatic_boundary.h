// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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

        /**
         * This variable is read from the parameter file through a parameter called 'Isotherm temperature'.
         */
        double isotherm_temperature;
        /**
         * This variable is read from the parameter file through a parameter called 'Surface temperature'.
         */
        double surface_temperature;
        /**
         * This variable is read from the parameter file through a parameter called 'Adiabatic temperature gradient'.
         */
        double temperature_gradient;
        types::boundary_id surface_boundary_id;

    };
  }
}
#endif

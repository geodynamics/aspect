// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2014-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------



#ifndef _aspect_initial_temperature_continental_geotherm_h
#define _aspect_initial_temperature_continental_geotherm_h

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
    class ContinentalGeotherm : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        ContinentalGeotherm ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * Return the initial temperature as a function of depth,
         * based on the solution of the steady-state heat conduction
         * differential equation.
         */
        double initial_temperature (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Surface temperature
         *
         * This variable is read from the parameter file through a parameter called 'Surface temperature'.
         */
        double T0;

        /**
         * Value of the isotherm defining the
         * Lithosphere-Asthenosphere boundary.
         * This variable is read from the parameter file through a parameter called 'Lithosphere-Asthenosphere boundary isotherm'.
         */
        double LAB_isotherm;

        /**
         * Vector for the thicknesses of the compositional fields
         * representing the layers 'upper_crust', 'lower_crust' and 'lithospheric_mantle'.
         *
         * This variable is read from the parameter file through a parameter called 'Layer thicknesses'.
         */
        std::vector<double> thicknesses;

        /**
         * Vector for the heat production rates of the different
         * compositional fields, read from parameter file.
         */
        std::vector<double> heat_productivities;

        /**
         * Vector for the thermal conductivities of the different
         * compositional fields, read from parameter file.
         */
        std::vector<double> conductivities;

        /**
         * Vector for the densities of the different compositional
         * fields, read from parameter file.
         */
        std::vector<double> densities;
    };
  }
}


#endif

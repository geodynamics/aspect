// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_initial_temperature_ascii_profile_h
#define _aspect_initial_temperature_ascii_profile_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/point.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that implements a prescribed temperature field determined from
     * a AsciiDataProfile input file.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class AsciiProfile : public Utilities::AsciiDataProfile<dim>, public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        AsciiProfile ();

        /**
         * Initialization function.
         */
        void initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataProfile<dim>::initialize;

        /**
         * Return the temperature at a given point of the domain.
         */
        double initial_temperature (const Point<dim> &p) const override;

        /**
         * Declare the parameters for the input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * The column index of the temperature in the data file.
         */
        unsigned int temperature_index;
    };
  }
}


#endif

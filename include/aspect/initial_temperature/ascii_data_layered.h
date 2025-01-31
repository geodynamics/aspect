/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_ascii_data_layered_h
#define _aspect_initial_temperature_ascii_data_layered_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that implements an initial temperature field determined from
     * AsciiData input files. Each file defines an isotherm as a grid of points.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class AsciiDataLayered : public Utilities::AsciiDataLayered<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiDataLayered ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataLayered<dim>::initialize;

        /**
         * Return the initial temperature as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_temperature (const Point<dim> &position) const override;

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
    };
  }
}


#endif

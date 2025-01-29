/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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

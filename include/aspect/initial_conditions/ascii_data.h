/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__initial_conditions_ascii_data_h
#define __aspect__initial_conditions_ascii_data_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that implements a prescribed temperature field
     * determined from a AsciiData input file.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        /**
         * Return the boundary temperature as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_temperature (const Point<dim> &position) const;

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
        parse_parameters (ParameterHandler &prm);
    };
  }
}


#endif
